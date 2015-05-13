#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#include "./lib/matrix.h"   //kroosh tridiagonal solver
#include "./lib/omp.h"      //kroosh tridiagonal parallel solver

#define A_DIFF 1.0      //kroosh
#define N_FUNC 1.0
#define C1_FUNC 1.0
#define C2_FUNC 1.0

#define X_MIN 0.1
#define X_MAX 0.7
#define T_MIN 0.1
#define T_MAX 0.5

#define X_POINTS_AMOUNT 15

#define DEBUG_FILE_NAME "data/debug.dat"
#define FILE_NAME "data/result.dat"

#define NEWTON_METHOD_TOLERANCE 0.0001  //condition to end iteration in Newton method: 
					// max(vector_x_i - vector_X_i-1) < NEWTON_METHOD_TOLERANCE

#define MAX_SIGMA 0.1  //should be less than 0.5

// defined in matrix.h
// #define allocate(type, size) (type*)malloc(sizeof(type) * size)
// #define fill_array(arr, size, default_value) for (int _iqw_ = 0; _iqw_ < size; ++_iqw_) {arr[_iqw_] = default_value;}

#define square(x) (x * x)
#define cube(x) (x * x * x)

#define X_STEP ((X_MAX - X_MIN) / (double) X_POINTS_AMOUNT)
#define T_STEP (MAX_SIGMA * X_STEP * X_STEP)
#define T_POINTS_AMOUNT (int) ((T_MAX - T_MIN) / T_STEP)

#define PSI (T_STEP / (4 * X_STEP * X_STEP))

double** create_matrix(int N, int M);
double exact_solution_func(double x, double t);
void init_boundaries(double** matrix);
void solve_pde(double** matrix);

double approx_t_first_deriv(double** matrix, int j, int i);   //i for x axis, j for t axis
double approx_x_first_deriv(double** matrix, int j, int i);
double approx_x_second_deriv(double** matrix, int j, int i);

int is_finish_condition(double* v1, double* v2, int size);
void add_to_first_vector(double* v1, double* v2, int size);

double finite_difference_function(double** matrix, int i, int j);
double previous_partial_derivative(double** matrix, int i, int j);
double current_partial_derivative(double** matrix, int i, int j);
double next_partial_derivative(double** matrix, int i, int j);

double* encode_Matrix(Matrix* m);
Matrix* decode_Matrix(double* encoded, int matrix_size);
int get_encoded_matrix_size(Matrix* m);
void tridiagonal_mpi_solve_second_process();
double* tridiagonal_mpi_solve_first_process(Matrix* mtr);

void write_matrix_to_file(double** matrix);

void copy_arr(double* src_arr, double* dest_arr, int size);

int main(int argc, char const *argv[])
{
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_size != 2) {
    fprintf(stderr, "This taks can be executed in 2 processes only %s\n");
    MPI_Abort(MPI_COMM_WORLD, 1); 

    return NULL;
  }

  if(world_rank == 0)
  {
    double** matrix = create_matrix(T_POINTS_AMOUNT, X_POINTS_AMOUNT);

    printf("before init_boundaries\n");
    init_boundaries(matrix);
    
    printf("solve_pde\n");
    solve_pde(matrix);

    printf("write matrix\n");
    write_matrix_to_file(matrix);

    printf("T_POINTS_AMOUNT: %d\n", T_POINTS_AMOUNT);
  }

  if(world_rank == 1)
  {
    int stop_signal;

    while(1)
    {
      MPI_Recv(&stop_signal, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if(stop_signal)
      {
        break;
      }

      tridiagonal_mpi_solve_second_process();
    }
  }

  MPI_Finalize();

  return 0;
}

double** create_matrix(int N, int M)
{
  double** matrix_pointer = allocate(double*, N);
  for (int i = 0; i < N; ++i)
  {
    matrix_pointer[i] = allocate(double, M);
    fill_array(matrix_pointer[i], M, 0.0);
  }
  return matrix_pointer;
}

double exact_solution_func(double x, double t)
{
  return 1.0 / (N_FUNC * x - A_DIFF * C1_FUNC * C1_FUNC * t + C2_FUNC);   //kroosh
}

void init_boundaries(double** matrix)
{
  //init T_MIN line
  for (int i = 0; i < X_POINTS_AMOUNT; ++i)
  {
    matrix[0][i] = exact_solution_func(X_MIN + X_STEP * i, T_MIN);
  }

  //init X_MIN, X_MAX
  for(int j = 1; j < T_POINTS_AMOUNT; ++j)
  {
    matrix[j][0] = exact_solution_func(X_MIN, T_MIN + T_STEP * j);
    matrix[j][X_POINTS_AMOUNT - 1] = exact_solution_func(X_MAX, T_MIN + T_STEP * j);
  }
}

void solve_pde(double** matrix)   //solve using implicit method
{
  int stop_signal = 0;
	// at each iteration we are calculating values for j row
  for (int j = 1; j < T_POINTS_AMOUNT - 1; ++j)  //i for x axis, j for t axis
  {
    for (int i = 1; i < X_POINTS_AMOUNT - 1; ++i)   //starting Newton method
    {
    	matrix[j][i] = matrix[j-1][i];
    }

    int matrix_size = X_POINTS_AMOUNT - 2;

    double* prev_values = allocate(double, matrix_size);
    double* delta_x;
    double* current_values = &(matrix[j][1]);
    do
    {
    	Matrix* mx = allocate(Matrix, 1);
    	mx->size = matrix_size;
    	mx->b = allocate(double, matrix_size);

    	mx->A = allocate(double*, matrix_size);
    	for (int i = 0; i < matrix_size; ++i)
    	{
    		mx->A[i] = allocate(double, matrix_size);
    		fill_array(mx->A[i], mx->size, 0.0);
    	}

    	//initilize vector b in mx
    	#pragma omp parallel for
    	for (int i = 0; i < matrix_size; ++i)
    	{
    		mx->b[i] = -finite_difference_function(matrix, i+1, j);
    	}

      mx->A[0][0] = current_partial_derivative(matrix, 1, j);
      mx->A[0][1] = next_partial_derivative(matrix, 1, j);

    	#pragma omp parallel for
      	for(int i = 1; i < matrix_size - 1; ++i)
      	{
        	mx->A[i][i-1] = previous_partial_derivative(matrix, i + 1, j);
       		mx->A[i][i]   = current_partial_derivative(matrix, i + 1, j);
        	mx->A[i][i+1] = next_partial_derivative(matrix, i + 1, j);
     	}

      mx->A[mx->size - 1][mx->size - 2] = previous_partial_derivative(matrix, X_POINTS_AMOUNT - 2, j);
      	mx->A[mx->size - 1][mx->size - 2] = previous_partial_derivative(matrix, X_POINTS_AMOUNT - 2, j);
    	mx->A[mx->size - 1][mx->size - 1] = current_partial_derivative(matrix, X_POINTS_AMOUNT - 2, j);

      MPI_Send(&stop_signal, 1, MPI_INT, 1, 5, MPI_COMM_WORLD);
      delta_x = tridiagonal_mpi_solve_first_process(mx);
    	copy_arr(current_values, prev_values, matrix_size);
    	add_to_first_vector(current_values, delta_x, matrix_size);
    	delete_Matrix(mx);
	} while (!is_finish_condition(current_values, prev_values, matrix_size));
  }

  stop_signal = 1;
  MPI_Send(&stop_signal, 1, MPI_INT, 1, 5, MPI_COMM_WORLD);
}

double approx_x_first_deriv(double** matrix, int j, int i)
{
  return (matrix[j][i + 1] - matrix[j][i - 1]) / (2.0 * X_STEP);
}

double approx_x_second_deriv(double** matrix, int j, int i)
{
  return (matrix[j][i - 1] - 2.0 * matrix[j][i] + matrix[j][i + 1]) / (X_STEP * X_STEP);
}

void write_matrix_to_file(double** matrix)
{
  FILE* f = fopen(FILE_NAME, "w");
  FILE* debug = fopen(DEBUG_FILE_NAME, "w");

  //write meta info
  fprintf(f, "%lf %d %lf\n", X_MIN, X_POINTS_AMOUNT, X_STEP);
  fprintf(f, "%lf %d %lf\n", T_MIN, T_POINTS_AMOUNT, T_STEP);

  fprintf(debug, "%lf %d %lf\n", X_MIN, X_POINTS_AMOUNT, X_STEP);
  fprintf(debug, "%lf %d %lf\n", T_MIN, T_POINTS_AMOUNT, T_STEP);

  //write matrix
  for (int j = 0; j < T_POINTS_AMOUNT; ++j)
  {
    for (int i = 0; i < X_POINTS_AMOUNT; ++i)
    {
      fprintf(debug, "[x: %lf t: %lf u: %lf]\n", X_MIN + X_STEP*i, T_MIN + T_STEP*j, matrix[j][i]);
      fprintf(f, "%lf ", matrix[j][i]);
    }
    fprintf(f, "\n");
    fprintf(debug, "\n");
  }
}

double previous_partial_derivative(double** matrix, int i, int j)
{
  return 2 * (T_STEP / (4.0 * X_STEP * X_STEP)) * (3 * matrix[j][i+1] - matrix[j][i-1]);
}

double current_partial_derivative(double** matrix, int i, int j)
{
  return -3 * matrix[j][i]*matrix[j][i] + 2 * matrix[j-1][i] * matrix[j][i] - 8 * (T_STEP / (4.0 * X_STEP * X_STEP)) * matrix[j][i+1];
}

double next_partial_derivative(double** matrix, int i, int j)
{
  return (T_STEP / (4.0 * X_STEP * X_STEP)) * (6 * matrix[j][i+1] + 6 * matrix[j][i-1] - 8 * matrix[j][i]);
}

double finite_difference_function(double** matrix, int i, int j)
{
  return - (matrix[j][i]) * (matrix[j][i]) * (matrix[j][i])
         + matrix[j-1][i] * matrix[j][i] * matrix[j][i]
         - PSI * (matrix[j][i + 1]) * (matrix[j][i + 1])
         + 2 * PSI * matrix[j][i+1] * matrix[j][i-1]
         - PSI * matrix[j][i-1] * matrix[j][i-1] 
         + 4 * PSI * matrix[j][i+1] * matrix[j][i+1]
         - 8 * PSI * matrix[j][i] * matrix[j][i+1]
         + 4 * PSI * matrix[j][i+1] * matrix[j][i-1];
}

int is_finish_condition(double* v1, double* v2, int size)
{
	for (int i = 0; i < size; ++i)
	{
		if (fabs(v1[i] - v2[i]) > NEWTON_METHOD_TOLERANCE )
		{
			return 0;
		}
	}
	return 1;
}

void copy_arr(double* src_arr, double* dest_arr, int size)
{
	for (int i = 0; i < size; ++i)
	{
		dest_arr[i] = src_arr[i];
	}
}

void add_to_first_vector(double* v1, double* v2, int size)
{
  for(int i = 0; i < size; ++i)
  {
    v1[i] += v2[i];
  }
}

double* encode_Matrix(Matrix* m)
{
  double* encoded = allocate(double, get_encoded_matrix_size(m));
  int iters = m->size - 1;

  int upper_diag_shift = 0;
  int middle_diag_shift = m->size - 1;
  int lower_diag_shift = middle_diag_shift + m->size;
  int b_vector_shift = lower_diag_shift + m->size - 1;

  for (int i = 0; i < iters; ++i)
  {
    encoded[upper_diag_shift + i] = m->A[i][i + 1];
    encoded[middle_diag_shift + i] = m->A[i][i];
    encoded[lower_diag_shift + i] = m->A[i + 1][i];
    encoded[b_vector_shift + i] = m->b[i];
  }
  encoded[middle_diag_shift + iters] = m->A[iters][iters];
  encoded[b_vector_shift + iters] = m->b[iters];

  return encoded;
}

Matrix* decode_Matrix(double* encoded, int matrix_size)
{
  // Matrix* m = allocate(Matrix, 1);
  // m->size = matrix_size;
  // m->b = allocate(double, matrix_size);
  // m->A = allocate(double*, matrix_size);
  // for (int i = 0; i < count; ++i)
  // {
  //   

  Matrix* m = create_Matrix(matrix_size);

  int upper_diag_shift = 0;
  int middle_diag_shift = m->size - 1;
  int lower_diag_shift = middle_diag_shift + m->size;
  int b_vector_shift = lower_diag_shift + m->size - 1;

  int iters = m->size - 1;
  for (int i = 0; i < iters; ++i)
  {
    m->A[i][i + 1] = encoded[upper_diag_shift + i];
    m->A[i][i] = encoded[middle_diag_shift + i];
    m->A[i + 1][i] = encoded[lower_diag_shift + i];
    m->b[i] = encoded[b_vector_shift + i];
  }

  m->A[iters][iters] = encoded[middle_diag_shift + iters];
  m->b[iters] = encoded[b_vector_shift + iters];

  return m;
}

int get_encoded_matrix_size(Matrix* m)
{
	return m->size * 4 - 2;
}

double* tridiagonal_mpi_solve_first_process(Matrix* mtr)
{
  int p = mtr->size / 2 + (mtr->size % 2);
  double* alphas = allocate(double, mtr->size);
  double* betas  = allocate(double, mtr->size);
  double* xies   = allocate(double, mtr->size);
  double* etas   = allocate(double, mtr->size);
  double* xs     = allocate(double, mtr->size);
  double xie_p, eta_p, alpha_p, beta_p;

  int encoded_size = get_encoded_matrix_size(mtr);
  MPI_Send(&encoded_size, 1, MPI_INT, 1, 9, MPI_COMM_WORLD);
  MPI_Send(encode_Matrix(mtr), encoded_size, MPI_DOUBLE, 1, 10, MPI_COMM_WORLD);
   
  calculate_alphas_and_betas(mtr, alphas, betas, p);
  
  MPI_Send(&(alphas[p+1]), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
  MPI_Send(&(betas[p+1]), 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

  alpha_p = alphas[p+1];
  beta_p = betas[p+1];

  MPI_Recv(&xie_p, 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&eta_p, 1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  xs[p] = (alpha_p*eta_p + beta_p) / (1 - alpha_p*xie_p);


  for(int i = p - 1; i >= 0; --i)
  {
    xs[i] = alphas[i+1] * xs[i+1] + betas[i+1];
  }

  double* xs_temp = allocate(double, mtr->size);
  MPI_Recv(xs_temp, mtr->size, MPI_DOUBLE, 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  for(int i = p; i < mtr->size; ++i)
  {
    xs[i] = xs_temp[i];
  }

  return xs;
}

void tridiagonal_mpi_solve_second_process()
{
  Matrix *mtr = allocate(Matrix, 1);

  int encoded_size;
  MPI_Recv(&encoded_size, 1, MPI_INT, 0, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  double* mtr_raw = allocate(double, encoded_size);
  MPI_Recv(mtr_raw, encoded_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  mtr = decode_Matrix(mtr_raw, (encoded_size + 2) / 4);

  int p = mtr->size / 2 + (mtr->size % 2);
  double* alphas = allocate(double, mtr->size);
  double* betas  = allocate(double, mtr->size);
  double* xies   = allocate(double, mtr->size);
  double* etas   = allocate(double, mtr->size);
  double* xs     = allocate(double, mtr->size);
  double xie_p, eta_p, alpha_p, beta_p;
  
  calculate_xies_and_etas(mtr, xies, etas, p);

  MPI_Recv(&alpha_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&beta_p, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Send(&(xies[p+1]), 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
  MPI_Send(&(etas[p+1]), 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

  xie_p = xies[p+1];
  eta_p = etas[p+1];

  xs[p] = (alpha_p*eta_p + beta_p) / (1 - alpha_p*xie_p);

  for(int i = p; i < mtr->size - 1; ++i)
  {
    xs[i + 1] = xies[i+1] * xs[i] + etas[i+1];
  }

  MPI_Send(xs, mtr->size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
}