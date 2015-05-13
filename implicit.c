#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "./lib/matrix.h"   //kroosh tridiagonal solver
// #include "./lib/omp.h"      //kroosh tridiagonal parallel solver

#define A_DIFF 1.0       //vmykh   
#define A_FUNC 1.0
#define B_FUNC 1.0

#define X_MIN 0.1
#define X_MAX 0.7
#define T_MIN 0.1
#define T_MAX 0.5

#define X_POINTS_AMOUNT 5

#define DEBUG_FILE_NAME "data/debug.dat"
#define FILE_NAME "data/result.dat"

#define NEWTON_METHOD_TOLERANCE 0.0001  //condition to end iteration in Newton method: 
					// max(vector_x_i - vector_X_i-1) < NEWTON_METHOD_TOLERANCE

#define MAX_SIGMA 0.1  //should be less than 0.5

#define allocate(type, size) (type*)malloc(sizeof(type) * size)
#define fill_array(arr, size, default_value) for (int _iqw_ = 0; _iqw_ < size; ++_iqw_) {arr[_iqw_] = default_value;}
#define square(x) (x * x)
#define cube(x) (x * x * x)

#define X_STEP ((X_MAX - X_MIN) / (double) X_POINTS_AMOUNT)
#define T_STEP (MAX_SIGMA * X_STEP * X_STEP)
#define T_POINTS_AMOUNT (int) ((T_MAX - T_MIN) / T_STEP)

#define ALPHA (A_DIFF / (X_STEP * X_STEP))

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

void write_matrix_to_file(double** matrix);

void copy_arr(double* src_arr, double* dest_arr, int size);

int main(int argc, char const *argv[])
{
  
  MPI_Init (&argc, &argv);
  	
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int PROCS_AMOUNT;
  MPI_Comm_size(MPI_COMM_WORLD, &PROCS_AMOUNT);



  if (rank == 0)
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
  else
  {
  	
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
	return sqrt(pow(x - A_FUNC, 2.0) / (4 * A_DIFF * (B_FUNC - t)));
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
	// at each iteration we are calculating values for j row
  for (int j = 1; j < T_POINTS_AMOUNT - 1; ++j)  //i for x axis, j for t axis
  {

  	#pragma omp parallel for
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
    	mx->A[mx->size - 1][mx->size - 1] = current_partial_derivative(matrix, X_POINTS_AMOUNT - 2, j);


      	delta_x = tridiagonalmatrix_right_solve(mx);


      	copy_arr(current_values, prev_values, matrix_size);
      	add_to_first_vector(current_values, delta_x, matrix_size);
      	delete_Matrix(mx);
	} while (!is_finish_condition(current_values, prev_values, matrix_size));

  }
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
	return ALPHA * (square(matrix[j][i]) - matrix[j][i] * matrix[j][i+1] + matrix[j][i] * matrix[j][i - 1]);
}

double current_partial_derivative(double** matrix, int i, int j)
{
	return 2 * ALPHA * matrix[j][i] * matrix[j][i - 1] - 6.0 * ALPHA * square(matrix[j][i]) +
			2.0 * ALPHA * matrix[j][i] * matrix[j][i+1] + 

			0.5 * ALPHA * square(matrix[j][i+1]) - ALPHA * matrix[j][i+1] * matrix[j][i-1] + 
			0.5 * ALPHA * square(matrix[j][i-1])   

			- 1 / T_STEP;
}

double next_partial_derivative(double** matrix, int i, int j)
{
	return ALPHA * square(matrix[j][i]) + ALPHA * matrix[j][i] * matrix[j][i+1] - ALPHA * matrix[j][i] * matrix[j][i-1];
}

double finite_difference_function(double** matrix, int i, int j)
{
	return ALPHA * ( square(matrix[j][i]) * matrix[j][i-1] - 2.0 * cube(matrix[j][i]) + square(matrix[j][i]) * matrix[j][i+1] )


	+ 0.5 * ALPHA * (  matrix[j][i] * square(matrix[j][i+1]) - 2.0 * matrix[j][i] * matrix[j][i+1] * matrix[j][i-1] 
		+ matrix[j][i] * square(matrix[j][i-1])  )

	-( matrix[j][i]/ T_STEP ) + matrix[j-1][i] / T_STEP;
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