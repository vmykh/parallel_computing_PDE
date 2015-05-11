#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A_DIFF 1.0
#define N_FUNC 1.0
#define C1_FUNC 1.0
#define C2_FUNC 1.0

#define X_MIN 0.1
#define X_MAX 0.7
#define T_MIN 0.1
#define T_MAX 0.5

#define X_POINTS_AMOUNT 30

#define DEBUG_FILE_NAME "data/debug.dat"
#define FILE_NAME "data/result.dat"

#define MAX_SIGMA 0.05  //should be less than 0.5

#define allocate(type, size) (type*)malloc(sizeof(type) * size)
#define fill_array(arr, size, default_value) for (int _iqw_ = 0; _iqw_ < size; ++_iqw_) {arr[_iqw_] = default_value;}

#define X_STEP (X_MAX - X_MIN) / (double) X_POINTS_AMOUNT
#define T_STEP (MAX_SIGMA * X_STEP * X_STEP)
#define T_POINTS_AMOUNT (int) ((T_MAX - T_MIN) / T_STEP)

double** create_matrix(int N, int M);
double exact_solution_func(double x, double t);
void init_boundaries(double** matrix);
void solve_pde(double** matrix);
double calculate_next_layer_point(double** matrix, int j, int i);  //i for x axis, j for t axis

double approx_t_first_deriv(double** matrix, int j, int i);   //i for x axis, j for t axis
double approx_x_first_deriv(double** matrix, int j, int i);
double approx_x_second_deriv(double** matrix, int j, int i);

void write_matrix_to_file(double** matrix);

int main(int argc, char const *argv[])
{
	double** matrix = create_matrix(T_POINTS_AMOUNT, X_POINTS_AMOUNT);
	init_boundaries(matrix);
	solve_pde(matrix);
	write_matrix_to_file(matrix);

	printf("T_POINTS_AMOUNT: %d\n", T_POINTS_AMOUNT);

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
	return 1.0 / (N_FUNC * x - A_DIFF * C1_FUNC * C1_FUNC * t + C2_FUNC);
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

void solve_pde(double** matrix)
{
	for (int j = 0; j < T_POINTS_AMOUNT - 1; ++j)  //i for x axis, j for t axis
	{
		for (int i = 1; i < X_POINTS_AMOUNT - 1; ++i)
		{
			matrix[j+1][i] = calculate_next_layer_point(matrix, j, i);
		}
	}
}

double calculate_next_layer_point(double** matrix, int j, int i)   //  i for x axis, j for t axis
{
	return matrix[j][i] + A_DIFF * T_STEP * (
		- 1.0 / (matrix[j][i]	* matrix[j][i]) * approx_x_first_deriv(matrix, j, i) * approx_x_first_deriv(matrix, j, i)
		+ 1.0 / matrix[j][i] * approx_x_second_deriv(matrix, j, i)
	);
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