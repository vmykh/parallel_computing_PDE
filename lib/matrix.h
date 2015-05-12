#pragma once

#include <stdlib.h>
#include <stdio.h>

#define allocate(type, size) (type*)malloc(sizeof(type)*size)
#define print_vector(vector, size, pattern) for(int i = 0; i < size; ++i) printf(pattern, vector[i]);

typedef struct _matrix{
  int      size;
  double*  b;
  double** A;
} matrix;

/**
* read_matrix
* @param  char*   filename - name of file with matrix
* @return matrix* pointer at the type matrix
* @brief reads matrix, the vector of free members and its size
*/
matrix* read_matrix(char * filename)
{
  FILE * file = fopen(filename, "r");

  if(file == NULL)
  {
    return NULL;
  }

  matrix* m = allocate(matrix, 1);
  
  fscanf(file, "%d", &(m->size));

  if(m->size < 2)
  {
    return NULL;
  }

  m->A = allocate(double*, m->size);

  for(int i = 0; i < m->size; ++i)
  {
    m->A[i] = allocate(double, m->size);

    for(int j = 0; j < m->size; ++j)
    {
      fscanf(file, "%lf", &(m->A[i][j]));
    }
  }

  m->b = allocate(double, m->size);
  for(int i = 0; i < m->size; ++i)
  {
    fscanf(file, "%lf", &(m->b[i]));
  }

  return m;
}

/**
* write_result
* @param char* filename - name of destination file
* @param double* result - vector of found results
* @param int size       - size of vector
* @return void
* @brief writes vector of results to file
*/
void write_result(char * filename, double* result, int size)
{
  FILE * file = fopen(filename, "w");

  for(int i = 0; i < size; ++i)
  {
    fprintf(file, "%lf\n", result[i]);
  }
}

double* tridiagonalmatrix_right_solve(matrix* mtr)
{
  double** A    = mtr->A;
  double*  b    = mtr->b;
  int      size = mtr->size;

  double* alphas = allocate(double, size);
  double* betas  = allocate(double, size);
  double* xs     = allocate(double, size);

  // calculate initial alpha and beta
  alphas[0] = 0;
  betas[0]  = 0;
  alphas[1] = - A[0][1] / A[0][0];
  betas[1]  = b[0] / A[0][0];

  double denominator;
  for(int i = 1; i < size - 1; ++i)
  {
    denominator = (A[i][i] + alphas[i]*A[i][i-1]);
    alphas[i + 1] = (-A[i][i+1]) / denominator;
    betas[i + 1]  = (b[i] - A[i][i-1]*betas[i]) / denominator;
  }

  xs[size-1] = (b[size-1] - A[size-1][size-2]*betas[size-1]) / (A[size-1][size-2]*alphas[size-1] + A[size-1][size-1]);

  for(int i = size - 2; i >= 0; --i)
  {
    xs[i] = alphas[i+1] * xs[i+1] + betas[i+1];
  }

  return xs;
}

double* tridiagonalmatrix_left_solve(matrix* mtr)
{
  double** A    = mtr->A;
  double*  b    = mtr->b;
  int      size = mtr->size;

  double* xies = allocate(double, size);
  double* etas = allocate(double, size);
  double* xs   = allocate(double, size);

  // calculate initial alpha and beta
  xies[0] = 0;
  etas[0] = 0;
  xies[size-1] = - A[size-1][size-2] / A[size-1][size-1];
  etas[size-1] = b[size-1] / A[size-1][size-1];

  double denominator;
  for(int i = size - 2; i > 0; --i)
  {
    denominator = (A[i][i] + xies[i+1]*A[i][i+1]);
    xies[i] = (-A[i][i-1]) / denominator;
    etas[i] = (b[i] - A[i][i+1]*etas[i+1]) / denominator;
  }

  xs[0] = (b[0] - A[0][1]*etas[1]) / (A[0][1]*xies[1] + A[0][0]);

  for(int i = 0; i < size - 1; ++i)
  {
    xs[i + 1] = xies[i+1] * xs[i] + etas[i+1];
  }

  return xs;
}

void calculate_alphas_and_betas(matrix* mtr, double* alphas, double* betas, int p)
{
  double** A    = mtr->A;
  double*  b    = mtr->b;
  int      size = mtr->size;

  // calculate initial alpha and beta
  alphas[0] = 0;
  betas[0]  = 0;
  alphas[1] = - A[0][1] / A[0][0];
  betas[1]  = b[0] / A[0][0];

  double denominator;
  for(int i = 1; i <= p; ++i)
  {
    denominator = (A[i][i] + alphas[i]*A[i][i-1]);
    alphas[i + 1] = (-A[i][i+1]) / denominator;
    betas[i + 1]  = (b[i] - A[i][i-1]*betas[i]) / denominator;
  }
}

void calculate_xies_and_etas(matrix* mtr, double* xies, double* etas, int p)
{
  double** A    = mtr->A;
  double*  b    = mtr->b;
  int      size = mtr->size;

  // calculate initial alpha and beta
  xies[0] = 0;
  etas[0]  = 0;
  xies[size-1] = - A[size-1][size-2] / A[size-1][size-1];
  etas[size-1] = b[size-1] / A[size-1][size-1];

  double denominator;
  for(int i = size - 2; i >= p; --i)
  {
    denominator = (A[i][i] + xies[i+1]*A[i][i+1]);
    xies[i] = (-A[i][i-1]) / denominator;
    etas[i] = (b[i] - A[i][i+1]*etas[i+1]) / denominator;
  }
}