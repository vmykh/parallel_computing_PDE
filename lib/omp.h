#pragma once

#include <omp.h>
#include "./matrix.h"

double* tridiagonalmatrix_parallel_solve(matrix* mtr)
{
  int p = mtr->size / 2 + (mtr->size % 2);
  double* alphas = allocate(double, mtr->size);
  double* betas  = allocate(double, mtr->size);
  double* xies   = allocate(double, mtr->size);
  double* etas   = allocate(double, mtr->size);
  double* xs     = allocate(double, mtr->size);

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      calculate_alphas_and_betas(mtr, alphas, betas, p);
    }
    #pragma omp section
    {
      calculate_xies_and_etas(mtr, xies, etas, p);
    }
  }

  xs[p] = (alphas[p+1]*etas[p+1] + betas[p+1]) / (1 - alphas[p+1]*xies[p+1]);

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      for(int i = p - 1; i >= 0; --i)
      {
        xs[i] = alphas[i+1] * xs[i+1] + betas[i+1];
      }
    }
    #pragma omp section
    {
      for(int i = p; i < mtr->size - 1; ++i)
      {
        xs[i + 1] = xies[i+1] * xs[i] + etas[i+1];
      }
    }
  }

  return xs;
}