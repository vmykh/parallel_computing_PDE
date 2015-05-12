##include "./matrix.h"

double* tridiagonal_mpi_solve(Matrix* mtr)
{
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_size != 2) {
    fprintf(stderr, "This taks can be executed in 2 processes only %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1); 

    return NULL;
  }

  int p = mtr->size / 2 + (mtr->size % 2);
  double* alphas = allocate(double, mtr->size);
  double* betas  = allocate(double, mtr->size);
  double* xies   = allocate(double, mtr->size);
  double* etas   = allocate(double, mtr->size);
  double* xs     = allocate(double, mtr->size);
  double xie_p, eta_p, alpha_p, beta_p;
   
  if(world_rank == 0)
  {
    calculate_alphas_and_betas(mtr, alphas, betas, p);
    
    MPI_Send(&(alphas[p+1]), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    MPI_Send(&(betas[p+1]), 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

    alpha_p = alphas[p+1];
    beta_p = betas[p+1];

    MPI_Recv(&xie_p, 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&eta_p, 1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  {
    calculate_xies_and_etas(mtr, xies, etas, p);

    MPI_Recv(&alpha_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&beta_p, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Send(&(xies[p+1]), 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    MPI_Send(&(etas[p+1]), 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

    xie_p = xies[p+1];
    eta_p = etas[p+1];
  }

  xs[p] = (alpha_p*eta_p + beta_p) / (1 - alpha_p*xie_p);


  if(world_rank == 0)
  {
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

    printf("result:\n");
    print_vector(xs, mtr->size, "%lf ");
    printf("\n");
  }
  else
  {
    for(int i = p; i < mtr->size - 1; ++i)
    {
      xs[i + 1] = xies[i+1] * xs[i] + etas[i+1];
    }

    MPI_Send(xs, mtr->size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
  }

  return xs;
}