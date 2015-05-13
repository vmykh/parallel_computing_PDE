# Parallel computing PDE
====================================

## MPI

To compile mpi, run:
```
mpicc explicit.c -o ./bin/mpiversion -std=c99 -lm -fopenmp
```

To run mpi:
```
mpirun -np 2 ./bin/mpiversion
```

##plotter.py
---------
Builds 3dPlot using matplotlib and numpy.

For installing:

```python
sudo apt-get install python-matplotlib
sudo apt-get install python-numpy
```

Data format:

```
x_min x_n step_x
t_min t_n step_t
Ui[t][x]
```
