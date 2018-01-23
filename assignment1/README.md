# deqn - a simple 2D heat diffusion simulation

`deqn` is a program that solves the 2D diffusion equation on a structured mesh.
Three different schemes are provided:

1. A simple explicit update scheme,
2. An iterative (matrix-free) Jacobi scheme, and
3. A scheme that uses the HYPRE library to solve the system of linear equations.

## Building

- To build the first two schemes, just type `make`.
- To build `deqn` with HYPRE, set the variable `HYPRE_DIR` in the Makefile, and
  ensure that the `CXXFLAGS` variable defines `HAVE_HYPRE`.

## Running

Using your favourite MPI flavour:

    mpirun -n 4 ./deqn <input-file>

## Input Files

The file `test/square.in` demonstrates the supported input parameters, most importantly:

- `scheme <scheme>` can be used to select which scheme to use (`explicit`, `jacobi`, or `hypre`).
- `vis_frequency <n>` controls how often visualisation files are written out.
- `subregion <xmin> <ymin> <xmax> <ymax>` specifies the region of the problem domain that will be initially heated.
