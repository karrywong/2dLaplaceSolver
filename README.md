# 2D Laplace Solver

This is my implementation in Fortran 90 for solving homogeneous Laplace equation with Dirichlet boundary condition and the boundary is a parametrized smooth closed curve. This set of code are ongoing work and to be further developed to solve the same boundary value problem where the boundary curve has corners.

Instructions for Linux:
  1. Make sure that there are libraries BLAS and LAPACK. For example, install [OpenBLAS](https://www.openblas.net/).
  
  2. Put all files **legendre.f90**, **disc.f90**, **matrix.f90** in the same directory.

  3. Compile on terminal using `gfortran -lblas -llapack matrix.f90 disc.f90 legendre.f90` and run executable **a.out**.
 
We can select different test cases in the file **matrix.f90** by (un)commenting the corresponding lines of code.
