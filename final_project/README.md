# Background

`EZ_PARALLEL` is a double-precision `MPI` module for Fortran created to help simplify the processes of parallelizing serial 2-D finite difference and spectral codes, by including pre-made subroutines for grid decomposition, adjacent sub-grid communication, and more. We cannot guarantee that this library will yield the most optimal results, but our hope is that people without the time to learn and/or experience with MPI can use this library to improve the scope of their work without many changes to their initial serial code.

`EZ_PARALLEL` was developed by Jason Turner under the guidance of Samuel Stechmann at the University of Wisconsin-Madison. The [distribution](https://github.com/jasonlturner/EZ_PARALLEL_project) contains [`DFFTPACK`](https://www.netlib.org/fftpack/dp.tgz), a double-precision clone of the popular fast Fourier transform (FFT) library [`FFTPACK`](https://www.netlib.org/fftpack/), on which `EZ_PARALLEL` depends. 

# Introduction

Similarly to the `plan`s of [`FFTW`](http://www.fftw.org/), `EZ_PARALLEL` is centered around a data structure called a `scheme` which contains information about the grid, grid-decomposition, FFT initialization, and more. A `scheme` can be created for a grid and used repeadtly through a code for multiple grids through use of the `CREATE_SCHEME` subroutine. For documentation detailing the various interfaces and subroutines, see the `doxygen_bin` subdirectory of the distribution.

For the purpose of explaining the module in a concise manner, we will use a prototypical example of a serial finite difference code with the following structure:

1. Initialize simulations parameters, e.g., the size of the array containing the domain (with allocating memory to that array), setting physical constants, etc.
2. Initialize the domain grid. This includes allocating memory to the array containing the domain and setting the initial condition.
3. Step forward in time. This includes calculating each subsequent time step and outputting the domain grid to a file.

This is not to say that every code must be in this format to use `EZ_PARALLEL`, but it is the basic form that is used in the examples in the [distribution](https://github.com/jasonlturner/EZ_PARALLEL_project).

# Scheme Creation

To use `EZ_PARALLEL`, the user must include the lines `USE EZ_PARALLEL_STRUCTS` and `USE EZ_PARALLEL` before the `IMPLICIT NONE` statment, in addition to initializing `MPI` before any calls to `EZ_PARALLEL` subroutines.

`scheme` initialization comes in three parts, each of which requires the former is completed. There is an additional `scheme` creation subroutine, which may be considered as a fourth part:
1. `create_scheme` is the general `scheme` creation subroutine, enabling grid decomposition and sub-grid boundary communication.
2. `create_scheme_fft` enables a `scheme` to be used in performing FFTs.
3. `create_scheme_spec_drv` enables a `scheme` to be used in performing spectral derivatives.
4. `create_scheme_zero_pad` creates a new `scheme` based on a previous `scheme` for a zero-padded version of the grid using the [3/2s-rule](http://www.math.jhu.edu/~feilu/notes/DealiasingFFT.pdf) for de-aliasing FFTs.

The first three of these are used to create a `scheme` in its entirity, if necessary for the application. The fourth creates a `scheme` for a zero-padded version of a grid associated with the original `scheme`.

A user can call the various members of a `scheme`, which is useful when parallelizing their serial code. The base members of a `scheme` (intialized with `create_scheme`) are:
	- `gridSize(0:1)`: The size in each dimension of the grid.
	- `colSpc`: The physical spacing between columns of the grid.
	- `comm`: The `MPI` communicator on which the grid resides.
	- `commSize`: The number of processes in `comm`.
	- `datatype`: The `MPI` datatype of the grid, one of `MPI_DOUBLE_PRECISION` or `MPI_DOUBLE_COMPLEX`.
	- `ovlp`: The overlap parameter of the finite difference method, i.e., the number of points on one side of a grid point needed to execute a method.
	- `rowDcmpSizes(0:commSize-1)`: The number of rows in each sub-grid of the horizontal-slab decomposition, excluding overlap.
	- `colDcmpSizes(0:commSize-1)`: The number of columns in each sub-grid of the vertical-slab decomposition, excluding overlap.
	- `colDcmpSizesOvlp(0:commSize-1)`: The number of columns in each sub-grid of the vertical-slab decomposition, including overlap.
  - `procID`: Rank of the local process in `comm`.
  - `hSlabSize(0:1)`: The size of the local sub-grid in the horizontal-slab decomposition, excluding overlap.
  - `vSlabSize(0:1)`: The size of the local sub-grid in the vertical-slab decomposition, excluding overlap.
  - `vSlabSizeOvlp(0:1)`: The size of the local sub-grid in the vertical-slab decomposition, including overlap.
  - `vSlabInt(0:1)`: The column indices of the interior of the local sub-grid in the vertical-slab decomposition.
  - `colRef`: The physical position of the local sub-grid reference point along the row.
  - `SEND_BOUNDARIES(0:1)`: `MPI` derived datatype corresponding to the interior parts of the local sub-grid that would be sent to neighboring sub-grids during a sub-grid boundary sharing communication.
  - `RECV_BOUNDARIES(0:1)`: `MPI` derived datatype corresponding to the boundary of the local sub-grid that would be recieved from neighboring sub-grids during a sub-grid boundary sharing communication.
  
The members of a `scheme` used for FFTs (initialized with `create_scheme_fft`) are:
	- `WSAVE1(gridSize(0))`: The array used in FFTs along the first dimension of the grid.
	- `WSAVE2(gridSize(1))`: The array used in FFTs along the second dimension of the grid.
	- `norm_1d_1`: The normalization coefficient for 1-D FFTs along the first dimension.
	- `norm_1d_2`: The normalization coefficient for 1-D FFTs along the second dimension.
	- `norm_2d`: The normalization coefficient for 2-D FFTs.
	- `SUBARRAYS(commSize)`: `MPI` `subarray`s used for the global redistribution of the grid in 2-D FFTs, as described [here](https://arxiv.org/pdf/1804.09536.pdf).
	- `counts(commSize)`: An array used in `MPI_ALLTOALLW` call in the global redistribution.
	- `displs(commSize)`: An array used in `MPI_ALLTOALLW` call in the global redistribution.

The members of a `scheme` used for spectral derivatives (initialized using `create_scheme_spec_drv`) are:
	- `wvNmbr1(vSlabSizeOvlp(0))`: The wavenumbers along the first dimension of the grid, e.g., `0, ..., gridSize(0)/2, -gridSize(0)/2, ..., -1` for `gridSize(0)` even and `0, ..., gridSize(0)/2, -(gridSize(0)-1)/2, ..., -1` for `gridSize(0)` odd.
	- `wvNmbr2(vSlabSizeOvlp(1))`: The wavenumbers along the second dimension of the grid corresponding to the local sub-grid.
	
# Sub-Grid Boundary Sharing and Scheme Execution Subroutines

Each local sub-grid consists of a vertical slab of the original grid, including an additional boundary region which contains the points from the interior of neighboring sub-grids needed to perform a time-step. To share sub-grid boundaries is to exchange information with neighboring sub-grids as to update the local sub-grid's boundary region(s). To share sub-grid boundaries, one simply needs to make a call to the `SHARE_SUBGRID_BDRY` subroutine with appropriate arguments.

There are several `scheme` execution subroutines, which perform FFTs, spectral derivatives, and zero-padding:
	- `EXECUTE_SCHEME_FFT` and `EXECUTE_SCHEME_IFFT`: The subroutines for performing FFTs and inverse FFTs on both `DOUBLE PRECISION` and `DOUBLE COMPLEX` grids. FFTs and not normalized, while inverse FFTs are normalized. There are three kinds of FFTs implemented, specified by the `kind` argument (note, the `kind` argument works for both FFTs and inverse FFTs):
	  + `FFT_1D_1`: 1-D FFT along the first dimension of the grid.
	  + `FFT_1D_2`: 1-D FFT along the second dimension of the grid.
	  + `FFT_2D`: 2-D FFT.
	- `EXECUTE_SCHEME_SPEC_DRV` and `EXECUTE_SCHEME_ISPEC_DRV`: The subroutines for performing spectral derivatives and inverse spectral derivatives on `DOUBLE COMPLEX` grids (the functionality for `DOUBLE PRECISION` grids has yet to be implemented). Note that the `subGrid` argument should already be in spectral space when calling this subroutine. There are two kinds of spectral derivatives implemented, specified by the `kind` argument (note, the `kind` argument works for both spectral derivatives and inverse spectral derivatives):
	  + `SPEC_DRV_1D_1`: 1-D spectral derivative along the first dimension.
	  + `SPEC_DRV_1D_2`: 1-D spectral derivative along the second dimension.
	Alternatively, one may use the `CREATE_SPEC_DRV` to generate an array to perform the appropriate spectral derivative through standard multiplication. This is recommended if there are numerous spectral derivatives needed, if spectral derivatives are needed along both dimensions, or if the spectral derivative operator needs to be changed for the specific purpose.
	- `EXECUTE_SCHEME_ZERO_PAD` and `EXECUTE_SCHEME_IZERO_PAD`: The zero-padding subroutines for `DOUBLE PRECISION` and `DOUBLE COMPLEX` grids. These require a `scheme` for both the original grid and the zero-padded grid.

# Destroying Schemes

A `scheme` may be destoryed by calling the `DESTROY_SCHEME` subroutine, which deallocates all arrays and sets many of the other members of the `scheme` to null values. Note that it does not change the `procID` member so that the user may use it following the destruction to specify which processes should execute `PRINT` statements (among other uses).

# Additional Subroutines

`EZ_PARALLEL` also includes subroutines for finding the maximum and minimum values of a `DOUBLE PRECISION` grid (`MAX_VAL` and `MIN_VAL`, respectively), which may be useful in adaptive time-stepping methods (among other uses).

# Notes on Usage

## Initializing Sub-Grids

The `CREATE_SCHEME` subroutine returns the column count and column reference position (assuming a constant column spacing) of the sub-grid, easing the update from serial code to parallel code. For example, if the serial code initializes the grid using
```
ALLOCATE(grid(0:rowCount-1,0:colCount-1))
DO j = 0, colCount-1
   DO i = 0, rowCount-1
      grid(i,j) = INITIAL_CONDITION(rowRef + i*rowSpc,colRef + j*colSpc)
   END DO
END DO
```
then initializing the sub-grid using `EZ_PARALLEL` would require
```
CALL CREATE_SCHEME(rowCount,colCount,colSpc,colRef,comm,mpiDatatype,ovlp,sch)
ALLOCATE(grid(0:rowCount-1,0:colCount-1))
DO j = 0, colCount-1
   DO i = 0, rowCount-1
      grid(i,j) = INITIAL_CONDITION(rowRef + i*rowSpc,colRef + j*colSpc)
   END DO
END DO
```
If the serial code uses Dirichlet boundary conditions such as
```
grid(0,:) = 0.0
grid(rowCount-1,:) = 0.0
grid(:,0) = 0.0
grid(:,colCount-1) = 0.0
```
then the parallel code would simply require one call the `SHARE_SUBGRID_BDRY` afterward
```
grid(0,:) = 0.0
grid(rowCount-1,:) = 0.0
grid(:,0) = 0.0
grid(:,colCount-1) = 0.0
CALL SHARE_SUBGRID_BDRY(grid,sch)
```
Note that the time-stepping routine would require an overlap of at least 1.

## Writing Output Files

When writing output files, one must take caution to give each processor a unique output file name. For example, if the serial code chooses and output file name using
```
WRITE(fileName,'(A,I0.8,A)') 'out_', timestep, '.csv'
```
then each processor will attempt to write to the same file. This may be resolved by using the `procID` member of the `scheme`, which is unique to each processor
```
WRITE(fileName,'(A,I0.8,A,I0.3,A)') 'out_', step, '_', sch%procID, '.csv'
```

## Usage in 1-D and 3-D  Codes

Although we have not used `EZ_PARALLEL` for anything besides 2-D codes, it is theoretically possible to use `EZ_PARALLEL` for 1-D and 3-D codes by allocating 1-D arrays as 2-D arrays with a second dimension of size 1 or by calling 2-D slices of 3-D arrays, as long as one only needs 1-D or 2-D FFTs. 

For example, the two-level quasigeostrophic (QG) example code utilizes `EZ_PARALLEL` on a 3-D array, consisting of two levels. Since we never require FFTs or spectral derivatives between levels, we simply call the `EZ_PARALLEL` subroutines with 2-D slices of the 3-D array.

# Contact Information

If you have any comments, questions, concerns, or ideas about this module, please contact Jason Turner at jlturner5@wisc.edu.

# License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
