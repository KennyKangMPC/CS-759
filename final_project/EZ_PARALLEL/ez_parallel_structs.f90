!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STRUCTS
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
! This file is part of EZ_PARALLEL.
!
! EZ_PARALLEL is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EZ_PARALLEL is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EZ_PARALLEL.  If not, see <https://www.gnu.org/licenses/>.
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> \brief The <tt>EZ_PARALLEL</tt> structures module.
!> Contains all <tt>EZ_PARALLEL</tt> derived datatypes and global parameters.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  PRIVATE

  !> @class SCHEME.
  !! The <tt>SCHEME</tt> derived datatype contains all information about
  !! the grid and the sub-grid decomposition of the grid, as well as information
  !! needed for FFTs.
  TYPE, PUBLIC :: SCHEME
     ! Grid/Global variables.
     INTEGER :: gridSize(0:1) !< Size in each dimension of the grid.
     DOUBLE PRECISION :: colSpc !< Physical spacing between columns of the grid.
     INTEGER :: comm !< Communicator that holds the grid.
     INTEGER :: commSize !< Number of processes in the MPI communicator.
     INTEGER :: datatype !< Datatype of the array.
     INTEGER :: ovlp !< Number of extra columns needed by each sub-grid to
     !! successfully step forward in time.
     INTEGER, ALLOCATABLE :: rowDcmpSizes(:) !< Number of rows in each sub-grid
     !! of the horizontal-slab decomposition, excluding overlap.
     INTEGER, ALLOCATABLE :: colDcmpSizes(:) !< Number of columns in each
     !! sub-grid of the vertical-slab decomposition, excluding overlap.
     INTEGER, ALLOCATABLE :: colDcmpSizesOvlp(:) !< Number of columns in each
     !! sub-grid of the vertical slab decomposition, including overlap.

     ! Sub-grid/Local variables.
     INTEGER :: procID !< Processor ID.
     INTEGER :: hSlabSize(0:1) !< Sizes in each dimension of the sub-grid in
     !! the horizontal-slab decomposition of the global array.
     INTEGER :: vSlabSize(0:1) !< Sizes in each dimension of the sub-grid in
     !! the vertical-slab decomposition of the global array, excluding overlap.
     INTEGER :: vSlabSizeOvlp(0:1) !< Sizes in each dimension of the sub-grid in
     !! the vertical-slab decomposition of the global array, including overlap.
     INTEGER :: vSlabInt(0:1) !< Column indices of the interior of the vertical
     !! slab.
     DOUBLE PRECISION :: colRef !< The physical position of the reference point
     !! in the dimension along the rows (corresponding to a column).
     INTEGER :: SEND_BOUNDARIES(0:1) !< MPI derived datatype for sending
     !! sub-grid boundaries to neightboring sub-grids (0 = left, 1 = right).
     INTEGER :: RECV_BOUNDARIES(0:1) !< MPI derived datatype for recieving
     !! sub-grid boundaries from neightboring sub-grids (0 = left, 1 = right).

     ! Additional variables for FFTs.
     DOUBLE PRECISION, ALLOCATABLE :: WSAVE1(:) !< Holds initialization info
     !! for DFFTPACK 1-D FFTS along first dimension.
     DOUBLE PRECISION, ALLOCATABLE :: WSAVE2(:) !< Holds initialization info
     !! for DFFTPACK 1-D FFTS along second dimension.
     DOUBLE PRECISION :: norm_1D_1 !< The normalization coefficient for possible
     !! 1-D FFTs along the first dimension.
     DOUBLE PRECISION :: norm_1D_2 !< The normalization coefficient for possible
     !! 1-D FFTs along the second dimension.
     DOUBLE PRECISION :: norm_2D !< The normalization coefficient for
     !! possible 2D FFTs.
     INTEGER, ALLOCATABLE :: SUBARRAYS(:,:) !< Holds the datatypes necessary to
     !! perform the transposition.
     INTEGER, ALLOCATABLE :: counts(:), displs(:) !< Arrays for use in global
     !! redistributiuon.

     ! Additional variables for spectral derivatives.
     DOUBLE COMPLEX, ALLOCATABLE :: wvNmbr1(:) !< Holds coefficients for spectral
     !! derivative along the first dimension.
     DOUBLE COMPLEX, ALLOCATABLE :: wvNmbr2(:) !< Holds coefficients for spectral
     !! derivative along the second dimension.

     ! Error handling variables.
     LOGICAL :: initScheme = .FALSE. !< Checks if <tt>SCHEME</tt> was created
     !! already.
     LOGICAL :: initFFT = .FALSE. !< Checks if FFTs for <tt>SCHEME</tt> were
     !! initialized already.
     LOGICAL :: initSpecDrv = .FALSE. !< Checks if spectral derivates for
     !! <tt>SCHEME</tt> were initialized already.
  END TYPE SCHEME

  PUBLIC :: FFT_1D_1
  INTEGER, PARAMETER :: FFT_1D_1 = 51 !< Flag to mark execution of 1D FFTs along
  !! the first dimension.
  PUBLIC :: FFT_1D_2
  INTEGER, PARAMETER :: FFT_1D_2 = 46 !< Flag to mark execution of 1D FFTs along
  !! the second dimension.
  PUBLIC :: FFT_2D
  INTEGER, PARAMETER :: FFT_2D = 95 !< Flag to mark execution of 2D FFTs.

  PUBLIC :: SPEC_DRV_1D_1
  INTEGER, PARAMETER :: SPEC_DRV_1D_1 = 23 !< Flag to mark execution of 1D
  !! spectral derivatives along the first dimension.
  PUBLIC :: SPEC_DRV_1D_2
  INTEGER, PARAMETER :: SPEC_DRV_1D_2 = 78 !< Flag to mark execution of 1D
  !! spectral derivatives along the second dimension.

  
END MODULE EZ_PARALLEL_STRUCTS
