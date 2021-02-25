!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The <tt>SCHEME</tt> FFT creation subroutine.
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
!> @param[inout] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE CREATE_SCHEME_FFT_SBR(sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  TYPE(SCHEME), INTENT(INOUT) :: sch

  ! Check for errors in user input.
  CALL CREATE_SCHEME_FFT_EH(sch)

  ! Create arrays needed for 1-D FFTs.
  ALLOCATE(sch%WSAVE1(4*sch%gridSize(0)+15))
  ALLOCATE(sch%WSAVE2(4*sch%gridSize(1)+15))
  IF (sch%datatype .EQ. MPI_DOUBLE_COMPLEX) THEN
     CALL ZFFTI(sch%gridSize(0), sch%WSAVE1)
     CALL ZFFTI(sch%gridSize(1), sch%WSAVE2)
  ELSEIF (sch%datatype .EQ. MPI_DOUBLE_PRECISION) THEN
     CALL DFFTI(sch%gridSize(0), sch%WSAVE1)
     CALL DFFTI(sch%gridSize(1), sch%WSAVE2)
  END IF

  ! Calculate normalization coefficients.
  sch%norm_1D_1 = DBLE(1/DBLE(sch%gridSize(0)))
  sch%norm_1D_2 = DBLE(1/DBLE(sch%gridSize(1)))
  sch%norm_2D = DBLE(1/DBLE(PRODUCT(sch%gridSize)))

  ! Create subarrays
  ALLOCATE(sch%SUBARRAYS(sch%commSize,0:1))
  CALL SUBARRAY(sch%vSlabSize, 0, sch%commSize, sch%datatype, sch%SUBARRAYS(:,0))
  CALL SUBARRAY(sch%hSlabSize, 1, sch%commSize, sch%datatype, sch%SUBARRAYS(:,1))

  ! Allocate arrays needed for ALLTOALLW communication.
  ALLOCATE(sch%counts(sch%commSize))
  sch%counts = 1
  ALLOCATE(sch%displs(sch%commSize))
  sch%displs = 0

  ! Set scheme as initialized for FFTs.
  sch%initFFT = .TRUE.

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine for the <tt>SCHEME</tt> creation subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    <li> is already initialized for FFTs
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE CREATE_SCHEME_FFT_EH(sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE


    TYPE(SCHEME), INTENT(IN) :: sch
    INTEGER :: error_code !< Integer argument for error flag for
    !! <tt>MPI_ABORT</tt>.
    INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
    !! <tt>MPI</tt> subroutine calls.
    LOGICAL :: error_flag = .FALSE. !< Error flag to stop all processors.

    ! Check if sch is not initialized.
    IF (.NOT. sch%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is already initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if sch is already initialized for FFTs.
    IF (sch%initFFT) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is already initialized for FFTs."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE CREATE_SCHEME_FFT_EH

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Creates a decomposition of a single dimension of the 2-D array.
  !> @param[IN] nelems Number of elements along the dimension of the array, >=0.
  !> @param[IN] nparts Number of parts to divide dimension into, >0.
  !> @param[IN] pidx Part index, >=0 and <nparts.
  !> @param[OUT] nelemspart Number of elements in part.
  !> @param[OUT] sidx Start index of part.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE DECOMPOSE(nelems, nparts, pidx, nelemspart, sidx)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nelems
    INTEGER, INTENT(IN) :: nparts
    INTEGER, INTENT(IN) :: pidx
    INTEGER, INTENT(OUT) :: nelemspart
    INTEGER, INTENT(OUT) :: sidx
    INTEGER :: q, r !< Temporary integer variables.

    q = nelems / nparts
    r = MOD(nelems, nparts)
    IF (r > pidx) THEN
       nelemspart = q + 1
       sidx = nelemspart * pidx
    ELSE
       nelemspart = q
       sidx = nelemspart * pidx + r
    END IF

  END SUBROUTINE DECOMPOSE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Creates the subarray datatypes for the local array.
  !> @param[IN] sizes Local sizes of array.
  !> @param[IN] axis Axis to partition, 0 <= axis < 2.
  !> @param[IN] nparts Number of parts, nparts > 0.
  !> @param[IN] datatype MPI datatype descriptor.
  !> @param[OUT] subarrays Subarray datatype descriptors.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE SUBARRAY(sizes, axis, nparts, datatype, subarrays)

    USE MPI

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: sizes(0:1)
    INTEGER, INTENT(IN) :: axis
    INTEGER, INTENT(IN) :: nparts
    INTEGER, INTENT(IN) :: datatype
    INTEGER, INTENT(OUT) :: subarrays(0:nparts-1)
    INTEGER :: subsizes(0:1) !< Sizes of the subarrays.
    INTEGER :: substarts(0:1) !< Start indices of the subarrays.
    INTEGER :: n, s, p !< Temporary integer variables.
    INTEGER :: ierror !< Additional integer for all Fortran MPI calls.

    subsizes = sizes
    substarts = 0

    DO p = 0, nparts-1
       CALL DECOMPOSE(sizes(axis), nparts, p, n, s)
       subsizes(axis) = n
       substarts(axis) = s
       CALL MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, substarts, &
            MPI_ORDER_FORTRAN, datatype, subarrays(p), ierror)
       CALL MPI_TYPE_COMMIT(subarrays(p), ierror)
    END DO

  END SUBROUTINE SUBARRAY
 
END SUBROUTINE CREATE_SCHEME_FFT_SBR




