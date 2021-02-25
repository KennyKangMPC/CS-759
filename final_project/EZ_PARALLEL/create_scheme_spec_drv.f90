!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The <tt>SCHEME</tt> zero-padding initialization subroutine (DOUBLE PRECISION).
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
!> @param[inout] sch <tt>SCHEME</tt> that holds information for grid
!! decomposition, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE CREATE_SCHEME_SPEC_DRV_SBR(sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  TYPE(SCHEME), INTENT(INOUT) :: sch
  INTEGER :: i !< Counter for DO loops.
  INTEGER :: minWvNmbr !< Minimum wavenumber.
  INTEGER :: stIdx !< Start index of second-dimension wavenumbers needed by
  !! sub-grid.
  DOUBLE COMPLEX, ALLOCATABLE :: wvNmbr2All(:) !< Array to store all wavenumbers
  !! along the second dimension.

  ! Check for errors in user input.
  CALL CREATE_SCHEME_SPEC_DRV_EH(sch)

  ! Create array for wavenumbers along the first dimension.
  ALLOCATE(sch%wvNmbr1(0:sch%gridSize(0)-1))
  DO i = 0, sch%gridSize(0)/2
     sch%wvNmbr1(i) = DCMPLX(0.0, i)
  END DO
  minWvNmbr = -1 * (sch%gridSize(0)-1)/2
  DO i = 0, (sch%gridSize(0)-3)/2
     sch%wvNmbr1(sch%gridSize(0)/2+1+i) = DCMPLX(0.0, minWvNmbr + i)
  END DO

  ! Create array for wavenumbers along the second dimension.
  ! We actually calculate all of the wavenumbers, then just store the values we
  ! need.
  ALLOCATE(wvNmbr2All(0:sch%gridSize(1)-1))
  DO i = 0, sch%gridSize(1)/2
     wvNmbr2All(i) = DCMPLX(0.0, i)
  END DO
  minWvNmbr = -1 * (sch%gridSize(1)-1)/2
  DO i = 0, (sch%gridSize(1)-3)/2
     wvNmbr2All(sch%gridSize(1)/2+1+i) = DCMPLX(0.0, minWvNmbr + i)
  END DO
  
  IF (sch%procID .EQ. 0) THEN
     stIdx = 0
  ELSE
     stIdx = SUM(sch%colDcmpSizes(0:sch%procID)) - sch%colDcmpSizes(sch%procID) &
          - sch%ovlp
  END IF

  ALLOCATE(sch%wvNmbr2(0:sch%vSlabSizeOvlp(1)-1))
  DO i = 0, sch%vSlabSizeOvlp(1)-1
     sch%wvNmbr2(i) = wvNmbr2All(stIdx+i)
  END DO
  
  ! Set scheme as initialized for FFTs.
  sch%initSpecDrv = .TRUE.

  ! Deallocate locally allocated arrays.
  DEALLOCATE(wvNmbr2All)

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
  !!    <li> is not initialized for FFTs
  !!    <li> is already initialized for spectral derivatives
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE CREATE_SCHEME_SPEC_DRV_EH(sch)

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

    ! Check if sch is not initialized for FFTs.
    IF (.NOT. sch%initFFT) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is already initialized for FFTs."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if sch is already initialized for spectral derivatives.
    IF (sch%initSpecDrv) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is already initialized for spectral derivatives."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE CREATE_SCHEME_SPEC_DRV_EH
 
END SUBROUTINE CREATE_SCHEME_SPEC_DRV_SBR

