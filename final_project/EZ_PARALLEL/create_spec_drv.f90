!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Generates a desired spectral derivative matrix.
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
!> @param[in] order1 Order of the derivative in the first direction.
!> @param[in] order2 Order of the derivative in the second dimension.
!> @param[inout] specDrv Allocatable array to hold the spectral derivative
!! matrix.
!> @param[in] sch <tt>SCHEME</tt> that holds information for grid
!! decomposition, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE CREATE_SPEC_DRV_SBR(order1, order2, specDrv, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: order1
  INTEGER, INTENT(IN) :: order2
  DOUBLE COMPLEX, ALLOCATABLE, INTENT(INOUT) :: specDrv(:,:)
  TYPE(SCHEME), INTENT(IN) :: sch
  INTEGER :: i, j !< Counter for DO loops.
  DOUBLE COMPLEX, ALLOCATABLE :: wvNmbrs1(:) !< Wavenumbers along the first
  !! dimension for the spectral derivative.
  DOUBLE COMPLEX, ALLOCATABLE :: wvNmbrs2(:) !< Wavenumbers along the second
  !! dimension for the spectral derivative.

  ! Check for errors in user input.
  CALL CREATE_SPEC_DRV_EH(order1, order2, sch)

  ! Get the wavenumbers along the first dimension.
  ALLOCATE(wvNmbrs1(0:sch%vSlabSizeOvlp(0)-1))
  wvNmbrs1 = (0.0, 0.0)
  IF (order1 .NE. 0) THEN
     wvNmbrs1 = sch%wvNmbr1

     ! If gridSize(0) even and derivative order odd, zero out highest wavenumber.
     IF ((MOD(sch%gridSize(0),2) .EQ. 0) .AND. (MOD(order1,2) .EQ. 1)) THEN
        wvNmbrs1(sch%gridSize(0)/2) = 0.0
     END IF

     wvNmbrs1 = wvNmbrs1**(DBLE(order1))

     ! Get rid of floating point artifacts.
     IF (MOD(order1, 2) .EQ. 0) THEN
        wvNmbrs1 = DBLE(wvNmbrs1)
     ELSE
        wvNmbrs1 = DCMPLX(0.0, DBLE((0.0, -1.0)*wvNmbrs1))
     END IF
     
  END IF

  ! Get the wavenumbers along the second dimension.
  ALLOCATE(wvNmbrs2(0:sch%vSlabSizeOvlp(1)-1))
  wvNmbrs2 = (0.0, 0.0)
  IF (order2 .NE. 0) THEN
     wvNmbrs2 = sch%wvNmbr2

     ! If gridSize(1) even and derivative order odd, zero out highest wavenumber.
     IF ((MOD(sch%gridSize(1),2) .EQ. 0) .AND. (MOD(order2,2) .EQ. 1)) THEN
        DO i = 0, sch%vSlabSizeOvlp(1)-1
           IF (ABS(wvNmbrs2(i)) .EQ. sch%gridSize(1)/2) THEN
              wvNmbrs2(i) = (0.0, 0.0)
           END IF
        END DO
     END IF

     wvNmbrs2 = wvNmbrs2**(DBLE(order2))

     ! Get rid of floating point artifacts.
     IF (MOD(order2, 2) .EQ. 0) THEN
        wvNmbrs2 = DBLE(wvNmbrs2)
     ELSE
        wvNmbrs2 = DCMPLX(0.0, DBLE((0.0, -1.0)*wvNmbrs2))
     END IF
  END IF

  ! Allocate and fill in the specDrv array.
  ALLOCATE(specDrv(0:sch%vSlabSizeOvlp(0)-1,0:sch%vSlabSizeOvlp(1)-1))
  DO j = 0, sch%vSlabSizeOvlp(1)-1
     DO i = 0, sch%vSlabSizeOvlp(0)-1
        specDrv(i,j) = wvNmbrs1(i) + wvNmbrs2(j)
     END DO
  END DO

  ! Deallocate locally allocated arrays.
  DEALLOCATE(wvNmbrs1)
  DEALLOCATE(wvNmbrs2)

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    <li> is not initialized for FFTs
  !!    <li> is not initialized for spectral derivatives
  !!    </ol>
  !! <li> order1...
  !!    <ol>
  !!    <li> is negative
  !!    </ol>
  !! <li> order 2...
  !!    <ol>
  !!    <li> is negative
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE CREATE_SPEC_DRV_EH(order1, order2, sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: order1
    INTEGER, INTENT(IN) :: order2
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

    ! Check if sch is not initialized for spectral derivatives.
    IF (.NOT. sch%initSpecDrv) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is not initialized for spectral derivatives."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if order1 is negative.
    IF (order1 .LT. 0) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Spectral derivative order1 is negative."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if order2 is negative.
    IF (order2 .LT. 0) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Spectral derivative order2 is negative."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE CREATE_SPEC_DRV_EH
 
END SUBROUTINE CREATE_SPEC_DRV_SBR

