!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The inverse spectral derivative execution subroutine (DOUBLE PRECISION).
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
!> @param[inout] subGrid The local sub-grid whose boundary will be shared.
!> @param[in] kind Type of spectral derivative to execute, one of SPEC_DRV_1D_1,
!! SPEC_DRV_1D_2.
!> @param[in] order Order of the spectral derivative.
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DBLE_SBR(subGrid, kind, order, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(INOUT) :: subGrid(0:,0:)
  INTEGER, INTENT(IN) :: kind
  INTEGER, INTENT(IN) :: order
  TYPE(SCHEME), INTENT(IN) :: sch
  INTEGER :: error_code !< Integer argument for error flag for
  !! <tt>MPI_ABORT</tt>.
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.

  ! Check for errors in user input.
  CALL EXECUTE_SCHEME_ISPEC_DRV_DBLE_EH(kind, order, sch)

  IF (sch%procID .EQ. 0) THEN
     PRINT *, "ERROR: Inverse spectral derivatives for DOUBLE PRECISION ", &
          "grids is not yet supported."
  END IF

  CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
  

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> kind...
  !!    <ol> is not supported
  !! <li> order...
  !!    <ol> is negative
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    <li> is not initialized for FFTs
  !!    <li> is not initialized for spectral derivatives
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DBLE_EH(kind, order, sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kind
    INTEGER, INTENT(IN) :: order
    TYPE(SCHEME), INTENT(IN) :: sch
    INTEGER :: error_code !< Integer argument for error flag for
    !! <tt>MPI_ABORT</tt>.
    INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
    !! <tt>MPI</tt> subroutine calls.
    LOGICAL :: error_flag = .FALSE. !< Error flag to stop all processors.

    ! Check if kind is not supported
    IF ((kind .NE. SPEC_DRV_1D_1) .AND. (kind .NE. SPEC_DRV_1D_2)) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Requested spectral derivative type is not supported."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if order is negative
    IF (order .LT. 0) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Requested order of spectral derivative is negative."
       END IF
       error_flag = .TRUE.
    END IF
    
    ! Check if sch is not initialized.
    IF (.NOT. sch%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is not initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if sch is not initialized for FFTs.
    IF (.NOT. sch%initFFT) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is not initialized for FFTs."
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

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DBLE_EH
  
END SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DBLE_SBR


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The spectral derivative execution subroutine (DOUBLE COMPLEX).
!
!> @param[inout] subGrid The local sub-grid whose boundary will be shared.
!> @param[in] kind Type of spectral derivative to execute, one of SPEC_DRV_1D_1,
!! SPEC_DRV_1D_2.
!> @param[in] order Order of the spectral derivative.
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DCMPX_SBR(subGrid, kind, order, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(0:,0:)
  INTEGER, INTENT(IN) :: kind
  INTEGER, INTENT(IN) :: order
  TYPE(SCHEME), INTENT(IN) :: sch
  DOUBLE COMPLEX, ALLOCATABLE :: wvNmbrs(:) !< Array to store wave numbers.
  INTEGER :: i, j !< Counters for DO loops.
  INTEGER :: error_code !< Integer argument for error flag for
  !! <tt>MPI_ABORT</tt>.
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.

  ! Check for errors in user input.
  CALL EXECUTE_SCHEME_ISPEC_DRV_DCMPX_EH(kind, order, sch)

  SELECT CASE (kind)

  CASE (SPEC_DRV_1D_1)

     ALLOCATE(wvNmbrs(0:sch%gridSize(0)-1))
     wvNmbrs = sch%wvNmbr1

     ! If gridSize(0) even and derivative order odd, zero out highest wavenumber.
     IF ((MOD(sch%gridSize(0),2) .EQ. 0) .AND. (MOD(order,2) .EQ. 1)) THEN
        wvNmbrs(sch%gridSize(0)/2) = 0.0
     END IF

     wvNmbrs = DCMPLX(1.0)/wvNmbrs
     DO i = 0, sch%gridSize(0)-1
        IF ((wvNmbrs(i) .EQ. wvNmbrs(i)-(1.0,0)) .OR. & ! check if get infinity
             (wvNmbrs(i) .NE. wvNmbrs(i))) THEN ! or NAN when dividing 1 by wavenumber
           wvNmbrs(i) = 0.0
        END IF
     END DO

     wvNmbrs = wvNmbrs**(DBLE(order))

     ! Get rid of floating point artifacts.
     IF (MOD(order, 2) .EQ. 0) THEN
        wvNmbrs = DBLE(wvNmbrs)
     ELSE
        wvNmbrs = DCMPLX(0.0, DBLE((0.0, -1.0)*wvNmbrs))
     END IF
        
     DO j = sch%vSlabInt(0), sch%vSlabInt(1)
        subGrid(:,j) = subGrid(:,j) * wvNmbrs
     END DO

     DEALLOCATE(wvNmbrs)

  CASE (SPEC_DRV_1D_2)

     ALLOCATE(wvNmbrs(0:sch%vSlabSizeOvlp(1)-1))
     wvNmbrs = sch%wvNmbr2
     ! If gridSize(0) even and derivative order odd, zero out highest wavenumber.
     IF ((MOD(sch%gridSize(0),2) .EQ. 0) .AND. (MOD(order,2) .EQ. 1)) THEN
        DO i = 0, sch%vSlabSizeOvlp(1)-1
           IF (ABS(wvNmbrs(i)) .EQ. sch%gridSize(1)/2) THEN
              wvNmbrs(i) = 0.0
           END IF
        END DO
     END IF

     wvNmbrs = DCMPLX(1.0)/wvNmbrs
     DO i = 0, sch%vSlabSizeOvlp(1)-1
        IF ((wvNmbrs(i) .EQ. wvNmbrs(i)-(1.0,0)) .OR. & ! check if get infinity
             (wvNmbrs(i) .NE. wvNmbrs(i))) THEN ! or NAN when dividing 1 by wavenumber
           wvNmbrs(i) = 0.0
        END IF
     END DO

     wvNmbrs = wvNmbrs**(DBLE(order))

     ! Get rid of floating point artifacts.
     IF (MOD(order, 2) .EQ. 0) THEN
        wvNmbrs = DBLE(wvNmbrs)
     ELSE
        wvNmbrs = DCMPLX(0.0, DBLE((0.0, -1.0)*wvNmbrs))
     END IF
     
        
     DO j = sch%vSlabInt(0), sch%vSlabInt(1)
        subGrid(:,j) = subGrid(:,j) * wvNmbrs(j)
     END DO

     DEALLOCATE(wvNmbrs)

  CASE DEFAULT

     IF (sch%procID .EQ. 0) THEN
        PRINT *, "ERROR: Requested spectral derivative type is not supported."
     END IF

     CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)

  END SELECT
  
CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> kind...
  !!    <ol> is not supported
  !! <li> order...
  !!    <ol> is negative
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    <li> is not initialized for FFTs
  !!    <li> is not initialized for spectral derivatives
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DCMPX_EH(kind, order, sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kind
    INTEGER, INTENT(IN) :: order
    TYPE(SCHEME), INTENT(IN) :: sch
    INTEGER :: error_code !< Integer argument for error flag for
    !! <tt>MPI_ABORT</tt>.
    INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
    !! <tt>MPI</tt> subroutine calls.
    LOGICAL :: error_flag = .FALSE. !< Error flag to stop all processors.

    ! Check if kind is not supported
    IF ((kind .NE. SPEC_DRV_1D_1) .AND. (kind .NE. SPEC_DRV_1D_2)) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Requested spectral derivative type is not supported."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if order is negative
    IF (order .LT. 0) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Requested order of spectral derivative is negative."
       END IF
       error_flag = .TRUE.
    END IF
    
    ! Check if sch is not initialized.
    IF (.NOT. sch%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is not initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if sch is not initialized for FFTs.
    IF (.NOT. sch%initFFT) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is not initialized for FFTs."
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

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DCMPX_EH
  
END SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DCMPX_SBR
