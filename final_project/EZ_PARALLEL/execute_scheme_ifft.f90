!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The inverse FFT execution subroutine (DOUBLE PRECISION).
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
!> @param[inout] subGrid The local sub-grid.
!> @param[in] kind Type of FFT to execute, one of FFT_1D_1, FFT_1D_2, or FFT_2D.
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE EXECUTE_SCHEME_IFFT_DBLE_SBR(subGrid, kind, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(INOUT) :: subGrid(0:,0:)
  INTEGER, INTENT(IN) :: kind
  TYPE(SCHEME), INTENT(IN) :: sch
  DOUBLE PRECISION, ALLOCATABLE :: subGridInt(:,:) !< Interior of the sub-grid.
  DOUBLE PRECISION, ALLOCATABLE :: subGridInt_T(:,:) !< Transpose of sub-grid
  !! interior.
  INTEGER :: i, j !< Counters for DO loops.
  INTEGER :: error_code !< Integer argument for error flag for
  !! <tt>MPI_ABORT</tt>.
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.

  ! Check for errors in user input.
  CALL EXECUTE_SCHEME_IFFT_DBLE_EH(kind, sch)

  SELECT CASE (kind)

  CASE(FFT_1D_1)

     ! FFT the columns.
     DO j = sch%vSlabInt(0), sch%vSlabInt(1)
        CALL DFFTB(sch%vSlabSizeOvlp(0), subGrid(:,j), sch%WSAVE1)
     END DO

     subGrid = sch%norm_1D_1 * subGrid

  CASE(FFT_1D_2)

     ! Transpose a la Dalcin et. al. 2019.
     ALLOCATE(subGridInt(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
     ALLOCATE(subGridInt_T(0:sch%hSlabSize(0)-1, 0:sch%hSlabSize(1)-1))
     subGridInt = subGrid(:,sch%vSlabInt(0):sch%vSlabInt(1))
     CALL MPI_ALLTOALLW(subGridInt, sch%counts, sch%displs, sch%SUBARRAYS(:,0), &
          subGridInt_T, sch%counts, sch%displs, sch%SUBARRAYS(:,1), sch%comm, &
          ierror)

     ! FFT the rows.
     DO i = 0, sch%hslabSize(0)-1
        CALL DFFTB(sch%hSlabSize(1), subGridInt_T(i,:), sch%WSAVE2)
     END DO

     ! Transpose and copy back.
     CALL MPI_ALLTOALLW(subGridInt_T, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,1), subGridInt, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,0), sch%comm, ierror)

     subGridInt = sch%norm_1D_2 * subGridInt
     
     subGrid(:,sch%vSlabInt(0): sch%vSlabInt(1)) = subGridInt
     

  CASE(FFT_2D)

     ! FFT the columns.
     DO j = sch%vSlabInt(0), sch%vSlabInt(1)
        CALL DFFTB(sch%vSlabSizeOvlp(0), subGrid(:,j), sch%WSAVE1)
     END DO

     ! Transpose a la Dalcin et. al. 2019.
     ALLOCATE(subGridInt(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
     ALLOCATE(subGridInt_T(0:sch%hSlabSize(0)-1, 0:sch%hSlabSize(1)-1))
     subGridInt = subGrid(:,sch%vSlabInt(0):sch%vSlabInt(1))
     
     CALL MPI_ALLTOALLW(subGridInt, sch%counts, sch%displs, sch%SUBARRAYS(:,0), &
          subGridInt_T, sch%counts, sch%displs, sch%SUBARRAYS(:,1), sch%comm, &
          ierror)

     ! FFT the rows.
     DO i = 0, sch%hslabSize(0)-1
        CALL DFFTB(sch%hSlabSize(1), subGridInt_T(i,:), sch%WSAVE2)
     END DO

     ! Transpose and copy back.
     CALL MPI_ALLTOALLW(subGridInt_T, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,1), subGridInt, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,0), sch%comm, ierror)

     subGridInt = sch%norm_2D * subGridInt

     subGrid(:,sch%vSlabInt(0): sch%vSlabInt(1)) = subGridInt

  CASE DEFAULT ! Should never be hit, should be caught by error handling.
     IF (sch%procID .EQ. 0) THEN
        PRINT *, "ERROR: FFT kind not supported."
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
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    <li> is not initialized for FFTs
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE EXECUTE_SCHEME_IFFT_DBLE_EH(kind, sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kind
    TYPE(SCHEME), INTENT(IN) :: sch
    INTEGER :: error_code !< Integer argument for error flag for
    !! <tt>MPI_ABORT</tt>.
    INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
    !! <tt>MPI</tt> subroutine calls.
    LOGICAL :: error_flag = .FALSE. !< Error flag to stop all processors.

    ! Check if kind is not supported
    IF ((kind .NE. FFT_1D_1) .AND. (kind .NE. FFT_1D_2) .AND. &
         (kind .NE. FFT_2D)) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Request FFT type is not supported."
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

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE EXECUTE_SCHEME_IFFT_DBLE_EH
  
END SUBROUTINE EXECUTE_SCHEME_IFFT_DBLE_SBR


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The inverse FFT execution subroutine (DOUBLE COMPLEX).
!
!> @param[inout] subGrid The local sub-grid whose boundary will be shared.
!> @param[in] kind Type of FFT to execute, one of FFT_1D_1, FFT_1D_2, or FFT_2D.
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE EXECUTE_SCHEME_IFFT_DCMPX_SBR(subGrid, kind, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(0:,0:)
  INTEGER, INTENT(IN) :: kind
  TYPE(SCHEME), INTENT(IN) :: sch
  DOUBLE COMPLEX, ALLOCATABLE :: subGridInt(:,:) !< Interior of the sub-grid.
  DOUBLE COMPLEX, ALLOCATABLE :: subGridInt_T(:,:) !< Transpose of sub-grid
  !! interior.
  INTEGER :: i, j !< Counters for DO loops.
  INTEGER :: error_code !< Integer argument for error flag for
  !! <tt>MPI_ABORT</tt>.
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.

  ! Check for errors in user input.
  CALL EXECUTE_SCHEME_IFFT_DCMPX_EH(kind, sch)

  SELECT CASE (kind)

  CASE(FFT_1D_1)

     ! FFT the columns.
     DO j = sch%vSlabInt(0), sch%vSlabInt(1)
        CALL ZFFTB(sch%vSlabSizeOvlp(0), subGrid(:,j), sch%WSAVE1)
     END DO

     subGrid = sch%norm_1D_1 * subGrid

  CASE(FFT_1D_2)

     ! Transpose a la Dalcin et. al. 2019.
     ALLOCATE(subGridInt(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
     ALLOCATE(subGridInt_T(0:sch%hSlabSize(0)-1, 0:sch%hSlabSize(1)-1))
     subGridInt = subGrid(:,sch%vSlabInt(0):sch%vSlabInt(1))
     CALL MPI_ALLTOALLW(subGridInt, sch%counts, sch%displs, sch%SUBARRAYS(:,0), &
          subGridInt_T, sch%counts, sch%displs, sch%SUBARRAYS(:,1), sch%comm, &
          ierror)

     ! FFT the rows.
     DO i = 0, sch%hslabSize(0)-1
        CALL ZFFTB(sch%hSlabSize(1), subGridInt_T(i,:), sch%WSAVE2)
     END DO

     ! Transpose and copy back.
     CALL MPI_ALLTOALLW(subGridInt_T, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,1), subGridInt, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,0), sch%comm, ierror)

     subGridInt = sch%norm_1D_2 * subGridInt

     subGrid(:,sch%vSlabInt(0): sch%vSlabInt(1)) = subGridInt
     

  CASE(FFT_2D)

     ! FFT the columns.
     DO j = sch%vSlabInt(0), sch%vSlabInt(1)
        CALL ZFFTB(sch%vSlabSizeOvlp(0), subGrid(:,j), sch%WSAVE1)
     END DO

     ! Transpose a la Dalcin et. al. 2019.
     ALLOCATE(subGridInt(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
     ALLOCATE(subGridInt_T(0:sch%hSlabSize(0)-1, 0:sch%hSlabSize(1)-1))
     subGridInt = subGrid(:,sch%vSlabInt(0):sch%vSlabInt(1))
     
     CALL MPI_ALLTOALLW(subGridInt, sch%counts, sch%displs, sch%SUBARRAYS(:,0), &
          subGridInt_T, sch%counts, sch%displs, sch%SUBARRAYS(:,1), sch%comm, &
          ierror)

     ! FFT the rows.
     DO i = 0, sch%hslabSize(0)-1
        CALL ZFFTB(sch%hSlabSize(1), subGridInt_T(i,:), sch%WSAVE2)
     END DO

     ! Transpose and copy back.
     CALL MPI_ALLTOALLW(subGridInt_T, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,1), subGridInt, sch%counts, sch%displs, &
          sch%SUBARRAYS(:,0), sch%comm, ierror)

     subGridInt = sch%norm_2D * subGridInt

     subGrid(:,sch%vSlabInt(0): sch%vSlabInt(1)) = subGridInt

  CASE DEFAULT ! Should never be hit, should be caught by error handling.
     IF (sch%procID .EQ. 0) THEN
        PRINT *, "ERROR: FFT kind not supported."
     END IF
     CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
  END SELECT
  
  

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine for the sub-grid boundary sharing subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> kind...
  !!    <ol> is not supported
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    <li> is not initialized for FFTs
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE EXECUTE_SCHEME_IFFT_DCMPX_EH(kind, sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kind
    TYPE(SCHEME), INTENT(IN) :: sch
    INTEGER :: error_code !< Integer argument for error flag for
    !! <tt>MPI_ABORT</tt>.
    INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
    !! <tt>MPI</tt> subroutine calls.
    LOGICAL :: error_flag = .FALSE. !< Error flag to stop all processors.

    ! Check if kind is not supported
    IF ((kind .NE. FFT_1D_1) .AND. (kind .NE. FFT_1D_2) .AND. &
         (kind .NE. FFT_2D)) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Request FFT type is not supported."
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

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE EXECUTE_SCHEME_IFFT_DCMPX_EH
  
END SUBROUTINE EXECUTE_SCHEME_IFFT_DCMPX_SBR



