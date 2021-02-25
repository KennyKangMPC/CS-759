!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The <tt>SCHEME</tt> zero-padding execution subroutine (DOUBLE PRECISION).
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
!> @param[in] arr Sub-grid that will be zero padded.
!> @param[in] sch <tt>SCHEME</tt> that is used for arr.
!> @param[inout] arrZP Zero-padded sub-grid.
!> @param[inout] schZP <tt>SCHEME</tt> that is made for arrZP.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DBLE_SBR(arr, sch, arrZP, schZP)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: arr(0:,0:)
  TYPE(SCHEME), INTENT(IN) :: sch
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arrZP(:,:)
  TYPE(SCHEME), INTENT(IN) :: schZP
  DOUBLE PRECISION, ALLOCATABLE :: arrZP1(:,:) !< arr padded along the first
  !! dimension.
  DOUBLE PRECISION, ALLOCATABLE :: arrZP1T(:,:) !< Globally redistributed
  !! arrZP1, used for padding along the second dimension.
  DOUBLE PRECISION, ALLOCATABLE :: arrZPT(:,:) !< arr padded along the first
  !! and second dimensions, transposed.
  INTEGER :: ierror !< Additional integer argument for Fortran MPI calls.
  INTEGER :: i !< Counter for DO loops.
  INTEGER, ALLOCATABLE :: SUBARRAYS(:,:) !< Subarrays for global data
  !! redistribution.

  ! Check for errors in user input.
  CALL EXECUTE_SCHEME_ZERO_PAD_DBLE_EH(sch, schZP)

  ! To create the zero-padded array, we will pad along the first dimension and
  ! then the second, by performing a global redistribution such that the dimension
  ! we are padding along is contiguous in memory.
  ! We will work without the overlap, then copy it to an array that has overlap.
  ALLOCATE(arrZP1(0:schZP%vSlabSize(0)-1,0:sch%vSlabSize(1)-1))
  arrZP1 = 0.0
  ! Copy the positive wavenumber entries to the padded array.
  DO i = 0, sch%vSlabSize(1)-1
     arrZP1(0:sch%vSlabSize(0)/2,i) = arr(0:sch%vSlabSize(0)/2,sch%vSlabInt(0)+i)
  END DO
  ! Copy the negative wavenumber entries to the padded array.
  DO i = 0, sch%vSlabSize(1)-1
     arrZP1(schZP%vSlabSize(0)-(sch%vSlabSize(0)-1)/2:schZP%vSlabSize(0)-1,i) &
          = arr(sch%vSlabSize(0)/2+1:sch%vSlabSize(0)-1,sch%vSlabInt(0)+i)
  END DO

  ! Now, we will redistribute the arrZP1 sub-grids such that each processor
  ! contains the entirity of the second dimension.
  ! Set up and perform global redistribution.
  ALLOCATE(SUBARRAYS(0:schZP%commSize-1,0:1))
  CALL SUBARRAY((/schZP%vSlabSize(0),sch%vSlabSize(1)/), 0, schZP%commSize, &
       schZP%datatype, SUBARRAYS(:,0))
  CALL SUBARRAY((/schZP%hSlabSize(0),sch%hSlabSize(1)/), 1, schZP%commSize, &
       schZP%datatype, SUBARRAYS(:,1))

  ALLOCATE(arrZP1T(0:schZP%hSlabSize(0)-1,0:sch%hSlabSize(1)-1))

  CALL MPI_ALLTOALLW(arrZP1, schZP%counts, schZP%displs, SUBARRAYS(:,0), &
       arrZP1T, schZP%counts, schZP%displs, SUBARRAYS(:,1), schZP%comm, &
       ierror) 

  ! Now, zero pad along the second dimension.
  ALLOCATE(arrZPT(0:schZP%hSlabSize(0)-1,0:schZP%hSlabSize(1)-1))
  arrZPT = 0.0
  ! Copy the positive wave number entries to the padded array.
  DO i = 0, sch%hSlabSize(1)/2
     arrZPT(:,i) = arrZP1T(:,i)
  END DO
  ! Copy the negative wavenumber entries to the padded array.
  DO i = 0, (sch%hSlabSize(1)-3)/2
     arrZPT(:,schZP%hSlabSize(1)-1-i) = arrZP1T(:,sch%hSlabSize(1)-1-i)
  END DO

  ! Transpose back into the final zero-padded array.
  ALLOCATE(arrZP(0:schZP%vSlabSize(0)-1,0:schZP%vSlabSize(1)-1))
  arrZP = 0.0
  ! Re-create the subarray datatypes for the reditribution.
  DO i = 0, schZP%commSize-1
     CALL MPI_TYPE_FREE(SUBARRAYS(i,0), ierror)
  END DO
  DO i = 0, schZP%commSize-1
     CALL MPI_TYPE_FREE(SUBARRAYS(i,1), ierror)
  END DO
  CALL SUBARRAY(schZP%hSlabSize, 1, schZP%commSize, schZP%datatype, &
       SUBARRAYS(:,0))
  CALL SUBARRAY(schZP%vSlabSize, 0, schZP%commSize, schZP%datatype, &
       SUBARRAYS(:,1))
  CALL MPI_ALLTOALLW(arrZPT, schZP%counts, schZP%displs, SUBARRAYS(:,0), &
       arrZP, schZP%counts, schZP%displs, SUBARRAYS(:,1), schZP%comm, &
       ierror)
    
  ! Deallocate locally allocated arrays.
  DEALLOCATE(arrZP1)
  DEALLOCATE(arrZP1T)

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
  !!    <li> is not initialized for spectral derivatives
  !!    </ol>
  !! <li> schZP...
  !!    <ol>
  !!    <li> is not initialized
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DBLE_EH(sch, schZP)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE


    TYPE(SCHEME), INTENT(IN) :: sch
    TYPE(SCHEME), INTENT(IN) :: schZP
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

    ! Check if schZP is not initialized.
    IF (.NOT. schZP%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme for zero-padding is not initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF
    

  END SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DBLE_EH
  
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
 
END SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DBLE_SBR


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The <tt>SCHEME</tt> zero-padding execution subroutine (DOUBLE COMPLEX).
!
!> @param[in] arr Sub-grid that will be zero padded.
!> @param[in] sch <tt>SCHEME</tt> that is used for arr.
!> @param[inout] arrZP Zero-padded sub-grid.
!> @param[inout] schZP <tt>SCHEME</tt> that is made for arrZP.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DCMPX_SBR(arr, sch, arrZP, schZP)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE COMPLEX, INTENT(IN) :: arr(0:,0:)
  TYPE(SCHEME), INTENT(IN) :: sch
  DOUBLE COMPLEX, ALLOCATABLE, INTENT(INOUT) :: arrZP(:,:)
  TYPE(SCHEME), INTENT(IN) :: schZP
  DOUBLE COMPLEX, ALLOCATABLE :: arrZP1(:,:) !< arr padded along the first
  !! dimension.
  DOUBLE COMPLEX, ALLOCATABLE :: arrZP1T(:,:) !< Globally redistributed
  !! arrZP1, used for padding along the second dimension.
  DOUBLE COMPLEX, ALLOCATABLE :: arrZPT(:,:) !< arr padded along the first
  !! and second dimensions, transposed.
  INTEGER :: ierror !< Additional integer argument for Fortran MPI calls.
  INTEGER :: i !< Counter for DO loops.
  INTEGER, ALLOCATABLE :: SUBARRAYS(:,:) !< Subarrays for global data
  !! redistribution.

  ! Check for errors in user input.
  CALL EXECUTE_SCHEME_ZERO_PAD_DCMPX_EH(sch, schZP)

  ! To create the zero-padded array, we will pad along the first dimension and
  ! then the second, by performing a global redistribution such that the dimension
  ! we are padding along is contiguous in memory.
  ! We will work without the overlap, then copy it to an array that has overlap.
  ALLOCATE(arrZP1(0:schZP%vSlabSize(0)-1,0:sch%vSlabSize(1)-1))
  arrZP1 = 0.0
  ! Copy the positive wavenumber entries to the padded array.
  DO i = 0, sch%vSlabSize(1)-1
     arrZP1(0:sch%vSlabSize(0)/2,i) = arr(0:sch%vSlabSize(0)/2,sch%vSlabInt(0)+i)
  END DO
  ! Copy the negative wavenumber entries to the padded array.
  DO i = 0, sch%vSlabSize(1)-1
     arrZP1(schZP%vSlabSize(0)-(sch%vSlabSize(0)-1)/2:schZP%vSlabSize(0)-1,i) &
          = arr(sch%vSlabSize(0)/2+1:sch%vSlabSize(0)-1,sch%vSlabInt(0)+i)
  END DO

  ! Now, we will redistribute the arrZP1 sub-grids such that each processor
  ! contains the entirity of the second dimension.
  ! Set up and perform global redistribution.
  ALLOCATE(SUBARRAYS(0:schZP%commSize-1,0:1))
  CALL SUBARRAY((/schZP%vSlabSize(0),sch%vSlabSize(1)/), 0, schZP%commSize, &
       schZP%datatype, SUBARRAYS(:,0))
  CALL SUBARRAY((/schZP%hSlabSize(0),sch%hSlabSize(1)/), 1, schZP%commSize, &
       schZP%datatype, SUBARRAYS(:,1))

  ALLOCATE(arrZP1T(0:schZP%hSlabSize(0)-1,0:sch%hSlabSize(1)-1))

  CALL MPI_ALLTOALLW(arrZP1, schZP%counts, schZP%displs, SUBARRAYS(:,0), &
       arrZP1T, schZP%counts, schZP%displs, SUBARRAYS(:,1), schZP%comm, &
       ierror) 

  ! Now, zero pad along the second dimension.
  ALLOCATE(arrZPT(0:schZP%hSlabSize(0)-1,0:schZP%hSlabSize(1)-1))
  arrZPT = 0.0
  ! Copy the positive wave number entries to the padded array.
  DO i = 0, sch%hSlabSize(1)/2
     arrZPT(:,i) = arrZP1T(:,i)
  END DO
  ! Copy the negative wavenumber entries to the padded array.
  DO i = 0, (sch%hSlabSize(1)-3)/2
     arrZPT(:,schZP%hSlabSize(1)-1-i) = arrZP1T(:,sch%hSlabSize(1)-1-i)
  END DO

  ! Transpose back into the final zero-padded array.
  ALLOCATE(arrZP(0:schZP%vSlabSize(0)-1,0:schZP%vSlabSize(1)-1))
  arrZP = 0.0
  ! Re-create the subarray datatypes for the reditribution.
  DO i = 0, schZP%commSize-1
     CALL MPI_TYPE_FREE(SUBARRAYS(i,0), ierror)
  END DO
  DO i = 0, schZP%commSize-1
     CALL MPI_TYPE_FREE(SUBARRAYS(i,1), ierror)
  END DO
  CALL SUBARRAY(schZP%hSlabSize, 1, schZP%commSize, schZP%datatype, &
       SUBARRAYS(:,0))
  CALL SUBARRAY(schZP%vSlabSize, 0, schZP%commSize, schZP%datatype, &
       SUBARRAYS(:,1))
  CALL MPI_ALLTOALLW(arrZPT, schZP%counts, schZP%displs, SUBARRAYS(:,0), &
       arrZP, schZP%counts, schZP%displs, SUBARRAYS(:,1), schZP%comm, &
       ierror)
    
  ! Deallocate locally allocated arrays.
  DEALLOCATE(arrZP1)
  DEALLOCATE(arrZP1T)

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
  !!    <li> is not initialized for spectral derivatives
  !!    </ol>
  !! <li> schZP...
  !!    <ol>
  !!    <li> is not initialized
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DCMPX_EH(sch, schZP)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE


    TYPE(SCHEME), INTENT(IN) :: sch
    TYPE(SCHEME), INTENT(IN) :: schZP
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

    ! Check if schZP is not initialized.
    IF (.NOT. schZP%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme for zero-padding is not initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF
    

  END SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DCMPX_EH
  
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
 
END SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DCMPX_SBR


