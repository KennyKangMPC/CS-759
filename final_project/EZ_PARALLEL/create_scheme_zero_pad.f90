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
!> @param[in] sch <tt>SCHEME</tt> that is used for arr.
!> @param[inout] schZP <tt>SCHEME</tt> that is made for arrZP.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE CREATE_SCHEME_ZERO_PAD_DBLE_SBR(sch, schZP)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  TYPE(SCHEME), INTENT(IN) :: sch
  TYPE(SCHEME), INTENT(INOUT) :: schZP
  INTEGER :: ierror !< Additional integer argument for Fortran MPI calls.
  INTEGER :: stdCount !< Integer for holding standard row/col count in
  !! decompositions.
  INTEGER :: i !< Counter for DO loops.
  INTEGER :: minWvNmbr !< Minimum wavenumber.
  INTEGER :: stIdx !< Start index of second-dimension wavenumbers needed by
  !! sub-grid.
  DOUBLE COMPLEX, ALLOCATABLE :: wvNmbr2All(:) !< Array to store all wavenumbers
  !! along the second dimension.

  ! Check for errors in user input.
  CALL CREATE_SCHEME_ZERO_PAD_DBLE_EH(sch, schZP)

  ! Create a new, fully initialized scheme for the zero-padded grid.
  ! Primary initialization.
  schZP%gridSize = 3*(sch%gridSize/2)
  schZP%colSpc = sch%colSpc
  schZP%comm = sch%comm
  schZP%commSize = sch%commSize
  schZP%datatype = sch%datatype
  schZP%ovlp = sch%ovlp

  ! Set the overlap to zero, although we will add compatibility for it later.
  schZP%ovlp = 0

  schZP%procID = sch%procID
  schZP%colRef = -1 ! arrZP takes values from arr, so colRef is not calculated.

  ! Fill in decomposition arrays for schZP.
  ALLOCATE(schZP%rowDcmpSizes(0:schZP%commSize-1))
  stdCount = schZP%gridSize(0)/schZP%commSize
  DO i = 0, schZP%CommSize-1
     IF (i .LT. MOD(schZP%gridSize(0), schZP%commSize)) THEN
        schZP%rowDcmpSizes(i) = stdCount + 1
     ELSE
        schZP%rowDcmpSizes(i) = stdCount
     END IF
  END DO

  ALLOCATE(schZP%colDcmpSizes(0:schZP%commSize-1))
  stdCount = schZP%gridSize(1)/schZP%commSize
  DO i = 0, schZP%CommSize-1
     IF (i .LT. MOD(schZP%gridSize(1), schZP%commSize)) THEN
        schZP%colDcmpSizes(i) = stdCount + 1
     ELSE
        schZP%colDcmpSizes(i) = stdCount
     END IF
  END DO

  ALLOCATE(schZP%colDcmpSizesOvlp(0:schZP%commSize-1))
  schZP%colDcmpSizesOvlp(0) = schZP%colDcmpSizes(0) + schZP%ovlp
  DO i = 1, schZP%commSize-2
     schZP%colDcmpSizesOvlp(i) = schZP%colDcmpSizes(i) + 2*schZP%ovlp
  END DO
  schZP%colDcmpSizesOvlp(schZP%commSize-1) = schZP%colDcmpSizes(schZP%commSize-1) &
       + schZP%ovlp

  ! Set slab sizes.
  schZP%hSlabSize(0) = schZP%rowDcmpSizes(schZP%procID)
  schZP%hSlabSize(1) = schZP%gridSize(1)

  schZP%vSlabSize(0) = schZP%gridSize(0)
  schZP%vSlabSize(1) = schZP%colDcmpSizes(schZP%procID)

  schZP%vSlabSizeOvlp(0) = schZP%gridSize(0)
  schZP%vSlabSizeOvlp(1) = schZP%colDcmpSizesOvlp(schZP%procID)

  ! Store indices of interior of vertical slab for later.
  IF (schZP%procID .EQ. 0) THEN
     schZP%vSlabInt(0) = 0
     schZP%vSlabInt(1) = schZP%vSlabSizeOvlp(1) - schZP%ovlp - 1
  ELSE IF (schZP%procId .EQ. schZP%commSize-1) THEN
     schZP%vSlabInt(0) = schZP%ovlp
     schZP%vSlabInt(1) = schZP%vSlabSizeOvlp(1) - 1
  ELSE
     schZP%vSlabInt(0) = schZP%ovlp
     schZP%vSlabInt(1) = schZP%vSlabSizeOvlp(1) - schZP%ovlp - 1
  END IF

  ! Set SEND boundary datatypes.
  IF ((schZP%commSize .GT. 1) .AND. (schZP%ovlp .GT. 0)) THEN
     CALL MPI_TYPE_CREATE_SUBARRAY(2, schZP%vSlabSizeOvlp, &
          (/schZP%vSlabSizeOvlp(0), schZP%ovlp/), (/0, schZP%ovlp/), &
          MPI_ORDER_FORTRAN, schZP%datatype, schZP%SEND_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_CREATE_SUBARRAY(2, schZP%vSlabSizeOvlp, &
          (/schZP%vSlabSizeOvlp(0), schZP%ovlp/), &
          (/0, schZP%vSlabSizeOvlp(1)-2*schZP%ovlp/), MPI_ORDER_FORTRAN, &
          schZP%datatype, schZP%SEND_BOUNDARIES(1), ierror)
     CALL MPI_TYPE_COMMIT(schZP%SEND_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_COMMIT(schZP%SEND_BOUNDARIES(1), ierror)

     ! Set RECV boundary datatypes.
     CALL MPI_TYPE_CREATE_SUBARRAY(2, schZP%vSlabSizeOvlp, &
          (/schZP%vSlabSizeOvlp(0), schZP%ovlp/), (/0, 0/), &
          MPI_ORDER_FORTRAN, schZP%datatype, schZP%RECV_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_CREATE_SUBARRAY(2, schZP%vSlabSizeOvlp, &
          (/schZP%vSlabSizeOvlp(0), schZP%ovlp/), &
          (/0, schZP%vSlabSizeOvlp(1)-schZP%ovlp/), MPI_ORDER_FORTRAN, &
          schZP%datatype, schZP%RECV_BOUNDARIES(1), ierror)
     CALL MPI_TYPE_COMMIT(schZP%RECV_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_COMMIT(schZP%RECV_BOUNDARIES(1), ierror)
  END IF

  ! Set schZP as created.
  schZP%initScheme = .TRUE.

  ! FFT initialization.
  ! Create arrays needed for 1-D FFTs.
  ALLOCATE(schZP%WSAVE1(4*schZP%gridSize(0)+15))
  ALLOCATE(schZP%WSAVE2(4*schZP%gridSize(1)+15))
  IF (schZP%datatype .EQ. MPI_DOUBLE_COMPLEX) THEN
     CALL ZFFTI(schZP%gridSize(0), schZP%WSAVE1)
     CALL ZFFTI(schZP%gridSize(1), schZP%WSAVE2)
  ELSEIF (schZP%datatype .EQ. MPI_DOUBLE_PRECISION) THEN
     CALL DFFTI(schZP%gridSize(0), schZP%WSAVE1)
     CALL DFFTI(schZP%gridSize(1), schZP%WSAVE2)
  END IF

  ! Calculate normalization coefficients.
  schZP%norm_1D_1 = DBLE(1/DBLE(schZP%gridSize(0)))
  schZP%norm_1D_2 = DBLE(1/DBLE(schZP%gridSize(1)))
  schZP%norm_2D = DBLE(1/DBLE(PRODUCT(schZP%gridSize)))

  ! Create subarrays
  ALLOCATE(schZP%SUBARRAYS(schZP%commSize,0:1))
  CALL SUBARRAY(schZP%vSlabSize, 0, schZP%commSize, schZP%datatype, &
       schZP%SUBARRAYS(:,0))
  CALL SUBARRAY(schZP%hSlabSize, 1, schZP%commSize, schZP%datatype, &
       schZP%SUBARRAYS(:,1))

  ! Allocate arrays needed for ALLTOALLW communication.
  ALLOCATE(schZP%counts(schZP%commSize))
  schZP%counts = 1
  ALLOCATE(schZP%displs(schZP%commSize))
  schZP%displs = 0

  ! Set schZP as initialized for FFTs.
  schZP%initFFT = .TRUE.

  ! Create array for wavenumbers along the first dimension.
  ALLOCATE(schZP%wvNmbr1(0:schZP%gridSize(0)-1))
  DO i = 0, schZP%gridSize(0)/2
     schZP%wvNmbr1(i) = DCMPLX(0.0, i)
  END DO
  minWvNmbr = -1 * (schZP%gridSize(0)-1)/2
  DO i = 0, (schZP%gridSize(0)-3)/2
     schZP%wvNmbr1(schZP%gridSize(0)/2+1+i) = DCMPLX(0.0, minWvNmbr + i)
  END DO

  ! Create array for wavenumbers along the second dimension.
  ! We actually calculate all of the wavenumbers, then just store the values we
  ! need.
  ALLOCATE(wvNmbr2All(0:schZP%gridSize(1)-1))
  DO i = 0, schZP%gridSize(1)/2
     wvNmbr2All(i) = DCMPLX(0.0, i)
  END DO
  minWvNmbr = -1 * (schZP%gridSize(1)-1)/2
  DO i = 0, (schZP%gridSize(1)-3)/2
     wvNmbr2All(schZP%gridSize(1)/2+1+i) = DCMPLX(0.0, minWvNmbr + i)
  END DO
  
  IF (schZP%procID .EQ. 0) THEN
     stIdx = 0
  ELSE
     stIdx = SUM(schZP%colDcmpSizes(0:schZP%procID)) &
          - schZP%colDcmpSizes(schZP%procID) - schZP%ovlp
  END IF

  ALLOCATE(schZP%wvNmbr2(0:schZP%vSlabSizeOvlp(1)-1))
  DO i = 0, schZP%vSlabSizeOvlp(1)-1
     schZP%wvNmbr2(i) = wvNmbr2All(stIdx+i)
  END DO
  
  ! Set schZP as initialized for FFTs.
  schZP%initSpecDrv = .TRUE.
    
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
  !!    <li> is not initialized for spectral derivatives
  !!    </ol>
  !! <li> schZP...
  !!    <ol>
  !!    <li> is already initialized
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE CREATE_SCHEME_ZERO_PAD_DBLE_EH(sch, schZP)

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

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

    ! Check if schZP is already initialized.
    IF (schZP%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme for zero-padding is already initialized."
       END IF
       error_flag = .TRUE.
    END IF
    

  END SUBROUTINE CREATE_SCHEME_ZERO_PAD_DBLE_EH
  
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
 
END SUBROUTINE CREATE_SCHEME_ZERO_PAD_DBLE_SBR


