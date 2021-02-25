!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The <tt>SCHEME</tt> creation subroutine.
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
!> @param[in] rowCount Number of rows in the grid.
!> @param[inout] colCount Number of columns in the grid. Returns the number
!! of columns in the sub-grid, including overlap.
!> @param[in] colSpc The spacing between columns of the grid, for reference
!! point identification.
!> @param[inout] colRef The physical position of the reference point in the
!! dimension along the rows (corresponding to a column).
!> @param[in] comm <tt>MPI</tt> communicator that host the processors that
!! contain the sub-grids.
!> @param[in] mpiDatatype <tt>MPI</tt> datatype corresponding to the
!! datatype of the grid.
!> @param[in] ovlp Number of extra columns needed by each sub-grid to
!! successfully step forward in time.
!> @param[inout] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE CREATE_SCHEME_SBR(rowCount, colCount, colSpc, colRef, comm, &
     mpiDatatype, ovlp, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: rowCount
  INTEGER, INTENT(INOUT) :: colCount
  DOUBLE PRECISION, INTENT(IN) :: colSpc
  DOUBLE PRECISION, INTENT(INOUT) :: colRef
  INTEGER, INTENT(IN) :: comm
  INTEGER, INTENT(IN) :: mpiDatatype
  INTEGER, INTENT(IN) :: ovlp
  TYPE(SCHEME), INTENT(INOUT) :: sch
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  INTEGER :: stdCount !< Integer for holding standard row/col count in
  !! decompositions.
  INTEGER :: i !< Iterator for DO loops.

  CALL MPI_COMM_RANK(comm, sch%procID, ierror)
  CALL MPI_COMM_SIZE(comm, sch%commSize, ierror)

  ! Check for errors in user input.
  CALL CREATE_SCHEME_EH(rowCount, colCount, ovlp, sch)

  ! Set variables that are equal to or easily calculable from inputs.
  sch%gridSize(0) = rowCount
  sch%gridSize(1) = colCount
  sch%colSpc = colSpc
  sch%comm = comm
  sch%datatype = mpiDatatype
  sch%ovlp = ovlp

  ! If there is only one processor in the communicator, we set the overlap to 0.
  IF (sch%commSize .EQ. 1) THEN
     sch%ovlp = 0
  END IF

  ! Fill in decomposition arrays.
  ALLOCATE(sch%rowDcmpSizes(0:sch%commSize-1))
  stdCount = rowCount/sch%commSize
  DO i = 0, sch%CommSize-1
     IF (i .LT. MOD(rowCount, sch%commSize)) THEN
        sch%rowDcmpSizes(i) = stdCount + 1
     ELSE
        sch%rowDcmpSizes(i) = stdCount
     END IF
  END DO

  ALLOCATE(sch%colDcmpSizes(0:sch%commSize-1))
  stdCount = colCount/sch%commSize
  DO i = 0, sch%CommSize-1
     IF (i .LT. MOD(colCount, sch%commSize)) THEN
        sch%colDcmpSizes(i) = stdCount + 1
     ELSE
        sch%colDcmpSizes(i) = stdCount
     END IF
  END DO

  ALLOCATE(sch%colDcmpSizesOvlp(0:sch%commSize-1))
  sch%colDcmpSizesOvlp(0) = sch%colDcmpSizes(0) + sch%ovlp
  DO i = 1, sch%commSize-2
     sch%colDcmpSizesOvlp(i) = sch%colDcmpSizes(i) + 2*sch%ovlp
  END DO
  sch%colDcmpSizesOvlp(sch%commSize-1) = sch%colDcmpSizes(sch%commSize-1) &
       + sch%ovlp

  ! Set slab sizes.
  sch%hSlabSize(0) = sch%rowDcmpSizes(sch%procID)
  sch%hSlabSize(1) = sch%gridSize(1)

  sch%vSlabSize(0) = sch%gridSize(0)
  sch%vSlabSize(1) = sch%colDcmpSizes(sch%procID)

  sch%vSlabSizeOvlp(0) = sch%gridSize(0)
  sch%vSlabSizeOvlp(1) = sch%colDcmpSizesOvlp(sch%procID)

  ! Store indices of interior of vertical slab for later.
  IF (sch%procID .EQ. 0) THEN
     sch%vSlabInt(0) = 0
     sch%vSlabInt(1) = sch%vSlabSizeOvlp(1) - sch%ovlp - 1
  ELSE IF (sch%procId .EQ. sch%commSize-1) THEN
     sch%vSlabInt(0) = sch%ovlp
     sch%vSlabInt(1) = sch%vSlabSizeOvlp(1) - 1
  ELSE
     sch%vSlabInt(0) = sch%ovlp
     sch%vSlabInt(1) = sch%vSlabSizeOvlp(1) - sch%ovlp - 1
  END IF
  
  ! Adjust colCount for return parameter.
  colCount = sch%vSlabSizeOvlp(1)

  ! Calculate the subgrid reference point.
  IF ((sch%commSize .EQ. 1) .OR. (sch%procID .EQ. 0)) THEN ! If only one
     ! processor or if processor 0, no calculation needed.
     sch%colRef = colRef
  ELSE
     sch%colRef = colRef + colSpc * &
          (SUM(sch%colDcmpSizes(0:sch%procID)) - sch%colDcmpSizes(sch%procID) &
          - sch%ovlp)
     colRef = sch%colRef
  END IF

  ! If only one processor and ovlp = 0, we don't create boundary datatypes.
  IF ((sch%commSize .GT. 1) .AND. (sch%ovlp .GT. 0)) THEN
     ! Set SEND boundary datatypes.
     CALL MPI_TYPE_CREATE_SUBARRAY(2, sch%vSlabSizeOvlp, &
          (/sch%vSlabSizeOvlp(0), sch%ovlp/), (/0, sch%ovlp/), &
          MPI_ORDER_FORTRAN, sch%datatype, sch%SEND_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_CREATE_SUBARRAY(2, sch%vSlabSizeOvlp, &
          (/sch%vSlabSizeOvlp(0), sch%ovlp/), &
          (/0, sch%vSlabSizeOvlp(1)-2*sch%ovlp/), MPI_ORDER_FORTRAN, &
          sch%datatype, sch%SEND_BOUNDARIES(1), ierror)
     CALL MPI_TYPE_COMMIT(sch%SEND_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_COMMIT(sch%SEND_BOUNDARIES(1), ierror)

     ! Set RECV boundary datatypes.
     CALL MPI_TYPE_CREATE_SUBARRAY(2, sch%vSlabSizeOvlp, &
          (/sch%vSlabSizeOvlp(0), sch%ovlp/), (/0, 0/), &
          MPI_ORDER_FORTRAN, sch%datatype, sch%RECV_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_CREATE_SUBARRAY(2, sch%vSlabSizeOvlp, &
          (/sch%vSlabSizeOvlp(0), sch%ovlp/), &
          (/0, sch%vSlabSizeOvlp(1)-sch%ovlp/), MPI_ORDER_FORTRAN, &
          sch%datatype, sch%RECV_BOUNDARIES(1), ierror)
     CALL MPI_TYPE_COMMIT(sch%RECV_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_COMMIT(sch%RECV_BOUNDARIES(1), ierror)
  END IF


  ! Set sch as created.
  sch%initScheme = .TRUE.

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine for the <tt>SCHEME</tt> creation subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> rowCount...
  !!    <ol>
  !!    <li> is not positive
  !!    </ol>
  !! <li> colCount...
  !!    <ol>
  !!    <li> is not positive
  !!    <li> is too small to be divided among the processors, ignoring the
  !!    overlap.
  !!    <li> is too small to be divided among the processors, including the
  !!    overlap.
  !!    </ol>
  !! <li> ovlp...
  !!    <ol>
  !!    <li> is negative
  !!    </ol>
  !! <li> sch...
  !!    <ol>
  !!    <li> is already initialized
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE CREATE_SCHEME_EH(rowCount, colCount, ovlp, sch)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: rowCount
    INTEGER, INTENT(IN) :: colCount
    INTEGER, INTENT(IN) :: ovlp
    TYPE(SCHEME), INTENT(IN) :: sch
    INTEGER :: error_code !< Integer argument for error flag for
    !! <tt>MPI_ABORT</tt>.
    INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
    !! <tt>MPI</tt> subroutine calls.
    LOGICAL :: error_flag = .FALSE. !< Error flag to stop all processors.

    ! Check if rowCount is not positive.
    IF (rowCount .LT. 1) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Non-positive grid row count. rowCount = ", rowCount
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if colCount is not positive.
    IF (colCount .LT. 1) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Non-positive grid column count. colCount = ", colCount
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if colCount is too small to be divided among processors, excluding
    ! overlap.
    IF (colCount .LT. sch%commSize) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Grid column count too small to be divided among ", &
               "the processors in the communicator, excluding overlap. ", &
               "colCount = ", colCount, " commSize = ", sch%commSize 
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if colCount is too small to be divided among processors, including
    ! overlap.
    ! If only one processor, overlap will be set to zero anyway.
    IF ((colCount .LT. ovlp*sch%commSize) .AND. (sch%commSize .NE. 1)) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Grid column count too small to be divided among ", &
               "the processors in the communicator, including overlap. ", &
               "colCount = ", colCount, " ovlp = ", ovlp
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if ovlp is negative.
    IF (ovlp .LT. 0) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Overlap is negative."
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if sch is already initialized.
    IF (sch%initScheme) THEN
       IF (sch%procID .EQ. 0) THEN
          PRINT *, "ERROR: Scheme is already initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE CREATE_SCHEME_EH
 
END SUBROUTINE CREATE_SCHEME_SBR


