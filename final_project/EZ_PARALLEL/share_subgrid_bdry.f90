!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The sub-grid boundary sharing subroutine (DOUBLE PRECISION).
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
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE SHARE_SUBGRID_BDRY_DBLE_SBR(subGrid, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(INOUT) :: subGrid(:,:)
  TYPE(SCHEME), INTENT(IN) :: sch
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  INTEGER :: status(MPI_STATUS_SIZE) !< Status flag for <tt>MPI_RECV</tt>.

  ! Check for errors in user input.
  CALL SHARE_SUBGRID_BDRY_DBLE_EH(sch)

  ! If only one processor or ovlp 0, we skip all communication.
  IF ((sch%commSize .EQ. 1) .OR. (sch%ovlp .EQ. 0)) THEN
     GOTO 1001
  END IF

  ! First sub-grid only SEND/RECV to/from one neighbor.
    IF (sch%procID .EQ. 0) THEN
       CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(1), &
            sch%procID+1, sch%procID+1, sch%comm, ierror)
       
       CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(1), &
            sch%procID+1, sch%procID, sch%comm, status, ierror)
       
    ! Last sub-grid only SEND/RECV to/from one neighbor.
    ELSE IF (sch%procID .EQ. sch%commSize-1) THEN
       ! If procID  even, SEND then RECV.
       IF (MOD(sch%procID,2) .EQ. 0) THEN
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)
          
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          
       ! If procID odd, RECV then SEND.
       ELSE
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)

       END IF
       
    ! Interior sub-grids SEND/RECV to/from two neighbors.
    ELSE
       ! If proc_id even, SEND then RECV.
       IF (MOD(sch%procID,2) .EQ. 0) THEN
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(1), &
               sch%procID+1, sch%procID+1, sch%comm, ierror)
          
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(1), &
               sch%procID+1, sch%procID, sch%comm, status, ierror)
          
       ! If proc_id odd, RECV then SEND
       ELSE
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(1), &
               sch%procID+1, sch%procID, sch%comm, status, ierror)
          
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(1), &
               sch%procID+1, sch%procID+1, sch%comm, ierror)

       END IF
    END IF

1001 CONTINUE

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine for the sub-grid boundary sharing subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SHARE_SUBGRID_BDRY_DBLE_EH(sch)

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
          PRINT *, "ERROR: Scheme is not initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE SHARE_SUBGRID_BDRY_DBLE_EH
  
END SUBROUTINE SHARE_SUBGRID_BDRY_DBLE_SBR


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The sub-grid boundary sharing subroutine (DOUBLE COMPLEX).
!
!> @param[inout] subGrid The local sub-grid whose boundary will be shared.
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!
! This is just a copy-paste of what's above, just for DOUBLE COMPLEX sub-grids.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE SHARE_SUBGRID_BDRY_DCMPX_SBR(subGrid, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(:,:)
  TYPE(SCHEME), INTENT(IN) :: sch
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  INTEGER :: status(MPI_STATUS_SIZE) !< Status flag for <tt>MPI_RECV</tt>.

  ! Check for errors in user input.
  CALL SHARE_SUBGRID_BDRY_DCMPX_EH(sch)

  ! If only one processor, we skip all communication.
  IF (sch%commSize .EQ. 1) THEN
     GOTO 100
  END IF

  ! First sub-grid only SEND/RECV to/from one neighbor.
    IF (sch%procID .EQ. 0) THEN
       CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(1), &
            sch%procID+1, sch%procID+1, sch%comm, ierror)
       
       CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(1), &
            sch%procID+1, sch%procID, sch%comm, status, ierror)
       
    ! Last sub-grid only SEND/RECV to/from one neighbor.
    ELSE IF (sch%procID .EQ. sch%commSize-1) THEN
       ! If procID  even, SEND then RECV.
       IF (MOD(sch%procID,2) .EQ. 0) THEN
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)
          
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          
       ! If procID odd, RECV then SEND.
       ELSE
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)

       END IF
       
    ! Interior sub-grids SEND/RECV to/from two neighbors.
    ELSE
       ! If proc_id even, SEND then RECV.
       IF (MOD(sch%procID,2) .EQ. 0) THEN
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(1), &
               sch%procID+1, sch%procID+1, sch%comm, ierror)
          
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(1), &
               sch%procID+1, sch%procID, sch%comm, status, ierror)
          
       ! If proc_id odd, RECV then SEND
       ELSE
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(0), &
               sch%procID-1, sch%procID, sch%comm, status, ierror)
          CALL MPI_RECV(subGrid, 1, sch%RECV_BOUNDARIES(1), &
               sch%procID+1, sch%procID, sch%comm, status, ierror)
          
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(0), &
               sch%procID-1, sch%procID-1, sch%comm, ierror)
          CALL MPI_SEND(subGrid, 1, sch%SEND_BOUNDARIES(1), &
               sch%procID+1, sch%procID+1, sch%comm, ierror)

       END IF
    END IF

    100 CONTINUE

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine for the sub-grid boundary sharing subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !! <li> sch...
  !!    <ol>
  !!    <li> is not initialized
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SHARE_SUBGRID_BDRY_DCMPX_EH(sch)

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
          PRINT *, "ERROR: Scheme is not initialized."
       END IF
       error_flag = .TRUE.
    END IF

    ! If an error is found, abort everything.
    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE SHARE_SUBGRID_BDRY_DCMPX_EH
  
END SUBROUTINE SHARE_SUBGRID_BDRY_DCMPX_SBR

