!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The max value subroutine.
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
!> @param[in] subGrid The local sub-grid.
!> @param[inout] maxValue The variable to store the maximum value across all
!! sub-grids.
!> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE MAX_VAL_SBR(subGrid, maxValue, sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: subGrid(0:,0:)
  DOUBLE PRECISION, INTENT(INOUT) :: maxValue
  DOUBLE PRECISION :: maxValueTemp
  TYPE(SCHEME), INTENT(IN) :: sch
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.

  ! Check for errors in user input.
  CALL MAX_VAL_EH(sch)

  maxValueTemp = MAXVAL(subGrid)

  CALL MPI_ALLREDUCE(maxValueTemp, maxValue, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       sch%comm, ierror)


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
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE MAX_VAL_EH(sch)

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

  END SUBROUTINE MAX_VAL_EH
 
END SUBROUTINE MAX_VAL_SBR


