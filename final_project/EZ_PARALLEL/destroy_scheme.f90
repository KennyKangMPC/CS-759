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
!> @param[inout] sch <tt>SCHEME</tt> that will be used for the grid
!! decomposition, parallel FFTs, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE DESTROY_SCHEME_SBR(sch)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  TYPE(SCHEME), INTENT(INOUT) :: sch
  INTEGER :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.

  ! Check for errors in user input.
  CALL DESTROY_SCHEME_EH(sch)
  
  ! Set variables to dummy values
  sch%gridSize = -1
  sch%colSpc = -1
  sch%datatype = -1
  sch%norm_1D_1 = -1
  sch%norm_1D_2 = -1
  sch%norm_2D = -1
  sch%colRef = -1

  ! We keep procID and comm for printing error messages.
  
  ! Deallocate allocated arrays.
  DEALLOCATE(sch%rowDcmpSizes)
  DEALLOCATE(sch%colDcmpSizes)
  DEALLOCATE(sch%colDcmpSizesOvlp)

  ! Set slab sizes to dummy values.
  sch%hSlabSize = -1
  sch%vSlabSize = -1
  sch%vSlabSizeOvlp = -1

  ! Destroy the boundary communication datatypes if they were created
  ! (non-zero overlap, more than one processor).
  IF ((sch%commSize .GT. 1) .AND. (sch%ovlp .GT. 0)) THEN
     CALL MPI_TYPE_FREE(sch%SEND_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_FREE(sch%SEND_BOUNDARIES(1), ierror)
     CALL MPI_TYPE_FREE(sch%RECV_BOUNDARIES(0), ierror)
     CALL MPI_TYPE_FREE(sch%RECV_BOUNDARIES(1), ierror)
  END IF

  ! If FFTs were initialized, destroy those as well.
  IF (sch%initFFT) THEN
     DEALLOCATE(sch%WSAVE1)
     DEALLOCATE(sch%WSAVE2)
     DEALLOCATE(sch%SUBARRAYS)
     DEALLOCATE(sch%counts)
     DEALLOCATE(sch%displs)
  END IF
  sch%initFFT = .FALSE.

  ! If spectral derivatives were initialized, destroy those as well.
  IF (sch%initSpecDrv) THEN
     DEALLOCATE(sch%wvNmbr1)
     DEALLOCATE(sch%wvNmbr2)
  END IF
  sch%initSpecDrv = .FALSE.

  ! Set last variables to dummy values.
  sch%commSize = -1
  sch%ovlp = -1
  
  ! Set sch as destroyed.
  sch%initScheme = .FALSE.
  

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The error handling subroutine for the <tt>SCHEME</tt> creation subroutine.
  !> Aborts the program if any of the input parameters satisfy any of the
  !! following conditions:
  !! <ul>
  !!    <ol>
  !!    <li> is not initialized
  !!    </ol>
  !! </ul>
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE DESTROY_SCHEME_EH(sch)

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

  END SUBROUTINE DESTROY_SCHEME_EH
  
END SUBROUTINE DESTROY_SCHEME_SBR



