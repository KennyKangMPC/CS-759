!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STOCH_HEAT_PARALLEL OUTPUT
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the data output subroutines needed for the example
!! parallel stochastic heat equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE OUTPUT

  IMPLICIT NONE

  PRIVATE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  PUBLIC :: WRITE_OUTPUT

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Writes the output to a .csv file in the output_data subdirectory.
  !
  !> @param[in] step Time-step number.
  !> @param[in] stepParity Parity of the time-step.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE WRITE_OUTPUT(step, stepParity)

    USE INITIALIZE
    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: step
    INTEGER(qb), INTENT(IN) :: stepParity
    INTEGER(qb) :: i, j !< Counters for DO loops.
    CHARACTER(LEN=34) :: fileName !< Filename for output data. CHANGED FOR
    !! PARALLEL.

    WRITE(fileName,'(A,I0.8,A,I0.3,A)') './output_data/out_', step, '_', &
         sch%procID, '.csv' !< CHANGED FOR PARALLEL.

    OPEN(100,file=fileName,form='formatted')

    DO j = 0, yLen-1
       DO i = 0, xLen-1
          WRITE(100,'(E32.16,A,1x)',ADVANCE='NO') tempGrid(i,j,stepParity), ','
       END DO
       WRITE(100,'(1x)')
    END DO
    CLOSE(100)

    PRINT *, 'Wrote grid to ', fileName, '.'

  END SUBROUTINE WRITE_OUTPUT

END MODULE

