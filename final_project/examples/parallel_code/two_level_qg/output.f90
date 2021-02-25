!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : TWO_LEVEL_QG_SERIAL OUTPUT
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the data output subroutines needed for the example
!! serial two-level QG equation code.
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
  !> @param[in] time Current time in simulation.
  !> @param[in] dt Current time-step size.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE WRITE_OUTPUT(step, time, dt)

    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: step
    REAL(dp), INTENT(IN) :: time
    REAL(dp), INTENT(IN) :: dt
    INTEGER(qb) :: i, j !< Counters for DO loops.
    CHARACTER(LEN = 37) :: level1FileName !< Output file name for level 1.
    CHARACTER(LEN = 37) :: level2FileName !< Output file name for level 2.
    CHARACTER(LEN = 43) :: timestepInfoFileName !< Output file name for
    !! time-step information.

    WRITE(level1FileName,'(A,I0.8,A,I0.3,A)') './output_data/level1_', step, &
         '_', sch%procID, '.csv'
    WRITE(level2FileName,'(A,I0.8,A,I0.3,A)') './output_data/level2_', step, &
         '_', sch%procID, '.csv'
    WRITE(timestepInfoFileName,'(A,I0.8,A,I0.3,A)') './output_data/out_', step, &
         '_', sch%procID, '_info.txt'

    OPEN(1001, file = level1FileName, form = 'formatted')
    OPEN(1002, file = level2FileName, form = 'formatted')

    DO j = 0 , yLen-1
       DO i = 0, xLen-1
          WRITE(1001, '(E32.16)', ADVANCE = 'NO') &
               REAL(physPotVortGrid(i,j,1), dp)
          WRITE(1001, '(A)', ADVANCE = 'NO') ','

          WRITE(1002, '(E32.16)', ADVANCE = 'NO') &
               REAL(physPotVortGrid(i,j,2), dp)
          WRITE(1002, '(A)', ADVANCE = 'NO') ','
       END DO
       WRITE(1001, '(1x)')
       WRITE(1002, '(1x)')
    END DO

    CLOSE(1001)
    CLOSE(1002)


    OPEN(1005, file = timestepInfoFileName, form = 'formatted')
    WRITE(1005,'(A,I0.8,1x)') 'step = ', step
    WRITE(1005,'(A,E32.16,1x)') 'time = ', time
    WRITE(1005,'(A,E32.16,1x)') 'dt = ', dt
    CLOSE(1005)

    WRITE(*,'(A,I0.8,A)') 'Wrote outputs for step ', step, '.'
  END SUBROUTINE WRITE_OUTPUT

END MODULE OUTPUT
