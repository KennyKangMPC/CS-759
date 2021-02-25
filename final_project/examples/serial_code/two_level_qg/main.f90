!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : TWO_LEVEL_QG_SERIAL
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief An example serial two-level quasigeostrophic (QG) equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM TWO_LEVEL_QG_SOLVER_SERIAL

  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  REAL(dp) :: startTime
  REAL(dp) :: endTime

  CALL CPU_TIME(startTime)
  CALL MAIN
  CALL CPU_TIME(endTime)

  WRITE(*,'(A,F16.8,A)') 'Execution time: ', endTime - startTime, '.'
  WRITE(*,*) 'TWO_LEVEL_QG_SOLVER_SERIAL execution complete. Normal ', &
       & 'termination...'

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> The main program driver, calls all necessary steps of the simulator.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE MAIN

    USE INITIALIZE
    USE TIME_STEPPER

    IMPLICIT NONE

    CALL INITIALIZE_PARAMETERS

    CALL INITIALIZE_GRID

    CALL TIME_STEP

    WRITE(*,*) 'Number of attempted time-steps: ', timestepCount, '.'

    RETURN

  END SUBROUTINE MAIN

END PROGRAM TWO_LEVEL_QG_SOLVER_SERIAL
