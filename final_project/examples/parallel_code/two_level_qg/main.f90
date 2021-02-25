!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : TWO_LEVEL_QG_PARALLEL
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief An example parallel two-level quasigeostrophic (QG) equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM TWO_LEVEL_QG_SOLVER_PARALLEL

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
  WRITE(*,*) 'TWO_LEVEL_QG_SOLVER_PARALLEL execution complete. Normal ', &
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
    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    INTEGER(qb) :: ierror !< Integer parameter for Fortran MPI calls.
    !! ADDED TO PARALLEL.

    CALL MPI_INIT(ierror) !< ADDED TO PARALLEL.

    CALL INITIALIZE_PARAMETERS

    CALL INITIALIZE_GRID

    CALL TIME_STEP

    !> ADDED TO PARALLEL.
    IF (sch%procID .NE. 0_qb) THEN
       CALL DESTROY_SCHEME(sch)
       CALL DESTROY_SCHEME(schZP)
       CALL MPI_FINALIZE(ierror)
       STOP
    END IF
    CALL DESTROY_SCHEME(sch)
    CALL DESTROY_SCHEME(schZP)
    CALL MPI_FINALIZE(ierror)

    WRITE(*,*) 'Number of attempted time-steps: ', timestepCount, '.'

    RETURN

  END SUBROUTINE MAIN

END PROGRAM TWO_LEVEL_QG_SOLVER_PARALLEL
