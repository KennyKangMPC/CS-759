!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STOCH_HEAT_PARALLEL TIME_STEPPER
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the time-stepping subroutine needed for the example
!! parallel stochastic heat equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

  IMPLICIT NONE

  PRIVATE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  PUBLIC :: TIME_STEP

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Driver for the time-stepping. First ensures that the CFL condition is met
  !! and that the outputFreq is not too large.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE TIME_STEP

    USE INITIALIZE
    USE OUTPUT

    IMPLICIT NONE

    INTEGER(qb) :: step !< Current time-step number.
    INTEGER(qb) :: maxStep !< Maximum number of time-steps based on time-step
    !! size and end time.
    INTEGER(qb) :: stepParity !< Parity of the time-step, to determine if the
    !! step is stored in tempGrid(:,:,0) or tempGrid(:,:,1).
    REAL(dp) :: time !< Current time in the simulation.
    REAL(dp) :: cflCondition !< CFL condition parameter.

    ! Stability constraint derived from our time step method, which is Forward
    ! Euler in time, and a 1st order centered difference for each spatial
    ! derivative (using Von Neumann analysis).
    cflCondition = 1.0_dp/(2.0_dp * (xDiffus/(dx**2.0_dp) &
         + yDiffus/(dy**2.0_dp)))

    IF ((dt .GT. cflCondition)) THEN
       PRINT *, "WARNING: Time step size ", dt, " exceeds the CFL parameter ", &
            cflCondition, ". Simulation may be unstable."
    END IF

    stepParity = 0_qb
    step = 0_qb
    time = initTime
    maxStep = INT((finTime - initTime)/dt,qb)

    IF (maxStep .LT. outputFreq) THEN
       ERROR STOP 'Output frequency too large. Please reduce.'
    END IF

    ! Output the initial condition.
    IF (outputFreq .NE. 0) THEN
       CALL WRITE_OUTPUT(step, stepParity)
    END IF

    ! Step simulation forward in time.
    DO WHILE (time .LE. finTime)
       stepParity = 2_qb - (stepParity + 1_qb)
       CALL TIME_STEP_SCHEME(stepParity)
       time = time + dt
       step = step + 1_qb

       ! Write output at desired frequency.
       IF (outputFreq .NE. 0) THEN
          IF (MOD(step, outputFreq) .EQ. 0) THEN
             CALL WRITE_OUTPUT(step, stepParity)
          END IF
       END IF
    END DO

  END SUBROUTINE TIME_STEP

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Executes the numerical scheme for the time-step, using a Forward Euler
  !! scheme in time and a first-order center difference in space. We also
  !! include a deterministic forcing, a relaxation term, and white noise
  !! (see Equation (2) of Hottovy, S., Stechmann, S. (2015)).
  !
  !> @param[in] stepParity Parity of the time-step.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE TIME_STEP_SCHEME(stepParity)

    USE INITIALIZE
    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: stepParity
    INTEGER(qb) :: prevStepParity !< Parity of the previous time-step.
    REAL(dp) :: whiteNoise(1:xLen-2,1:yLen-2) !< White noise for the
    !! stochasticity of the time-step.

    ! Get a random number for the white noise.
    CALL RANDOM_NUMBER(whiteNoise)
    whiteNoise = 2.0_dp * (whiteNoise - 0.5_dp)

    ! Update the step parity.
    prevStepParity = 2_qb - (stepParity + 1_qb)

    ! Step forward in time.
    tempGrid(1:xLen-2, 1:yLen-2, stepParity) = &
         tempGrid(1:xLen-2, 1:yLen-2, prevStepParity) &
         + (dt * xDiffus / (dx**2.0_dp)) &
         * (tempGrid(2:xLen-1, 1:yLen-2, prevStepParity) &
         - 2.0_dp * tempGrid(1:xLen-2, 1:yLen-2, prevStepParity) &
         + tempGrid(0:xLen-3, 1:yLen-2, prevStepParity)) &
         + (dt * yDiffus / (dy**2.0_dp)) &
         * (tempGrid(1:xLen-2, 2:yLen-1, prevStepParity) &
         - 2.0_dp * tempGrid(1:xLen-2, 1:yLen-2, prevStepParity) &
         + tempGrid(1:xLen-2, 0:yLen-3, prevStepParity)) &
         + deterForce - (1.0_dp/dampCoeff) &
         * (tempGrid(1:xLen-2, 1:yLen-2, stepParity) - relaxTrgt) &
         + stochMag * whiteNoise

    CALL SHARE_SUBGRID_BDRY(tempGrid(:,:,stepParity), sch) !< ADDED TO PARALLEL.

  END SUBROUTINE TIME_STEP_SCHEME

END MODULE TIME_STEPPER
