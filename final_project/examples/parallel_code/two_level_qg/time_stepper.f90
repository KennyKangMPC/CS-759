!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : TWO_LEVEL_QG_PARALLEL TIME_STEPPER
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the time-stepping subroutine needed for the example
!! parallel two-level QG equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

  IMPLICIT NONE

  PRIVATE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  INTEGER(qb), PUBLIC :: timestepCount !< Counts the total number of attempted
  !! time-steps.

  PUBLIC :: TIME_STEP

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Driver for time-stepping. Checks if outputFreq is too large, and FFTs the
  !! initial condition.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE TIME_STEP

    USE INITIALIZE
    USE JACOBIAN_EKMAN_SHEAR_SOLVE
    USE OUTPUT
    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.
    
    IMPLICIT NONE

    INTEGER(qb) :: step !< Current time-step number.
    REAL(dp) :: time !< Current time in simulation.
    REAL(dp) :: dt !< Current time-step size.
    REAL(dp) :: errorToler_0 !< Current time-step error tolerance.
    COMPLEX(dp) :: esdirkCoeff(5,5) !< Coefficients for the explicit singly
    !! diagonally implicit Runge-Kutta (ESDIRK) portion of the additive
    !! Runge-Kutta (ARK) method, used for the 'fast' phenomena.
    COMPLEX(dp) :: erkCoeff(6,7) !< The coefficients for the explicit
    !! Runge-Kutta (ERK) portion of the ARK method, used for the 'slow'
    !! phenomena.
    COMPLEX(dp) :: specXDeriv(0:xLen-1,0:yLen-1) !< Array for storing the
    !! wavenumbers along the x-dimension in spectral space.
    COMPLEX(dp) :: specYDeriv(0:xLen-1,0:yLen-1) !< Array for storing the
    !! wavenumbers along the y-dimension in spectral space.
    COMPLEX(dp) :: specBiharm(0:xLen-1,0:yLen-1,2)
    COMPLEX(dp) :: freqPotVortGrid(0:xLen-1,0:yLen-1,2)
    COMPLEX(dp), ALLOCATABLE :: arrTemp(:,:) !< ADDED FOR PARALLEL. Array for
    !! storing values temporarily.

    ! Abort if outputFreq is too big.
    IF (numTimesteps .LT. outputFreq) THEN
       ERROR STOP 'Output frequency too large. Please reduce.'
    END IF

    ! Output the initial condition.
    time = initTime
    dt = initDt
    step = 0_qb
    IF (outputFreq .NE. 0_qb) THEN
       CALL WRITE_OUTPUT(step, time, dt)
    END IF

    ! FFT the physical potential vorticity to frequency space for timestepping.
    freqPotVortGrid = physPotVortGrid
    !> CHANGED FOR PARALLEL.
    !CALL CFFT2DF(xLen, yLen, freqPotVortGrid(:,:,1))
    !CALL CFFT2DF(xLen, yLen, freqPotVortGrid(:,:,2))
    CALL EXECUTE_SCHEME_FFT(freqPotVortGrid(:,:,1), FFT_2D, sch)
    CALL EXECUTE_SCHEME_FFT(freqPotVortGrid(:,:,2), FFT_2D, sch)

    ! Set up hyperviscosity for time-stepping.
    !> CHANGED FOR PARALLEL.
    !CALL SPECTRAL_X_DERIVATIVE(xLen, yLen, specXDeriv, 2_qb * biharmOrder)
    !CALL SPECTRAL_Y_DERIVATIVE(xLen, yLen, specYDeriv, 2_qb * biharmOrder)
    CALL CREATE_SPEC_DRV(2_qb * biharmOrder, 2_qb * biharmOrder, arrTemp, sch)
    specBiharm(:,:,1) = (-1.0_dp, 0.0_dp)**(REAL(biharmOrder + 1_qb, dp)) &
         * biharmViscCoeff * arrTemp
    DEALLOCATE(arrTemp)
    specBiharm(:,:,2) = specBiharm(:,:,1)

    ! Calculate the ARK coefficients.
    CALL CALCULATE_ARK_COEFF(erkCoeff, esdirkCoeff)

    timestepCount = 0_qb
    
    ! Step simulation forward in time.
    DO step = 1, numTimesteps
       errorToler_0 = 0.8_dp * errorToler
       CALL TIME_STEP_SCHEME(freqPotVortGrid, specBiharm, time, &
            & errorToler_0, dt, erkCoeff, esdirkCoeff)

       ! Write output at desired frequency.
       IF (outputFreq .NE. 0_qb) THEN
          IF (MOD(step, outputFreq) .EQ. 0_qb) THEN
             ! Output the current step.
             CALL WRITE_OUTPUT(step, time, dt)
          END IF
       END IF
    END DO

  END SUBROUTINE TIME_STEP

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Executes the numerical scheme for the time-step, using an ARK method
  !! described by Kennedy, C. A., et. al in "Additive Runge-Kutta Schemes for
  !! Convection-Diffusion-Reaction Equations" (July 2001). We use
  !! ARK4(3)6L[2]SA - ERK for the Jacobian, Ekamn friction, and vertical shear
  !! terms (which act on slow time scales) and ARK4(3)6L[2]SA - ESDIRK for
  !! the hyperviscosity term (which acts on fast time scales). Since the latter
  !! is diagonally implicit, the equation for each stage has been rearranged to
  !! isolate the new value of potential vorticity at that stage. The coefficient
  !! that arises from this rearrangement is store in stageCoeff.
  !
  !> @param[inout] freqPotVortGrid The potential vorticity grid in frequency
  !! space.
  !> @param[in] specBiharm The spectral biharmonic operator used in the
  !! hyperviscosity term of the QG equations.
  !> @param[inout] time Current time in simulation.
  !> @param[inout] errorToler_0 Current one-step error tolerance for adaptive
  !! time-stepping.
  !> @param[inout] dt Current time-step size.
  !> @param[in] erkCoeff The coefficients for the ERK portion of the ARK method,
  !! used for the 'slow' phenomena.
  !> @param[in] esdirkCoeff Coefficients for the ESDIRK portion of the ARK
  !! method, used for the 'fast' phenomena.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE TIME_STEP_SCHEME(freqPotVortGrid, specBiharm, time, &
       & errorToler_0, dt, erkCoeff, esdirkCoeff)

    USE INITIALIZE
    USE JACOBIAN_EKMAN_SHEAR_SOLVE
    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    COMPLEX(dp), INTENT(INOUT) :: freqPotVortGrid(0:xLen-1,0:yLen-1,2)
    COMPLEX(dp), INTENT(IN) :: specBiharm(0:xLen-1,0:yLen-1,2)
    REAL(dp), INTENT(INOUT) :: time
    REAL(dp), INTENT(INOUT) :: errorToler_0
    REAL(dp), INTENT(INOUT) :: dt
    COMPLEX(dp), INTENT(IN) :: erkCoeff(6,7)
    COMPLEX(dp), INTENT(IN) :: esdirkCoeff(5,5)
    REAL(dp) :: errorToler_1 !< One-step error tolerance for next step for
    !! adaptive time-stepping.
    REAL(dp) :: errorToler_temp !< ADDED TO PARALLEL. Temporary variable for
    !! one-step error tolerance.
    REAL(dp) :: maxPotVort !< The maximum potential vorticity in the current
    !! time-step.  
    COMPLEX(dp) :: stageCoeff(0:xLen-1,0:yLen-1,2) !< The multiplicative
    !! coefficient common across all stages of the ARK method.
    COMPLEX(dp) :: jacobianEkmanShear_0(0:xLen-1,0:yLen-1,2) !< The sum of the
    !! Jacobian, Ekman friction, and background shear terms at stage 0.
    COMPLEX(dp) :: biharmVisc_0(0:xLen-1,0:yLen-1,2) !< The hyperviscosity term
    !! at stage 0.
    COMPLEX(dp) :: freqPotVort_1(0:xLen-1,0:yLen-1,2) !< The potential vorticity
    !! in frequency space at stage 1.
    COMPLEX(dp) :: jacobianEkmanShear_1(0:xLen-1,0:yLen-1,2) !< The sum of the
    !! Jacobian, Ekman friction, and background shear terms at stage 1.
    COMPLEX(dp) :: biharmVisc_1(0:xLen-1,0:yLen-1,2) !< The hyperviscosity term
    !! at stage 1.
    COMPLEX(dp) :: freqPotVort_2(0:xLen-1,0:yLen-1,2) !< The potential vorticity
    !! in frequency space at stage 2.
    COMPLEX(dp) :: jacobianEkmanShear_2(0:xLen-1,0:yLen-1,2) !< The sum of the
    !! Jacobian, Ekman friction, and background shear terms at stage 2.
    COMPLEX(dp) :: biharmVisc_2(0:xLen-1,0:yLen-1,2) !< The hyperviscosity term
    !! at stage 2.
    COMPLEX(dp) :: freqPotVort_3(0:xLen-1,0:yLen-1,2) !< The potential vorticity
    !! in frequency space at stage 3.
    COMPLEX(dp) :: jacobianEkmanShear_3(0:xLen-1,0:yLen-1,2) !< The sum of the
    !! Jacobian, Ekman friction, and background shear terms at stage 3.
    COMPLEX(dp) :: biharmVisc_3(0:xLen-1,0:yLen-1,2) !< The hyperviscosity term
    !! at stage 3.
    COMPLEX(dp) :: freqPotVort_4(0:xLen-1,0:yLen-1,2) !< The potential vorticity
    !! in frequency space at stage 4.
    COMPLEX(dp) :: jacobianEkmanShear_4(0:xLen-1,0:yLen-1,2) !< The sum of the
    !! Jacobian, Ekman friction, and background shear terms at stage 4.
    COMPLEX(dp) :: biharmVisc_4(0:xLen-1,0:yLen-1,2) !< The hyperviscosity term
    !! at stage 4.
    COMPLEX(dp) :: freqPotVort_5(0:xLen-1,0:yLen-1,2) !< The potential vorticity
    !! in frequency space at stage 5.
    COMPLEX(dp) :: jacobianEkmanShear_5(0:xLen-1,0:yLen-1,2) !< The sum of the
    !! Jacobian, Ekman friction, and background shear terms at stage 5.
    COMPLEX(dp) :: biharmVisc_5(0:xLen-1,0:yLen-1,2) !< The hyperviscosity term
    !! at stage 5.
    COMPLEX(dp) :: errorControl(0:xLen-1,0:yLen-1,2) !< Stores one-step error at
    !! each grid point.
    INTEGER(qb) :: ierror !< ADDED TO PARALLEL. Integer argument for Fortran MPI
    !! calls.
    INTEGER(qb) :: errorCode !< ADDED TO PARALLEL. Integer argument in case
    !! of MPI_ABORT.
    REAL(dp) :: maxPotVort1 !< ADDED FOR PARALLEL.


1000 CONTINUE

    ! Iterate the time-step counter.
    timestepCount = timestepCount + 1_qb 

    stageCoeff = 1.0_dp/(1.0_dp - 0.25_dp*dt*specBiharm)

    ! First stage of additive RK method.
    ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
    ! at this stage.
    jacobianEkmanShear_0 = (0.0_dp, 0.0_dp)
    biharmVisc_0 = (0.0_dp, 0.0_dp)
    CALL JACOBIAN_EKMAN_SHEAR(freqPotVortGrid, xLen, yLen, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         jacobianEkmanShear_0)
    
    !> EDITTED FOR PARALLEL.
    biharmVisc_0 = specBiharm * freqPotVortGrid
    
    ! Calculate potential vorticity in frequency space at this stage.
    freqPotVort_1 = (0.0_dp, 0.0_dp)
    freqPotVort_1 = stageCoeff * (freqPotVortGrid + CMPLX(dt, 0.0_dp, dp) &
         * (erkCoeff(1,1) * jacobianEkmanShear_0 &
         + esdirkCoeff(1,1) * biharmVisc_0))

    ! Second stage of additive RK method.
    ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
    ! at this stage.
    jacobianEkmanShear_1 = (0.0_dp, 0.0_dp)
    biharmVisc_1 = (0.0_dp, 0.0_dp)
    CALL JACOBIAN_EKMAN_SHEAR(freqPotVort_1, xLen, yLen, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         jacobianEkmanShear_1)

    !> EDITTED FOR PARALLEL.
    biharmVisc_1 = specBiharm * freqPotVort_1
    
    ! Calculate potential vorticity in frequency space at this stage.
    freqPotVort_2 = (0.0_dp, 0.0_dp)
    freqPotVort_2 = stageCoeff * (freqPotVortGrid + CMPLX(dt, 0.0_dp, dp) &
         * (erkCoeff(1,2) * jacobianEkmanShear_0 &
         + erkCoeff(2,2) * jacobianEkmanShear_1 &
         + esdirkCoeff(1,2) * biharmVisc_0 &
         + esdirkCoeff(2,2) * biharmVisc_1))

    ! Third stage of additive RK method.
    ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
    ! at this stage.
    jacobianEkmanShear_2 = (0.0_dp, 0.0_dp)
    biharmVisc_2 = (0.0_dp, 0.0_dp)
    CALL JACOBIAN_EKMAN_SHEAR(freqPotVort_2, xLen, yLen, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         jacobianEkmanShear_2)

    !> EDITTED FOR PARALLEL.
    biharmVisc_2 = specBiharm * freqPotVort_2
    
    ! Calculate potential vorticity in frequency space at this stage.
    freqPotVort_3 = (0.0_dp, 0.0_dp)
    freqPotVort_3 = stageCoeff * (freqPotVortGrid + CMPLX(dt, 0.0_dp, dp) &
         * (erkCoeff(1,3) * jacobianEkmanShear_0 &
         + erkCoeff(2,3) * jacobianEkmanShear_1 &
         + erkCoeff(3,3) * jacobianEkmanShear_2 &
         + esdirkCoeff(1,3) * biharmVisc_0 &
         + esdirkCoeff(2,3) * biharmVisc_1 &
         + esdirkCoeff(3,3) * biharmVisc_2))

    ! Fourth stage of additive RK method.
    ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
    ! at this stage.
    jacobianEkmanShear_3 = (0.0_dp, 0.0_dp)
    biharmVisc_3 = (0.0_dp, 0.0_dp)
    CALL JACOBIAN_EKMAN_SHEAR(freqPotVort_3, xLen, yLen, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         jacobianEkmanShear_3)
    
    !> EDITTED FOR PARALLEL.
    biharmVisc_3 = specBiharm * freqPotVort_3
    
    ! Calculate potential vorticity in frequency space at this stage.
    freqPotVort_4 = (0.0_dp, 0.0_dp)
    freqPotVort_4 = stageCoeff * (freqPotVortGrid + CMPLX(dt, 0.0_dp, dp) &
         * (erkCoeff(1,4) * jacobianEkmanShear_0 &
         + erkCoeff(2,4) * jacobianEkmanShear_1 &
         + erkCoeff(3,4) * jacobianEkmanShear_2 &
         + erkCoeff(4,4) * jacobianEkmanShear_3 &
         + esdirkCoeff(1,4) * biharmVisc_0 &
         + esdirkCoeff(2,4) * biharmVisc_1 &
         + esdirkCoeff(3,4) * biharmVisc_2 &
         + esdirkCoeff(4,4) * biharmVisc_3))

    ! Fifth stage of additive RK method.
    ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
    ! at this stage.
    jacobianEkmanShear_4 = (0.0_dp, 0.0_dp)
    biharmVisc_4 = (0.0_dp, 0.0_dp)
    CALL JACOBIAN_EKMAN_SHEAR(freqPotVort_4, xLen, yLen, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         jacobianEkmanShear_4)

    !> EDITTED FOR PARALLEL.
    biharmVisc_4 = specBiharm * freqPotVort_4
    
    ! Calculate potential vorticity in frequency space at this stage.
    freqPotVort_5 = (0.0_dp, 0.0_dp)
    freqPotVort_5 = stageCoeff * (freqPotVortGrid + CMPLX(dt, 0.0_dp, dp) &
         * (erkCoeff(1,5) * jacobianEkmanShear_0 &
         + erkCoeff(2,5) * jacobianEkmanShear_1 &
         + erkCoeff(3,5) * jacobianEkmanShear_2 &
         + erkCoeff(4,5) * jacobianEkmanShear_3 &
         + erkCoeff(5,5) * jacobianEkmanShear_4 &
         + esdirkCoeff(1,5) * biharmVisc_0 &
         + esdirkCoeff(3,5) * biharmVisc_2 &
         + esdirkCoeff(4,5) * biharmVisc_3 &
         + esdirkCoeff(5,5) * biharmVisc_4))

    ! Sixth stage of additive RK method.
    ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
    ! at this stage.
    jacobianEkmanShear_5 = (0.0_dp, 0.0_dp)
    biharmVisc_5 = (0.0_dp, 0.0_dp)
    CALL JACOBIAN_EKMAN_SHEAR(freqPotVort_5, xLen, yLen, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         jacobianEkmanShear_5)

    !> EDITTED FOR PARALLEL.
    biharmVisc_5 = specBiharm * freqPotVort_5
    
    ! Error control, see the final two rows of the Butcher Tableau for each RK
    ! method.
    ! Calculate the error for adaptive time stepping.
    errorControl = (0.0_dp, 0.0_dp)
    errorControl = erkCoeff(1,6) * (jacobianEkmanShear_0 + biharmVisc_0) &
         + erkCoeff(3,6) * (jacobianEkmanShear_2 + biharmVisc_2) &
         + erkCoeff(4,6) * (jacobianEkmanShear_3 + biharmVisc_3) &
         + erkCoeff(5,6) * (jacobianEkmanShear_4 + biharmVisc_4) &
         + erkCoeff(6,6) * (jacobianEkmanShear_5 + biharmVisc_5)
    
    ! If the one-step error exceeds the tolerated value, decrease dt by 3/4 and
    ! try again.
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DB(xLen, yLen, errorControl(:,:,1))
    !CALL CFFT2DB(xLen, yLen, errorControl(:,:,2))
    !errorToler_1 = dt * MAXVAL(ABS(errorControl))
    CALL EXECUTE_SCHEME_IFFT(errorControl(:,:,1), FFT_2D, sch)
    CALL EXECUTE_SCHEME_IFFT(errorControl(:,:,2), FFT_2D, sch)
    CALL MAX_VAL(ABS(errorControl(:,:,1)), errorToler_1, sch)
    CALL MAX_VAL(ABS(errorControl(:,:,2)), errorToler_temp, sch)
    errorToler_1 = dt * MAX(errorToler_1, errorToler_temp)
    
    IF (errorToler_1 .GT. errorToler) THEN
       dt = 0.75_dp * dt
       !IF (outputFreq .NE. 0_qb) THEN
       !   PRINT *, 'Step error too large, retrying...'
       !END IF
       GO TO 1000
    END IF

    ! Successful step, proceed to evaluation.
    time = time + dt

    ! Store RK coefficients for calculating the next time step. Since the ERK
    ! and ESDIRK coefficients match, we only need to use the ERK ones.
    
    freqPotVortGrid = freqPotVortGrid + CMPLX(dt, 0.0_dp, dp) &
         * (erkCoeff(1,7) * (jacobianEkmanShear_0 + biharmVisc_0) &
         + erkCoeff(3,7) * (jacobianEkmanShear_2 + biharmVisc_2) &
         + erkCoeff(4,7) * (jacobianEkmanShear_3 + biharmVisc_3) &
         + erkCoeff(5,7) * (jacobianEkmanShear_4 + biharmVisc_4) &
         + erkCoeff(6,7) * (jacobianEkmanShear_5 + biharmVisc_5))
    
    ! Tranform solution to physical space to check against upper bound of
    ! physical potential vorticity.
    physPotVortGrid = freqPotVortGrid
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DB(xLen, yLen, physPotVortGrid(:,:,1))
    !CALL CFFT2DB(xLen, yLen, physPotVortGrid(:,:,2))
    CALL EXECUTE_SCHEME_IFFT(physPotVortGrid(:,:,1), FFT_2D, sch)
    CALL EXECUTE_SCHEME_IFFT(physPotVortGrid(:,:,2), FFT_2D, sch)
    
    ! Zero-out complex part artifacts leftover from inverse FFT.
    physPotVortGrid = REAL(physPotVortGrid)
    ! Transform back to frequency space now that complex artifacts are gone.
    freqPotVortGrid = physPotVortGrid
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DF(xLen, yLen, freqPotVortGrid(:,:,1))
    !CALL CFFT2DF(xLen, yLen, freqPotVortGrid(:,:,2))
    CALL EXECUTE_SCHEME_FFT(freqPotVortGrid(:,:,1), FFT_2D, sch)
    CALL EXECUTE_SCHEME_FFT(freqPotVortGrid(:,:,2), FFT_2D, sch)
    
    
    ! Step size adjustment: EPS, PI.3.4.
    dt = ((0.75_dp * errorToler)/errorToler_1)**(0.075_dp) &
         * (errorToler_0/errorToler_1)**(0.1_dp) * dt
    errorToler_0 = errorToler_1

    !> EDITTED FOR PARALLEL.
    !maxPotVort = MAXVAL(REAL(ABS(physPotVortGrid), dp))
    CALL MAX_VAL(REAL(ABS(physPotVortGrid(:,:,1)), dp), maxPotVort, sch)
    CALL MAX_VAL(REAL(ABS(physPotVortGrid(:,:,2)), dp), maxPotVort1, sch)
    maxPotVort = MAX(maxPotVort, maxPotVort1)
    IF (maxPotVort .GT. potVortBound) THEN
       !> EDITTED FOR PARALLEL.
       !WRITE(*,'(A,F10.3,A,F10.3,A)') 'ERROR: Max physPotVortGrid = ', &
       !     & maxPotVort, ' exceeds potVortBound = ', potVortBound, '.'
       !ERROR STOP
       IF (sch%procId .NE. 0) THEN
          CALL MPI_ABORT(MPI_COMM_WORLD, errorCode, ierror)
       END IF
       WRITE(*,'(A,F10.3,A,F10.3,A)') 'ERROR: Max physPotVortGrid = ', &
            & maxPotVort, ' exceeds potVortBound = ', potVortBound, '.'
       CALL MPI_ABORT(MPI_COMM_WORLD, errorCode, ierror)
    END IF
    
  END SUBROUTINE TIME_STEP_SCHEME

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the coefficients for an additive Runge-Kutta method
  !! described by Kennedy, C. A., et. al in "Additive Runge-Kutta Schemes for
  !! Convection-Diffusion-Reaction Equations" (July 2001)
  !
  !> @param[inout] erkCoeff The coefficients for the ERK portion of the ARK
  !! method, used for the 'slow' phenomena.
  !> @param[inout] esdirkCoeff Coefficients for the ESDIRK portion of the ARK
  !! method, used for the 'fast' phenomena.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE CALCULATE_ARK_COEFF(erkCoeff, esdirkCoeff)

    IMPLICIT NONE
    
    COMPLEX(dp) :: esdirkCoeff(5,5)
    COMPLEX(dp) :: erkCoeff(6,7)
    
    ! Zero-out the arrays before filling them with the proper coefficients.
    erkCoeff(:,:) = CMPLX(0.0_dp, 0.0_dp, dp)
    esdirkCoeff(:,:) = CMPLX(0.0_dp, 0.0_dp, dp)
    
    ! ARK coefficients for calculating the first stage.
    erkCoeff(1,1) = CMPLX((1.0_dp/2.0_dp), 0.0_dp, dp)
    esdirkCoeff(1,1) = CMPLX((1.0_dp/4.0_dp), 0.0_dp, dp)
    
    ! ARK coefficients for calculating the second stage.
    erkCoeff(1,2) = CMPLX((13861.0_dp/62500.0_dp), 0.0_dp, dp)
    erkCoeff(2,2) = CMPLX((6889.0_dp/62500.0_dp), 0.0_dp, dp)
    esdirkCoeff(1,2) = CMPLX((8611.0_dp/62500.0_dp), 0.0_dp, dp)
    esdirkCoeff(2,2) = CMPLX((-1743.0_dp/31250.0_dp), 0.0_dp, dp)
    
    ! ARK coefficients for calculating the third stage.
    erkCoeff(1,3) = CMPLX((-116923316275.0_dp/2393684061468.0_dp), 0.0_dp, dp)
    erkCoeff(2,3) = CMPLX((-2731218467317.0_dp/15368042101831.0_dp), 0.0_dp, dp)
    erkCoeff(3,3) = CMPLX((9408046702089.0_dp/11113171139209.0_dp), 0.0_dp, dp)
    esdirkCoeff(1,3) = CMPLX((5012029.0_dp/34652500.0_dp), 0.0_dp, dp)
    esdirkCoeff(2,3) = CMPLX((-654441.0_dp/2922500.0_dp), 0.0_dp, dp)
    esdirkCoeff(3,3) = CMPLX((174375.0_dp/388108.0_dp), 0.0_dp, dp)
    
    ! ARK coefficients for calculating the fourth stage.
    erkCoeff(1,4) = CMPLX((-451086348788.0_dp/2902428689909.0_dp), 0.0_dp, dp)
    erkCoeff(2,4) = CMPLX((-2682348792572.0_dp/7519795681897.0_dp), 0.0_dp, dp)
    erkCoeff(3,4) = CMPLX((12662868775082.0_dp/11960479115383.0_dp), 0.0_dp, dp)
    erkCoeff(4,4) = CMPLX((3355817975965.0_dp/11060851509271.0_dp), 0.0_dp, dp)
    esdirkCoeff(1,4) = CMPLX((15267082809.0_dp/155376265600.0_dp), 0.0_dp, dp)
    esdirkCoeff(2,4) = CMPLX((-71443401.0_dp/120774400.0_dp), 0.0_dp, dp)
    esdirkCoeff(3,4) = CMPLX((730878875.0_dp/902184768.0_dp), 0.0_dp, dp)
    esdirkCoeff(4,4) = CMPLX((2285395.0_dp/8070912.0_dp), 0.0_dp, dp)
    
    ! ARK coefficients for calculating the fifth stage.
    erkCoeff(1,5) = CMPLX((647845179188.0_dp/3216320057751.0_dp), 0.0_dp, dp)
    erkCoeff(2,5) = CMPLX((73281519250.0_dp/8382639484533.0_dp), 0.0_dp, dp)
    erkCoeff(3,5) = CMPLX((552539513391.0_dp/3454668386233.0_dp), 0.0_dp, dp)
    erkCoeff(4,5) = CMPLX((3354512671639.0_dp/8306763924573.0_dp), 0.0_dp, dp)
    erkCoeff(5,5) = CMPLX((4040.0_dp/17871.0_dp), 0.0_dp, dp)
    esdirkCoeff(1,5) = CMPLX((82889.0_dp/524892.0_dp), 0.0_dp, dp)
    esdirkCoeff(3,5) = CMPLX((15625.0_dp/83664.0_dp), 0.0_dp, dp)
    esdirkCoeff(4,5) = CMPLX((69875.0_dp/102672.0_dp), 0.0_dp, dp)
    esdirkCoeff(5,5) = CMPLX((-2260.0_dp/8211.0_dp), 0.0_dp, dp)
    
    ! ARK coefficients for calculating the error control.
    erkCoeff(1,6) = &
         CMPLX((82889.0_dp/524892.0_dp) - (4586570599.0_dp/29645900160.0_dp), &
         0.0_dp, dp)
    erkCoeff(3,6) = &
         CMPLX((15625.0_dp/83664.0_dp) - (178811875.0_dp/945068544.0_dp), &
         0.0_dp, dp)
    erkCoeff(4,6) = &
         CMPLX((69875.0_dp/102672.0_dp) - (814220225.0_dp/1159782912.0_dp), &
         0.0_dp, dp)
    erkCoeff(5,6) = &
         CMPLX((-2260.0_dp/8211.0_dp) - (-3700637.0_dp/11593932.0_dp), &
         0.0_dp, dp)
    erkCoeff(6,6) = &
         CMPLX((1.0_dp/4.0_dp) - (61727.0_dp/225920.0_dp), 0.0_dp, dp)
    
    ! ARK coefficients for calculating the next time step.
    erkCoeff(1,7) = CMPLX((82889.0_dp/524892.0_dp), 0.0_dp, dp)
    erkCoeff(3,7) = CMPLX((15625.0_dp/83664.0_dp), 0.0_dp, dp)
    erkCoeff(4,7) = CMPLX((69875.0_dp/102672.0_dp), 0.0_dp, dp)
    erkCoeff(5,7) = CMPLX((-2260.0_dp/8211.0_dp), 0.0_dp, dp)
    erkCoeff(6,7) = CMPLX((1.0_dp/4.0_dp), 0.0_dp, dp)
    
  END SUBROUTINE CALCULATE_ARK_COEFF
  
END MODULE TIME_STEPPER
