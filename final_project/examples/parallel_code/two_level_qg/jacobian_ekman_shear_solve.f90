!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : TWO_LEVEL_QG_PARALLEL JACOBIAN_EKMAN_SHEAR_SOLVE
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the solver for the Jacobian, Ekman friction, and
!! vertical shear terms subroutine needed for the example parallel two-level QG
!! equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE JACOBIAN_EKMAN_SHEAR_SOLVE
  
  IMPLICIT NONE

  PRIVATE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'


  COMPLEX(dp), ALLOCATABLE :: specLaplacian(:,:) !< The discrete Laplacian
  !! operator in spectral space.
  COMPLEX(dp), ALLOCATABLE :: specInvBarotropic(:,:) !< The discrete inverse
  !! barotropic operator in spectral space (see two-layer QG equations in Di, 
  !! Q., & Majda, A. (2015)).
  COMPLEX(dp), ALLOCATABLE :: specInvBaroclinic(:,:)!< The discrete inverse
  !! baroclinic operator in spectral space (see two-layer QG equations in Di, 
  !! Q., & Majda, A. (2015)).
  COMPLEX(dp), ALLOCATABLE :: specXDeriv(:,:,:) !< The wavenumber array for
  !! calculating the x-derivative in spectral space.
  COMPLEX(dp), ALLOCATABLE :: specYDeriv(:,:,:) !< The wavenumber array for
  !! calculating the y-derivative in spectral space.
  COMPLEX(dp), ALLOCATABLE :: scaledSpecXDeriv(:,:,:) !< The wavenumber array
  !! for calculating the x-derivative in spectral space for the zero-padded
  !! grid.
  COMPLEX(dp), ALLOCATABLE :: scaledSpecYDeriv(:,:,:) !< The wavenumber array
  !! for calculating the y-derivative in spectral space for the zero-padded
  !! grid.
  LOGICAL :: ran = .FALSE. !< Checks if the discrete differential spectral
  !! operators needed for the solver.

  PUBLIC :: SPECTRAL_X_DERIVATIVE
  PUBLIC :: SPECTRAL_Y_DERIVATIVE
  PUBLIC :: JACOBIAN_EKMAN_SHEAR

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the Jacobian term, the Ekman friction term, and the vertical
  !! shear term in the QG equation for use in progressing the simulation
  !! forward in time.
  !
  !> If the subroutine hasn't been called before, calculate all of the
  !! differential operators needed. Calculate the contributions from the mean
  !! vertical shear and rotation deformation in levels 1 and 2, and the Ekman
  !! friction in layer 2. Rescale the grid to de-alias it and calculate the
  !! derivatives in spectral space. Transform into physical space to perform
  !! the multiplications to calculate the Jacobian term, and transform back
  !! into frequency space to return the final result.
  !
  !> @param[in] freqPotVortGrid Potential vorticity grid in frequency space.
  !> @param[in] xLenSub Size of the grid along the first dimension, local to
  !! subroutine.
  !> @param[in] yLenSub Size of the grid along the second dimension, local to
  !! subroutine.
  !> @param[in] deformWavenumSub The baroclinic deformation wavenumber
  !! corresponding to the Rossby radius of deformation, local to subroutine.
  !> @param[in] rotateWavenumSub The wavenumber corresponding to the
  !! rotation coefficient, which controls the advection of streamfunctions,
  !! local to subroutine.
  !> @param[in] vertShearSub Large-scale vertical shear, opposite direction
  !! in each level to induce baroclinic instability, local to subroutine.
  !> @param[in] ekamnFricCoeffSub The coefficient for Ekman friction in level 2,
  !! local to subroutine.
  !> @param[inout] jacobianEkmanShearGrid Contains sum of the Jacobian, Ekman
  !! friction, and mean shear terms.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE JACOBIAN_EKMAN_SHEAR(freqPotVortGrid, xLenSub, yLenSub, &
       & deformWavenumSub, rotatWavenumSub, vertShearSub, ekmanFricCoeffSub, &
       & jacobianEkmanShearGrid)

    !> ADDED TO PARALLEL.
    USE MPI
    USE EZ_PARALLEL
    USE EZ_PARALLEL_STRUCTS
    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: xLenSub 
    INTEGER(qb), INTENT(IN) :: yLenSub
    COMPLEX(dp), INTENT(IN) :: freqPotVortGrid(xLenSub,yLenSub,2)
    REAL(dp), INTENT(IN) :: deformWavenumSub
    REAL(dp), INTENT(IN) :: rotatWavenumSub
    REAL(dp), INTENT(IN) :: vertShearSub
    REAL(dp), INTENT(IN) :: ekmanFricCoeffSub
    COMPLEX(dp), INTENT(INOUT) :: jacobianEkmanShearGrid(xLenSub,yLenSub,2)
    INTEGER(qb) :: scaledXLen !< Length of the zero-padded grid along the first
    !! dimension (xLenSub/2*3).
    INTEGER(qb) :: scaledYLen !< Length of the zero-padded grid along the second
    !! dimension (yLenSub/2*3).
    INTEGER(qb) :: i, j !< Counters for DO loops.
    COMPLEX(dp) :: barotropicPotVort(xLenSub,yLenSub) !< Barotropic mode of
    !! potential vorticity.
    COMPLEX(dp) :: baroclinicPotVort(xLenSub,yLenSub) !< Baroclinic mode of
    !! potential vorticity.
    COMPLEX(dp) :: barotropicPotVortStrmfunc(xLenSub,yLenSub) !< Barotropic mode
    !! of the potential vorticity streamfunction.
    COMPLEX(dp) :: baroclinicPotVortStrmfunc(xLenSub,yLenSub) !< Baroclinic mode
    !! of the potential vorticity streamfunction.
    COMPLEX(dp) :: freqStrmfunc(xLenSub,yLenSub,2) !< Potential vorticity
    !! streamfunction in frequency space.
    COMPLEX(dp) :: freqJacobian(xLenSub,yLenSub,2) !< Jacobian term in frequency
    !! space.
    COMPLEX(dp), ALLOCATABLE :: scaledFreqPotVortGrid(:,:,:) !< Zero-padded
    !! potential vorticity grid in frequency space.
    COMPLEX(dp), ALLOCATABLE :: scaledFreqStrmfunc(:,:,:) !< Zero-padded
    !! streamfunction grid in frequency space.
    COMPLEX(dp), ALLOCATABLE :: scaledStrmfuncXDeriv(:,:,:) !< Zero-padded
    !! x-derivative of the streamfunction in frequency space.
    COMPLEX(dp), ALLOCATABLE :: scaledStrmfuncYDeriv(:,:,:) !< Zero-padded
    !! y-derivative of the streamfunction in frequency space.
    COMPLEX(dp), ALLOCATABLE :: scaledPotVortXDeriv(:,:,:) !< Zero-padded
    !! x-derivative of the potential vorticity in frequency space.
    COMPLEX(dp), ALLOCATABLE :: scaledPotVortYDeriv(:,:,:) !< Zero-padded
    !! y-derivative of the potential vorticity in frequency space.
    COMPLEX(dp), ALLOCATABLE :: scaledFreqJacobian(:,:,:) !< Zero-padded
    !! Jacobian term in frequency space.
    COMPLEX(dp), ALLOCATABLE :: arrTemp(:,:) !< ADDED TO PARALLEL. Temporary
    !! array for storing values.

    ! If it hasn't been called before, calculate all differential operators
    ! needed.
    IF (.NOT. ran) THEN
       ! Get the spectral Laplacian operator.
       !> EDITTED FOR PARALLEL.
       !ALLOCATE(specXDeriv(xLenSub, yLenSub, 1))
       !ALLOCATE(specYDeriv(xLenSub, yLenSub, 1))
       !CALL SPECTRAL_X_DERIVATIVE(xLenSub, yLenSub, specXDeriv, 2_qb)
       !CALL SPECTRAL_Y_DERIVATIVE(xLenSub, yLenSub, specYDeriv, 2_qb)
       !ALLOCATE(specLaplacian(xLenSub, yLenSub))
       !specLaplacian = specXDeriv(:,:,1) + specYDeriv(:,:,1)
       !DEALLOCATE(specXDeriv)
       !DEALLOCATE(specYDeriv)
       CALL CREATE_SPEC_DRV(2,2,specLaplacian,sch)

       ! Get the inverse spectral barotropic and inverse baroclinic operators.
       ALLOCATE(specInvBarotropic(xLenSub, yLenSub))
       ALLOCATE(specInvBaroclinic(xLenSub, yLenSub))
       specInvBarotropic = 1.0_dp/specLaplacian
       specInvBaroclinic = 1.0_dp/(specLaplacian - deformWavenumSub**(2.0_dp))
       DO j = 1, yLen
          DO i = 1, xLen
             IF (specInvBarotropic(i,j) .NE. specInvBarotropic(i,j)) THEN
                specInvBarotropic(i,j) = 0.0
             END IF
             IF (specInvBaroclinic(i,j) .NE. specInvBaroclinic(i,j)) THEN
                specInvBaroclinic(i,j) = 0.0
             END IF
          END DO
       END DO

       ! Get the first-order differential operators.
       ALLOCATE(specXDeriv(xLenSub, yLenSub, 2))
       ALLOCATE(specYDeriv(xLenSub, yLenSub, 2))
       !> EDITTED FOR PARALLEL.
       !CALL SPECTRAL_X_DERIVATIVE(xLenSub, yLenSub, specXDeriv(:,:,1), 1_qb)
       !specXDeriv(:,:,2) = specXDeriv(:,:,1)
       !CALL SPECTRAL_Y_DERIVATIVE(xLenSub, yLenSub, specYDeriv(:,:,1), 1_qb)
       !specYDeriv(:,:,2) = specYDeriv(:,:,1)
       CALL CREATE_SPEC_DRV(1,0,arrTemp,sch)
       specXDeriv(:,:,1) = arrTemp
       specXDeriv(:,:,2) = arrTemp
       DEALLOCATE(arrTemp)
       CALL CREATE_SPEC_DRV(0,1,arrTemp,sch)
       specYDeriv(:,:,1) = arrTemp
       specYDeriv(:,:,2) = arrTemp
       DEALLOCATE(arrTemp)

       !> EDITTED FOR PARALLEL.
       ALLOCATE(scaledSpecXDeriv(schZP%vSlabSizeOvlp(0), schZP%vSlabSizeOvlp(1), 2))
       ALLOCATE(scaledSpecYDeriv(schZP%vSlabSizeOvlp(0), schZP%vSlabSizeOvlp(1), 2))
       !CALL SPECTRAL_X_DERIVATIVE(3_qb * xLenSub/2_qb, 3_qb * yLenSub/2_qb, &
       !     & scaledSpecXDeriv(:,:,1), 1_qb)
       !scaledSpecXDeriv(:,:,2) = scaledSpecXDeriv(:,:,1)
       !CALL SPECTRAL_Y_DERIVATIVE(3_qb * xLenSub/2_qb, 3_qb * yLenSub/2_qb, &
       !     & scaledSpecYDeriv(:,:,1), 1_qb)
       !scaledSpecYDeriv(:,:,2) = scaledSpecYDeriv(:,:,1)
       CALL CREATE_SPEC_DRV(1,0,arrTemp,schZP)
       scaledSpecXDeriv(:,:,1) = arrTemp
       scaledSpecXDeriv(:,:,2) = arrTemp
       DEALLOCATE(arrTemp)
       CALL CREATE_SPEC_DRV(0,1,arrTemp,schZP)
       scaledSpecYDeriv(:,:,1) = arrTemp
       scaledSpecYDeriv(:,:,2) = arrTemp
       DEALLOCATE(arrTemp)
       
       ran = .TRUE.
    END IF

    ! Calculate the spectral baroclinic and barotropic potential vorticities.
    barotropicPotVort = 0.5_dp * (freqPotVortGrid(:,:,1) &
         + freqPotVortGrid(:,:,2))
    baroclinicPotVort = 0.5_dp * (freqPotVortGrid(:,:,1) &
         - freqPotVortGrid(:,:,2))

    ! Calculate the streamfunctions for the spectral baroclinic and barotropic
    ! potential vorticities.
    barotropicPotVortStrmfunc = specInvBarotropic * barotropicPotVort
    baroclinicPotVortStrmfunc = specInvBaroclinic * baroclinicPotVort

    ! Calculate the strmfunction for the spectral potential vort.
    freqStrmfunc(:,:,1) = barotropicPotVortStrmfunc &
         + baroclinicPotVortStrmfunc
    freqStrmfunc(:,:,2) = barotropicPotVortStrmfunc &
         - baroclinicPotVortStrmfunc

    ! Calculate the contributions from the mean shear, rotation deformation, and
    ! the Ekman friction (the last in layer 2 only).
    jacobianEkmanShearGrid = (0.0_dp, 0.0_dp)
    ! Set layer 1.
    jacobianEkmanShearGrid(:,:,1) = CMPLX((-1.0_dp), 0.0_dp, dp) &
         * CMPLX(vertShearSub, 0.0_dp, dp) * specXDeriv(:,:,1) &
         * freqPotVortGrid(:,:,1) - CMPLX((rotatWavenumSub**(2.0_dp) &
         + vertShearSub * deformWavenumSub**(2.0_dp)), 0.0_dp, dp) &
         * specXDeriv(:,:,1) * freqStrmfunc(:,:,1)
    ! Set layer 2.
    jacobianEkmanShearGrid(:,:,2) = CMPLX((1.0_dp), 0.0_dp, dp) &
         * CMPLX(vertShearSub, 0.0_dp, dp) * specXDeriv(:,:,1) &
         * freqPotVortGrid(:,:,2) - CMPLX((rotatWavenumSub**(2.0_dp) &
         - vertShearSub * deformWavenumSub**(2.0_dp)), 0.0_dp, dp) &
         * specXDeriv(:,:,1) * freqStrmfunc(:,:,2) - CMPLX(ekmanFricCoeffSub, &
         0.0_dp, dp) * specLaplacian * freqStrmfunc(:,:,2)

    ! Must rescale the potential vort and the streamfunction grids to dealias the
    ! Jacobian, see Orszag, S. "On the Elimination of Aliasing in Finite-Difference 
    ! Schemes by Filtering High-Wavenumber Components" (1971).
    !> EDITTED FOR PARALLEL.
    !scaledXLen = xLenSub/2_qb * 3_qb
    !scaledYLen = yLenSub/2_qb * 3_qb
    scaledXLen = schZP%vSlabSize(0)
    scaledYLen = schZP%vSlabSize(1)
    ! Rescale the frequency potential vorticity.
    ALLOCATE(scaledFreqPotVortGrid(scaledXLen, scaledYLen, 2))
    !> EDITTED FOR PARALLEL.
    !scaledFreqPotVortGrid = (0.0_dp, 0.0_dp)
    !CALL ZERO_PADDING(xLenSub, yLenSub, freqPotVortGrid(:,:,1), &
    !     scaledXLen, scaledYLen, scaledFreqPotVortGrid(:,:,1))
    !CALL ZERO_PADDING(xLenSub, yLenSub, freqPotVortGrid(:,:,2), &
    !     scaledXLen, scaledYLen, scaledFreqPotVortGrid(:,:,2))
    CALL EXECUTE_SCHEME_ZERO_PAD(freqPotVortGrid(:,:,1), sch, arrTemp, schZP)
    scaledFreqPotVortGrid(:,:,1) = arrTemp
    DEALLOCATE(arrTemp)
    CALL EXECUTE_SCHEME_ZERO_PAD(freqPotVortGrid(:,:,2), sch, arrTemp, schZP)
    scaledFreqPotVortGrid(:,:,2) = arrTemp
    DEALLOCATE(arrTemp)
    scaledFreqPotVortGrid = (2.25_dp, 0.0_dp) * scaledFreqPotVortGrid

    
    ! Rescale the potential vorticity streamfunction.
    ALLOCATE(scaledFreqStrmfunc(scaledXLen, scaledYLen, 2))
    !> EDITTED FOR PARALLEL.
    !scaledFreqStrmfunc = (0.0_dp, 0.0_dp)
    !CALL ZERO_PADDING(xLenSub, yLenSub, freqStrmfunc(:,:,1), scaledXLen, &
    !     scaledYLen, scaledFreqStrmfunc(:,:,1))
    !CALL ZERO_PADDING(xLenSub, yLenSub, freqStrmfunc(:,:,2), scaledXLen, &
    !     scaledYLen, scaledFreqStrmfunc(:,:,2))
    CALL EXECUTE_SCHEME_ZERO_PAD(freqStrmfunc(:,:,1), sch, arrTemp, schZP)
    scaledFreqStrmfunc(:,:,1) = arrTemp
    DEALLOCATE(arrTemp)
    CALL EXECUTE_SCHEME_ZERO_PAD(freqStrmfunc(:,:,2), sch, arrTemp, schZP)
    scaledFreqStrmfunc(:,:,2) = arrTemp
    DEALLOCATE(arrTemp)
    scaledFreqStrmfunc = (2.25_dp, 0.0_dp) * scaledFreqStrmfunc

    ! To avoid convolution, we will calculate the Jacobian in physical space, and
    ! then transform it back to freq space. We want
    ! J(streamfunction, potential vorticity).
    ! Calculate x-derivative of the potential vorticity.
    ALLOCATE(scaledStrmfuncXDeriv(scaledXLen, scaledYLen, 2))
    scaledStrmfuncXDeriv = scaledSpecXDeriv * scaledFreqStrmfunc
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DB(scaledXLen, scaledYLen, scaledStrmfuncXDeriv(:,:,1))
    !CALL CFFT2DB(scaledXLen, scaledYLen, scaledStrmfuncXDeriv(:,:,2))
    CALL EXECUTE_SCHEME_IFFT(scaledStrmfuncXDeriv(:,:,1), FFT_2D, schZP)
    CALL EXECUTE_SCHEME_IFFT(scaledStrmfuncXDeriv(:,:,2), FFT_2D, schZP)
    ! Take only the real part to get rid of machine-epsilon errors from the
    ! inverse FFT.
    scaledStrmfuncXDeriv = REAL(scaledStrmfuncXDeriv)
    ! Calculate y-derivative of the streamfunction.
    ALLOCATE(scaledStrmfuncYDeriv(scaledXLen, scaledYLen, 2))
    scaledStrmfuncYDeriv = scaledSpecYDeriv * scaledFreqStrmfunc
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DB(scaledXLen, scaledXLen, scaledStrmfuncYDeriv(:,:,1))
    !CALL CFFT2DB(scaledYLen, scaledYLen, scaledStrmfuncYDeriv(:,:,2))
    CALL EXECUTE_SCHEME_IFFT(scaledStrmfuncYDeriv(:,:,1), FFT_2D, schZP)
    CALL EXECUTE_SCHEME_IFFT(scaledStrmfuncYDeriv(:,:,2), FFT_2D, schZP)
    ! Take only the real part to get rid of machine-epsilon errors from the
    ! inverse FFT.
    scaledStrmfuncYDeriv = REAL(scaledStrmfuncYDeriv)
    ! Calculate x-derivative of the streamfunction.
    ALLOCATE(scaledPotVortXDeriv(scaledXLen, scaledYLen, 2))
    scaledPotVortXDeriv = scaledSpecXDeriv * scaledFreqPotVortGrid
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DB(scaledXLen, scaledYLen, scaledPotVortXDeriv(:,:,1))
    !CALL CFFT2DB(scaledXLen, scaledYLen, scaledPotVortXDeriv(:,:,2))
    CALL EXECUTE_SCHEME_IFFT(scaledPotVortXDeriv(:,:,1), FFT_2D, schZP)
    CALL EXECUTE_SCHEME_IFFT(scaledPotVortXDeriv(:,:,2), FFT_2D, schZP)
    ! Take only the real part to get rid of machine-epsilon errors from the
    ! inverse FFT.
    scaledPotVortXDeriv = REAL(scaledPotVortXDeriv)
    ! Calculate y-derivative of the potential vorticity.
    ALLOCATE(scaledPotVortYDeriv(scaledXLen, scaledYLen, 2))
    scaledPotVortYDeriv = scaledSpecYDeriv * scaledFreqPotVortGrid
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DB(scaledXLen, scaledYLen, scaledPotVortYDeriv(:,:,1))
    !CALL CFFT2DB(scaledXLen, scaledYLen, scaledPotVortYDeriv(:,:,2))
    CALL EXECUTE_SCHEME_IFFT(scaledPotVortYDeriv(:,:,1), FFT_2D, schZP)
    CALL EXECUTE_SCHEME_IFFT(scaledPotVortYDeriv(:,:,2), FFT_2D, schZP)
    ! Take only the real part to get rid of machine-epsilon errors from the
    ! inverse FFT.
    scaledPotVortYDeriv = REAL(scaledPotVortYDeriv)
    ! Calculate the actual Jacobian.
    ALLOCATE(scaledFreqJacobian(scaledXLen, scaledYLen, 2))
    scaledFreqJacobian = scaledStrmfuncXDeriv * scaledPotVortYDeriv &
         & - scaledStrmfuncYDeriv * scaledPotVortXDeriv
    !> EDITTED FOR PARALLEL.
    !CALL CFFT2DF(scaledXLen, scaledYLen, scaledFreqJacobian(:,:,1))
    !CALL CFFT2DF(scaledXLen, scaledYLen, scaledFreqJacobian(:,:,2))
    CALL EXECUTE_SCHEME_FFT(scaledFreqJacobian(:,:,1), FFT_2D, schZP)
    CALL EXECUTE_SCHEME_FFT(scaledFreqJacobian(:,:,2), FFT_2D, schZP)
    ! The larger grid size means the entries of the scaled Jacobian are 9/4 times
    ! larger than they should be, so we must correct that.
    scaledFreqJacobian = CMPLX((4.0_dp/9.0_dp), 0.0_dp, dp) * scaledFreqJacobian

    ! Reduce the Jacobian to the original grid size.
    !> EDITTED FOR PARALLEL.
    !CALL ZERO_PADDING_INV(xLenSub, yLenSub, freqJacobian(:,:,1), &
    !     scaledXLen, scaledYLen, scaledFreqJacobian(:,:,1))
    !CALL ZERO_PADDING_INV(xLenSub, yLenSub, freqJacobian(:,:,2), &
    !     scaledXLen, scaledYLen, scaledFreqJacobian(:,:,2))
    CALL EXECUTE_SCHEME_IZERO_PAD(freqJacobian(:,:,1), sch, &
         scaledFreqJacobian(:,:,1), schZP)
    CALL EXECUTE_SCHEME_IZERO_PAD(freqJacobian(:,:,2), sch, &
         scaledFreqJacobian(:,:,2), schZP)

    ! Deallocate all scaled arrays.
    DEALLOCATE(scaledFreqPotVortGrid)
    DEALLOCATE(scaledFreqStrmfunc)
    DEALLOCATE(scaledStrmfuncXDeriv)
    DEALLOCATE(scaledStrmfuncYDeriv)
    DEALLOCATE(scaledPotVortXDeriv)
    DEALLOCATE(scaledPotVortYDeriv)
    DEALLOCATE(scaledFreqJacobian)

    ! Add in the Jacobian to the output matrix.
    jacobianEkmanShearGrid = jacobianEkmanShearGrid - freqJacobian

  END SUBROUTINE JACOBIAN_EKMAN_SHEAR

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the wavenumber array needed to take the derivative in the
  !! x-direction in frequency space.
  !
  !> Calculates the wavenumbers in the x-direction. If xLenSub even and order
  !! odd, zero out the xLenSub/2 wavenumber. Populate specXDeriv array with
  !! (sqrt(-1)*x wavenumber)**order.
  !
  !> @param[in] xLenSub Size of the grid along the first dimension, local to
  !! subroutine.
  !> @param[in] yLenSub Size of the grid along the second dimension, local to
  !! subroutine.
  !> @param[inout] specXDeriv Array to store coefficients needed for numerical
  !! derivative in spectral space.
  !> @param[in] order Order of the derivative needed.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SPECTRAL_X_DERIVATIVE(xLenSub, yLenSub, specXDeriv, order)

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: xLenSub
    INTEGER(qb), INTENT(IN) :: yLenSub
    INTEGER(qb), INTENT(IN) :: order
    COMPLEX(dp), INTENT(INOUT) :: specXDeriv(xLenSub,yLenSub)
    INTEGER(qb) ::  i !< Counter for DO loops.
    REAL(dp) :: xWavenums(xLenSub) !< Array used to store wavenumbers
    !! in the x-direction.

    ! Calculate the x-direction wavenumbers.
    xWavenums = 0.0_dp
    xWavenums(2) = 1.0_dp
    xWavenums(xLenSub) = -1.0_dp
    DO i = 1_qb, xLenSub/2_qb - 1_qb
       xWavenums(xLenSub - i) = REAL(-i - 1_qb, dp)
       xWavenums(i + 2_qb) = REAL(i + 1_qb, dp)
    END DO
    ! If xLenSub even and order odd, we have to zero out highest wavenumber for
    ! derivative.
    IF ((MOD(xLenSub, 2_qb) .EQ. 0_qb) .AND. (MOD(order, 2_qb) .EQ. 1_qb)) THEN
       xWavenums(xLenSub/2_qb + 1_qb) = 0.0_dp
    END IF

    DO i = 1, yLenSub
       specXDeriv(:,i) = (CMPLX(0.0, xWavenums(:), dp))**(REAL(order, dp))
    END DO

  END SUBROUTINE SPECTRAL_X_DERIVATIVE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the wavenumber array needed to take the derivative in the
  !! y-direction in frequency space.
  !
  !> Calculates the wavenumbers in the y-direction. If yLenSub even and order
  !! odd, zero out the xLenSub/2 wavenumber. Populate specYDeriv array with
  !! (sqrt(-1)*y wavenumber)**order.
  !
  !> @param[in] xLenSub Size of the grid along the first dimension, local to
  !! subroutine.
  !> @param[in] yLenSub Size of the grid along the second dimension, local to
  !! subroutine.
  !> @param[inout] specYDeriv Array to store coefficients needed for numerical
  !! derivative in spectral space.
  !> @param[in] order Order of the derivative needed.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE SPECTRAL_Y_DERIVATIVE(xLenSub, yLenSub, specYDeriv, order)

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: xLenSub
    INTEGER(qb), INTENT(IN) :: yLenSub
    INTEGER(qb), INTENT(IN) :: order
    COMPLEX(dp), INTENT(INOUT) :: specYDeriv(xLenSub,yLenSub)
    INTEGER(qb) ::  i !< Counter for DO loops.
    REAL(dp) :: yWavenums(xLenSub) !< Array used to store wavenumbers
    !! in the x-direction.

    ! Calculate the y-direction wavenumbers.
    yWavenums = 0.0_dp
    yWavenums(2) = 1.0_dp
    yWavenums(xLenSub) = -1.0_dp
    DO i = 1_qb, yLenSub/2_qb - 1_qb
       yWavenums(yLenSub - i) = REAL(-i - 1_qb, dp)
       yWavenums(i + 2_qb) = REAL(i + 1_qb, dp)
    END DO
    ! If xLenSub even and order odd, we have to zero out highest wavenumber for
    ! derivative.
    IF ((MOD(yLenSub, 2_qb) .EQ. 0_qb) .AND. (MOD(order, 2_qb) .EQ. 1_qb)) THEN
       yWavenums(yLenSub/2_qb + 1_qb) = 0.0_dp
    END IF

    DO i = 1, xLenSub
       specYDeriv(i,:) = (CMPLX(0.0, yWavenums(:), dp))**(REAL(order, dp))
    END DO

  END SUBROUTINE SPECTRAL_Y_DERIVATIVE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Zero-pads array going by 3/2s-rule.
  !
  !> Calculates the index in the scaled matrix of the largest negative
  !! wavenumber in both dim1 and dim2. Fills in the non-zero entries of the
  !! scaled matrix (which correspond to the wavenumbers of the original matrix)
  !! with the entries of the original matrix.
  !
  !> @param[in] dim1Len Length along first dimension of array to zero-pad.
  !> @param[in] dim2Len Length along second dimension of array to zero-pad.
  !> @param[in] matrix Array to zero-pad.
  !> @param[in] sclDim1Len Length along first dimension of zero-padded
  !! array.
  !> @param[in] sclDim2Len Length along second dimension of zero-padded
  !! array.
  !> @param[inout] scaledMatrix Zero-padded array.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE ZERO_PADDING(dim1Len, dim2Len, matrix, sclDim1Len, sclDim2Len, &
       scaledMatrix)

    IMPLICIT NONE
    
    INTEGER(qb), INTENT(IN) :: dim1Len
    INTEGER(qb), INTENT(IN) :: dim2Len
    INTEGER(qb), INTENT(IN) :: sclDim1Len
    INTEGER(qb), INTENT(IN) :: sclDim2Len
    COMPLEX(dp), INTENT(IN) :: matrix(dim1Len,dim2Len)
    COMPLEX(dp), INTENT(INOUT) :: scaledMatrix(sclDim1Len,sclDim2Len)
    INTEGER(qb) :: scaledDim1NegWavenumLowIndex !< Lower index indicative the
    !! lowest (most negative) wavenumber along the first dimension for the
    !! zero-padded array.
    INTEGER(qb) :: scaledDim2NegWavenumLowIndex !< Lower index indicative the
    !! lowest (most negative) wavenumber along the second dimension for the
    !! zero-padded array.

    scaledMatrix = (0.0_dp, 0.0_dp)
    scaledDim1NegWavenumLowIndex = 0_qb
    IF (MOD(dim1Len, 2_qb) .EQ. 0_qb) THEN
       scaledDim1NegWavenumLowIndex = dim1Len + 2_qb
    ELSE IF (MOD(dim1Len, 2_qb) .EQ. 1_qb) THEN
       scaledDim1NegWavenumLowIndex = dim1Len + 1_qb
    END IF

    scaledDim2NegWavenumLowIndex = 0_qb
    IF (MOD(dim2Len, 2_qb) .EQ. 0_qb) THEN
       scaledDim2NegWavenumLowIndex = dim2Len + 2_qb
    ELSE IF (MOD(dim2Len, 2_qb) .EQ. 1_qb) THEN
       scaledDim2NegWavenumLowIndex = dim2Len + 1_qb
    END IF

    ! This matrix is scale by 3/2 in each dimension, with the indices corresponding
    ! to the largest wavenumbers (in magnitude) zeroed out.
    scaledMatrix(1:dim1Len/2+1, 1:dim2Len/2+1) = &
         matrix(1:dim1Len/2+1, 1:dim2Len/2+1)
    scaledMatrix(1:dim1Len/2+1, scaledDim2NegWavenumLowIndex:sclDim2Len) = &
         matrix(1:dim1Len/2+1, dim2Len/2+2:dim2Len)
    scaledMatrix(scaledDim1NegWavenumLowIndex:sclDim1Len, 1:dim2Len/2+1) = &
         matrix(dim1Len/2+2:dim1Len, 1:dim2Len/2+1)
    scaledMatrix(scaledDim1NegWavenumLowIndex:sclDim1Len, &
         scaledDim2NegWavenumLowIndex:sclDim2Len) = &
         matrix(dim1Len/2+2:dim1Len, dim2Len/2+2:dim2Len)

  END SUBROUTINE ZERO_PADDING

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Un-zero-pads array going by 3/2s-rule.
  !
  !> Calculates the index in the scaled matrix of the largest negative
  !! wavenumber in both dim1 and dim2. Fills in the non-zero entries of the
  !! scaled matrix (which correspond to the wavenumbers of the original matrix)
  !! with the entries of the original matrix.
  !
  !> @param[in] dim1Len Length along first dimension of array to zero-pad.
  !> @param[in] dim2Len Length along second dimension of array to zero-pad.
  !> @param[inout] matrix Array to zero-pad.
  !> @param[in] sclDim1Len Length along first dimension of zero-padded
  !! array.
  !> @param[in] sclDim2Len Length along second dimension of zero-padded
  !! array.
  !> @param[in] scaledMatrix Zero-padded array.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE ZERO_PADDING_INV(dim1Len, dim2Len, matrix, sclDim1Len, &
       sclDim2Len, scaledMatrix)


    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: dim1Len
    INTEGER(qb), INTENT(IN) :: dim2Len
    INTEGER(qb), INTENT(IN) :: sclDim1Len
    INTEGER(qb), INTENT(IN) :: sclDim2Len
    COMPLEX(dp), INTENT(INOUT) :: matrix(dim1Len,dim2Len)
    COMPLEX(dp), INTENT(IN) :: scaledMatrix(sclDim1Len,sclDim2Len)
    INTEGER(qb) :: scaledDim1NegWavenumLowIndex !< Lower index indicative the
    !! lowest (most negative) wavenumber along the first dimension for the
    !! zero-padded array.
    INTEGER(qb) :: scaledDim2NegWavenumLowIndex !< Lower index indicative the
    !! lowest (most negative) wavenumber along the second dimension for the
    !! zero-padded array.

    matrix = (0.0_dp, 0.0_dp)
    scaledDim1NegWavenumLowIndex = 0_qb
    IF (MOD(dim1Len, 2_qb) .EQ. 0_qb) THEN
       scaledDim1NegWavenumLowIndex = dim1Len + 2_qb
    ELSE IF (MOD(dim1Len, 2_qb) .EQ. 1_qb) THEN
       scaledDim1NegWavenumLowIndex = dim1Len + 1_qb
    END IF

    scaledDim2NegWavenumLowIndex = 0_qb
    IF (MOD(dim2Len, 2_qb) .EQ. 0_qb) THEN
       scaledDim2NegWavenumLowIndex = dim2Len + 2_qb
    ELSE IF (MOD(dim2Len, 2_qb) .EQ. 1_qb) THEN
       scaledDim2NegWavenumLowIndex = dim2Len + 1_qb
    END IF

    ! This matrix is scale by 3/2 in each dimension, with the indices corresponding
    ! to the largest wavenumbers (in magnitude) zeroed out.
    matrix(1:dim1Len/2+1, 1:dim2Len/2+1) = &
         & scaledMatrix(1:dim1Len/2+1, 1:dim2Len/2+1)

    matrix(1:dim1Len/2+1, dim2Len/2+2:dim2Len) = &
         & scaledMatrix(1:dim1Len/2+1, scaledDim2NegWavenumLowIndex:sclDim2Len)

    matrix(dim1Len/2+2:dim1Len, 1:dim2Len/2+1) = &
         & scaledMatrix(scaledDim1NegWavenumLowIndex:sclDim1Len, &
         & 1:dim2Len/2+1)

    matrix(dim1Len/2+2:dim1Len, dim2Len/2+2:dim2Len) = &
         & scaledMatrix(scaledDim1NegWavenumLowIndex:sclDim1Len, &
         & scaledDim2NegWavenumLowIndex:sclDim2Len)

  END SUBROUTINE ZERO_PADDING_INV

END MODULE JACOBIAN_EKMAN_SHEAR_SOLVE
