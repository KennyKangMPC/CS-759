!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : TWO_LEVEL_QG_PARALLEL INITIALIZE
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the initialization subroutines needed for the
!! example parallel two-level QG equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

  USE MPI !< ADDED TO PARALLEL.
  USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
  USE EZ_PARALLEL !< ADDED TO PARALLEL.

  IMPLICIT NONE

  PRIVATE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  INTEGER(qb), PUBLIC :: xLen !< Size of physPotVortGrid along the first
  !! dimension.
  INTEGER(qb), PUBLIC :: yLen !< Size of physPotVortGrid along the second
  !! dimension.
  INTEGER(qb), PUBLIC :: numTimesteps !< Number of time-steps to execute.
  INTEGER(qb), PUBLIC :: biharmOrder !< Order of the biharmonic operator for
  !! the hyperviscosity term.
  INTEGER(qb), PUBLIC :: outputFreq !< Output frequency for the simulation.
  !! 1 = output files each time-steps. 0 = no output.
  REAL(dp), PUBLIC :: xRef !< Reference position along the x-dimension for the
  !! initial condition.
  REAL(dp), PUBLIC :: yRef !< Reference position along the y-dimension for the
  !! initial condition.
  REAL(dp), PUBLIC :: dx !< Physical spacing of grid points along the
  !! x-dimension.
  REAL(dp), PUBLIC :: dy !< Physical spacing of grid points along the
  !! y-dimension.
  REAL(dp), PUBLIC :: initDt !< Initial time-step size.
  REAL(dp), PUBLIC :: potVortBound !< Program aborts if potential vorticity
  !! exceeds this bound.
  REAL(dp), PUBLIC :: deformWavenum !< The baroclinic deformation wavenumber
  !! corresponding to the Rossby radius of deformation.
  REAL(dp), PUBLIC :: rotatWavenum !< The wavenumber corresponding to the
  !! rotation coefficient, which controls the advection of streamfunctions.
  REAL(dp), PUBLIC :: vertShear !< Large-scale vertical shear, opposite
  !! direction in each level to induce baroclinic instability.
  REAL(dp), PUBLIC :: ekmanFricCoeff !< The coefficient for Ekman firction in
  !! level 2.
  REAL(dp), PUBLIC :: biharmViscCoeff !< The coefficient for hyperviscosity.
  REAL(dp), PUBLIC :: initTime !< The starting time for the simulation.
  REAL(dp), PUBLIC :: errorToler !< Tolerance parameter for adaptive
  !! time-stepping.
  COMPLEX(dp), ALLOCATABLE, PUBLIC :: physPotVortGrid(:,:,:) !< The potential
  !! voprticity grid in physical space.

  TYPE(SCHEME), PUBLIC :: sch !< <tt>SCHEME</tt> for use of EZ_PARALLEL. ADDED
  !! TO PARALLEL.
  TYPE(SCHEME), PUBLIC :: schZP !< Zero-padded <tt>SCHEME</tt> for use of
  !! EZ_PARALLEL. ADDED TO PARALLEL.
  
  PUBLIC :: INITIALIZE_PARAMETERS
  PUBLIC :: INITIALIZE_GRID

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Reads the NAMELIST for input parameters.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE INITIALIZE_PARAMETERS
    
    IMPLICIT NONE
    
    NAMELIST /model/ xLen, yLen, dx, dy, numTimesteps, initDt, potVortBound, &
         deformWavenum, rotatWavenum, vertShear, ekmanFricCoeff, &
         biharmViscCoeff, biharmOrder, xRef, yRef, initTime, errorToler, &
         outputFreq
    
    OPEN(1000, file = "NAMELIST")
    READ(1000, nml = model)
    CLOSE(1000)

  END SUBROUTINE INITIALIZE_PARAMETERS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Allocates memory to physPotVortGrid and fills in the initial condition.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  SUBROUTINE INITIALIZE_GRID

    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    INTEGER(qb) :: i, j !< Counter for DO loops.

    !> ADDED TO PARALLEL.
    CALL CREATE_SCHEME(xLen, yLen, dy, yRef, MPI_COMM_WORLD, &
         MPI_DOUBLE_COMPLEX, 0, sch)
    CALL CREATE_SCHEME_FFT(sch)
    CALL CREATE_SCHEME_SPEC_DRV(sch)
    CALL CREATE_SCHEME_ZERO_PAD(sch, schZP)

    ! The physical potential vorticity grid contains layer 1 (:,:,1) and
    ! layer 2 (:, :, 2).
    ALLOCATE(physPotVortGrid(0:xLen-1,0:yLen-1,2))
    physPotVortGrid(:,:,:) = (0.0_dp, 0.0_dp)

    ! Fill in the initial condition for the physical potential vorticity grid.
    DO j = 0, yLen-1
       DO i = 0, xLen-1
          physPotVortGrid(i,j,1) = &
               & CMPLX(INITIAL_CONDITION_1(xRef + i*dx, yRef + j*dy), &
               & 0.0_dp, dp)
          physPotVortGrid(i,j,2) = &
               & CMPLX(INITIAL_CONDITION_2(xRef + i*dx, yRef + j*dy), &
               & 0.0_dp, dp)
       END DO
    END DO
    
  END SUBROUTINE INITIALIZE_GRID

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> The initial condition function for level 1.
  !
  !> @param[in] xPos The position along the x-direction.
  !> @param[in] yPos The position along the y-axis.
  !> @param[out] output The value of the initial condition at (xPos, yPos).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  REAL(dp) FUNCTION INITIAL_CONDITION_1(xPos, yPos) RESULT(output)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: xPos
    REAL(dp), INTENT(IN) :: yPos
    
    ! Random initial condition.
    CALL RANDOM_NUMBER(output)
    
  END FUNCTION INITIAL_CONDITION_1

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> The initial condition function for level 2.
  !
  !> @param[in] xPos The position along the x-direction.
  !> @param[in] yPos The position along the y-axis.
  !> @param[out] output The value of the initial condition at (xPos, yPos).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  REAL(dp) FUNCTION INITIAL_CONDITION_2(xPos, yPos) RESULT(output)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: xPos
    REAL(dp), INTENT(IN) :: yPos
    
    ! Random initial condition.
    CALL RANDOM_NUMBER(output)
    
  END FUNCTION INITIAL_CONDITION_2

END MODULE INITIALIZE
