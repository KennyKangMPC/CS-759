!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STOCH_HEAT_PARALLEL INITIALIZE
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the initialization subroutines needed for the
!! example parallel stochastic heat equation code.
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

  INTEGER(qb), PUBLIC :: xLen !< Size of tempGrid along the first dimension.
  INTEGER(qb), PUBLIC :: yLen !< Size of tempGrid along the second dimension.
  INTEGER(qb), PUBLIC :: outputFreq !< Output frequency for the simulation.
  !! 1 = output a file each time-step. 0 = no output.
  REAL(dp), PUBLIC :: xRef !< Reference position along the x-dimension for the
  !! initial condition.
  REAL(dp), PUBLIC :: yRef !< Reference position along the y-dimension for the
  !! initial condition.
  REAL(dp), PUBLIC :: dx !< Physical spacing of grid points along the
  !! x-dimension.
  REAL(dp), PUBLIC :: dy !< Physical spacing of grid points along the
  !! y-dimension.
  REAL(dp), PUBLIC :: dt !< Time-step size.
  REAL(dp), PUBLIC :: xDiffus !< Heat diffusivity along the x-dimension.
  REAL(dp), PUBLIC :: yDiffus !< Heat diffusivity along the y-dimension.
  REAL(dp), PUBLIC :: deterForce !< Magnitude of the deterministic forcing.
  REAL(dp), PUBLIC :: dampCoeff !< Magnitude of the dampening coefficient.
  REAL(dp), PUBLIC :: stochMag !< Magnitude of the stochastic noise.
  REAL(dp), PUBLIC :: relaxTrgt !< Relaxation target.
  REAL(dp), PUBLIC :: dirBdyVal !< Boundary value for the Dirichlet boundary
  !! condition.
  REAL(dp), PUBLIC :: initTime !< The starting time for the simulation.
  REAL(dp), PUBLIC :: finTime !< The ending time for the simulation.
  REAL(dp), ALLOCATABLE, PUBLIC :: tempGrid(:,:,:) !< Array for storing the
  !! temperature at each grid point. The last dimension is used index two
  !! time-steps.

  TYPE(SCHEME), PUBLIC :: sch !< <tt>SCHEME</tt> For use of EZ_PARALLEL. ADDED
  !! TO PARALLEL.
  
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

    NAMELIST /model/ xLen, yLen, dx, dy, dt, xDiffus, yDiffus, deterForce, &
         dampCoeff, stochMag, relaxTrgt, dirBdyVal, xRef, yRef, initTime, &
         finTime, outputFreq

    OPEN(1000, file = "NAMELIST")
    READ(1000, nml = model)
    CLOSE(1000)

  END SUBROUTINE INITIALIZE_PARAMETERS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Allocates memory to tempGrid and fills in the initial condition.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_GRID

    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    INTEGER(qb) :: i, j !< Counters for DO loops.

    CALL CREATE_SCHEME(xLen, yLen, dy, yRef, MPI_COMM_WORLD, &
         MPI_DOUBLE_PRECISION, 1, sch) !< ADDED TO PARALLEL.
    
    ! tempGrid contains two time steps (third index).
    ALLOCATE(tempGrid(0:xLen-1,0:yLen-1,0:1))
    tempGrid(:,:,0) = 0.0_dp
    
    ! Fill in the initial condition, based on the INITIAL_CONDITION function
    ! defined below.
    DO j = 1, yLen-2
       DO i = 1, xLen-2
          tempGrid(i,j,0) = INITIAL_CONDITION(xRef + i*dx, yRef + j*dy)
       END DO
    END DO
    
    ! Set boundary value to a fixed value.
    ! Faster to assign value to array which is contiguous in memory all at once.
    tempGrid(:,0,0) = dirBdyVal
    tempGrid(:,yLen-1,0) = dirBdyVal
    ! Faster to assign value to array which is not contiguous in memory with DO.
    DO j = 0, yLen-1
       tempGrid(0,j,0) = dirBdyVal
       tempGrid(xLen-1,j,0) = dirBdyVal
    END DO

    CALL SHARE_SUBGRID_BDRY(tempGrid(:,:,0), sch) !< ADDED TO PARALLEL.
    
    tempGrid(:,:,1) = tempGrid(:,:,0)

  END SUBROUTINE INITIALIZE_GRID

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> The initial condition function.
  !
  !> @param[in] xPos The position along the x-direction.
  !> @param[in] yPos The position along the y-axis.
  !> @param[out] output The value of the initial condition at (xPos, yPos).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  REAL(dp) FUNCTION INITIAL_CONDITION(xPos, yPos) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: xPos
    REAL(dp), INTENT(IN) :: yPos

    ! Random initial condition.
    CALL RANDOM_NUMBER(output)
    output = relaxTrgt + 5.0_dp * output

  END FUNCTION INITIAL_CONDITION

END MODULE INITIALIZE
