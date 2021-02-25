!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : MAIN
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
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
!> @author
!> Jason Turner
!
!> \brief The <tt>EZ_PARALLEL</tt> module.
!> Contains <tt>EZ_PARALLEL</tt> subroutines and their interfaces.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE EZ_PARALLEL
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CREATE_SCHEME
  PUBLIC :: CREATE_SCHEME_FFT
  PUBLIC :: CREATE_SCHEME_SPEC_DRV
  PUBLIC :: CREATE_SCHEME_ZERO_PAD
  PUBLIC :: DESTROY_SCHEME
  PUBLIC :: SHARE_SUBGRID_BDRY
  PUBLIC :: EXECUTE_SCHEME_FFT
  PUBLIC :: EXECUTE_SCHEME_IFFT
  PUBLIC :: EXECUTE_SCHEME_SPEC_DRV
  PUBLIC :: EXECUTE_SCHEME_ISPEC_DRV
  PUBLIC :: EXECUTE_SCHEME_ZERO_PAD
  PUBLIC :: EXECUTE_SCHEME_IZERO_PAD
  PUBLIC :: MAX_VAL
  PUBLIC :: MIN_VAL
  PUBLIC :: CREATE_SPEC_DRV


  ! Interface blocks for each of the EZ_PARALLEL subroutines.

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> creation subroutine.
  !
  !> The <tt>CREATE_SCHEME</tt> subroutine initializes a <tt>SCHEME</tt> which
  !! holds information about the grid and grid decomposition, as well as
  !! variables for use in FFTs.
  !> For greater detail on the <tt>SCHEME</tt> creation subroutine, see
  !! create_scheme.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE CREATE_SCHEME
     SUBROUTINE CREATE_SCHEME_SBR(rowCount, colCount, colSpc, colRef, comm, &
          mpiDatatype, ovlp, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       INTEGER, INTENT(IN) :: rowCount
       INTEGER, INTENT(INOUT) :: colCount
       DOUBLE PRECISION, INTENT(IN) :: colSpc
       DOUBLE PRECISION, INTENT(INOUT) :: colRef
       INTEGER, INTENT(IN) :: comm
       INTEGER, INTENT(IN) :: mpiDatatype
       INTEGER, INTENT(IN) :: ovlp
       TYPE(SCHEME), INTENT(INOUT) :: sch
     END SUBROUTINE CREATE_SCHEME_SBR
  END INTERFACE CREATE_SCHEME

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> FFT initialization subroutine.
  !
  !> The <tt>CREATE_SCHEME_FFT</tt> subroutine initializes a <tt>SCHEME</tt>
  !! for future FFTs.
  !> For greater detail on the <tt>SCHEME</tt> FFT initialization subroutine, see
  !! create_scheme_fft.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE CREATE_SCHEME_FFT
     !> @param[inout] sch <tt>SCHEME</tt> to initialize for FFTs.
     SUBROUTINE CREATE_SCHEME_FFT_SBR(sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       TYPE(SCHEME), INTENT(INOUT) :: sch
     END SUBROUTINE CREATE_SCHEME_FFT_SBR
  END INTERFACE CREATE_SCHEME_FFT

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> spectral derivative
  !! initialization subroutine.
  !
  !> The <tt>CREATE_SCHEME_SPEC_DERV</tt> subroutine initializes a
  !! <tt>SCHEME</tt> for future spectral derivates.
  !> For greater detail on the <tt>SCHEME</tt> spectral derivative
  !! initialization subroutine, see create_scheme_spec_drv.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE CREATE_SCHEME_SPEC_DRV
     !> @param[inout] sch <tt>SCHEME</tt> to initialize for spectral
     !! derivatives.
     SUBROUTINE CREATE_SCHEME_SPEC_DRV_SBR(sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       TYPE(SCHEME), INTENT(INOUT) :: sch
     END SUBROUTINE CREATE_SCHEME_SPEC_DRV_SBR
  END INTERFACE CREATE_SCHEME_SPEC_DRV

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> zero-padding initializtion
  !! subroutine.
  !
  !> The <tt>CREATE_SCHEME_ZERO_PAD</tt> subroutine initializes a
  !! <tt>SCHEME</tt> to handle a zero-padded version of the global array.
  !> For greater detail on the <tt>SCHEME</tt> spectral derivative initialization
  !! subroutine, see create_scheme_zero_pad.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE CREATE_SCHEME_ZERO_PAD
     !> @param[in] sch Scheme associated with arr.
     !> @param[inout] schZP <tt>SCHEME</tt> associated with arrZP.
     SUBROUTINE CREATE_SCHEME_ZERO_PAD_DBLE_SBR(sch, schZP)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       TYPE(SCHEME), INTENT(IN) :: sch
       TYPE(SCHEME), INTENT(INOUT) :: schZP
     END SUBROUTINE CREATE_SCHEME_ZERO_PAD_DBLE_SBR
  END INTERFACE CREATE_SCHEME_ZERO_PAD

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> destruction subroutine.
  !
  !> The <tt>DESTROY_SCHEME</tt> subroutine deallocates a <tt>SCHEME</tt>.
  !> For greater detail on the <tt>SCHEME</tt> destruction subroutine, see
  !! destroy_scheme.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE DESTROY_SCHEME
     !> @param[in] sch <tt>SCHEME</tt> to be destroyed.
     SUBROUTINE DESTROY_SCHEME_SBR(sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       TYPE(SCHEME), INTENT(INOUT) :: sch
     END SUBROUTINE DESTROY_SCHEME_SBR
  END INTERFACE DESTROY_SCHEME

  !> \public
  !! @brief The sub-grid boundary communication interface for the
  !! <tt>EZ_PARALLEL</tt> module.
  !
  !> Communicates the sub-grid boundary to neighboring sub-grids.
  !!
  !> For greater detail on the sub-grid boundary communication subroutine,
  !! see share_subgrid_bdry.f90.
  INTERFACE SHARE_SUBGRID_BDRY
     !> @param[inout] subGrid The sub-grid belonging to the processor.
     !> @param[in] sch <tt>SCHEME</tt> that holds grid decomposition
     !! information, etc.
     SUBROUTINE SHARE_SUBGRID_BDRY_DBLE_SBR(subGrid, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(INOUT) :: subGrid(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE SHARE_SUBGRID_BDRY_DBLE_SBR

     !> @param[inout] subGrid The sub-grid belonging to the processor.
     !> @param[in] sch <tt>SCHEME</tt> that holds grid decomposition
     !! information, etc.
     SUBROUTINE SHARE_SUBGRID_BDRY_DCMPX_SBR(subGrid, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE SHARE_SUBGRID_BDRY_DCMPX_SBR
  END INTERFACE SHARE_SUBGRID_BDRY

  !> \public
  !! @brief The FFT execution interface for the <tt>EZ_PARALLEL</tt> module.
  !
  !> For greater detail on the FFT execution subroutine, see
  !! execute_scheme_fft.f90.
  INTERFACE EXECUTE_SCHEME_FFT
     !> @param[inout] subGrid The local sub-grid whose boundary will be shared.
     !> @param[in] kind Type of FFT to execute, one of FFT_1D_1, FFT_1D_2,
     !! or FFT_2D.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_FFT_DBLE_SBR(subGrid, kind, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_FFT_DBLE_SBR
     !> @param[inout] subGrid The local sub-grid whose boundary will be shared.
     !> @param[in] kind Type of FFT to execute, one of FFT_1D_1, FFT_1D_2,
     !! or FFT_2D.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_FFT_DCMPX_SBR(subGrid, kind, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_FFT_DCMPX_SBR
  END INTERFACE EXECUTE_SCHEME_FFT

  !> \public
  !! @brief The inverse FFT execution interface for the <tt>EZ_PARALLEL</tt>
  !! module.
  !
  !> For greater detail on the inverse FFT execution subroutine, see
  !! execute_scheme_ifft.f90.
  INTERFACE EXECUTE_SCHEME_IFFT
     !> @param[inout] subGrid The local sub-grid whose boundary will be shared.
     !> @param[in] kind Type of FFT to execute, one of FFT_1D_1, FFT_1D_2,
     !! or FFT_2D.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_IFFT_DBLE_SBR(subGrid, kind, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_IFFT_DBLE_SBR
     !> @param[inout] subGrid The local sub-grid whose boundary will be shared.
     !> @param[in] kind Type of FFT to execute, one of FFT_1D_1, FFT_1D_2,
     !! or FFT_2D.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_IFFT_DCMPX_SBR(subGrid, kind, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_IFFT_DCMPX_SBR
  END INTERFACE EXECUTE_SCHEME_IFFT

  !> \public
  !! @brief The spectral derivative execution interface for the
  !! <tt>EZ_PARALLEL</tt> module.
  !
  !> For greater detail on the spectral derivative execution subroutine, see
  !! execute_scheme_spec_drv.f90.
  INTERFACE EXECUTE_SCHEME_SPEC_DRV
     !> @param[inout] subGrid The local sub-grid whose boundary will
     !! be shared.
     !> @param[in] kind Type of spectral derivative to execute, one of
     !! SPEC_DRV_1D_1, SPEC_DERV_1D_2.
     !> @param[in] order Order of the spectral derivative.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_SPEC_DRV_DBLE_SBR(subGrid, kind, order, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       INTEGER, INTENT(IN) :: order
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_SPEC_DRV_DBLE_SBR
     !> @param[inout] subGrid The local sub-grid whose boundary will
     !! be shared.
     !> @param[in] kind Type of spectral derivative to execute, one of
     !! SPEC_DRV_1D_1, SPEC_DERV_1D_2.
     !> @param[in] order Order of the spectral derivative.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_SPEC_DRV_DCMPX_SBR(subGrid, kind, order, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       INTEGER, INTENT(IN) :: order
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_SPEC_DRV_DCMPX_SBR
  END INTERFACE EXECUTE_SCHEME_SPEC_DRV

  !> \public
  !! @brief The inverse spectral derivative execution interface for the
  !! <tt>EZ_PARALLEL</tt> module.
  !
  !> For greater detail on the inverse spectral derivative execution subroutine,
  !! see execute_scheme_ispec_drv.f90.
  INTERFACE EXECUTE_SCHEME_ISPEC_DRV
     !> @param[inout] subGrid The local sub-grid whose boundary will
     !! be shared.
     !> @param[in] kind Type of spectral derivative to execute, one of
     !! SPEC_DRV_1D_1, SPEC_DERV_1D_2.
     !> @param[in] order Order of the spectral derivative.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DBLE_SBR(subGrid, kind, order, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       INTEGER, INTENT(IN) :: order
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DBLE_SBR
     !> @param[inout] subGrid The local sub-grid whose boundary will
     !! be shared.
     !> @param[in] kind Type of spectral derivative to execute, one of
     !! SPEC_DRV_1D_1, SPEC_DERV_1D_2.
     !> @param[in] order Order of the spectral derivative.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DCMPX_SBR(subGrid, kind, order, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(:,:)
       INTEGER, INTENT(IN) :: kind
       INTEGER, INTENT(IN) :: order
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE EXECUTE_SCHEME_ISPEC_DRV_DCMPX_SBR
  END INTERFACE EXECUTE_SCHEME_ISPEC_DRV

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> zero-padding execution
  !! subroutine.
  !
  !> The <tt>EXECUTE_SCHEME_ZERO_PAD</tt> subroutine zero-pads the global array.
  !> For greater detail on the <tt>SCHEME</tt> zero-padding
  !! subroutine, see execute_scheme_zero_pad.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE EXECUTE_SCHEME_ZERO_PAD
     !> @param[in] arr Array to be zero-padded.
     !> @param[in] sch Scheme associated with arr.
     !> @param[inout] arrZP Allocatable array that will hold the zero-padded arr.
     !> @param[in] schZP <tt>SCHEME</tt> associated with arrZP.
     SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DBLE_SBR(arr, sch, arrZP, schZP)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(IN) :: arr(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
       DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arrZP(:,:)
       TYPE(SCHEME), INTENT(IN) :: schZP
     END SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DBLE_SBR
     !> @param[in] arr Array to be zero-padded.
     !> @param[in] sch Scheme associated with arr.
     !> @param[inout] arrZP Allocatable array that will hold the zero-padded arr.
     !> @param[in] schZP <tt>SCHEME</tt> associated with arrZP.
     SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DCMPX_SBR(arr, sch, arrZP, schZP)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(IN) :: arr(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
       DOUBLE COMPLEX, ALLOCATABLE, INTENT(INOUT) :: arrZP(:,:)
       TYPE(SCHEME), INTENT(IN) :: schZP
     END SUBROUTINE EXECUTE_SCHEME_ZERO_PAD_DCMPX_SBR
  END INTERFACE EXECUTE_SCHEME_ZERO_PAD

  !> \public
  !> @brief The interface for the <tt>SCHEME</tt> zero-padding removal
  !! subroutine.
  !
  !> The <tt>EXECUTE_SCHEME_IZERO_PAD</tt> subroutine returns a non-zero-padded
  !! version of the array.
  !> For greater detail on the <tt>SCHEME</tt> zero-padding removal subroutine,
  !! see execute_scheme_izero_pad.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE EXECUTE_SCHEME_IZERO_PAD
     !> @param[inout] arr Array to be un-zero-padded.
     !> @param[in] sch Scheme associated with arr.
     !> @param[in] arrZP Allocatable array that will hold the zero-padded arr.
     !> @param[in] schZP <tt>SCHEME</tt> associated with arrZP.
     SUBROUTINE EXECUTE_SCHEME_IZERO_PAD_DBLE_SBR(arr, sch, arrZP, schZP)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(INOUT) :: arr(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
       DOUBLE PRECISION, INTENT(IN) :: arrZP(:,:)
       TYPE(SCHEME), INTENT(IN) :: schZP
     END SUBROUTINE EXECUTE_SCHEME_IZERO_PAD_DBLE_SBR
     !> @param[inout] arr Array to be un-zero-padded.
     !> @param[in] sch Scheme associated with arr.
     !> @param[in] arrZP Allocatable array that will hold the zero-padded arr.
     !> @param[in] schZP <tt>SCHEME</tt> associated with arrZP.
     SUBROUTINE EXECUTE_SCHEME_IZERO_PAD_DCMPX_SBR(arr, sch, arrZP, schZP)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE COMPLEX, INTENT(INOUT) :: arr(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
       DOUBLE COMPLEX, INTENT(IN) :: arrZP(:,:)
       TYPE(SCHEME), INTENT(IN) :: schZP
     END SUBROUTINE EXECUTE_SCHEME_IZERO_PAD_DCMPX_SBR
  END INTERFACE EXECUTE_SCHEME_IZERO_PAD

  !> \public
  !! @brief The maximum value interface for the <tt>EZ_PARALLEL</tt> module.
  !
  !> For greater detail on the maximum value subroutine, see max_val.f90.
  INTERFACE MAX_VAL
     !> @param[in] subGrid The local sub-grid.
     !> @param[inout] maxVal The variable to store the maximum value across all
     !! sub-grids.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE MAX_VAL_SBR(subGrid, maxValue, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(IN) :: subGrid(:,:)
       DOUBLE PRECISION, INTENT(INOUT) :: maxValue
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE MAX_VAL_SBR
  END INTERFACE MAX_VAL

  !> \public
  !! @brief The minimum value interface for the <tt>EZ_PARALLEL</tt> module.
  !
  !> For greater detail on the minimum value subroutine, see min_val.f90.
  INTERFACE MIN_VAL
     !> @param[in] subGrid The local sub-grid.
     !> @param[inout] minVal The variable to store the minimum value across all
     !! sub-grids.
     !> @param[in] sch <tt>SCHEME</tt> that will be used for the grid
     !! decomposition, parallel FFTs, etc.
     SUBROUTINE MIN_VAL_SBR(subGrid, minValue, sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       DOUBLE PRECISION, INTENT(IN) :: subGrid(:,:)
       DOUBLE PRECISION, INTENT(INOUT) :: minValue
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE MIN_VAL_SBR
  END INTERFACE MIN_VAL

  !> \public
  !> @brief The interface for the spectral derivative array creation
  !! subroutine.
  !
  !> The <tt>CREATE_SPEC_DERV</tt> subroutine fills in an array for future
  !! spectral derivatives.
  !> For greater detail on the spectral derivative array creation
  !! subroutine, see create_spec_drv.f90.
  !> For greater detail on the <tt>SCHEME</tt> datatype, see
  !! ez_parallel_structs.f90.
  INTERFACE CREATE_SPEC_DRV
     !> @param[in] order1 Order of the derivative in the first direction.
     !> @param[in] order2 Order of the derivative in the second dimension.
     !> @param[inout] sch <tt>SCHEME</tt> to initialize for spectral
     !! derivatives.
     !> @param[in] sch <tt>SCHEME</tt> that holds information for grid
     !! decomposition, etc.
     SUBROUTINE CREATE_SPEC_DRV_SBR(order1,order2,specDrv,sch)
       USE MPI
       USE EZ_PARALLEL_STRUCTS
       INTEGER, INTENT(IN) :: order1
       INTEGER, INTENT(IN) :: order2
       DOUBLE COMPLEX, ALLOCATABLE, INTENT(INOUT) :: specDrv(:,:)
       TYPE(SCHEME), INTENT(IN) :: sch
     END SUBROUTINE CREATE_SPEC_DRV_SBR
  END INTERFACE CREATE_SPEC_DRV
  
END MODULE EZ_PARALLEL

INCLUDE 'create_scheme.f90'
INCLUDE 'create_scheme_fft.f90'
INCLUDE 'create_scheme_spec_drv.f90'
INCLUDE 'create_scheme_zero_pad.f90'
INCLUDE 'destroy_scheme.f90'
INCLUDE 'share_subgrid_bdry.f90'
INCLUDE 'execute_scheme_fft.f90'
INCLUDE 'execute_scheme_ifft.f90'
INCLUDE 'execute_scheme_spec_drv.f90'
INCLUDE 'execute_scheme_ispec_drv.f90'
INCLUDE 'execute_scheme_zero_pad.f90'
INCLUDE 'execute_scheme_izero_pad.f90'
INCLUDE 'max_val.f90'
INCLUDE 'min_val.f90'
INCLUDE 'create_spec_drv.f90'
