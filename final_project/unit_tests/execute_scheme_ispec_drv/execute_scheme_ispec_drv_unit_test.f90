!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : EXECUTE_SCHEME_ISPEC_DRV UNIT TEST
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> @brief The unit test for the <tt>EXECUTE_SCHEME_ISPEC_DRV</tt> subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM EXECUTE_SCHEME_ISPEC_DRV_UNIT_TEST

  USE MPI
  USE EZ_PARALLEL
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  REAL(dp), PARAMETER :: pi_dp = 4.0_dp * ATAN(1.0_dp)
  INTEGER(qb) :: rowCount = 4
  INTEGER(qb) :: colCount = 4
  REAL(dp) :: colSpc
  REAL(dp) :: colRef = 0.0
  REAL(dp) :: rowSpc
  REAL(dp) :: rowRef = 0.0
  INTEGER(qb) :: comm = MPI_COMM_WORLD
  INTEGER(qb) :: mpiDatatype = MPI_DOUBLE_COMPLEX
  INTEGER(qb) :: ovlp = 1
  COMPLEX(dp), ALLOCATABLE :: arr(:,:)
  TYPE(SCHEME) :: sch
  INTEGER(qb) :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  INTEGER(qb) :: i, j !< Counters for DO loops.

  CALL MPI_INIT(ierror)

  colSpc = 2.0_dp*pi_dp/REAL(colCount, dp)
  rowSpc = 2.0_dp*pi_dp/REAL(rowCount, dp)
  

  CALL CREATE_SCHEME(rowCount, colCount, colSpc, colRef, comm, mpiDatatype, &
       ovlp, sch)

  ALLOCATE(arr(0:rowCount-1,0:colCount-1))
  arr = 0.0
  DO j = 0, colCount-1
     DO i = 0, rowCount-1
        arr(i,j) = SIN(colRef + j*colSpc)
     END DO
  END DO

  ! Print the sub-grid before boundary communication.
  DO i = 0, sch%commSize-1
     IF (i .EQ. sch%procID) THEN
        PRINT *, " "
        PRINT *, "procID: ", sch%procID
        PRINT *, arr
     END IF
     CALL MPI_BARRIER(sch%comm, ierror)
  END DO

  CALL CREATE_SCHEME_FFT(sch)
  CALL EXECUTE_SCHEME_FFT(arr, FFT_1D_2, sch)
  CALL CREATE_SCHEME_SPEC_DRV(sch)
  CALL EXECUTE_SCHEME_ISPEC_DRV(arr, SPEC_DRV_1D_2, 1, sch)
  CALL EXECUTE_SCHEME_IFFT(arr, FFT_1D_2, sch)

  ! Print the sub-grid after boundary communication.
  PRINT *, " "
  DO i = 0, sch%commSize-1
     IF (i .EQ. sch%procID) THEN
        PRINT *, " "
        PRINT *, "procID: ", sch%procID
        PRINT *, arr
     END IF
     CALL MPI_BARRIER(sch%comm, ierror)
  END DO

  CALL DESTROY_SCHEME(sch)

  CALL MPI_FINALIZE(ierror)
  
END PROGRAM EXECUTE_SCHEME_ISPEC_DRV_UNIT_TEST


