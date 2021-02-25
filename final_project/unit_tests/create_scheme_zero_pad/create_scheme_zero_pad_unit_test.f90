!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : CREATE_SCHEME_ZERO_PAD UNIT TEST
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> @brief The unit test for the <tt>CREATE_SCHEME_ZERO_PAD</tt> subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM CREATE_SCHEME_ZERO_PAD_UNIT_TEST

  USE MPI
  USE EZ_PARALLEL
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  REAL(dp), PARAMETER :: pi_dp = 4.0_dp * ATAN(1.0_dp)
  INTEGER(qb) :: rowCount = 5
  INTEGER(qb) :: colCount = 4
  REAL(dp) :: colSpc
  REAL(dp) :: colRef = 0.0
  REAL(dp) :: rowSpc
  REAL(dp) :: rowRef = 0.0
  INTEGER(qb) :: comm = MPI_COMM_WORLD
  INTEGER(qb) :: mpiDatatype = MPI_DOUBLE_PRECISION
  INTEGER(qb) :: ovlp = 1
  INTEGER(qb) :: minWvNmbr
  REAL(dp), ALLOCATABLE :: arr(:,:)
  REAL(dp), ALLOCATABLE :: arrZP(:,:)
  TYPE(SCHEME) :: sch
  TYPE(SCHEME) :: schZP
  INTEGER(qb) :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  INTEGER(qb) :: i, j !< Counters for DO loops.

  CALL MPI_INIT(ierror)

  colSpc = 2.0_dp*pi_dp/REAL(colCount, dp)
  rowSpc = 2.0_dp*pi_dp/REAL(rowCount, dp)

  CALL CREATE_SCHEME(rowCount, colCount, colSpc, colRef, comm, mpiDatatype, &
       ovlp, sch)

  ALLOCATE(arr(0:rowCount-1,0:colCount-1))
  DO j = 0, colCount-1
     DO i = 0, rowCount-1
        arr(i,j) = 100*sch%procId + 10*j + 1*i
     END DO
  END DO

  ! Print the sub-grid before boundary communication.
  DO i = 0, sch%commSize-1
     IF (i .EQ. sch%procID) THEN
        PRINT *, " "
        PRINT *, "procID: ", sch%procID
        PRINT *, SHAPE(arr)
        PRINT *, TRANSPOSE(arr)
     END IF
     CALL MPI_BARRIER(sch%comm, ierror)
  END DO

  CALL CREATE_SCHEME_FFT(sch)
  CALL CREATE_SCHEME_SPEC_DRV(sch)
  CALL CREATE_SCHEME_ZERO_PAD(sch, schZP)

  CALL DESTROY_SCHEME(sch)
  CALL DESTROY_SCHEME(schZP)

  CALL MPI_FINALIZE(ierror)
  
END PROGRAM CREATE_SCHEME_ZERO_PAD_UNIT_TEST


