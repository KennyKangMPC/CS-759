!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : MIN_VAL UNIT TEST
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> @brief The unit test for the <tt>MIN_VAL</tt> subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM MIN_VAL_UNIT_TEST

  USE MPI
  USE EZ_PARALLEL
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  INTEGER(qb) :: rowCount = 3
  INTEGER(qb) :: colCount = 5
  REAL(dp) :: colSpc = 0.1
  REAL(dp) :: colRef = 0.0
  REAL(dp) :: rowSpc = 0.01
  REAL(dp) :: rowRef = 0.0
  REAL(dp) :: minValue
  INTEGER(qb) :: comm = MPI_COMM_WORLD
  INTEGER(qb) :: mpiDatatype = MPI_DOUBLE_PRECISION
  INTEGER(qb) :: ovlp = 2
  REAL(dp), ALLOCATABLE :: arr(:,:)
  TYPE(SCHEME) :: sch
  INTEGER(qb) :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  INTEGER(qb) :: i, j !< Counters for DO loops.

  CALL MPI_INIT(ierror)

  CALL CREATE_SCHEME(rowCount, colCount, colSpc, colRef, comm, mpiDatatype, &
       ovlp, sch)

  ALLOCATE(arr(0:rowCount-1,0:colCount-1))
  arr = 0.0
  DO j = 0, colCount-1
     DO i = 0, rowCount-1
        arr(i,j) = -100*sch%procID - 10*j + 1*i
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

  CALL MIN_VAL(arr, minValue, sch)

  ! Print the minimum value.
  PRINT *, " "
  DO i = 0, sch%commSize-1
     IF (i .EQ. sch%procID) THEN
        PRINT *, " "
        PRINT *, "procID: ", sch%procID, " minValue: ", minValue
     END IF
     CALL MPI_BARRIER(sch%comm, ierror)
  END DO

  CALL MPI_FINALIZE(ierror)
  
END PROGRAM MIN_VAL_UNIT_TEST


