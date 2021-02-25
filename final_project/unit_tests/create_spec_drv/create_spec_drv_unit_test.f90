!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : CREATE_SPEC_DRV UNIT TEST
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> @brief The unit test for the <tt>CREATE_SPEC_DRV</tt> subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM CREATE_SPEC_DRV_UNIT_TEST

  USE MPI
  USE EZ_PARALLEL
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  INTEGER(qb) :: rowCount = 5
  INTEGER(qb) :: colCount = 8
  REAL(dp) :: colSpc = 0.1
  REAL(dp) :: colRef = 0
  INTEGER(qb) :: comm = MPI_COMM_WORLD
  INTEGER(qb) :: mpiDatatype = MPI_DOUBLE_COMPLEX
  INTEGER(qb) :: ovlp = 1
  TYPE(SCHEME) :: sch
  INTEGER(qb) :: i !< Counter for DO loops.
  INTEGER(qb) :: ierror !< Additional integer argument for <tt>FORTRAN</tt>
  !! <tt>MPI</tt> subroutine calls.
  COMPLEX(dp), ALLOCATABLE :: specDrv(:,:) !< Array to hold spectral derivative
  !! operator.

  CALL MPI_INIT(ierror)

  CALL CREATE_SCHEME(rowCount, colCount, colSpc, colRef, comm, mpiDatatype, &
       ovlp, sch)

  CALL CREATE_SCHEME_FFT(sch)

  CALL CREATE_SCHEME_SPEC_DRV(sch)

  CALL CREATE_SPEC_DRV(1,1,specDrv,sch)

  DO i = 0, sch%commSize-1
     IF (i .EQ. sch%procID) THEN
        PRINT *, "procID: ", sch%procID
        PRINT *, SHAPE(specDrv)
        PRINT *, TRANSPOSE(specDrv)
     END IF
     CALL MPI_BARRIER(comm, ierror)
  END DO

  CALL DESTROY_SCHEME(sch)

  CALL MPI_FINALIZE(ierror)
  
END PROGRAM CREATE_SPEC_DRV_UNIT_TEST


