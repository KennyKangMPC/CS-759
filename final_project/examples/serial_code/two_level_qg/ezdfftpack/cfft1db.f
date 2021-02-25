!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     This subroutine is based on the FFT2DC subroutine included in the
!     UCLA LES t12 code, and the name Bjorn Stevens is included as the
!     the author of the documentation.
!
!     Further, this subroutine is built on the FFTPACK Fortran library,
!     developed by Paul N. Schwartztrauber (SW). The double precision
!     version, used here, was converted from the original by Hugh C.
!     Pumphrey, and is available at
!          https://person.hst.aau.dk/magnus/pkgsrc-kolga/math/dfftpack/
!
!     Written By: Jason Turner
!     Last Updated: August 24, 2019
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE CFFT1DB(dim1, matrix)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     This routine computes the backward 1D Fourier transform of a
!     double precision complex dim1 array.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     VARIABLE KEY:
!     - matrix: Input array, of size dim1.
!     - norm: Normalization coefficient.
!     - dim1: Dimension of input array matrix.
!     - WSAVE: Work array to store the prime factorization of
!              dim1 and tabulate trigonometric functions. See
!              documentation on SW's FFT routines.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE

      INTEGER :: dim1
      DOUBLE PRECISION :: norm
      DOUBLE PRECISION, DIMENSION(4*dim1+15) :: WSAVE
      DOUBLE COMPLEX, DIMENSION(dim1) :: matrix

!     SW complex backward 1D FFT (BFFT) initialization, fills WSAVE.
      WSAVE = 0.
      CALL ZFFTI(dim1, WSAVE)

!     Obtain normalization coefficient.
      norm = 1./(DBLE(dim1))

!     Perform 1D BFFT on matrix.
      CALL ZFFTB(dim1, matrix(:), WSAVE) ! SW BFFT routine.

!     Perform normalization.
      matrix = NORM * matrix

      END SUBROUTINE CFFT1DB
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
