!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     This subroutine is based on the FFT2DC subroutine included in the
!     UCLE LES t12 code, and the name Bjorn Stevens is included as the
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
      SUBROUTINE CFFT2DF(dim1, dim2, matrix)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     This routine computes the forward 2D Fourier transform of a
!     double precision complex dim1-by-dim2 array.
!
!     We begin by calculating the Fourier transform of the rows of
!     matrix, transpose matrix to obtain matrix_transpose, then
!     calculate the Fourier transform of the rows of matrix_transpose,
!     and transposing back.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     VARIABLE KEY:
!     - matrix: Input array, of size dim1-by-dim2.
!     - dim1, dim2: Dimensions of input array matrix.
!     - WSAVE1, WSAVE2: Work array to store the prime factorization of
!                       dim1, dim2 and tabulate trigonometric
!                       functions. See documentation on SW's FFT
!                       routines.
!     - ii: Counter for DO loops.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE

      INTEGER :: dim1, dim2, ii
      DOUBLE PRECISION, DIMENSION(4*dim1+15) :: WSAVE1
      DOUBLE PRECISION, DIMENSION(4*dim2+15) :: WSAVE2
      DOUBLE COMPLEX, DIMENSION(dim1,dim2) :: matrix
      DOUBLE COMPLEX, DIMENSION(dim2,dim1) :: matrix_transpose

!     SW complex forward 1D FFT (FFFT) initialization, fills WSAVE1.
      WSAVE1 = 0.
      CALL ZFFTI(dim1, WSAVE1)

!     Perform 1D FFFT on the rows of matrix.
      DO ii = 1, dim2
        CALL ZFFTF(dim1, matrix(:,ii), WSAVE1) ! SW FFFT routine.
      END DO
      matrix_transpose = TRANSPOSE(matrix)

      WSAVE2 = 0.
      CALL ZFFTI(dim2, WSAVE2)
!     SW 1D (FFFT) initialization, fills WSAVE2.

!     Perform 1D FFFT on the columns of matrix_transpose.
      DO ii = 1, dim1
        CALL ZFFTF(dim2, matrix_transpose(:,ii), WSAVE2)
      END DO
      matrix = TRANSPOSE(matrix_transpose)

      END SUBROUTINE CFFT2DF
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
