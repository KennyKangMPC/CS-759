! Integer precision types: single byte, double bytes, quad bytes, oct bytes.
! E.g., sb can hold up to 10^2, db can hold up to 10^4.
INTEGER, PARAMETER :: sb = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: db = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: qb = SELECTED_INT_KIND(8)
INTEGER, PARAMETER :: ob = SELECTED_INT_KIND(16)
