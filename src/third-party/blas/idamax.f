      integer function idamax( n, x, incx )
      implicit double precision (a-h,o-z)
      integer         n, incx
      dimension       x( * )

c  idamax returns the smallest value of i such that
c
c     abs( x( i ) ) = max( abs( x( j ) ) )
c                      j
c
c  nag fortran 77 version of the blas routine idamax.
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 31-may-1983.
c     sven hammarling, nag central office.

      intrinsic           abs
      integer             i     , imax  , ix

      if( n.lt.1 )then
         idamax = 0
         return
      end if

      imax = 1
      if( n.gt.1 )then
         xmax = abs( x( 1 ) )
         ix   = 1
         do 10, i = 2, n
            ix = ix + incx
            if( xmax.lt.abs( x( ix ) ) )then
               xmax = abs( x( ix ) )
               imax = i
            end if
   10    continue
      end if

      idamax = imax
      return

*     end of idamax.

      end function idamax
