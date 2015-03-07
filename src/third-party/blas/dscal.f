      subroutine dscal ( n, alpha, x, incx )
      implicit double precision (a-h,o-z)
      integer          n, incx
      dimension        x( * )

c  dscal  performs the operation
c
c     x := alpha*x
c
c
c  modified nag fortran 77 version of the blas routine dscal .
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            ix
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

      if( n.ge.1 )then
         if( alpha.eq.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
   10       continue
         else if( alpha.eq.( -one ) )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = -x( ix )
   20       continue
         else if( alpha.ne.one )then
            do 30, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = alpha*x( ix )
   30       continue
         end if
      end if

      return

*     end of dscal .

      end subroutine dscal
