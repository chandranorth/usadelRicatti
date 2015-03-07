      double precision function ddot  ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer           n, incx, incy
      dimension         x( * ), y( * )

c  ddot   returns the value
c
c     ddot   = x'y
c
c
c  modified nag fortran 77 version of the blas routine ddot  .
c
c  -- written on 21-september-1982.
c     sven hammarling, nag central office.

      integer             i     , ix    , iy
      parameter         ( zero  = 0.0d+0 )

      sum = zero
      if( n.ge.1 )then
         if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               sum = sum + x( ix )*y( ix )
   10       continue
         else
            if( incy.ge.0 )then
               iy = 1
            else
               iy = 1 - ( n - 1 )*incy
            end if
            if( incx.gt.0 )then
               do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                  sum = sum + x( ix )*y( iy )
                  iy  = iy  + incy
   20          continue
            else
               ix = 1 - ( n - 1 )*incx
               do 30, i = 1, n
                  sum = sum + x( ix )*y( iy )
                  ix  = ix  + incx
                  iy  = iy  + incy
   30          continue
            end if
         end if
      end if

      ddot   = sum
      return

*     end of ddot  .

      end function ddot
