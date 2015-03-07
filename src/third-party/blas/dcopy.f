      subroutine dcopy ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer            n, incx, incy
      dimension          x( * ), y( * )

c  dcopy  performs the operation
c
c     y := x
c
c  nag fortran 77 version of the blas routine dcopy .
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy

      if( n.lt.1 )return

      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy
            y( iy ) = x( iy )
   10    continue
      else
         if( incx.ge.0 )then
            ix = 1
         else
            ix = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( ix )
               ix      = ix + incx
   20       continue
         else
            iy = 1 - ( n - 1 )*incy
            do 30, i = 1, n
               y( iy ) = x( ix )
               iy      = iy + incy
               ix      = ix + incx
   30       continue
         end if
      end if

      return

*     end of dcopy .

      end subroutine dcopy
