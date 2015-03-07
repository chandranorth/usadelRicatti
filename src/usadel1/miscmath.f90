! -*-f90-*-

! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
! Miscellaneous mathematical functions
!
MODULE MISCMATH
  IMPLICIT NONE

  double complex, PARAMETER :: imagu = (0d0, 1d0)
  double precision, PARAMETER :: PI = 3.14159265353d0
CONTAINS
  
  double precision FUNCTION intervalpoint(i,N,xmin,xmax)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: N, i
    double precision, INTENT(in) :: xmin, xmax
    
    IF (N == 1) THEN
       intervalpoint = xmin
    ELSE
       intervalpoint = xmin + REAL(i-1) * (xmax-xmin) / REAL(N-1)
    END IF
  END FUNCTION intervalpoint

  FUNCTION artanh(x)
    IMPLICIT NONE
    double complex :: x, artanh
    artanh = .5*LOG((1+x)/(1-x))
  END FUNCTION artanh

  FUNCTION ctanh(x)
    IMPLICIT NONE
    double complex :: x, ctanh
    ctanh = -imagu * sin(imagu*x)/cos(imagu*x)
  END FUNCTION ctanh

  FUNCTION csinh(x)
    IMPLICIT NONE
    double complex :: x, csinh
    csinh = -imagu * sin(imagu*x)
  END FUNCTION csinh

  FUNCTION ccosh(x)
    IMPLICIT NONE
    double complex :: x, ccosh
    ccosh = cos(imagu*x)
  END FUNCTION ccosh

  function ccsch(x)
    IMPLICIT NONE
    double complex :: x, ccsch
    ccsch = 1d0/csinh(x)
  end function ccsch

  function ccoth(x)
    IMPLICIT NONE
    double complex :: x, ccoth
    ccoth = 1d0/ctanh(x)
  end function ccoth

  function contains_nan_0(x) result(y)
    double precision, intent(in) :: x
    logical :: y
    y = .not. (x >= 0 .or. x < 0)
  end function contains_nan_0
  
  function contains_nan_1(x) result(y)
    double precision, intent(in), dimension(:) :: x
    logical :: y
    integer :: i1
    y = .false.
    do i1=1, size(x)
       if (contains_nan_0(x(i1))) y = .true.
    end do
  end function contains_nan_1

  function contains_nan_2(x) result(y)
    double precision, intent(in), dimension(:,:) :: x
    logical :: y
    integer :: i1, i2
    y = .false.
    do i1=1, size(x,1)
       do i2=1, size(x,2)
          if (contains_nan_0(x(i1,i2))) y = .true.
       end do
    end do
  end function contains_nan_2


  !! Calculate the intersection point of two arbitrary planes and the plane z=0
  !    x, y, z1    Three points on the first plane [in]
  !    x, y, z2    Three points on the second plane [in]
  !    ix, iy      The intersection point [out]
  SUBROUTINE intersection_of_planes_and_zero(x,y,z1,z2,ix,iy)
    IMPLICIT NONE
    double precision, DIMENSION(3), INTENT(in) :: x,y,z1,z2
    double precision, INTENT(out) :: ix,iy

    double precision, DIMENSION(2) :: a,b,c
    double precision, DIMENSION(3,3) :: coeff
    double precision :: detinv

    ! Calculate plane constants
    detinv = ((x(3)-x(2))*y(1) + (x(1)-x(3))*y(2) + (x(2)-x(1))*y(3))
    IF (detinv .EQ. 0) PAUSE '% Singular matrix in intersection_of_planes_and_zero'

    coeff = 1d0 / detinv * TRANSPOSE(RESHAPE((/ &
         y(2)-y(3), y(3)-y(1), y(1)-y(2), &
         x(3)-x(2), x(1)-y(3), x(2)-x(1), &
         x(2)*y(3)-x(3)*y(2), x(3)*y(1)-x(1)*y(3), x(1)*y(2)-x(2)*y(1) /), &
      (/ 3, 3 /)))

    a(1) = coeff(1,1)*z1(1) + coeff(1,2)*z1(2) + coeff(1,3)*z1(3)
    b(1) = coeff(2,1)*z1(1) + coeff(2,2)*z1(2) + coeff(2,3)*z1(3)
    c(1) = coeff(3,1)*z1(1) + coeff(3,2)*z1(2) + coeff(3,3)*z1(3)

    a(2) = coeff(1,1)*z2(1) + coeff(1,2)*z2(2) + coeff(1,3)*z2(3)
    b(2) = coeff(2,1)*z2(1) + coeff(2,2)*z2(2) + coeff(2,3)*z2(3)
    c(2) = coeff(3,1)*z2(1) + coeff(3,2)*z2(2) + coeff(3,3)*z2(3)

    ! Calculate intersection
    detinv = a(1)*b(2) - a(2)*b(1)
    IF (detinv .EQ. 0) PAUSE '% Parallel etc. planes in intersection_of_planes_and_zero'
    
    ix = ( b(1)*c(2) - b(2)*c(1) ) / detinv
    iy = ( a(2)*c(1) - a(1)*c(2) ) / detinv
    
  END SUBROUTINE intersection_of_planes_and_zero
END MODULE MISCMATH
