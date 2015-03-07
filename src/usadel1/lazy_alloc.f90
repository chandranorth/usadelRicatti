!! -*-F90-*-
module lazy_alloc
  implicit none
contains
  subroutine lazy_alloc_r1(a, d1)
    implicit none
    double precision, pointer :: a(:)
    integer, intent(in) :: d1
    if (.not. associated(a) .or. size(a,1) /= d1) then
       if (associated(a)) deallocate(a)
       allocate(a(d1))
    end if
  end subroutine lazy_alloc_r1

  subroutine lazy_alloc_r2(a, d1, d2)
    implicit none
    double precision, pointer :: a(:,:)
    integer, intent(in) :: d1, d2
    if (.not. associated(a) .or. size(a,1) /= d1 .or. size(a,2) /= d2) then
       if (associated(a)) deallocate(a)
       allocate(a(d1,d2))
    end if
  end subroutine lazy_alloc_r2

  subroutine lazy_alloc_r3(a, d1, d2, d3)
    implicit none
    double precision, pointer :: a(:,:,:)
    integer, intent(in) :: d1, d2, d3
    if (.not. associated(a) .or. size(a,1) /= d1 .or. size(a,2) /= d2 &
         .or. size(a,3) /= d3) then
       if (associated(a)) deallocate(a)
       allocate(a(d1,d2,d3))
    end if
  end subroutine lazy_alloc_r3

  subroutine lazy_alloc_c1(a, d1)
    implicit none
    double complex, pointer :: a(:)
    integer, intent(in) :: d1
    if (.not. associated(a) .or. size(a,1) /= d1) then
       if (associated(a)) deallocate(a)
       allocate(a(d1))
    end if
  end subroutine lazy_alloc_c1

  subroutine lazy_alloc_i1(a, d1)
    implicit none
    integer, pointer :: a(:)
    integer, intent(in) :: d1
    if (.not. associated(a) .or. size(a,1) /= d1) then
       if (associated(a)) deallocate(a)
       allocate(a(d1))
    end if
  end subroutine lazy_alloc_i1

  subroutine lazy_alloc_i2(a, d1, d2)
    implicit none
    integer, pointer :: a(:,:)
    integer, intent(in) :: d1, d2
    if (.not. associated(a) .or. size(a,1) /= d1 .or. size(a,2) /= d2) then
       if (associated(a)) deallocate(a)
       allocate(a(d1,d2))
    end if
  end subroutine lazy_alloc_i2

  subroutine lazy_alloc_i3(a, d1, d2, d3)
    implicit none
    integer, pointer :: a(:,:,:)
    integer, intent(in) :: d1, d2, d3
    if (.not. associated(a) .or. size(a,1) /= d1 .or. size(a,2) /= d2 &
         .or. size(a,3) /= d3) then
       if (associated(a)) deallocate(a)
       allocate(a(d1,d2,d3))
    end if
  end subroutine lazy_alloc_i3
end module lazy_alloc
