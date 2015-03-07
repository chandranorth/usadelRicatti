! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE INTERPOLATE
  USE MISCMATH

  TYPE Interpolation
     double precision, DIMENSION(:), POINTER :: X => NULL()
     double precision, DIMENSION(:,:), POINTER :: Y => NULL()
     INTEGER :: M, N, ILO
     LOGICAL :: LIMIT_TO_BOUNDARY = .TRUE.
  END TYPE Interpolation

  PRIVATE

  PUBLIC Interpolation, interpolation_new, interpolation_destroy, &
       interpolation_linear_sum, interpolation_get, &
       interpolation_free_array, lower_bound
  
CONTAINS

  SUBROUTINE interpolation_alloc(intp, M, N)
    IMPLICIT NONE
    TYPE(Interpolation), INTENT(INOUT) :: intp
    INTEGER, INTENT(IN) :: M, N

    IF (.NOT. ASSOCIATED(intp%X) .OR. .NOT. ASSOCIATED(intp%Y) &
                          .OR. intp%N /= N .OR. intp%M /= M) THEN
       IF (ASSOCIATED(intp%X)) THEN
          DEALLOCATE(intp%X)
          NULLIFY(intp%X)
       END IF
       IF (ASSOCIATED(intp%Y)) THEN
          DEALLOCATE(intp%Y)
          NULLIFY(intp%Y)
       END IF
       
       intp%M = M
       intp%N = N
       ALLOCATE( intp%X(N), intp%Y(M,N) )
    END IF
  END SUBROUTINE interpolation_alloc

  SUBROUTINE interpolation_new(intp, M, N, X, Y, limit_to_bnd)
    IMPLICIT NONE
    TYPE(Interpolation), INTENT(INOUT) :: intp
    INTEGER, INTENT(IN) :: M, N
    double precision, DIMENSION(M,N), INTENT(IN) :: Y
    double precision, DIMENSION(N), INTENT(IN) :: X
    LOGICAL, OPTIONAL :: limit_to_bnd

    CALL interpolation_alloc(intp, M, N)

    intp%X = X
    intp%Y = Y
    intp%ILO = 1
    IF (PRESENT(limit_to_bnd)) THEN
       intp%LIMIT_TO_BOUNDARY = limit_to_bnd
    ELSE
       intp%LIMIT_TO_BOUNDARY = .TRUE.
    END IF
  END SUBROUTINE interpolation_new

  SUBROUTINE interpolation_destroy(intp)
    IMPLICIT NONE
    TYPE(Interpolation), INTENT(INOUT) :: intp

    IF (ASSOCIATED(intp%X)) THEN
       DEALLOCATE(intp%X)
       NULLIFY(intp%X)
    END IF
    IF (ASSOCIATED(intp%Y)) THEN
       DEALLOCATE(intp%Y)
       NULLIFY(intp%Y)
    END IF
  END SUBROUTINE interpolation_destroy

  SUBROUTINE interpolation_linear_sum(A, B, TO, x)
    IMPLICIT NONE
    TYPE(Interpolation), INTENT(IN) :: A, B
    TYPE(Interpolation), INTENT(INOUT) :: TO
    double precision, INTENT(IN) :: x

    CALL interpolation_alloc(TO, A%M, A%N)

    to%x = a%x
    to%y = (1-x)*a%y + x*b%y
  END SUBROUTINE interpolation_linear_sum

  FUNCTION interpolation_get(intp, X) RESULT(Y)
    IMPLICIT NONE

    TYPE(Interpolation), INTENT(INOUT) :: intp
    double precision, INTENT(IN) :: X

    double precision, DIMENSION(intp%M) :: Y, t

    INTEGER :: ILO, ILEFT, MFLAG

    CALL dintrv(intp%x, intp%n, X, intp%ILO, ILEFT, MFLAG)

    IF (MFLAG == 1) THEN
       ILEFT = ILEFT - 1
    END IF
    IF (MFLAG == -1 .OR. (MFLAG == 1 .AND. X > intp%X(intp%n))) THEN
       IF (intp%LIMIT_TO_BOUNDARY) THEN
          IF (MFLAG == -1) THEN
             Y = intp%y(:,1)
          ELSE
             Y = intp%y(:,intp%n)
          END IF
       ELSE
          PAUSE 'out of interval in interpolation_get'
       END IF
    ELSE
       t = (x - intp%x(ILEFT)) / (intp%x(ILEFT+1) - intp%x(ILEFT))
       Y = (1-t)*intp%y(:,ILEFT) + t*intp%y(:,ILEFT+1)
    END IF
  END FUNCTION interpolation_get

  ! Binary search for i in and ordered array s.t. array(i) <= x < array(i+1)
  ! Returns i=0,   if x < array(i) for all i and
  !         i=N+1, if array(i) <= x for all i.
  SUBROUTINE lower_bound(N, array, x, i)
    USE MISCMATH
    IMPLICIT NONE
    INTEGER, INTENT(in) :: N
    double precision, DIMENSION(N), INTENT(in) :: array
    double precision, INTENT(in) :: x
    INTEGER, INTENT(out) :: i

    INTEGER :: j, k

    IF (x < array(1)) THEN
       i = 0
       RETURN
    ELSE IF (x >= array(N)) THEN
       i = N+1
       RETURN
    END IF

    i=1   
    j=N
    DO 
       k=(i+j)/2   
       IF (x < array(k)) THEN 
          j=k  
       ELSE
          i=k
       END IF
       IF (i+1 >= j) EXIT
    END DO

    ! For debugging...
    IF (.NOT. (array(i) <= x .AND. x < array(i+1))) THEN
       PAUSE '% Binary search error in lower_bound'
    END IF
  END SUBROUTINE lower_bound

END MODULE INTERPOLATE
