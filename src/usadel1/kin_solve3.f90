! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE KIN_SOLVE3
  USE PARAMS
  USE KIN_EQUATIONS
  use lazy_alloc
  
  PUBLIC SOLVE, INIT_SOLVER, GET_VARS

  PRIVATE

  INTEGER :: NRWTOP, NRWBOT, NOVRLP
  INTEGER :: NRWBLK, NCLBLK, NBLOKS
  INTEGER :: N, NX

  double precision, DIMENSION(:,:), POINTER :: TOPBLK, BOTBLK
  double precision, DIMENSION(:,:,:), POINTER :: ARRAY

  INTEGER, DIMENSION(:), POINTER :: PIVOT

  double precision, DIMENSION(:), POINTER :: B
  double precision, DIMENSION(:), POINTER :: XX

CONTAINS

  SUBROUTINE INIT_SOLVER(NREC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NREC

    NX = WIRE_coefs(1)%N
    !NX = 2

    NRWTOP = MSTAR - NREC
    NRWBOT = NREC

    NOVRLP = MSTAR

    NRWBLK = MSTAR
    NCLBLK = 2*MSTAR

    NBLOKS = NX - 1

    N = NBLOKS * NRWBLK + NOVRLP

    call lazy_alloc_r3(ARRAY, NRWBLK, NCLBLK, NBLOKS)
    call lazy_alloc_r2(TOPBLK, NRWTOP, NOVRLP)
    call lazy_alloc_r2(BOTBLK, NRWBOT, NOVRLP)
    call lazy_alloc_i1(PIVOT, N)
    call lazy_alloc_r1(B, N)
    call lazy_alloc_r1(XX, NX)

    XX = WIRE_coefs(1)%X

    ! Boundary condition positions
    call lazy_alloc_r1(ZETA, MSTAR)
    ZETA = 1d0
    ZETA(1:(MSTAR - NREC)) = 0d0
  END SUBROUTINE INIT_SOLVER

  SUBROUTINE SOLVE(CONVERGED)
    USE TWPBVPC_MOD
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: CONVERGED

    INTEGER :: K, M, IFLAG, JOB
    double precision :: h, x

    double precision, DIMENSION(MSTAR) :: Z

    CALL INTERPOLATE_KINETIC_FOR_ENERGY(real(CURRENT_ENERGY))

    !!! Construct the Jacobian in the COLROW form
    B = 0
    ARRAY = 0
    Z = 0

    ! left bc
    DO K = 1, NRWTOP
       CALL GSUB(K, Z, B(K))
       CALL DGSUB(K, Z, TOPBLK(K, :))
    END DO
    B(1:NRWTOP) = -B(1:NRWTOP)

    ! right bc
    DO K = 1, NRWBOT
       CALL GSUB(NRWTOP + K, Z, B(N - NRWBOT + K))
       CALL DGSUB(NRWTOP + K, Z, BOTBLK(K, :))
    END DO
    B((N-NRWBOT):N) = -B((N-NRWBOT):N)

    ! equations
    DO K = 1, NBLOKS
       x = .5*(XX(K+1) + XX(K))
       h = XX(K+1) - XX(K)
       CALL DFSUB(x, Z, ARRAY(1:NCOMP, 1:MSTAR, K))

       ARRAY(1:NCOMP, 1:MSTAR, K) = .5*h*ARRAY(1:NCOMP, 1:MSTAR, K)
       ARRAY(1:NCOMP, (MSTAR+1):(2*MSTAR), K) = ARRAY(1:NCOMP, 1:MSTAR, K)

       DO M = 1, NCOMP
          ! Second derivative
          ARRAY(M, 2*(M-1) + 2, K) = ARRAY(M, 2*(M-1) + 2, K) + 1
          ARRAY(M, 2*(M-1) + 2 + MSTAR, K) = ARRAY(M, 2*(M-1) + 2 + MSTAR, K) - 1

          ! First derivative
          ARRAY(M+NCOMP, 2*(M-1) + 1, K)         =  1
          ARRAY(M+NCOMP, 2*(M-1) + 1 + MSTAR, K) = -1

          ARRAY(M+NCOMP, 2*(M-1) + 2, K)         = .5*h
          ARRAY(M+NCOMP, 2*(M-1) + 2 + MSTAR, K) = .5*h
       END DO
    END DO

    ! Solve the equations
    JOB = 0
    CALL COLROW (N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK, &
         NBLOKS,BOTBLK,NRWBOT,PIVOT,B,IFLAG,JOB)
    
    IF (IFLAG .NE. 0) THEN
       WRITE(0,*) '%% Failed to solve: ', IFLAG
       CONVERGED = .FALSE.
    ELSE
       CONVERGED = .TRUE.
    END IF
  END SUBROUTINE SOLVE

  SUBROUTINE DUMPMAT(M, N, A)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M, N
    double precision, INTENT(IN), DIMENSION(*) :: A
    INTEGER :: I, J
    WRITE(*,*) '|---'
    DO I = 1, M
       DO J = 1, N
          WRITE(*,'(1X,G7.1,1X)',ADVANCE='NO') A(I + M*(J-1))
       END DO
       WRITE(*,'(1X)')
    END DO
  END SUBROUTINE DUMPMAT

  SUBROUTINE get_vars(X, fL, dfL, fT, dfT, jL, jT)
    USE INTERPOLATE
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(in) :: X
    double precision, DIMENSION(:,:), INTENT(out) :: fL, fT, dfL, dfT, &
         jL, jT
    INTEGER :: ix, wire
    TYPE(Interpolation), SAVE :: u_interp
    
    double precision, DIMENSION(MSTAR) :: z
    
    INTEGER :: jjfL, jjfT, jjdfL, jjdfT

    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT

    CALL interpolation_new(u_interp, MSTAR, NX, XX, &
         RESHAPE(B, (/ MSTAR, NX /)))

    DO ix = 1, SIZE(X)
       z = interpolation_get(u_interp, X(ix))
       CALL GET_KINETIC(X(ix),DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
       DO wire = 1, nwire
          CALL get_vars_from_Z(X(ix), z, wire, &
               fL(ix,wire), dfL(ix,wire), fT(ix,wire), dfT(ix,wire))

          ! Scale lengths properly
          dfL(ix,wire) = dfL(ix,wire) / WIRE_LENGTH(wire)
          dfT(ix,wire) = dfT(ix,wire) / WIRE_LENGTH(wire)

          ! Evaluate currents
          jS(wire) = jS(wire) / WIRE_LENGTH(wire)
          jT(ix,wire) = DT(wire)*dfT(ix,wire) + TT(wire)*dfL(ix,wire) + jS(wire)*fL(ix,wire)
          jL(ix,wire) = DL(wire)*dfL(ix,wire) - TT(wire)*dfT(ix,wire) + jS(wire)*fT(ix,wire)
       END DO
    END DO
  END SUBROUTINE get_vars
END MODULE KIN_SOLVE3
