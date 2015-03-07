! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE KIN_EQUATIONS
  USE PARAMS
  USE MISCMATH
  use lazy_alloc

  INTEGER :: NCOMP, MSTAR
  INTEGER, DIMENSION(:), POINTER :: MSTAROFFSET  ! [NWIRE]
  INTEGER, DIMENSION(:), POINTER :: NCOMPOFFSET  ! [NWIRE]
  double precision, DIMENSION(:), POINTER :: lratio
  
  INTEGER, DIMENSION(:), POINTER :: WIRE_TYPE     ! [NWIRE]
  INTEGER, DIMENSION(:),   POINTER :: BCTYPE      ! [MSTAR]
  double precision, DIMENSION(:,:),  POINTER :: BCDATA ! [MSTAR,16]
  INTEGER, DIMENSION(:,:), POINTER :: BCCONNECT   ! [MSTAR, NWIRE+1]

  double precision, DIMENSION(:), POINTER :: ZETA        ! [MSTAR]

  double precision, PARAMETER :: BIGNUM = 1d3, SMALLNUM = 1d-10

CONTAINS

  SUBROUTINE INIT_EQUATIONS(WIRE_TYPE_NEW)
    IMPLICIT NONE
    INTEGER, DIMENSION(NWIRE), INTENT(IN) :: WIRE_TYPE_NEW
    
    INTEGER :: I
    
    call lazy_alloc_i1(MSTAROFFSET, NWIRE)
    call lazy_alloc_i1(NCOMPOFFSET, NWIRE)
    call lazy_alloc_i1(WIRE_TYPE, NWIRE)

    WIRE_TYPE = WIRE_TYPE_NEW
    
    NCOMP = 0
    MSTAR = 0
    DO I = 1, NWIRE
       MSTAROFFSET(I) = MSTAR
       NCOMPOFFSET(I) = NCOMP
       NCOMP = NCOMP + 2
       MSTAR = MSTAR + 4
    END DO
    
    call lazy_alloc_i1(BCTYPE, MSTAR)
    call lazy_alloc_r2(BCDATA, MSTAR, 16)
    call lazy_alloc_i2(BCCONNECT, MSTAR, NWIRE+1)
    call lazy_alloc_r1(lratio, NWIRE)
    
    lratio = WIRE_LENGTH
  END SUBROUTINE INIT_EQUATIONS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! All the equations are directly pasted from a Mathematica notebook

subroutine indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
  implicit none
  integer, intent(in) :: wire
  integer, intent(out) :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL

  jd2fT = (1 + NCOMPOFFSET(wire))
  jd2fL = (2 + NCOMPOFFSET(wire))
  
  jfT =  (1 + MSTAROFFSET(wire))
  jfL =  (3 + MSTAROFFSET(wire))
  
  jdfT = (2 + MSTAROFFSET(wire))
  jdfL = (4 + MSTAROFFSET(wire))
end subroutine indices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_d(fL,fT,dfL,dfT,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT,&
     d2fL, d2fT)
  implicit none
  double precision, intent(in) :: &
       fL,fT,dfL,dfT,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT
  double precision, intent(out) :: d2fL, d2fT
     
  d2fL = &
       -((dDL*dfL*DT - dfT*DT*dTT + djS*DT*fT + dfT*DT*jS + dDT*dfT*TT &
       + dfL*dTT*TT - cTL*fL*TT + djS*fL*TT - cTT*fT*TT + dfL*jS*TT)&
       /(DL*DT + TT**2))

  d2fT = &
       -((dDT*dfT*DL + dfL*DL*dTT - cTL*DL*fL + djS*DL*fL - cTT*DL*fT &
       + dfL*DL*jS - dDL*dfL*TT + dfT*dTT*TT - djS*fT*TT - dfT*jS*TT) &
       /(DL*DT + TT**2))
end subroutine eval_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_dd(fL,fT,dfL,dfT,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT,&
     d2fL_fL, d2fL_fT, d2fL_dfL, d2fL_dfT, &
     d2fT_fL, d2fT_fT, d2fT_dfL, d2fT_dfT)
  implicit none
  double precision, intent(in) :: &
       fL,fT,dfL,dfT,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT
  double precision, intent(out) :: &
       d2fL_fL, d2fL_fT, d2fL_dfL, d2fL_dfT, &
       d2fT_fL, d2fT_fT, d2fT_dfL, d2fT_dfT
       
  d2fL_fL = &
       ((cTL - djS)*TT)/(DL*DT + TT**2)
  
  d2fL_fT = &
       (-(djS*DT) + cTT*TT)/(DL*DT + TT**2)
  
  d2fL_dfL = &
       -((dDL*DT + (dTT + jS)*TT)/(DL*DT + TT**2))
  
  d2fL_dfT = &
       (DT*dTT - DT*jS - dDT*TT)/(DL*DT + TT**2)
  
  !!----
  
  d2fT_fL = &
       ((cTL - djS)*DL)/(DL*DT + TT**2)
  
  d2fT_fT = &
       (cTT*DL + djS*TT)/(DL*DT + TT**2)
  
  d2fT_dfL = &
       (-(DL*(dTT + jS)) + dDL*TT)/(DL*DT + TT**2)
  
  d2fT_dfT = &
       -((dDT*DL + dTT*TT - jS*TT)/(DL*DT + TT**2))
end subroutine eval_dd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE FSUB(X, Z, F)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, DIMENSION(NCOMP), INTENT(OUT) :: F

    INTEGER :: wire
    integer :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL
    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT
    double precision :: d2fT, d2fL

    F = 0

    CALL GET_KINETIC(X,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)

    DO wire=1,NWIRE
       IF (WIRE_TYPE(wire) .EQ. WIRE_TYPE_NULL) CYCLE

       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF ((DL(wire)*DT(wire) + TT(wire)**2) < SMALLNUM) THEN
          F(jd2fL) = Z(jdfL) * BIGNUM  ! --> \nabla f_L = 0
          F(jd2fT) = Z(jfT) * BIGNUM   ! --> f_T = 0
          CYCLE
       END IF

       IF (WIRE_TYPE(wire) .EQ. WIRE_TYPE_S) THEN
          call eval_d(Z(jfL), Z(jfT), Z(jdfL), Z(jdfT), &
               DL(wire), DT(wire), TT(wire), jS(wire), &
               dDL(wire), dDT(wire), dTT(wire), djS(wire), &
               cTL(wire), cTT(wire), &
               d2fL, d2fT)
       ELSE
          call eval_d(Z(jfL), Z(jfT), Z(jdfL), Z(jdfT), &
               DL(wire), DT(wire), TT(wire), jS(wire), &
               dDL(wire), dDT(wire), dTT(wire), 0d0, &
               0d0, 0d0, &
               d2fL, d2fT)
       END IF

       F(jd2fT) = d2fT
       F(jd2fL) = d2fL
    END DO
  END SUBROUTINE FSUB

  SUBROUTINE DFSUB(X, Z, DF)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, DIMENSION(NCOMP, MSTAR), INTENT(OUT) :: DF

    INTEGER :: wire
    integer :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL
    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT
    double precision  :: &
         d2fL_fL, d2fL_fT, d2fL_dfL, d2fL_dfT, &
         d2fT_fL, d2fT_fT, d2fT_dfL, d2fT_dfT
    
    DF = 0

    CALL GET_KINETIC(X, DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)

    DO wire=1,NWIRE
       IF (WIRE_TYPE(wire) .EQ. WIRE_TYPE_NULL) CYCLE

       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF ((DL(wire)*DT(wire) + TT(wire)**2) < SMALLNUM) THEN
          DF(jd2fL, :) = 0
          DF(jd2fT, :) = 0

          DF(jd2fL, jdfL) = BIGNUM  ! --> \nabla f_L = 0
          DF(jd2fT, jfT)  = BIGNUM  ! --> f_T = 0
          CYCLE
       END IF

       IF (WIRE_TYPE(wire) .EQ. WIRE_TYPE_S) THEN
          call eval_dd(Z(jfL), Z(jfT), Z(jdfL), Z(jdfT), &
               DL(wire), DT(wire), TT(wire), jS(wire), &
               dDL(wire), dDT(wire), dTT(wire), djS(wire), &
               cTL(wire), cTT(wire), &
               d2fL_fL, d2fL_fT, d2fL_dfL, d2fL_dfT, &
               d2fT_fL, d2fT_fT, d2fT_dfL, d2fT_dfT)
       ELSE
          call eval_dd(Z(jfL), Z(jfT), Z(jdfL), Z(jdfT), &
               DL(wire), DT(wire), TT(wire), jS(wire), &
               dDL(wire), dDT(wire), dTT(wire), 0d0, &
               0d0, 0d0, &
               d2fL_fL, d2fL_fT, d2fL_dfL, d2fL_dfT, &
               d2fT_fL, d2fT_fT, d2fT_dfL, d2fT_dfT)
       END IF

       !! d2fT
       DF(jd2fT, jfL)  = d2fT_fL
       DF(jd2fT, jfT)  = d2fT_fT
       DF(jd2fT, jdfL) = d2fT_dfL
       DF(jd2fT, jdfT) = d2fT_dfT
       
       !! d2fL
       DF(jd2fL, jfL)  = d2fL_fL
       DF(jd2fL, jfT)  = d2fL_fT
       DF(jd2fL, jdfL) = d2fL_dfL
       DF(jd2fL, jdfT) = d2fL_dfT
    END DO
  END SUBROUTINE DFSUB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GSUB(I, Z, G)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, INTENT(OUT) :: G

    INTEGER :: J, BCI, BCT, BCS
    INTEGER :: wire, wire2, term

    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT
    integer :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL, &
         jfT2, jfL2, jdfT2, jdfL2

    double precision :: BCX, jL, jT

    BCX = ZETA(I)

    BCS = IAND(ISHFT(BCTYPE(I), -16), 255)
    BCT = IAND(ISHFT(BCTYPE(I), -8), 255)
    BCI = IAND(BCTYPE(I), 255)

    G = 0

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (BCT .EQ. BCTYPE_CLEAN_S_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF (BCI == 1) THEN
          IF (BCS .EQ. BCSUB_CHARACTERISTIC_ZERO) THEN
             G = Z(jfT) - 0
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_L) THEN
             G = Z(jfT) - 0
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_T) THEN
             G = Z(jfT) - 1
          ELSE
             G = Z(jfT) - 0
          END IF
       ELSE IF (BCI == 2) THEN
          IF (BCS .EQ. BCSUB_CHARACTERISTIC_ZERO) THEN
             G = Z(jfL) - 0
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_L) THEN
             G = Z(jfL) - 1
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_T) THEN
             G = Z(jfL) - 0
          ELSE
             G = Z(jfL) - bulk_fL(real(CURRENT_ENERGY), term)
          END IF
       ELSE IF (BCI == 3) THEN
          G = Z(jdfL) - 0
       ELSE
          WRITE(*,*) '%% Invalid BC 1: ', BCI
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_N_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF (BCI == 1) THEN
          IF (BCS .EQ. BCSUB_CHARACTERISTIC_ZERO) THEN
             G = Z(jfT) - 0
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_L) THEN
             G = Z(jfT) - 0
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_T) THEN
             G = Z(jfT) - 1
          ELSE
             G = Z(jfT) - bulk_fT(real(CURRENT_ENERGY), term)
          END IF
       ELSE IF (BCI == 2) THEN
          IF (BCS .EQ. BCSUB_CHARACTERISTIC_ZERO) THEN
             G = Z(jfL) - 0
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_L) THEN
             G = Z(jfL) - 1
          ELSE IF (BCS .EQ. BCSUB_CHARACTERISTIC_T) THEN
             G = Z(jfL) - 0
          ELSE
             G = Z(jfL) - bulk_fL(real(CURRENT_ENERGY), term)
          END IF
       ELSE
          WRITE(*,*) '%% Invalid BC 2'
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
       call indices(wire2, jd2fT, jd2fL, jfT2, jfL2, jdfT2, jdfL2)

       IF (BCI == 1) THEN
          G = Z(jfL) - Z(jfL2)
       ELSE IF (BCI == 2) THEN
          G = Z(jfT) - Z(jfT2)
       ELSE IF (BCI == 3) THEN
          CALL GET_KINETIC(BCX,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

             jT = DT(wire)*Z(jdfT) + TT(wire)*Z(jdfL) + jS(wire)*Z(jfL)
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * jT
          END DO
       ELSE IF (BCI == 4) THEN
          CALL GET_KINETIC(BCX,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

             jL = DL(wire)*Z(jdfL) - TT(wire)*Z(jdfT) + jS(wire)*Z(jfT)
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * jL
          END DO
       ELSE
          WRITE(*,*) '%% Invalid BC 3: BCI=', BCI
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_FREE_INTERFACE) THEN
       wire = BCCONNECT(I, 1)
       CALL GET_KINETIC(BCX,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF (BCI == 1) THEN
          jT = DT(wire)*Z(jdfT) + TT(wire)*Z(jdfL) + jS(wire)*Z(jfL)*0
          G = jT - 0
       ELSE IF (BCI == 2) THEN
          jL = DL(wire)*Z(jdfL) - TT(wire)*Z(jdfT) + jS(wire)*Z(jfT)*1
          G = jL - 0
       ELSE
          WRITE(*,*) '%% Invalid BC 5'
       END IF
    ELSE
       WRITE(*,*) '%% Invalid BC 4'
    END IF
  END SUBROUTINE GSUB

  SUBROUTINE DGSUB(I, Z, DG)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, DIMENSION(MSTAR), INTENT(OUT) :: DG

    INTEGER :: wire, wire2, J, term, BCS, BCI, BCT
    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT
    integer :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL, &
         jfT2, jfL2, jdfT2, jdfL2
    
    double precision :: BCX

    BCX = ZETA(I)

    BCS = IAND(ISHFT(BCTYPE(I), -16), 255)
    BCT = IAND(ISHFT(BCTYPE(I), -8), 255)
    BCI = IAND(BCTYPE(I), 255)


    DG = 0

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (BCT .EQ. BCTYPE_CLEAN_S_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF (BCI == 1) THEN
          DG(jfT) = 1
       ELSE IF (BCI == 2) THEN
          DG(jfL) = 1
       ELSE IF (BCI == 3) THEN
          DG(jdfL) = 1
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_N_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

       IF (BCI == 1) THEN
          DG(jfT) = 1
       ELSE IF (BCI == 2) THEN
          DG(jfL) = 1
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
       call indices(wire2, jd2fT, jd2fL, jfT2, jfL2, jdfT2, jdfL2)

       IF (BCI == 1) THEN
          DG(jfL) = 1
          DG(jfL2) = -1
       ELSE IF (BCI == 2) THEN
          DG(jfT) = 1
          DG(jfT2) = -1
       ELSE IF (BCI == 3) THEN
          CALL GET_KINETIC(BCX,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

             DG(jdfT) = WIRE_CONDUCT(wire)/lratio(wire) * DT(wire)
             DG(jdfL) = WIRE_CONDUCT(wire)/lratio(wire) * TT(wire)
             DG(jfL)  = WIRE_CONDUCT(wire)/lratio(wire) * jS(wire)
          END DO
       ELSE IF (BCI == 4) THEN
          CALL GET_KINETIC(BCX,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
             
             DG(jdfL) = WIRE_CONDUCT(wire)/lratio(wire) * DL(wire)
             DG(jdfT) = WIRE_CONDUCT(wire)/lratio(wire) * (-TT(wire))
             DG(jfT)  = WIRE_CONDUCT(wire)/lratio(wire) * jS(wire)
          END DO
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_FREE_INTERFACE) THEN
       wire = BCCONNECT(I, 1)
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
       CALL GET_KINETIC(BCX,DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)

       IF (BCI == 1) THEN
          DG(jdfT) = DT(wire)
          DG(jdfL) = TT(wire)
          DG(jfL)  = jS(wire)*0
       ELSE IF (BCI == 2) THEN
          DG(jdfL) = DL(wire)
          DG(jdfT) = -TT(wire)
          DG(jfT)  = jS(wire)*1
       END IF
    END IF
  END SUBROUTINE DGSUB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initial guess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GUESS(X, Z, DM)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(OUT) :: Z
    double precision, DIMENSION(NCOMP), INTENT(OUT) :: DM

    INTEGER :: wire
    integer :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL

    DM = 0

    DO wire = 1, NWIRE
       call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
       Z(jfL) = 0
       Z(jfT) = 0
       
       Z(jdfT) = 0
       Z(jdfL) = 0
    END DO
  END SUBROUTINE GUESS


  SUBROUTINE get_vars_from_Z(X, Z, wire, fL, dfL, fT, dfT)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    INTEGER, INTENT(IN) :: wire
    double precision, INTENT(OUT) :: fL, fT, dfL, dfT

    integer :: jd2fT, jd2fL, jfT, jfL, jdfT, jdfL
    
    call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)

    fL  = Z(jfL)
    fT  = Z(jfT)
    dfT = Z(jdfT)
    dfL = Z(jdfL)
  END SUBROUTINE get_vars_from_Z

  SUBROUTINE get_var_indices(wire, jfL, jfT, jdfL, jdfT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: wire
    INTEGER, INTENT(OUT) :: jfL, jfT, jdfL, jdfT
    integer :: jd2fT, jd2fL

    call indices(wire, jd2fT, jd2fL, jfT, jfL, jdfT, jdfL)
  END SUBROUTINE get_var_indices
END MODULE KIN_EQUATIONS
