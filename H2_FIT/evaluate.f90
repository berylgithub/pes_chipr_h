!######################################################################################################
! caller (evaluator):
!% For diatomics:
!CALL CHIPR_DIAT(R,POT,DER,DVDR) 
!% For three-body energies
!CALL CHIPR_TRIAT(R,POT,DER,DVDR) 
!% For four-body energies
!CALL CHIPR_TETRA(R,POT,DER,DVDR)

program evaluator
    IMPLICIT NONE  
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R_test_array, V_test_array, V_pred_array
    DOUBLE PRECISION :: DVDR = 0, rmse_var
    LOGICAL :: DER = .false.
    INTEGER :: i, nlines

    !allocate vars:
    nlines = 0
    open (1, file = 'test_data.txt')
    do
        read (1,*, END=10)
        nlines = nlines + 1
    end do
    10 close (1)
    ALLOCATE(R_test_array(nlines))
    ALLOCATE(V_test_array(nlines))
    ALLOCATE(V_pred_array(nlines))

    !load test data:
    open (1, file = 'test_data.txt')
    do i=1, nlines
        read(1,*) R_test_array(i), V_test_array(i)
    end do
    close (1)

    !evaluate test data points:
    do i=1, nlines
        call CHIPR_DIAT(R_test_array(i), V_pred_array(i), DER, DVDR)
    end do

    !compute RMSE:
    rmse_var = RMSE(V_test_array, V_pred_array)

    !save evaluated points to file:
    open(1, file = "eval_data.txt")
    do i=1, nlines
        write(1, *) R_test_array(i), V_pred_array(i)
    end do
    write(1, *) "#RMSE:"
    write(1, *) rmse_var
    close(1)


    contains
        FUNCTION RMSE(Y, Y_pred) RESULT(rmse_var)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:) :: Y, Y_pred
            DOUBLE PRECISION :: rmse_var, N

            N = size(Y)
            rmse_var = SQRT(SUM((Y-Y_pred)**2)/N)
        END FUNCTION RMSE

end program evaluator

!######################################################################################################

      SUBROUTINE CHIPR_DIAT(R,POT,DER,DVDR)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=31
      INTEGER :: I,J
      INTEGER :: BSORDER
      INTEGER :: POLORDER
      INTEGER :: NCBAS,NCPOL
      DOUBLE PRECISION, DIMENSION(2) :: Z
      DOUBLE PRECISION, DIMENSION(NC) :: C
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BS
      DOUBLE PRECISION :: Y   
      DOUBLE PRECISION :: R,POT
      LOGICAL :: DER
      DOUBLE PRECISION :: DYDR
      DOUBLE PRECISION :: DVDRPART1
      DOUBLE PRECISION :: DVDRPART2
      DOUBLE PRECISION :: DVDR

      Z(  1)=  0.10000000000000000000000000000000000D+01
      Z(  2)=  0.10000000000000000000000000000000000D+01

      C(  1)=  0.11050783889712387786374847564729862D+01
      C(  2)=  0.65771041832943666349819977767765522D+00
      C(  3)= -0.35232136145224129952779890118108597D+00
      C(  4)=  0.21388582632944883599179775046650320D+01
      C(  5)=  0.80933503975620113557454260444501415D+00
      C(  6)= -0.94417555656220442106274504112661816D+00
      C(  7)=  0.80162668074669785056585169513709843D+01
      C(  8)=  0.10533350316412120051268175302539021D+02
      C(  9)= -0.12036682619703291674539968880708329D+01
      C( 10)=  0.42384100718910202587608182511758059D+01
      C( 11)=  0.45073503370892469277464442711789161D+01
      C( 12)= -0.17837370886131186153988892328925431D+02
      C( 13)= -0.14117518989782993799053656402975321D+02
      C( 14)= -0.18481454464664265335827053604589310D+00
      C( 15)=  0.53068799183913129002831965408404358D+00
      C( 16)=  0.93656577710856925289562013858812861D+00
      C( 17)= -0.14462118474956997538072300812928006D+00
      C( 18)=  0.48313463053277350134351308952318504D+00
      C( 19)= -0.91513098501637391013474598366883583D+00
      C( 20)=  0.17902803665316935344264948071213439D+00
      C( 21)= -0.16016776574749681572029658127576113D+03
      C( 22)= -0.12166918754369115962532532648765482D+01
      C( 23)=  0.21795068532519525916768543538637459D+01
      C( 24)=  0.65673569974377754565608711345703341D+00
      C( 25)= -0.79323035562410293408674988313578069D+01
      C( 26)= -0.74240409795840225370966436457820237D+01
      C( 27)= -0.17656785429075322888577831426104581D-04
      C( 28)=  0.24882349188198209510858305293368176D+01
      C( 29)= -0.37504007621975169968209229409694672D+04
      C( 30)= -0.16848334068981456468350188515614718D+01
      C( 31)= -0.14210502161179198316043326144608727D-01

      BSORDER=8

      POLORDER=13

      NCBAS=2*BSORDER+2

      NCPOL=POLORDER

      ALLOCATE(BS(NCBAS))

      DO I=1,NCBAS
        BS(I)=0.00D+00
      END DO

      DO I=1,NCBAS
        J=NCPOL+I
        BS(I)=C(J)
      END DO
 
      Y=0.00D+00
      POT=0.00D+00

      CALL BASIS_CONTRACT(1,BSORDER,BS,R,Y)

      DO I=1,POLORDER
        POT=POT+(Z(1)*Z(2)/R)*C(I)*Y**(DBLE(I))
      END DO

!###### CALCULATING ANALYTIC DERIVATIVES IF DER=.TRUE.

      IF (DER) THEN 
        CALL DBASIS_CONTRACT(1,BSORDER,BS,R,DYDR)
        DVDRPART1=-POT/R
        DVDRPART2=0.00D+00
        DO I=1,POLORDER
          DVDRPART2=DVDRPART2+(Z(1)*Z(2)/R)*C(I)*DBLE(I)*Y**(DBLE(I-1))
        END DO
        DVDRPART2=DVDRPART2*DYDR
        DVDR=DVDRPART1+DVDRPART2
      ELSE
        DVDR=1000
      ENDIF

!######

      DEALLOCATE(BS)

      RETURN
      END

!####################################################################################

      SUBROUTINE CARTDERPES2BD(X,DVDX)
      IMPLICIT NONE 
      INTEGER :: I
      INTEGER, PARAMETER :: NATOM=2    
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: X
      DOUBLE PRECISION, DIMENSION(3*NATOM) :: DVDX,O
      DOUBLE PRECISION :: Y,R,DVDR
      DOUBLE PRECISION :: V

      CALL CART2INTER2BD(X,O,3*NATOM)

      R=O(1)

      CALL CHIPR_DIAT(R,V,.TRUE.,DVDR) 
      
      Y=DVDR/R

      DVDX(1)=-Y*(X(4)-X(1))
      DVDX(2)=-Y*(X(5)-X(2))
      DVDX(3)=-Y*(X(6)-X(3)) 
      DVDX(4)= Y*(X(4)-X(1))
      DVDX(5)= Y*(X(5)-X(2))
      DVDX(6)= Y*(X(6)-X(3))

      RETURN
      END

!####################################################################################

      SUBROUTINE CART2INTER2BD(X,R,NTA)
      IMPLICIT NONE
      INTEGER, PARAMETER :: NA=2,NT=3*NA
      INTEGER :: NTA,I,K,J
      DOUBLE PRECISION, DIMENSION (NT):: X,R
      K=0
      DO I=1,NTA,3
       DO J=I,NTA,3
        IF(I.NE.J)THEN
         K=K+1
         R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+&
     &   (X(I+2)-X(J+2))**2)
        END IF
       END DO
      END DO
      RETURN
      END

!######################################################################################################

      SUBROUTINE BASIS_CONTRACT(DEG,M,C,R,YVAL)
      IMPLICIT NONE
      INTEGER :: I,J,DEG
      INTEGER :: M
      DOUBLE PRECISION :: RREF0,ZETA
      DOUBLE PRECISION, DIMENSION(2*M+2) :: C
      DOUBLE PRECISION, DIMENSION(M) :: GAMA,VAL
      DOUBLE PRECISION :: R,YVAL

      DO I=1,M
        GAMA(I)=C(M+I)
      END DO

      RREF0=C(2*M+1)
      ZETA=C(2*M+2)

      DO I=1,M-1
        CALL PHISECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,R,VAL(I))
      END DO

      DO I=M,M
        CALL PHICSECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,6,R,VAL(I))
      END DO

      YVAL=0.00D+00

      DO J=1,M
        YVAL=YVAL+C(J)*VAL(J)
      END DO

      RETURN
      END

!######################################################################################################

      SUBROUTINE PHISECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,R,VAL)
      IMPLICIT NONE
      INTEGER :: DEG, IND, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, ORIG
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      VAL=(SECH(GAMA*RHO))**(DBLE(ETA))
      RETURN
      END

!######################################################################################################

      SUBROUTINE PHICSECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,LR,R,VAL)
      IMPLICIT NONE
      INTEGER :: DEG, IND, LR, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, BETA, FAC, ORIG
      BETA=1.00D+00/5.0D+00
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      FAC=(TANH(BETA*R)/R)**(DBLE(LR))
      VAL=FAC*(SECH(GAMA*RHO))**(DBLE(ETA))
      RETURN
      END

!######################################################################################################

      SUBROUTINE DBASIS_CONTRACT(DEG,M,C,R,DYDR)
      IMPLICIT NONE
      INTEGER :: I,J,DEG
      INTEGER :: M
      DOUBLE PRECISION :: RREF0,ZETA
      DOUBLE PRECISION, DIMENSION(2*M+2) :: C
      DOUBLE PRECISION, DIMENSION(M) :: GAMA,DPHIDR
      DOUBLE PRECISION :: R,DYDR

      DO I=1,M
        GAMA(I)=C(M+I)
      END DO

      RREF0=C(2*M+1)
      ZETA=C(2*M+2)

      DO I=1,M-1
        CALL DPHISECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,R,DPHIDR(I))
      END DO

      DO I=M,M
        CALL DPHICSECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,6,R,DPHIDR(I))
      END DO

      DYDR=0.00D+00

      DO J=1,M
        DYDR=DYDR+C(J)*DPHIDR(J)
      END DO

      RETURN
      END

!######################################################################################################

      SUBROUTINE DPHISECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,R,DPHIDR)
      IMPLICIT NONE
      INTEGER :: DEG, IND, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, ORIG
      DOUBLE PRECISION :: PART1,PART2
      DOUBLE PRECISION :: DPHIDR
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      PART1=(SECH(GAMA*RHO))**(DBLE(ETA))
      PART2=(TANH(GAMA*RHO))
      DPHIDR=-(DBLE(ETA))*GAMA*PART1*PART2
      RETURN
      END

!######################################################################################################

      SUBROUTINE DPHICSECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,LR,R,DPHIDR)
      IMPLICIT NONE
      INTEGER :: DEG, IND, LR, ETA
      DOUBLE PRECISION :: RREF0, ZETA, RREFIND 
      DOUBLE PRECISION :: GAMA, R, VAL, RHO
      DOUBLE PRECISION :: SECH, BETA, FAC, ORIG
      DOUBLE PRECISION :: FAC1PART1,FAC2PART1
      DOUBLE PRECISION :: FAC3PART1,FAC1PART2
      DOUBLE PRECISION :: FAC2PART2,FAC3PART2
      DOUBLE PRECISION :: PART1,PART2
      DOUBLE PRECISION :: DPHIDR
      BETA=1.00D+00/5.0D+00
      RREFIND=ORIG(IND,ZETA,RREF0)
      RHO=(R-RREFIND)
      FAC1PART1=(TANH(BETA*R)/R)**(DBLE(LR-1))
      FAC2PART1=(BETA*(SECH(BETA*R))**(2)/R)-(TANH(BETA*R)/R**2)
      FAC3PART1=(SECH(GAMA*RHO))**(DBLE(ETA))
      PART1=(DBLE(LR))*FAC1PART1*FAC2PART1*FAC3PART1
      FAC1PART2=(SECH(GAMA*RHO))**(DBLE(ETA))
      FAC2PART2=(TANH(GAMA*RHO))
      FAC3PART2=(TANH(BETA*R)/R)**(DBLE(LR))
      PART2=-(DBLE(ETA))*GAMA*FAC1PART2*FAC2PART2*FAC3PART2
      DPHIDR=PART1+PART2
      RETURN
      END

!######################################################################################################

      DOUBLE PRECISION FUNCTION ORIG(IND,ZETA,RREF0)
      IMPLICIT NONE
      INTEGER :: IND
      DOUBLE PRECISION :: ZETA,RREF0
      ORIG=ZETA*(RREF0)**(DBLE(IND)-1.0D+00)
      RETURN
      END

!######################################################################################################

      DOUBLE PRECISION FUNCTION SECH(X)
      IMPLICIT NONE
      DOUBLE PRECISION :: X
      SECH=1.00D+00/(COSH(X))      
      RETURN
      END
