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
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R_array, R_pred_array, V_array, V_pred_array
    DOUBLE PRECISION :: DVDR = 0, rmse_var
    LOGICAL :: DER = .false.
    INTEGER :: i, nlines

    !load R from file:
    nlines = 0
    open (1, file = 'abinitio_data.txt')
    do
        read (1,*, END=10)
        nlines = nlines + 1
    end do
    10 close (1)
    ALLOCATE(R_array(nlines))
    ALLOCATE(R_pred_array(nlines))
    ALLOCATE(V_array(nlines))
    ALLOCATE(V_pred_array(nlines))
    open (1, file = 'abinitio_data.txt')
    do i=1, nlines
        read(1,*) R_array(i), V_array(i)
    end do
    close (1)

    !load test data:
    !R_pred_array .....

    !evaluate test data points:
    do i=1, nlines
        call CHIPR_DIAT(R_array(i), V_pred_array(i), DER, DVDR)
    end do

    !compute RMSE:
    rmse_var = RMSE(V_array, V_pred_array)

    !save evaluated points to file:
    open(1, file = "data_evaluation.txt")
    do i=1, nlines
        write(1, *) R_array(i), V_pred_array(i)
    end do
    write(1, *) "#RMSE:"
    write(1, *) rmse_var
    close(1)


    contains
        FUNCTION RMSE(Y, Y_pred) RESULT(rmse_var)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:) :: Y, Y_pred
            DOUBLE PRECISION :: rmse_var

            rmse_var = SQRT(SUM((Y-Y_pred)**2))
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

      C(  1)=  0.11972558670992829998880324637866579D+01
      C(  2)= -0.15773156612563057898945828583237017D+00
      C(  3)= -0.54748741681060151265114654961507767D+01
      C(  4)= -0.49691983072670113941171621263492852D+01
      C(  5)= -0.93000115011802209075142400251934305D+00
      C(  6)= -0.83767962537636098119264715933240950D+01
      C(  7)= -0.25273808997492508865434501785784960D+02
      C(  8)= -0.20249977140037572098663076758384705D+02
      C(  9)=  0.12111046811709575976578889822121710D+02
      C( 10)=  0.17253113194749200687283519073389471D+02
      C( 11)= -0.18462712675332724643340043257921934D+02
      C( 12)= -0.33620096802461226559444185113534331D+02
      C( 13)= -0.13293426709624727166669799771625549D+02
      C( 14)= -0.50196819926494029573404986876994371D+00
      C( 15)=  0.87203348965986016150253590240026824D+00
      C( 16)=  0.51330396132145827969850415684049949D+00
      C( 17)= -0.10746494497134159024076538457848073D-01
      C( 18)= -0.11589697022802611381386839184415294D+00
      C( 19)= -0.96237350209106864440400386229157448D+00
      C( 20)=  0.53720000781503504660996384245663648D-01
      C( 21)= -0.23414650478207921224793608416803181D+01
      C( 22)= -0.89617493812937054631362343570799567D+00
      C( 23)=  0.46293840200470831591772480351210106D+00
      C( 24)= -0.14413069650176122404872103288653307D+01
      C( 25)=  0.28496702013404839881616226193727925D+01
      C( 26)= -0.46900526430132100097480929434823338D+00
      C( 27)= -0.58098802412764383173993110176525079D-01
      C( 28)= -0.12549190709485117167787393555045128D+02
      C( 29)=  0.96023584485018043110926555527839810D+00
      C( 30)=  0.26869992356681478007374153094133362D+01
      C( 31)=  0.63276051924359269507114333919162164D-01

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

