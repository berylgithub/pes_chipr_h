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
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R_array, POT_array
    DOUBLE PRECISION :: DVDR = 0
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
    !print*, nlines
    ALLOCATE(R_array(nlines))
    ALLOCATE(POT_array(nlines))
    !print*,shape(R_array)
    open (1, file = 'abinitio_data.txt')
    do i=1, nlines
        read(1,*) R_array(i), POT_array(i)
    end do
    close (1)

    !evaluate each point:
    do i=1, nlines
        call CHIPR_DIAT(R_array(i), POT_array(i), DER, DVDR)
    end do

    !save evaluated points to file:
    open(1, file = "h2_evaluation.txt")
    do i=1, nlines
        write(1, *) R_array(i), POT_array(i)
    end do
    close(1)

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

      C(  1)= -0.20244854223544681559587843366898596D+01
      C(  2)=  0.17435955846974378058433785554370843D+01
      C(  3)= -0.58116564368505563287214954470982775D+00
      C(  4)=  0.10722732579382783513199939307014574D+00
      C(  5)= -0.12544662979550811818252320506417163D-01
      C(  6)=  0.99309880037843148437026741248701001D-03
      C(  7)= -0.54955921722621863575574535643042395D-04
      C(  8)=  0.21527739708922375881042926559771900D-05
      C(  9)= -0.59475815914670122823631955326015297D-07
      C( 10)=  0.11337659640011901848842856322958222D-08
      C( 11)= -0.14193284000511143262418373344774410D-10
      C( 12)=  0.10501442237421527764275964460125086D-12
      C( 13)= -0.34800356676644887734694037831150740D-15
      C( 14)=  0.14890108571310017512856482824190607D-01
      C( 15)= -0.10892630024741943237098773522575357D-01
      C( 16)=  0.41903776810371553551703982520848513D+02
      C( 17)=  0.14966087118467594554993560507227812D-01
      C( 18)= -0.10827121326005622117816251659405680D-01
      C( 19)=  0.41921232150178731501455331454053521D+02
      C( 20)=  0.14764727948189199091544310249446426D-01
      C( 21)= -0.10945451229238045351421781958833890D-01
      C( 22)=  0.77344459896782902230683021116419695D+00
      C( 23)=  0.80282242682443238912526339845499024D+00
      C( 24)=  0.13562749309760493421350702192285098D+01
      C( 25)=  0.78285514468193673209839289484079927D+00
      C( 26)=  0.80595528189233489602116833339096047D+00
      C( 27)=  0.13562349503180692877890578529331833D+01
      C( 28)=  0.80842890481786067180536292653414421D+00
      C( 29)=  0.32360967508267611947303521446883678D+03
      C( 30)=  0.14402965687991096110920352657558396D+01
      C( 31)=  0.10469614339079305054269752872642130D+01

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
