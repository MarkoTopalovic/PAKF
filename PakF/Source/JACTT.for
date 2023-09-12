C==========================================================================
       SUBROUTINE JACTT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21,
     &DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,
     &FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVZK,ZVXOM,ZVYOM,ZVZOM,
     &ZVXV1,ZVXV2,ZVXV3,ZVYV1,ZVYV2,ZVYV3,ZVZV1,ZVZV2,ZVZV3,CORD,ID,
     &PRITISAK,PENALT,PRIT,IDPRIT,GUSM)
C
CE  Subroutine JACTT is used for integration 3D finite elements IN TURBULENCE MODE
C
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON /IZID/ IZID
       DIMENSION H(21),P(3,21),NEL(NDIM+1,*)
       DIMENSION XJ(3,3),LE(3),ME(3),CK(21,3),PJ(3,21),HP(8)
       DIMENSION TT21(56),V1(21),V2(21),V3(21),PP(8),XJJ(3,3),TEMP(21)
       DIMENSION VTANG(3,21),AN(3),VSHX(3),VSHY(3),VSHZ(3),SHEAR(*)
       DIMENSION IPERM(8),NOD9(13)
       DIMENSION ID(7,*),YPLUS(3),VPLUS(3),YP(3,21),CORD(3,*)
C       DIMENSION PRITISAK(1,*)
       DIMENSION PRIT(IDPRIT,*)
C--------------------------------------------------------------------
C     K-OMEGA TURB. MODEL
C--------------------------------------------------------------------
      DIMENSION AAKT(21)
      DIMENSION OMT(21)
      AKAPA=0.42
      PET=3.0
C--------------------------------------------------------------------
C       COMMON /VISKOZ/ AMI,INDAMI
C       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C       COMMON /SRPSKI/ ISRPS
C       COMMON /ULAZNI/ IULAZ,IIZLAZ
C       COMMON /AJVVV/ HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW
C       COMMON /UPWIND/ IUPWIN

C      H(21)           - Interolacione funkcije
C      P(3,21)         - Izvodi interpolacionih funkcija
C      NEL(NDIM,*)     - Koordinate cvorova elementa
C      XJ(3,3)         - Jakobijan
C      XJJ(3,3)        - Inverzni Jakobijan
C      LE(3)           - 
C      ME(3)           - 
C      CK(21,3)        - 
C      PJ(3,21)        - Matrica izvoda interpolacionih funkcija 
C                        po globalnim koordinatama
C      HP(8)           - 
C      TT21(*)         - 
C      V1(21)          - Brzina u X pravcu
C      V2(21)          - Brzina u Y pravcu
C      V3(21)          - Brzina u Z pravcu
C      PP(8)           - Pritisak
C      TEMP(21)        - Temperatura
C      VTANG(3,21)     - 
C      AN(3)           - 
C      VSHX(3)         - 
C      VSHY(3)         -  
C      VSHZ(3)         - 
C      SHEAR(*)        - 
C      IPERM(8)        - Oznacavanje cvorova od 1 do 8
C      NOD9(13)        - Oznacavanje cvorova od 9 do 21

       DATA IPERM/2,3,4,1,6,7,8,5/
       DATA NOD9/9,10,11,12,13,14,15,16,17,18,19,20,21/

       IELX=NDIM
       INDUP=0
       NND9=13
       CALL CLEAR(SHEAR,3)
C (zagrade kod interpolacionih funckija: P=+, M=-)
      RP=1.0+R
      SP=1.0+S
      TP=1.0+T
      RM=1.0-R
      TM=1.0-T
      SM=1.0-S
      RR=1.0-R*R
      SS=1.0-S*S
      TT=1.0-T*T
      DO 82 I=1,21
      H(I)=0.
      DO 82 J=1,3
   82 P(J,I)=0.
   
      DO I=1,NDIM
C Brzina I-tog cvora u X-pravcu
      V1(I)=TT21(I)
C Brzina I-tog cvora u Y-pravcu
      V2(I)=TT21(I+NDIM)
C Brzina I-tog cvora u Z-pravcu
      V3(I)=TT21(I+2*NDIM)
C Temperatura I-tog cvora
      TEMP(I)=TT21(I+3*NDIM+8)
      IF (I.LE.8) THEN
C Pritisak I-tog cvora za  IPERM cvorove
        PP(I)=TT21(I+3*NDIM)
      ENDIF
C--------------------------------------------------------------------
C DODATO ZA KT I OMT
      AAKT(I)=TT21(I+5*NDIM)
      OMT(I)=TT21(I+6*NDIM)
C--------------------------------------------------------------------
      ENDDO
      
C---------------------------------------------------
C     ODREDJIVANJE CVORA U ELEMENTU
C---------------------------------------------------
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)
      N5=NEL(5,NBREL)
      N6=NEL(6,NBREL)
      N7=NEL(7,NBREL)
      N8=NEL(8,NBREL)
      
C
C     INTERPOLACIJSKE FUNKCIJE I NJIHOVI IZVODI
CE    INTERPOLATION FUNCTIONS AND DERIVATIVES
C
C
CS    PRVIH 8 CVOROVA
CE    FIRST 4 NODES
C
      H(1)=0.125*RP*SP*TP
      H(2)=0.125*RM*SP*TP
      H(3)=0.125*RM*SM*TP
      H(4)=0.125*RP*SM*TP
      H(5)=0.125*RP*SP*TM
      H(6)=0.125*RM*SP*TM
      H(7)=0.125*RM*SM*TM
      H(8)=0.125*RP*SM*TM
      
 100  DO JJ=1,8
       HP(JJ)=H(JJ)
      ENDDO
C hi,1 za prvih 8 cvorova
      P(1,1)=0.125*SP*TP
      P(1,2)=-P(1,1)
      P(1,3)=-0.125*SM*TP
      P(1,4)=-P(1,3)
      P(1,5)=0.125*SP*TM
      P(1,6)=-P(1,5)
      P(1,7)=-0.125*SM*TM
      P(1,8)=-P(1,7)
C hi,2 za prvih 8 cvorova
      P(2,1)=0.125*RP*TP
      P(2,2)=0.125*RM*TP
      P(2,3)=-P(2,2)
      P(2,4)=-P(2,1)
      P(2,5)=0.125*RP*TM
      P(2,6)=0.125*RM*TM
      P(2,7)=-P(2,6)
      P(2,8)=-P(2,5)
C hi,3 za prvih 8 cvorova 
      P(3,1)=0.125*RP*SP
      P(3,2)=0.125*RM*SP
      P(3,3)=0.125*RM*SM
      P(3,4)=0.125*RP*SM
      P(3,5)=-P(3,1)
      P(3,6)=-P(3,2)
      P(3,7)=-P(3,3)
      P(3,8)=-P(3,4)
C
      IF (IELX.EQ.8) GO TO 50
C
CS    STEPENI SLOBODE ZA CVOROVE PREKO 8
CE    DEGREES OF FREEDOM FOR NODES OVER 8
C
      I=0
    2 I=I+1
      IF (I.GT.NND9) GO TO 30
      NN=NOD9(I)-8
      GO TO (9,10,11,12,13,14,15,16,17,18,19,20,21),NN
C
    9 H(9)=0.25*RR*SP*TP
      P(1,9)=-0.50*R*SP*TP
      P(2,9)=0.25*RR*TP
      P(3,9)=0.25*RR*SP
      GO TO 2
   10 H(10)=0.25*RM*SS*TP
      P(1,10)=-0.25*SS*TP
      P(2,10)=-0.50*RM*S*TP
      P(3,10)=0.25*RM*SS
      GO TO 2
   11 H(11)=0.25*RR*SM*TP
      P(1,11)=-0.50*R*SM*TP
      P(2,11)=-0.25*RR*TP
      P(3,11)=0.25*RR*SM
      GO TO 2
   12 H(12)=0.25*RP*SS*TP
      P(1,12)=0.25*SS*TP
      P(2,12)=-0.50*RP*S*TP
      P(3,12)=0.25*RP*SS
      GO TO 2
   13 H(13)=0.25*RR*SP*TM
      P(1,13)=-0.50*R*SP*TM
      P(2,13)=0.25*RR*TM
      P(3,13)=-0.25*RR*SP
      GO TO 2
   14 H(14)=0.25*RM*SS*TM
      P(1,14)=-0.25*SS*TM
      P(2,14)=-0.50*RM*S*TM
      P(3,14)=-0.25*RM*SS
      GO TO 2
   15 H(15)=0.25*RR*SM*TM
      P(1,15)=-0.50*R*SM*TM
      P(2,15)=-0.25*RR*TM
      P(3,15)=-0.25*RR*SM
      GO TO 2
   16 H(16)=0.25*RP*SS*TM
      P(1,16)=0.25*SS*TM
      P(2,16)=-0.50*RP*S*TM
      P(3,16)=-0.25*RP*SS
      GO TO 2
   17 H(17)=0.25*RP*SP*TT
      P(1,17)=0.25*SP*TT
      P(2,17)=0.25*RP*TT
      P(3,17)=-0.50*RP*SP*T
      GO TO 2
   18 H(18)=0.25*RM*SP*TT
      P(1,18)=-0.25*SP*TT
      P(2,18)=0.25*RM*TT
      P(3,18)=-0.50*RM*SP*T
      GO TO 2
   19 H(19)=0.25*RM*SM*TT
      P(1,19)=-0.25*SM*TT
      P(2,19)=-0.25*RM*TT
      P(3,19)=-0.50*RM*SM*T
      GO TO 2
   20 H(20)=0.25*RP*SM*TT
      P(1,20)=0.25*SM*TT
      P(2,20)=-0.25*RP*TT
      P(3,20)=-0.50*RP*SM*T
      GO TO 2
   21 H(21)=RR*SS*TT
      P(1,21)=-2.0*R*SS*TT
      P(2,21)=-2.0*S*RR*TT
      P(3,21)=-2.0*T*RR*SS
      GO TO 2
C
CS    KOREKCIJE PRVIH 20 FUNKCIJA AKO JE UPOTREBLJEN CVOR 21
CE    CORECTION OF FIRST 20 FUNCTIONS IF NODE 21 EXISTS
C
   30 IN=NOD9(NND9)
      IF(IN.NE.21) GO TO 40
      DO 36 I=1,8
      H(I)=H(I)-0.125*H(21)
      DO 36 J=1,3
   36 P(J,I)=P(J,I)-0.125*P(J,21)
      IF(NND9.EQ.1) GO TO 51
      DO 37 I=1,NND9-1
      IN=NOD9(I)
      H(IN)=H(IN)-0.25*H(21)
      DO 37 J=1,3
   37 P(J,IN)=P(J,IN)-0.25*P(J,21)
C
CS    KOREKCIJE PRVIH 8 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 8
CE    CORECTION OF FIRST 8 FUNCTIONS IF NODES OVER 8 EXISTS
C
   40 IH=0
   41 IH=IH+1
      IF(IH.GT.NND9) GO TO 50
      IN=NOD9(IH)
      IF(IN.GT.16) GO TO 46
      I1=IN-8
      I2=IPERM(I1)
C
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 45 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   45 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   46 IF(IN.EQ.21) GO TO 51
      I1=IN-16
      I2=I1+4
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 47 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   47 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   51 H(NND9+8)=H(21)
      DO 39 J=1,3
   39 P(J,NND9+8)=P(J,21)

      
   50 HH=0.D0
      DO I=1,NDIM
       HH=HH+H(I)
      ENDDO
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C
      DO I=1,3
      DO J=1,3
       XJ(I,J)=0.D0
       XJJ(I,J)=0.D0
        DO  KK=1,NDIM
C Matrica Jakobijana
         XJ(I,J)=XJ(I,J)+P(I,KK)*CK(KK,J)
C Inicijalizacija Inverzne Matrice Jakobijana
         XJJ(I,J)=XJJ(I,J)+P(I,KK)*CK(KK,J)
        ENDDO    
C      WRITE(IIZLAZ,*)'XJ=',XJ(I,J)
      ENDDO    
      ENDDO    

      HV1=DOT(H,V1,NDIM)
      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)

C Odredjivanje Inverzne Matrice Jakobijana
      CALL MINV(XJJ,3,DET1,LE,ME)

C      IF (DET1.LT.1.D-15) THEN
C       WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
C       WRITE(IIZLAZ,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
C       WRITE(IIZLAZ,*)'DETERMINANTE= ',DET1
C       WRITE(IIZLAZ,*)'NODES COORDINATES'
C       DO I=1,NDIM
C        WRITE(IIZLAZ,1000) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
C       ENDDO
C       STOP
C      ELSE
C      ENDIF
        
C 1000 FORMAT(I5,3(D13.5))

C     Matrica izvoda interpolacionih funkcija 
C     po globalnim koordinatama
      DO 85 I=1,3
      DO 85 JJ=1,NDIM
      PJ(I,JJ)=0.D0
      DO 85 K=1,3
      PJ(I,JJ)=PJ(I,JJ) + XJJ(I,K)*P(K,JJ)
   85 CONTINUE


      IF (IUPWIN.EQ.1.AND.INDUP.EQ.0) THEN
       CALL INTER1(CK,V1,V2,V3,H,PJ)
       INDUP=1
       GOTO 100
      ENDIF
      

      HXU=0.D0
      HYU=0.D0
      HZU=0.D0
      HXV=0.D0
      HYV=0.D0
      HZV=0.D0
      HXW=0.D0
      HYW=0.D0
      HZW=0.D0
      ZVXT=0.D0
      ZVYT=0.D0
      ZVZT=0.D0
      AKTX=0.D0
      AKTY=0.D0
      AKTZ=0.D0
      OMTX=0.D0
      OMTY=0.D0
      OMTZ=0.D0
      
      DO L =1,NDIM
        HXU=HXU+PJ(1,L)*V1(L)
        HYU=HYU+PJ(2,L)*V1(L)
        HZU=HZU+PJ(3,L)*V1(L)
        HXV=HXV+PJ(1,L)*V2(L)
        HYV=HYV+PJ(2,L)*V2(L)
        HZV=HZV+PJ(3,L)*V2(L)
        HXW=HXW+PJ(1,L)*V3(L)
        HYW=HYW+PJ(2,L)*V3(L)
        HZW=HZW+PJ(3,L)*V3(L)
        ZVXT=ZVXT+PJ(1,L)*TEMP(L)
        ZVYT=ZVYT+PJ(2,L)*TEMP(L)
        ZVZT=ZVZT+PJ(3,L)*TEMP(L)
C--------------------------------------------------------------------
C DODATO ZA KT I OMT
        AKTX=AKTX+PJ(1,L)*AAKT(L)
        AKTY=AKTY+PJ(2,L)*AAKT(L)
        AKTZ=AKTZ+PJ(3,L)*AAKT(L)
        OMTX=OMTX+PJ(1,L)*OMT(L)
        OMTY=OMTY+PJ(2,L)*OMT(L)
        OMTZ=OMTZ+PJ(3,L)*OMT(L)
        
C PROIZVOD IZVODA H I BRZINA
      ZVXV1=DOT(PJ(1,L),V1,NDIM)
      ZVXV2=DOT(PJ(1,L),V2,NDIM)
      ZVXV3=DOT(PJ(1,L),V3,NDIM)
      
      ZVYV1=DOT(PJ(2,L),V1,NDIM)
      ZVYV2=DOT(PJ(2,L),V2,NDIM)
      ZVYV3=DOT(PJ(2,L),V3,NDIM)
      
      ZVZV1=DOT(PJ(3,L),V1,NDIM)
      ZVZV2=DOT(PJ(3,L),V2,NDIM)
      ZVZV3=DOT(PJ(3,L),V3,NDIM)

C PROIZVOD IZVODA H I KT I OMT

      ZVXK=DOT(PJ(1,L),AAKT,NDIM)
      ZVYK=DOT(PJ(2,L),AAKT,NDIM)
      ZVZK=DOT(PJ(3,L),AAKT,NDIM)

      ZVXOM=DOT(PJ(1,L),OMT,NDIM)
      ZVYOM=DOT(PJ(2,L),OMT,NDIM)
      ZVZOM=DOT(PJ(3,L),OMT,NDIM)
      
C PROIZVOD KT, OMT I H

      HKT=DOT(H,AAKT,NDIM)
      HOM=DOT(H,OMT,NDIM)
      
      ENDDO
      
      IF(KFIX.GT.0) GO TO 70
      RETURN

C
CS     DETERMINATA POVRSINSKOG JAKOBIJANA
CE     SURFACE JACOBIAN DETERMINANT
C
  70   GO TO (71,72,73),KFIX
CS     KONSTANTNO KSI
CE     CONSTANT KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      ANX=R*(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))
      ANY=R*(XJ(3,1)*XJ(2,3)-XJ(3,3)*XJ(2,1))
      ANZ=R*(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ETA
CE     CONSTANT ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      ANX=-S*(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))
      ANY=S*(XJ(1,1)*XJ(3,3)-XJ(1,3)*XJ(3,1))
      ANZ=-S*(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ZETA
CE     CONSTANT ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
      ANX=T*(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))
      ANY=-T*(XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1))
      ANZ=T*(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))
   74 DET=DSQRT(DET)
      ANX=ANX/DET
      ANY=ANY/DET
      ANZ=ANZ/DET
      AN(1)=ANX
      AN(2)=ANY
      AN(3)=ANZ
C      WRITE(IIZLAZ,*)'ANX=',ANX
C      WRITE(IIZLAZ,*)'ANY=',ANY
C      WRITE(IIZLAZ,*)'ANZ=',ANZ
      CALL SHEARS(AN,TT21,VTANG,NDIM)
      
      V1X=0.D0
      V1Y=0.D0
      V1Z=0.D0
      V2X=0.D0
      V2Y=0.D0
      V2Z=0.D0
      V3X=0.D0
      V3Y=0.D0
      V3Z=0.D0
      VSHX(1)=0.D0
      VSHX(2)=0.D0
      VSHX(3)=0.D0
      VSHY(1)=0.D0
      VSHY(2)=0.D0
      VSHY(3)=0.D0
      VSHZ(1)=0.D0
      VSHZ(2)=0.D0
      VSHZ(3)=0.D0
      DO I=1,NDIM
       V1X=V1X+PJ(1,I)*V1(I)
       V1Y=V1Y+PJ(2,I)*V1(I)
       V1Z=V1Z+PJ(3,I)*V1(I)
       V2X=V2X+PJ(1,I)*V2(I)
       V2Y=V2Y+PJ(2,I)*V2(I)
       V2Z=V2Z+PJ(3,I)*V2(I)
       V3X=V3X+PJ(1,I)*V3(I)
       V3Y=V3Y+PJ(2,I)*V3(I)
       V3Z=V3Z+PJ(3,I)*V3(I)
       VSHX(1)=VSHX(1)+PJ(1,I)*VTANG(1,I)
       VSHX(2)=VSHX(2)+PJ(2,I)*VTANG(1,I)
       VSHX(3)=VSHX(3)+PJ(3,I)*VTANG(1,I)
       VSHY(1)=VSHY(1)+PJ(1,I)*VTANG(2,I)
       VSHY(2)=VSHY(2)+PJ(2,I)*VTANG(2,I)
       VSHY(3)=VSHY(3)+PJ(3,I)*VTANG(2,I)
       VSHZ(1)=VSHZ(1)+PJ(1,I)*VTANG(3,I)
       VSHZ(2)=VSHZ(2)+PJ(2,I)*VTANG(3,I)
       VSHZ(3)=VSHZ(3)+PJ(3,I)*VTANG(3,I)
      ENDDO
      SHEAR(1)=DOT(AN,VSHX,3)
      SHEAR(2)=DOT(AN,VSHY,3)
      SHEAR(3)=DOT(AN,VSHZ,3)
      PRIT1=DOT(HP,PP,8)
C      SF1=-PRIT*ANX
C      SF2=-PRIT*ANY
C      SF3=-PRIT*ANZ
      SF1=-PRIT1*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
      SF2=-PRIT1*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
      SF3=-PRIT1*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
      
C DODATO ZA PRITISAK

      DO I=1,NDIM
C        I1=I
C        I2=I+NDIM
C        I3=I+2*NDIM
      IF (PENALT.GT.1.D0) THEN
      
C      PRITISAK=-PENALT*(PJ(1,I)*TT21(I)+PJ(2,I)*TT21(I+NDIM)+
C     &PJ(3,I)*TT21(I+2*NDIM))
     
      ENDIF
      
C      WRITE(IIZLAZ,*)'PRITISAK= ',PRITISAK, NBREL

      ENDDO      

C------------------------------------------------------------
C DODATO ZA FSK1 I FSOMEGA1
      
      FSK1=(AKTX*ANX+AKTY*ANY+AKTZ*ANZ)
      FSOMEGA1=(OMTX*ANX+OMTY*ANY+OMTZ*ANZ)
      
C----------------------------------------------------- 
C----------------------------------------------------- 
C ZIDNE FUNKCIJE (WALL FUNCTIONS)
C----------------------------------------------------- 
C----------------------------------------------------- 
      IF (IZID.EQ.1) THEN
      
      DO 101 I=1,NDIM   
C AKO SU CVOROVI 5,6,7 I 8 = 0 
      IF (ID(1,N5).EQ.0.AND.ID(2,N5).EQ.0.AND.ID(3,N5).EQ.0.
     &AND.ID(1,N6).EQ.0.AND.ID(2,N6).EQ.0.AND.ID(3,N6).EQ.0.
     &AND.ID(1,N7).EQ.0.AND.ID(2,N7).EQ.0.AND.ID(3,N7).EQ.0.
     &AND.ID(1,N8).EQ.0.AND.ID(2,N8).EQ.0.AND.ID(3,N8).EQ.0.) THEN
     
      YP(1,I)=H(I)*CORD(1,I)
      YP(2,I)=H(I)*CORD(2,I)
      YP(3,I)=H(I)*CORD(3,I)
      
      YPLUS(1)=GUSM*YP(1,I)*VTANG(1,I)/AMI
      YPLUS(2)=GUSM*YP(2,I)*VTANG(2,I)/AMI
      YPLUS(3)=GUSM*YP(3,I)*VTANG(3,I)/AMI
     
      IF (YPLUS(1).EQ.0.OR.YPLUS(2).EQ.0.OR.YPLUS(3).EQ.0) THEN
      VPLUS(1)=(YPLUS(1))/AKAPA+PET
      VPLUS(2)=(YPLUS(2))/AKAPA+PET
      VPLUS(3)=(YPLUS(3))/AKAPA+PET
      ELSE
      VPLUS(1)=(LOG(YPLUS(1)))/AKAPA+PET
      VPLUS(2)=(LOG(YPLUS(2)))/AKAPA+PET
      VPLUS(3)=(LOG(YPLUS(3)))/AKAPA+PET
      ENDIF
      
      TT21(NDIM-3)=V1(NDIM-3)
      TT21(NDIM-2)=V1(NDIM-2)
      TT21(NDIM-1)=V1(NDIM-1)
      TT21(NDIM)=V1(NDIM)
      TT21(NDIM-7)=VPLUS(1)
      TT21(NDIM-6)=VPLUS(1)
      TT21(NDIM-5)=VPLUS(1)
      TT21(NDIM-4)=VPLUS(1)
      
      TT21(5+NDIM)=V2(NDIM-3)
      TT21(6+NDIM)=V2(NDIM-2)
      TT21(7+NDIM)=V2(NDIM-1)
      TT21(8+NDIM)=V2(NDIM)
      TT21(1+NDIM)=VPLUS(2)
      TT21(2+NDIM)=VPLUS(2)
      TT21(3+NDIM)=VPLUS(2)
      TT21(4+NDIM)=VPLUS(2)
      
c      TT21(3*NDIM-3)=V3(NDIM-3)
c      TT21(3*NDIM-2)=V3(NDIM-2)
c      TT21(3*NDIM-1)=V3(NDIM-1)
c      TT21(3*NDIM)=V3(NDIM)
c      TT21(3*NDIM-7)=VPLUS(3)
c      TT21(3*NDIM-6)=VPLUS(3)
c      TT21(3*NDIM-5)=VPLUS(3)
c      TT21(3*NDIM-4)=VPLUS(3)
      
      PRIT(IDPRIT,NBREL)=0.0
      
C AKO SU CVOROVI 2,6,7 I 3 = 0 

      ELSEIF (ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.
     &AND.ID(1,N6).EQ.0.AND.ID(2,N6).EQ.0.AND.ID(3,N6).EQ.0.
     &AND.ID(1,N7).EQ.0.AND.ID(2,N7).EQ.0.AND.ID(3,N7).EQ.0.
     &AND.ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.) THEN
     
      YP(1,I)=H(I)*CORD(1,I)
      YP(2,I)=H(I)*CORD(2,I)
      YP(3,I)=H(I)*CORD(3,I)
      
      YPLUS(1)=GUSM*YP(1,I)*VTANG(1,I)/AMI
      YPLUS(2)=GUSM*YP(2,I)*VTANG(2,I)/AMI
      YPLUS(3)=GUSM*YP(3,I)*VTANG(3,I)/AMI
      
      IF (YPLUS(1).EQ.0.OR.YPLUS(2).EQ.0.OR.YPLUS(3).EQ.0) THEN
      VPLUS(1)=(YPLUS(1))/AKAPA+PET
      VPLUS(2)=(YPLUS(2))/AKAPA+PET
      VPLUS(3)=(YPLUS(3))/AKAPA+PET
      ELSE
      VPLUS(1)=(LOG(YPLUS(1)))/AKAPA+PET
      VPLUS(2)=(LOG(YPLUS(2)))/AKAPA+PET
      VPLUS(3)=(LOG(YPLUS(3)))/AKAPA+PET
      ENDIF
      
      TT21(NDIM-6)=V1(NDIM-6)
      TT21(NDIM-2)=V1(NDIM-2)
      TT21(NDIM-1)=V1(NDIM-1)
      TT21(NDIM-5)=V1(NDIM-5)
      TT21(NDIM-7)=VPLUS(1)
      TT21(NDIM-3)=VPLUS(1)
      TT21(NDIM)=VPLUS(1)
      TT21(NDIM-4)=VPLUS(1)
      
      TT21(2+NDIM)=V2(NDIM-6)
      TT21(6+NDIM)=V2(NDIM-2)
      TT21(7+NDIM)=V2(NDIM-1)
      TT21(3+NDIM)=V2(NDIM-5)
      TT21(1+NDIM)=VPLUS(2)
      TT21(5+NDIM)=VPLUS(2)
      TT21(8+NDIM)=VPLUS(2)
      TT21(4+NDIM)=VPLUS(2)
      
c      TT21(3*NDIM-6)=V3(NDIM-6)
c      TT21(3*NDIM-2)=V3(NDIM-2)
c      TT21(3*NDIM-1)=V3(NDIM-1)
c      TT21(3*NDIM-5)=V3(NDIM-5)
c      TT21(3*NDIM-7)=VPLUS(3)
c      TT21(3*NDIM-3)=VPLUS(3)
c      TT21(3*NDIM)=VPLUS(3)
c      TT21(3*NDIM-4)=VPLUS(3)
      
      PRIT(IDPRIT,NBREL)=0.0
      
C AKO SU CVOROVI 2,6,5 I 1 = 0 

      ELSEIF (ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.
     &AND.ID(1,N6).EQ.0.AND.ID(2,N6).EQ.0.AND.ID(3,N6).EQ.0.
     &AND.ID(1,N5).EQ.0.AND.ID(2,N5).EQ.0.AND.ID(3,N5).EQ.0.
     &AND.ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.) THEN
     
      YP(1,I)=H(I)*CORD(1,I)
      YP(2,I)=H(I)*CORD(2,I)
      YP(3,I)=H(I)*CORD(3,I)
      
      YPLUS(1)=GUSM*YP(1,I)*VTANG(1,I)/AMI
      YPLUS(2)=GUSM*YP(2,I)*VTANG(2,I)/AMI
      YPLUS(3)=GUSM*YP(3,I)*VTANG(3,I)/AMI
      
      IF (YPLUS(1).EQ.0.OR.YPLUS(2).EQ.0.OR.YPLUS(3).EQ.0) THEN
      VPLUS(1)=(YPLUS(1))/AKAPA+PET
      VPLUS(2)=(YPLUS(2))/AKAPA+PET
      VPLUS(3)=(YPLUS(3))/AKAPA+PET
      ELSE
      VPLUS(1)=(LOG(YPLUS(1)))/AKAPA+PET
      VPLUS(2)=(LOG(YPLUS(2)))/AKAPA+PET
      VPLUS(3)=(LOG(YPLUS(3)))/AKAPA+PET
      ENDIF
            
      TT21(NDIM-6)=V1(NDIM-6)
      TT21(NDIM-2)=V1(NDIM-2)
      TT21(NDIM-3)=V1(NDIM-3)
      TT21(NDIM-7)=V1(NDIM-7)
      TT21(NDIM-5)=VPLUS(1)
      TT21(NDIM-1)=VPLUS(1)
      TT21(NDIM)=VPLUS(1)
      TT21(NDIM-4)=VPLUS(1)
      
      TT21(2+NDIM)=V2(NDIM-6)
      TT21(6+NDIM)=V2(NDIM-2)
      TT21(5+NDIM)=V2(NDIM-3)
      TT21(1+NDIM)=V2(NDIM-7)
      TT21(3+NDIM)=VPLUS(2)
      TT21(7+NDIM)=VPLUS(2)
      TT21(8+NDIM)=VPLUS(2)
      TT21(4+NDIM)=VPLUS(2)
      
c      TT21(3*NDIM-6)=V3(NDIM-6)
c      TT21(3*NDIM-2)=V3(NDIM-2)
c      TT21(3*NDIM-3)=V3(NDIM-3)
c      TT21(3*NDIM-7)=V3(NDIM-7)
c      TT21(3*NDIM-5)=VPLUS(3)
c      TT21(3*NDIM-1)=VPLUS(3)
c      TT21(3*NDIM)=VPLUS(3)
c      TT21(3*NDIM-4)=VPLUS(3)

      PRIT(IDPRIT,NBREL)=0.0
      
C AKO SU CVOROVI 1,2,3 I 4 = 0 

      ELSEIF (ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.
     &AND.ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.
     &AND.ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.
     &AND.ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0.) THEN
     
      YP(1,I)=H(I)*CORD(1,I)
      YP(2,I)=H(I)*CORD(2,I)
      YP(3,I)=H(I)*CORD(3,I)
      
      YPLUS(1)=GUSM*YP(1,I)*VTANG(1,I)/AMI
      YPLUS(2)=GUSM*YP(2,I)*VTANG(2,I)/AMI
      YPLUS(3)=GUSM*YP(3,I)*VTANG(3,I)/AMI
      
      IF (YPLUS(1).EQ.0.OR.YPLUS(2).EQ.0.OR.YPLUS(3).EQ.0) THEN
      VPLUS(1)=(YPLUS(1))/AKAPA+PET
      VPLUS(2)=(YPLUS(2))/AKAPA+PET
      VPLUS(3)=(YPLUS(3))/AKAPA+PET
      ELSE
      VPLUS(1)=(LOG(YPLUS(1)))/AKAPA+PET
      VPLUS(2)=(LOG(YPLUS(2)))/AKAPA+PET
      VPLUS(3)=(LOG(YPLUS(3)))/AKAPA+PET
      ENDIF
            
      TT21(NDIM-7)=V1(NDIM-7)
      TT21(NDIM-6)=V1(NDIM-6)
      TT21(NDIM-5)=V1(NDIM-5)
      TT21(NDIM-4)=V1(NDIM-4)
      TT21(NDIM-3)=VPLUS(1)
      TT21(NDIM-2)=VPLUS(1)
      TT21(NDIM-1)=VPLUS(1)
      TT21(NDIM)=VPLUS(1)
      
      TT21(1+NDIM)=V2(NDIM-7)
      TT21(2+NDIM)=V2(NDIM-6)
      TT21(3+NDIM)=V2(NDIM-5)
      TT21(4+NDIM)=V2(NDIM-4)
      TT21(5+NDIM)=VPLUS(2)
      TT21(6+NDIM)=VPLUS(2)
      TT21(7+NDIM)=VPLUS(2)
      TT21(8+NDIM)=VPLUS(2)
      
c      TT21(1+2*NDIM)=V3(NDIM-7)
c      TT21(2+2*NDIM)=V3(NDIM-6)
c      TT21(3+2*NDIM)=V3(NDIM-5)
c      TT21(4+2*NDIM)=V3(NDIM-4)
c      TT21(5+2*NDIM)=VPLUS(3)
c      TT21(6+2*NDIM)=VPLUS(3)
c      TT21(7+2*NDIM)=VPLUS(3)
c      TT21(8+2*NDIM)=VPLUS(3)
            
      PRIT(IDPRIT,NBREL)=0.0
      
C AKO SU CVOROVI 1,5,8 I 4 = 0 
      
      ELSEIF (ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.
     &AND.ID(1,N5).EQ.0.AND.ID(2,N5).EQ.0.AND.ID(3,N5).EQ.0.
     &AND.ID(1,N8).EQ.0.AND.ID(2,N8).EQ.0.AND.ID(3,N8).EQ.0.
     &AND.ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0.) THEN
     
      YP(1,I)=H(I)*CORD(1,I)
      YP(2,I)=H(I)*CORD(2,I)
      YP(3,I)=H(I)*CORD(3,I)
      
      YPLUS(1)=GUSM*YP(1,I)*VTANG(1,I)/AMI
      YPLUS(2)=GUSM*YP(2,I)*VTANG(2,I)/AMI
      YPLUS(3)=GUSM*YP(3,I)*VTANG(3,I)/AMI
      
      IF (YPLUS(1).EQ.0.OR.YPLUS(2).EQ.0.OR.YPLUS(3).EQ.0) THEN
      VPLUS(1)=(YPLUS(1))/AKAPA+PET
      VPLUS(2)=(YPLUS(2))/AKAPA+PET
      VPLUS(3)=(YPLUS(3))/AKAPA+PET
      ELSE
      VPLUS(1)=(LOG(YPLUS(1)))/AKAPA+PET
      VPLUS(2)=(LOG(YPLUS(2)))/AKAPA+PET
      VPLUS(3)=(LOG(YPLUS(3)))/AKAPA+PET
      ENDIF
      
      TT21(NDIM-7)=V1(NDIM-7)
      TT21(NDIM-3)=V1(NDIM-3)
      TT21(NDIM)=V1(NDIM)
      TT21(NDIM-4)=V1(NDIM-4)
      TT21(NDIM-6)=VPLUS(1)
      TT21(NDIM-2)=VPLUS(1)
      TT21(NDIM-1)=VPLUS(1)
      TT21(NDIM-5)=VPLUS(1)
      
      TT21(1+NDIM)=V2(NDIM-7)
      TT21(5+NDIM)=V2(NDIM-3)
      TT21(8+NDIM)=V2(NDIM)
      TT21(4+NDIM)=V2(NDIM-4)
      TT21(2+NDIM)=VPLUS(2)
      TT21(6+NDIM)=VPLUS(2)
      TT21(7+NDIM)=VPLUS(2)
      TT21(3+NDIM)=VPLUS(2)
      
C      TT21(1+2*NDIM)=V3(NDIM-7)
C      TT21(5+2*NDIM)=V3(NDIM-3)
C      TT21(8+2*NDIM)=V3(NDIM)
C      TT21(4+2*NDIM)=V3(NDIM-4)
C      TT21(2+2*NDIM)=VPLUS(3)
C      TT21(6+2*NDIM)=VPLUS(3)
C      TT21(7+2*NDIM)=VPLUS(3)
C      TT21(3+2*NDIM)=VPLUS(3)
            
      PRIT(IDPRIT,NBREL)=0.0
      
C AKO SU CVOROVI 3,7,8 I 4 = 0 

      ELSEIF (ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.
     &AND.ID(1,N7).EQ.0.AND.ID(2,N7).EQ.0.AND.ID(3,N7).EQ.0.
     &AND.ID(1,N8).EQ.0.AND.ID(2,N8).EQ.0.AND.ID(3,N8).EQ.0.
     &AND.ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0.) THEN
     
      YP(1,I)=H(I)*CORD(1,I)
      YP(2,I)=H(I)*CORD(2,I)
      YP(3,I)=H(I)*CORD(3,I)
      
      YPLUS(1)=GUSM*YP(1,I)*VTANG(1,I)/AMI
      YPLUS(2)=GUSM*YP(2,I)*VTANG(2,I)/AMI
      YPLUS(3)=GUSM*YP(3,I)*VTANG(3,I)/AMI
      
      IF (YPLUS(1).EQ.0.OR.YPLUS(2).EQ.0.OR.YPLUS(3).EQ.0) THEN
      VPLUS(1)=(YPLUS(1))/AKAPA+PET
      VPLUS(2)=(YPLUS(2))/AKAPA+PET
      VPLUS(3)=(YPLUS(3))/AKAPA+PET
      ELSE
      VPLUS(1)=(LOG(YPLUS(1)))/AKAPA+PET
      VPLUS(2)=(LOG(YPLUS(2)))/AKAPA+PET
      VPLUS(3)=(LOG(YPLUS(3)))/AKAPA+PET
      ENDIF
            
      TT21(NDIM-5)=V1(NDIM-5)
      TT21(NDIM-1)=V1(NDIM-1)
      TT21(NDIM)=V1(NDIM)
      TT21(NDIM-4)=V1(NDIM-4)
      TT21(NDIM-6)=VPLUS(1)
      TT21(NDIM-2)=VPLUS(1)
      TT21(NDIM-3)=VPLUS(1)
      TT21(NDIM-7)=VPLUS(1)
      
      TT21(3+NDIM)=V2(NDIM-5)
      TT21(7+NDIM)=V2(NDIM-1)
      TT21(8+NDIM)=V2(NDIM)
      TT21(4+NDIM)=V2(NDIM-4)
      TT21(2+NDIM)=VPLUS(2)
      TT21(6+NDIM)=VPLUS(2)
      TT21(5+NDIM)=VPLUS(2)
      TT21(1+NDIM)=VPLUS(2)
      
c      TT21(3+2*NDIM)=V3(NDIM-5)
c      TT21(7+2*NDIM)=V3(NDIM-1)
c      TT21(8+2*NDIM)=V3(NDIM)
c      TT21(4+2*NDIM)=V3(NDIM-4)
c      TT21(2+2*NDIM)=VPLUS(3)
c      TT21(6+2*NDIM)=VPLUS(3)
c      TT21(5+2*NDIM)=VPLUS(3)
c      TT21(1+2*NDIM)=VPLUS(3)
            
      PRIT(IDPRIT,NBREL)=0.0
      
      
      ENDIF
 101  CONTINUE 
      ENDIF
      
  
C------------------------------------------------------------   
      IF ( DET.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000) NBREL,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     *WRITE(3,6000) NBREL,KFIX,R,S,T,DET
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
C

      END
C=======================================================================