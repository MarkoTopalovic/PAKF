C=========================================================================
C============================================================================
C   SUBROUTINE EXIM2D
C              MASAM
C              ZADPRI
C              REZON1
C              LAGRAN
C              REMESH
C              POCET
C              ALEHV
C              ALEHV3
C              ALEHV4
C              IDALF
C              ADDST1
C              EXPAN
C              EXPAN1
C              BINARY
C              ASCIIF
C              BIFURC
C=========================================================================
C=========================================================================
      SUBROUTINE EXIM2D(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,AMASA,VMESH,CCORD,VELOC,IDALE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine RACU2D is used for 2D analysis
CE It is used global loop per elements
C

C SAMO PRIVREMENO UBACEN COMMON ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,DPRES
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)

      CHARACTER*250NASLOV
      DIMENSION GNODE(2,6,*),ALEVO(*),DESNO(*),SILE(*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION SKEF(NDES,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
      DIMENSION SPSIL(2,*),PRES(3,*)
      DIMENSION AMASA(*),VMESH(3,*),CCORD(3,*),VELOC(3,*)
      DIMENSION IDALE(3,*)

      DIMENSION X(9),Y(9),LM2(4*9),TT21(4*9),TT210(4*9),TT21A(2*9)
      DIMENSION R(3,3),S(3,3),W(3,3)
C      DIMENSION R1(2),S1(2),W1(2)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION AKVV2(9,9)
      DIMENSION AKMIV2(9,9)
      DIMENSION AKMIV3(9,9)
      DIMENSION A12(9,9)
      DIMENSION A21(9,9)
      DIMENSION AKK(9,9)
      DIMENSION AMV2(9,9)
      DIMENSION C(9,9)
      DIMENSION AJV2V2(9,9)
      DIMENSION AJV3V3(9,9)
      DIMENSION AJV2V3(9,9)
      DIMENSION AJV3V2(9,9)
      DIMENSION AKTV2(9,9)
      DIMENSION AKTV3(9,9)
      DIMENSION AKV2P(9,4)
      DIMENSION AKV3P(9,4),AKV2P1(9,4)
      DIMENSION AJK(9,9)
      DIMENSION RS2(9)
      DIMENSION RS3(9)
      DIMENSION RB2(9),RBP(9)
      DIMENSION RB3(9)
      DIMENSION F31(31)
      DIMENSION AKUAX(9,9),AKVAX(9,9),APAXIX(9,9),APAXIY(9,9)
      DIMENSION POISS(4,4)


      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      W(3,1)=0.55555555555556
      W(3,2)=0.88888888888889
      W(3,3)=0.55555555555556


      R(2,1)=0.57735026918963
      R(2,2)=-0.5773502691896
      S(2,1)=0.57735026918963
      S(2,2)=-0.5773502691896
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      W(1,1)=2.0


C      ISTEP=1
      ISTEP=2

 100  CALL CLEAR(ALEVO,NWK)
      CALL CLEAR(DESNO,NWK)
      CALL CLEAR(SILE,JEDN)
      CALL CLEAR(VELOC,3*NPT)


C GLAVNA PETLJA PO ELEMENTIMA
C GLOBAL LOOP PER ELEMENTS
C NBREL is counter of elements

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      VOLUM=0.D0
      DO 125 I=1,4*NDIM
      TT210(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0
C=========================================================================
      DO 130 KLM=1,NDIM
      X(KLM)=CCORD(1,NEL(KLM,NBREL))
      Y(KLM)=CCORD(2,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      IF (KLM.LE.4) LM2(KLM+2*NDIM)=ID(4,NEL(KLM,NBREL))
      LM2(KLM+2*NDIM+4)=ID(5,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
      DO I=1,NDIM
       NODE=NEL(I,NBREL)
       TT21A(I)=VMESH(1,NODE)
       TT21A(I+NDIM)=VMESH(2,NODE)
      ENDDO

      DO 140 KLM=1,NDIM
      DO 135 NR=1,2
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
      IF (KLM.LE.4) THEN
       TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
       TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
      TT21(KLM+2*NDIM+4)=GNODE(2,5,NEL(KLM,NBREL))
      TT210(KLM+2*NDIM+4)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE

C      DO 140 KLM=1,NDIM
C      DO 135 NR=1,2
C      NODE=NEL(KLM,NBREL)
C      JJ=ID(NR,NODE)
C      IF (JJ.NE.0) THEN
C       TT21(KLM+(NR-1)*NDIM)=TT1(JJ)
C       TT210(KLM+(NR-1)*NDIM)=TT10(JJ)
C      ENDIF
C 135  CONTINUE
C      IF (KLM.LE.4) THEN
C       JJ=ID(4,NODE)
C       IF (JJ.NE.0) THEN
C        TT21(KLM+2*NDIM)=TT1(JJ)
C        TT210(KLM+2*NDIM)=TT10(JJ)
C       ENDIF
C      ENDIF
C      JJ=ID(5,NODE)
C      IF (JJ.NE.0) THEN
C        TT21(KLM+2*NDIM+4)=TT1(JJ)
C        TT210(KLM+2*NDIM+4)=TT10(JJ)
C      ENDIF
C 140  CONTINUE
C
C=======================================================================
C=======================================================================
      DO I=1,3*NDIM+4
       IF(KKORAK.EQ.1.AND.ITER.EQ.1) TT210(I)=TT21(I)
      ENDDO
C=======================================================================
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
      IF (N.LE.4) THEN
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
       AKV2P1(K,N)=0.D0
       IF (K.LE.4) THEN
        POISS(K,N)=0.D0
       ENDIF
      ENDIF
      AKVV2(K,N)=0.D0
      APAXIX(K,N)=0.D0
      APAXIY(K,N)=0.D0
      AKUAX(K,N)=0.D0
      AKVAX(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKMIV2(K,N)=0.D0
      AKMIV3(K,N)=0.D0
      A12(K,N)=0.D0
      A21(K,N)=0.D0
      AMV2(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V2(K,N)=0.D0
      AKK(K,N)=0.D0
      AJK(K,N)=0.D0
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
       RS2(K)=0.D0
       RS3(K)=0.D0
       RB2(K)=0.D0
       RBP(K)=0.D0
       RB3(K)=0.D0
  163 CONTINUE


C===========================================================================

C INTEGRACIJA U GAUSOVIM TACKAMA
C integration per gauss points

      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT

C subroutine FNTERP return to us interpolation functions...

      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)

       CALL ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
C       CALL AXISYF(INDAX,DEBLJP,X,HP,4)

      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

      VOLUM=VOLUM+WDT
C       WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ*DEBLJ
C       WDTP=W(IBRGT,I)*W(IBRGT,J)*DETJ*DEBLJP

      IF (INDAMI.EQ.1) CALL NENJUT(ZVHX,ZVHY,TT210)
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM
      IF (INDAX.EQ.1) THEN
      AKUAX(K,N)=AKUAX(K,N)-WDT*H(K)*(ZVHX(N)/DEBLJ-H(N)/(DEBLJ**2))*AMI
      AKVAX(K,N)=AKVAX(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AMI
      ENDIF
      AKVV2(K,N)=AKVV2(K,N)+WDT*((H(K)*HV2*ZVHX(N)+
     1H(K)*HV3*ZVHY(N))*GUSM)
      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVXT
     1*GUSM*CC)
      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVYT
     1*GUSM*CC)
      AKMIV2(K,N)=AKMIV2(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI)
      AKK(K,N)=AKK(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AKT)
      IF (INDAX.EQ.1) THEN
       AKK(K,N)=AKK(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AKT
      ENDIF

      AKMIV3(K,N)=AKMIV2(K,N)
      AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM
C      AJV2V2(K,N)=AJV2V2(K,N)+WDT*(H(K)*ZVXV2*H(N)*GUSM)
C      AJV3V3(K,N)=AJV3V3(K,N)+WDT*(H(K)*ZVYV3*H(N)*GUSM)
C      AJV2V3(K,N)=AJV2V3(K,N)+WDT*(H(K)*ZVYV2*H(N)*GUSM)
C      AJV3V2(K,N)=AJV3V2(K,N)+WDT*(H(K)*ZVXV3*H(N)*GUSM)
       IF (N.LE.4.AND.K.LE.4) THEN
C      POISS(K,N)=POISS(K,N)+WDT*(ZVHXP(K)*ZVHXP(N)+ZVHYP(K)*ZVHYP(N))
      POISS(K,N)=POISS(K,N)+WDT*(ZVHX(K)*ZVHX(N)+ZVHY(K)*ZVHY(N))
       ENDIF
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P(K,N)=AKV2P(K,N)+WDT*(ZVHX(K)*HP(N))
      AKV3P(K,N)=AKV3P(K,N)+WDT*(ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
  164 CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
C       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(0.D0+BETA*TETAO)*WDT
C       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(0.D0+BETA*TETAO)*WDT
C       RBP(K)=RBP(K)+(ZVHXP(K)*FB2+ZVHYP(K)*FB3)*GUSM*WDT
       RBP(K)=RBP(K)+(ZVHX(K)*FB2+ZVHY(K)*FB3)*GUSM*WDT
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE

C===========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PENALTY FAKTOROM
C reduced integration if penalty function is defined

      IF (PENALT.GT.0.D0) THEN
      DO 200 I=1,IBRGT-1
      DO 195 J=1,IBRGT-1
      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
      DO 190 K=1,NDIM
      DO 185 N=1,NDIM
      AKMIV2(K,N)=AKMIV2(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT
      AKMIV3(K,N)=AKMIV3(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT
      A12(K,N)=A12(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      AJV2V3(K,N)=AJV2V3(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      A21(K,N)=A21(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
      AJV3V2(K,N)=AJV3V2(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
       IF (INDAX.EQ.1) THEN
       APAXIX(K,N)=APAXIX(K,N)+PENALT*WDT*(ZVHX(K)*H(N)/DEBLJ)
       APAXIY(K,N)=APAXIY(K,N)+PENALT*WDT*(ZVHY(K)*H(N)/DEBLJ)
       ENDIF
      AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
  185 CONTINUE
  190 CONTINUE
  195 CONTINUE
  200 CONTINUE
      ENDIF
C===========================================================================

C======================================================================= 
C POVRSINSKE SILE
C surface forces

       INDX=0
       INDY=0

      DO 250 JBRPS=1,MAXSIL
      IF (NBREL.EQ.NGPSIL(1,JBRPS)) THEN
C       WRITE(IIZLAZ,*)'IND,ELEM',NGPSIL(4,JBRPS),NBREL
        IF (PENALT.LE.1.D0) IBRGT=IBRGT+1
C        INDX=NGPSIL(6,JBRPS)
C        INDY=NGPSIL(7,JBRPS)

C        IF (NGPSIL(4,JBRPS).EQ.0) THEN
C        INDX=1
C        INDY=1
C        ELSE IF(NGPSIL(6,JBRPS).EQ.1) THEN
C        INDX=1
C        INDY=0
C        ELSE IF(NGPSIL(7,JBRPS).EQ.1) THEN
C        INDX=0
C        INDY=1
C        ENDIF

       TTAU=0.D0
       NPARAM=1
       IF(NGPSIL(4,JBRPS).EQ.3) NPARAM=2
       
      NODE1=NGPSIL(2,JBRPS)
      NODE2=NGPSIL(3,JBRPS)
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)
      DO 225 I=1,IBRGT-1
      IF ((NODE1.EQ.N1 .AND. NODE2.EQ.N4).OR.
     1(NODE1.EQ.N4 .AND. NODE2.EQ.N1)) THEN
        CALL FNTERP(1.D0,S(IBRGT-1,I),NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
        GOTO 205
      ENDIF 

      IF ((NODE1.EQ.N2 .AND. NODE2.EQ.N3).OR.
     1(NODE1.EQ.N3 .AND. NODE2.EQ.N2)) THEN
      CALL FNTERP(-1.D0,S(IBRGT-1,I),NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
        GOTO 205
      ENDIF 

      IF ((NODE1.EQ.N1 .AND. NODE2.EQ.N2).OR.
     1(NODE1.EQ.N2 .AND. NODE2.EQ.N1)) THEN
      CALL FNTERP(R(IBRGT-1,I),1.D0,NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
        GOTO 205
      ENDIF 


      IF ((NODE1.EQ.N3 .AND. NODE2.EQ.N4).OR.
     1(NODE1.EQ.N4 .AND. NODE2.EQ.N3)) THEN
      CALL FNTERP(R(IBRGT-1,I),-1.D0,NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
        GOTO 205
      ENDIF 


 205  CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
C      WDT=W(IBRGT-1,I)*DETJS*DEBLJ
      WDT=W(IBRGT-1,I)*DETJS
      IF(NGPSIL(4,JBRPS).EQ.3) GOTO 207
      DO 206 K=1,4
       FS2=
     &((TT21(K)-TT210(K))*ELX+(TT21(K+NDIM)-TT210(K+NDIM))*ELY)/TIME
       RS2(K)=RS2(K)+W(IBRGT-1,I)*(HP(K)*FS2*DETJS)

C       RS2(K)=RS2(K)+W(IBRGT-1,I)*(H(K)*FS2*DETJS)
C       RS3(K)=RS3(K)+WDT*(H(K)*FS3)*NGPSIL(5,JBRPS)
 206  CONTINUE
C 207   TTAU=TTAU+WDT*TAU
 207   WRITE(IIZLAZ,*)'NBREL= ',NBREL,'TAU= ',TAU
 225  CONTINUE

C       IF(NGPSIL(4,JBRPS).EQ.3) THEN
C        WRITE(IIZLAZ,*)'ELEM',NBREL,'TTAU',TTAU
C       ENDIF


      IF (PENALT.LE.1.D0) IBRGT=IBRGT-1
      ENDIF
 250  CONTINUE

C C matrix include heat conduction


      DO 255 I=1,NDIM
      DO 254 J=1,NDIM

      C(I,J)=AMV2(I,J)*CC              
C      AKK(I,J)=AKMIV2(I,J)*AKT/AMI
      AJK(I,J)=AKVV2(I,J)*CC

 254  CONTINUE
 255  CONTINUE

C=========================================================================
C INCIJALIZACIJA MATRICE SKEF I F31
      DO 260 I=1,NDES
      DO 258 J=1,NDES
       SKEF(I,J)=0.D0
 258  CONTINUE
       F31(I)=0.D0
 260  CONTINUE
C=========================================================================
C PRVI I CETVRTI KORAK:

	IF (ISTEP.EQ.1) THEN


       DO I=1,NDES
        LM2(I)=0
       ENDDO

 

       
      DO I=1,4
      LM2(I+2*NDIM)=ID(4,NEL(I,NBREL))
         K=I+2*NDIM
        F31(K)=F31(K)+RS2(I)
       DO J=1,NDIM
        F31(K)=F31(K)-
     &(1.D0/TIME)*(AKV2P(J,I)*TT210(J)+AKV3P(J,I)*TT210(J+NDIM))*GUSM
       ENDDO        
        F31(K)=F31(K)+RBP(I)
      ENDDO        

      DO I=1,4
       DO J=1,4
       SKEF(I+2*NDIM,J+2*NDIM)=POISS(I,J)
       ENDDO        
      ENDDO     

C      CALL WRRF(SKEF,NDES*NDES,'SKEF1',IIZLAZ)      
C      CALL WRRF(F31,NDES,'F31  ',IIZLAZ)      
C      CALL WRRF(TT210,NDES,'TT210',IIZLAZ)      
        
      ELSE IF (ISTEP.EQ.4) THEN

       DO I=1,NDES
        LM2(I)=0
       ENDDO

       
      DO I=1,4
      LM2(I+2*NDIM)=ID(4,NEL(I,NBREL))
         K=I+2*NDIM
        F31(K)=F31(K)+RS2(I)
       DO J=1,NDIM
        F31(K)=F31(K)-
     &(1.D0/TIME)*(AKV2P(J,I)*TT21(J)+AKV3P(J,I)*TT21(J+NDIM))*GUSM
       ENDDO        
        F31(K)=F31(K)+RBP(I)
      ENDDO        

      DO I=1,4
       DO J=1,4
       SKEF(I+2*NDIM,J+2*NDIM)=POISS(I,J)
       ENDDO        
      ENDDO     

      ELSE IF (ISTEP.EQ.2) THEN
C=========================================================================
C DRUGI  KORAK
      
       DO I=1,NDES
        LM2(I)=0
       ENDDO


      DO 270 I=1,NDIM
C      LM2(I)=ID(1,NEL(I,NBREL))
C      LM2(I+NDIM)=ID(2,NEL(I,NBREL))
C      LM2(I)=NEL(I,NBREL)
C      LM2(I+NDIM)=NEL(I,NBREL)+NPT
C      LM2(I+NDIM)=NEL(I,NBREL)
      DO 265 J=1,NDIM
      SKEF(I,J)=-AKMIV2(I,J)
      IF (J.LE.4) THEN
      SKEF(I,J+2*NDIM)=AKV2P(I,J)
      SKEF(I+NDIM,J+2*NDIM)=AKV3P(I,J)
      ENDIF
      SKEF(I+NDIM,J+NDIM)=-AKMIV3(I,J)
 265  CONTINUE
 270  CONTINUE


        DELI=1.D0

      DO  I=1,NDIM
        F31(I)=F31(I)+RB2(I)*TIME/GUSM
        F31(I+NDIM)=F31(I+NDIM)+RB3(I)*TIME/GUSM
      ENDDO

      DO  I=1,2*NDIM
       DO  J=1,2*NDIM+4
        F31(I)=F31(I)+SKEF(I,J)*TT21(J)*(TIME/GUSM)
       ENDDO
      ENDDO



C      CALL WRRF(F31,NDES,'F31  ',IIZLAZ)      
C      CALL WRRF(TT210,NDES,'TT210',IIZLAZ)      
C      CALL WRRF(TT21,NDES,'TT21 ',IIZLAZ)      
C      CALL WRRF(TT21A,2*9,'TT21A',IIZLAZ)      
C      CALL WRRF(AKVV2,9*9,'AKVV2',IIZLAZ)      


C      DO  I=1,NDES
C       DO  J=1,NDES
C        SKEF(I,J)=0.D0
C       ENDDO
C      ENDDO

      ELSE IF (ISTEP.EQ.5) THEN
C=========================================================================
C   PETI KORAK
      
       DO I=1,NDES
        LM2(I)=0
       ENDDO


      DO  I=1,NDIM
C      LM2(I)=ID(1,NEL(I,NBREL))
C      LM2(I+NDIM)=ID(2,NEL(I,NBREL))
      LM2(I)=NEL(I,NBREL)
C      LM2(I+NDIM)=NEL(I,NBREL)+NPT
      LM2(I+NDIM)=NEL(I,NBREL)
      DO  J=1,NDIM
      SKEF(I,J)=-AKMIV2(I,J)
      IF (J.LE.4) THEN
      SKEF(I,J+2*NDIM)=AKV2P(I,J)
      SKEF(I+NDIM,J+2*NDIM)=AKV3P(I,J)
      ENDIF
      SKEF(I+NDIM,J+NDIM)=-AKMIV3(I,J)
      ENDDO
      ENDDO

        DELI=1.D0

      DO  I=1,NDIM
        F31(I)=F31(I)+RB2(I)*TIME/GUSM
        F31(I+NDIM)=F31(I+NDIM)+RB3(I)*TIME/GUSM
      ENDDO

      DO  I=1,2*NDIM
       DO  J=1,2*NDIM+4
        F31(I)=F31(I)+SKEF(I,J)*TT21(J)*(TIME/GUSM)
       ENDDO
      ENDDO

C      DO  I=1,NDES
C       DO  J=1,NDES
C        SKEF(I,J)=0.D0
C       ENDDO
C      ENDDO

C=========================================================================
C FAZA 2 KOD ALE FORMULACIJE :

C TRECI KORAK:

	 ELSE IF(ISTEP.EQ.3) THEN

       DO I=1,NDES
        LM2(I)=0
       ENDDO


      DO  I=1,NDIM
      LM2(I)=NEL(I,NBREL)
      LM2(I+NDIM)=NEL(I,NBREL)
       DO  J=1,NDIM
        F31(I)=F31(I)+(-AKVV2(I,J)*TIME/GUSM)*TT210(J)
        F31(I+NDIM)=F31(I+NDIM)+(-AKVV2(I,J)*TIME/GUSM)*TT210(J+NDIM)
      ENDDO
      ENDDO


      ENDIF






     



C==========================================================================

C==========================================================================
C==========================================================================



C============================================================================
C  FOR FLUX CALCULATION
C       DO I=1,NDIM
C        N=NEL(I,NBREL)
C         POT=0.D0
C         DO J=1,NDES
C          POT=POT+SKEF(I+2*NDIM+4,J)*TT21(J)
C         ENDDO
C         PRES(1,N)=PRES(1,N)+POT
C       ENDDO
C============================================================================


C       IF (ISTEP.EQ.1) CALL MASAM(AMASA,VOLUM,NEL,NDIM,NBREL,2)
       IF (ISTEP.EQ.2) CALL MASAM(AMASA,VOLUM,NEL,NDIM,NBREL,2)


      INDSK=0
      IF(ISTEP.EQ.1.OR.ISTEP.EQ.4) THEN
       INDSK=1
       CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F31,MAXA,LM2,NDES,INDSK)
      ELSE
       CALL ADDST1(VELOC,F31,NEL,NBREL,NDIM,IDALE)
      ENDIF

C RACUNANJE SILA INTERAKCIJE KOJIMA FLUID DELUJE NA ZIDOVE
c calculation of interaction forces between fluid and solid
C	 IF (NSTAC.EQ.0) THEN
C        DO I=1,NDIM
C          II=I+NDIM
C         DO J=1,NDIM
C          JJ=J+NDIM
C          SKEF(I,J)=SKEF(I,J)-AMV2(I,J)*TT210(J)/TIME
C          SKEF(II,JJ)=SKEF(II,JJ)-AMV2(I,J)*TT210(JJ)/TIME
C         ENDDO
C        ENDDO
C       ENDIF

C       FACAXY=1.D0
C       IF(INDAX.EQ.1) FACAXY=8.D0*DATAN(1.D0)

C       DO I=1,NDIM
C       N=NEL(I,NBREL)
C        DO J=1,NDES
C         P1=SKEF(I,J)*TT21(J)
C         P2=SKEF(I+NDIM,J)*TT21(J)
C         SPSIL(1,N)=SPSIL(1,N)-P1/FACAXY
C         SPSIL(2,N)=SPSIL(2,N)-P2/FACAXY
C       ENDDO
C      ENDDO



C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C end of loop per elements
C=======================================================================
 400  CONTINUE


C       CALL WRRF(SILE,JEDN,'SILE1',IIZLAZ)

      
C ONLY FOR EXAMPLE SOLITARY WAVE PROPAGATION
C      IF (ISTEP.EQ.1.OR.ISTEP.EQ.4) THEN
C       CALL ZADPRI(ALEVO,ID,SILE,MAXA,CORD,CCORD,VVREME,ISTEP,NPT)
C      ENDIF

      IF (ISTEP.EQ.2.OR.ISTEP.EQ.3.OR.ISTEP.EQ.5) THEN
       IF (ISTEP.EQ.5) THEN
          CALL REZON1(GNODE,VELOC,AMASA,NPT,2)
       ELSE
          CALL REZON1(GNODE,VELOC,AMASA,NPT,1)
       ENDIF 
      ELSE
       CALL UACTCF(ALEVO,DESNO,SILE,MAXA,JEDN,1)
       CALL UACTCF(ALEVO,DESNO,SILE,MAXA,JEDN,2)
      ENDIF


      IF (ISTEP.EQ.1.OR.ISTEP.EQ.4) THEN
       DO I=1,NPT
           JJ=ID(4,I)
          IF (JJ.NE.0) THEN 
           GNODE(2,4,I)=SILE(JJ)
          ENDIF
       ENDDO
      ENDIF        
      

C ONLY FOR EXAMPLE SOLITARY WAVE PROPAGATION
C       IF (ISTEP.EQ.2) THEN
C        CALL LAGRAN(GNODE,NPT,CCORD,TIME)
C        CALL REMESH(CCORD,VMESH,NPT,TIME,CORD)
C       ENDIF


      WRITE(*,*)'ISTEP= ',ISTEP



       IF (ISTEP.LE.4) THEN 
        ISTEP=ISTEP+1
           GO TO 100
        ENDIF

       IF (ISTEP.EQ.5) RETURN



C End of subroutine
      END
C=======================================================================
C==========================================================================
      SUBROUTINE MASAM(AMASA,VOLUM,NEL,NDIM,NBREL,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DIMENSION AMASA(*),NEL(NDIM+1,*)

       IF (NDIM.EQ.8.AND.NETIP.EQ.3) THEN
       DO I=1,NDIM
       NODE=NEL(I,NBREL)
       AMASA(NODE)=AMASA(NODE)+VOLUM/8.D0
      ENDDO
      ENDIF

       IF (NDIM.EQ.4.AND.NETIP.EQ.2) THEN
      DO I=1,NDIM
       NODE=NEL(I,NBREL)
       AMASA(NODE)=AMASA(NODE)+VOLUM/4.D0
      ENDDO
      ENDIF

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ZADPRI(ALEVO,ID,SILE,MAXA,CORD,CCORD,VVREME,ITER,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
  
      DIMENSION ALEVO(*),SILE(*),CORD(3,*),CCORD(3,*)
      DIMENSION ID(6,*),MAXA(*)


       G=9.81D0
       RO=1000.D0
       H=2.D0
       T=VVREME
       D=10.D0

       C=DSQRT(G*D*(1.D0+0.5D0*H/D))
       Y1=DSQRT(3.D0*H/(4.D0*(D**3)))

       DO I=1,NPT
       KOL=I/6
         IF ((I-KOL*6).EQ.0) THEN
       JJP=ID(4,I)
       X=CORD(1,I)
       Y=CORD(2,I)
       Y2=Y1*(X-C*T)
       Y3=(1.D0/DCOSH(Y2))**2
       U=DSQRT(G*D)*(H/D)*Y3
       V1=DTANH(Y2)
       V2=DSQRT(3.D0*G*D)*(H/D)**(3.D0/2.D0)*(Y/D)
       V=V2*Y3*V1
       HH=D+H*Y3
C       P=RO*G*(HH-Y)
       P=0.D0
C        CCORD(2,I)=HH
        IF (JJP.NE.0.AND.(ITER.EQ.1.OR.ITER.EQ.4)) THEN
         ALEVO(MAXA(ID(4,I)))=1.0D35
         SILE(ID(4,I))=1.0D35*P
        ENDIF
         ENDIF
       ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE REZON1(GNODE,VELOC,AMASA,NPT,INDT0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION GNODE(2,6,*),VELOC(3,*),AMASA(*)

      DO NODE=1,NPT
        DO J=1,3
         GNODE(2,J,NODE)=GNODE(INDT0,J,NODE)+VELOC(J,NODE)/AMASA(NODE)
        ENDDO
      ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE LAGRAN(GNODE,NPT,CCORD,TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION GNODE(2,6,*),CCORD(3,*)
       

C SUBROUTINE FOR MOVING NODES-LAGRANGIAN METHOD
       DO 10 NODE=1,NPT
C        JX=ID(1,NODE)
C        JY=ID(2,NODE)
C      IF (JX.NE.0) CCORD(1,NODE)=CCORD(1,NODE)+
C     &0.5D0*(TT1(JX)+TT10(JX))*TIME  
C      IF (JY.NE.0) CCORD(2,NODE)=CCORD(2,NODE)+
C     &0.5D0*(TT1(JY)+TT10(JY))*TIME  

      CX=0.5D0*(GNODE(2,1,NODE)+GNODE(1,1,NODE))*TIME
      CY=0.5D0*(GNODE(2,2,NODE)+GNODE(1,2,NODE))*TIME
      CCORD(1,NODE)=CCORD(1,NODE)+CX
      CCORD(2,NODE)=CCORD(2,NODE)+CY

 10   CONTINUE
      
      END
C==========================================================================
C==========================================================================
      SUBROUTINE REMESH(CCORD,VMESH,NPT,TIME,CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CCORD(3,*),VMESH(3,*),CORD(3,*)
       
C      NH=11
      NH=6
C      NH=5

      DO J=0,NPT/NH-1
      IREF=NH+NH*J
      Y=CCORD(2,IREF)
C      X=CCORD(1,IREF)
      X=CORD(1,IREF)
      RNEWX=X
      ROLDX=CCORD(1,IREF)
      VMESH(1,IREF)=(RNEWX-ROLDX)/TIME
      VMESH(2,IREF)=0.D0
      CCORD(1,IREF)=RNEWX
      DO I=1,NH-1
       NODE=I+NH*J    
       OLDY=CCORD(2,NODE)   
       OLDX=CCORD(1,NODE)   
       ANEWY=(Y/(NH-1))*(I-1)
       ANEWX=X
       VMESH(1,NODE)=(ANEWX-OLDX)/TIME
       VMESH(2,NODE)=(ANEWY-OLDY)/TIME
       CCORD(1,NODE)=ANEWX       
       CCORD(2,NODE)=ANEWY
      ENDDO
      ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE POCET(GNODE,IDALE,CORD,CCORD,VVREME,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
  
      DIMENSION GNODE(2,6,*),IDALE(3,*),CORD(3,*),CCORD(3,*)

       G=9.81D0
       RO=1000.D0
       H=2.D0
       T=VVREME
       D=10.D0

       C=DSQRT(G*D*(1.D0+0.5D0*H/D))
       Y1=DSQRT(3.D0*H/(4.D0*(D**3)))

       DO I=1,NPT
       JJX=IDALE(1,I)
       JJY=IDALE(2,I)
C       JJP=ID(4,I)
       X=CORD(1,I)
       Y=CORD(2,I)
       Y2=Y1*(X-C*T)
       Y3=(1.D0/DCOSH(Y2))**2
       U=DSQRT(G*D)*(H/D)*Y3
       V1=DTANH(Y2)
       V2=DSQRT(3.D0*G*D)*(H/D)**(3.D0/2.D0)*(Y/D)
       V=V2*Y3*V1
       HH=D+H*Y3
       P=RO*G*(HH-Y)
       KOL=I/6
         IF ((I-KOL*6).EQ.0) THEN
C          CORD(2,I)=HH
          CCORD(2,I)=HH
C          WRITE(IIZLAZ,*)'CVOR: ',I
         ENDIF
       IF (JJX.NE.0) GNODE(2,1,I)=U
       IF (JJY.NE.0) GNODE(2,2,I)=V
C       IF (JJP.NE.0) TT1(JJP)=P
        GNODE(2,4,I)=P
       ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TT21(*),TT21A(*),H(*),V1(9),V2(9)

      DO I=1,NDIM
       V1(I)=TT21(I)-TT21A(I)
       V2(I)=TT21(I+NDIM)-TT21A(I+NDIM)
      ENDDO

      HV2=DOT(H,V1,NDIM)
      HV3=DOT(H,V2,NDIM)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ALEHV4(TT21,TT21A,HV2,HV3,NDIM,H,ZVXV2,ZVYV3,ZVYV2,
     &ZVXV3,ZVHX,ZVHY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TT21(*),TT21A(*),H(*),V1(9),V2(9),ZVHX(*),ZVHY(*)

      DO I=1,NDIM
       V1(I)=TT21(I)-TT21A(I)
       V2(I)=TT21(I+NDIM)-TT21A(I+NDIM)
      ENDDO

      HV2=DOT(H,V1,NDIM)
      HV3=DOT(H,V2,NDIM)

      ZVXV2=DOT(ZVHX,V1,NDIM)
      ZVYV3=DOT(ZVHY,V2,NDIM)
      ZVYV2=DOT(ZVHY,V1,NDIM)
      ZVXV3=DOT(ZVHX,V2,NDIM)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TT21(*),TT21A(*),H(*),V1(21),V2(21),V3(21)

      DO I=1,NDIM
       V1(I)=TT21(I)-TT21A(I)
       V2(I)=TT21(I+NDIM)-TT21A(I+NDIM)
       V3(I)=TT21(I+2*NDIM)-TT21A(I+2*NDIM)
      ENDDO

      HV1=DOT(H,V1,NDIM)
      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE IDALF(IDALE,ID,NPT,NEQF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IDALE(3,*),ID(6,*)

       KK=0
      DO I=1,NPT
       DO J=1,3
        IDALE(J,I)=ID(J,I)      
        ID(J,I)=0
       ENDDO
      DO JJ=1,5
       IF (ID(JJ,I).NE.0) THEN
        KK=KK+1
        ID(JJ,I)=KK
       ELSE
        ID(JJ,I)=0
       ENDIF
      ENDDO
      ENDDO
     
      NEQF=KK

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ADDST1(VELOC,F31,NEL,NBREL,NDIM,IDALE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VELOC(3,*),F31(*)
      DIMENSION NEL(NDIM+1,*),IDALE(3,*)
   
       DO 10 I=1,NDIM
        NODE=NEL(I,NBREL)
        JX=IDALE(1,NODE)
        JY=IDALE(2,NODE)
        JZ=IDALE(3,NODE)
        IF (JX.NE.0) VELOC(1,NODE)=VELOC(1,NODE)+F31(I) 
        IF (JY.NE.0) VELOC(2,NODE)=VELOC(2,NODE)+F31(I+NDIM) 
        IF (JZ.NE.0) VELOC(3,NODE)=VELOC(3,NODE)+F31(I+2*NDIM) 
 10    CONTINUE
      END
C==========================================================================
C==========================================================================
      SUBROUTINE JocaFluidOnly(CCORD,IZLAZF,NEL,NDIM,NPT,T,ID,GNODE,
     &CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CCORD(3,*), NEL(NDIM,*),CORD(3,*)
	DIMENSION ID(6,*)
      DIMENSION GNODE(2,6,*)  


      IU=NEL(4,1)-NEL(1,1)-1
	IW=NPTF/(IU+1)-1


      PI=4.D0*DATAN(1.D0)

      P=1.d0
      AN=2*PI/P
      C=0.25d0
      FI=(1.D0+C)**(1.D0/3.D0)
      AK=(FI-1.D0)/(FI+1.D0)
       
      N0=NPT-IU
      TT=T+P/4.D0
      XC=CORD(1,N0)
	YC=CORD(2,N0)

      DO I=1,NPT
        X=CORD(1,I)
        Y=CORD(2,I)
C       X0=CCORD(1,I)
C        Y0=CCORD(2,I)
        CCORD(1,I)=XC+(X-XC)*(1.D0+5.D0*AK*DSIN(AN*TT))
        CCORD(2,I)=YC+(Y-YC)*(1.D0+AK*DSIN(AN*TT))
C        VMESH(1,I)=(CCORD(1,I)-X0)/TIME
C        VMESH(2,I)=(CCORD(2,I)-Y0)/TIME
       IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0) THEN
        VX=5.D0*AN*AK*(X-XC)*DCOS(AN*TT)
        VY=AN*AK*(Y-YC)*DCOS(AN*TT)
        GNODE(2,1,I)=VX
        GNODE(2,2,I)=VY
       ENDIF
      ENDDO  

      END
C==========================================================================
C==========================================================================
      SUBROUTINE joca(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA



      do i=1,npt
	 if (id (1,i).eq.0.and.id (2,i).eq.0) then 
	   n2=i
	goto 100
       endif
	enddo


    
100   continue
      x=cord(1,n2)
	x2=0.4d0


      do j=0,npt,n2
      do i=1+j,n2+j
	  k=kkorak-(kkorak/100)*kkorak
	if (k.le.20)then
	  ccord(1,i)=cord(1,i)-cord(1,i)*(k/20)*x2
	  if (mod(i,n2).eq.0) gnode(2,1,i)=-(cord(1,i)*x2/20.d0)/time
        
      else if (k.gt.20.and.k.le.40) then
	  ccord(1,i)=cord(1,i)-cord(1,i)*x2+(k-20)/20*x2*cord(1,i)
	  if (mod(i,n2).eq.0) gnode(2,1,i)=(cord(1,i)*x2/20.d0)/time
      else if (k.gt.40) then
	  ccord(1,i)=cord(1,i)
	  if (mod(i,n2).eq.0) gnode(2,1,i)=-0.1*(cord(1,i)*x2/20.d0)/time
	endif
	enddo
	enddo



	end
C==========================================================================
      SUBROUTINE EXPAN(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA

      if (iakira.eq.11) then
      call joca(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK)
	return
	endif


C PRIVREMENO, ZA MILOSA, AUG 22, 2006
c      RETURN 
      TT=T-PERIO1/4.D0

      IF (IAKIRA.EQ.0) RETURN

      IF (IAKIRA.EQ.3) THEN
        CALL BIFUR1(GNODE,CORD,T,NPT,TIME,NSTAC,AMI,KKORAK)
        RETURN
      ENDIF
      
      IF (IAKIRA.EQ.4) RETURN

      IF (IAKIRA.EQ.5) THEN
        CALL EXPAN2(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME)
	RETURN
      ENDIF


      IF (IAKIRA.EQ.0.OR.NSTAC.EQ.1) RETURN
      IF (IAKIRA.EQ.2) THEN
      CALL EXPAN1(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,ZADVRE,AMI
     &,KKORAK)
	RETURN
	ENDIF

      PI=4.D0*DATAN(1.D0)

      MOVX=0
      MOVY=0
      IF (MOVEW.EQ.1) MOVX=1
      IF (MOVEW.EQ.2) MOVY=1
      IF (MOVEW.EQ.0) THEN
        MOVX=1
        MOVY=1
      ENDIF
 
      P=PERIO1
      AN=2.D0*PI/P
      C=SVEXC
      FI=(1.D0+C)**(1.D0/3.D0)
      AK=(FI-1.D0)/(FI+1.D0)
      CY=DCEXP
      CX=0.D0

      cx=-0.0215d0
	cy=0.d0


 
C=======================================================================
C=======================================================================
C SAMO ZA KOJICA SPECIJALNO DA PROPUSTI NEKE PRIMERE (10.11.2002):
C=======================================================================
C=======================================================================
      IF (NPRVEL.EQ.1) THEN
       DO I=1,MESH
        GNODE(2,2,I)=ZADVRE(I)*DSIN(AN*T)
        GNODE(2,2,NPT-(MESH+1)+I)=ZADVRE(I)*DSIN(AN*T)
       ENDDO
      RETURN
      ENDIF      
C=======================================================================
C=======================================================================
C=======================================================================

      DO 10 I=1,NPT
        VX=0.D0
        VY=0.D0
        X=CORD(1,I)
        Y=CORD(2,I)
        X0=CCORD(1,I)
        Y0=CCORD(2,I)
         IF (MOVY.EQ.1) THEN
           CCORD(2,I)=CY+(Y-CY)*(1.D0+AK*DSIN(AN*TT))
C           VY=(Y-CY)*AK*AN*DCOS(AN*T)
C Ovde izvrsena izmena za granican uslov skracenja duzine cevi
C odnosno zadavanje brzine na kracoj duzini prave cevi
           VY=(3.D0*ALEN-2.D0*Y)*AK*AN*DCOS(AN*TT)
         ENDIF
         IF (MOVX.EQ.1) THEN
           CCORD(1,I)=CX+(X-CX)*(1.D0+AK*DSIN(AN*TT))
           VX=(X-CX)*AK*AN*DCOS(AN*TT)
         ENDIF
        VMESH(1,I)=(CCORD(1,I)-X0)/TIME
        VMESH(2,I)=(CCORD(2,I)-Y0)/TIME
       IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0) THEN
        IF (NPRVEL.EQ.1.AND.I.LE.MESH) THEN
          GNODE(2,2,I)=ZADVRE(I)*(1.D0+AK*DSIN(AN*TT))
          GOTO 10
        ENDIF
c         GNODE(2,1,I)=VX
c         GNODE(2,2,I)=VY
       ENDIF
  10   CONTINUE
 

      END
C==========================================================================
C==========================================================================
      SUBROUTINE EXPAN2(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*)
      DIMENSION ID(6,*)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA

      PI=4.D0*DATAN(1.D0)

      P=PERIO1
      AN=2*PI/P
      C=SVEXC
      FI=(1.D0+C)**(1.D0/3.D0)
      AK=(FI-1.D0)/(FI+1.D0)
      CCY=DCEXP
      CCX=0.D0

C      RAA=0.00025D0
      RAA=RA
      ALZ=0.0125D0
C      C1=10.5500459D0
C      C2=4.345899419D-3
C      C3=41.38408529D0
C      C4=-8.650533063D-4
      OPEN (47,FILE='CONST.TXT')
      READ(47,*)C1
      READ(47,*)C2
      READ(47,*)C3
      READ(47,*)C4
	CLOSE (47)

C      C1=19.81913D0
C      C2=0.008164D0
C      C3=29.5755D0
C      C4=0.006515D0
      TOLER=1.0D-8

      DO 10 I=1,NPT
        IND=0
        VX=0.D0
        VY=0.D0
        X=CORD(1,I)
        Y=CORD(2,I)
        X0=CCORD(1,I)
        Y0=CCORD(2,I)




      DO J=1,NUMALV
      ALVY=DIS1AL+(J-1)*DISTAL

C Na uspravnim zidovima
      IF((DABS(DABS(Y-ALVY)-RAA).LT.TOLER).AND.(X>RD-ALZ+RAA-TOLER))
     *IND=2

      IF((DABS(DABS(Y-ALVY)-RAA).LT.TOLER).AND.(X>RD-ALZ+RAA-TOLER)
     *.AND.(X<RD-ALZ+2.D0*RAA+TOLER))IND=1

C Donji polukrug
      CX=RD-ALZ+RAA
      CY=ALVY
      D=DISTAN(X,Y,CX,CY)
      IF (DABS(D-RAA)<1000.D0*TOLER) IND=1
C Gornji levi polukrug

      CX=RD-RAA
      CY=ALVY-2.D0*RAA
      D=DISTAN(X,Y,CX,CY)
      IF((DABS(D-RAA)<1000.D0*TOLER).AND.(X.GT.CX-TOLER).AND.
     *(Y.GT.CY-TOLER)) IND=2

C Gornji desni polukrug

      CX=RD-RAA
      CY=ALVY+2.D0*RAA
      D=DISTAN(X,Y,CX,CY)
      IF((DABS(D-RAA)<1000.D0*TOLER).AND.(X.GT.CX-TOLER).AND.
     *(Y.LT.CY+TOLER)) IND=2

      ENDDO


100     IF (IND.EQ.1) THEN
        CCORD(1,I)=X+(C1*X+C2/X)*RD*(AK*DSIN(AN*T))
        VX=(C1*X+C2/X)*RD*AK*AN*DCOS(AN*T)
       ELSEIF(IND.EQ.2) THEN
        CCORD(1,I)=X+(C3*X+C4/X)*RD*(AK*DSIN(AN*T))
        VX=(C3*X+C4/X)*RD*AK*AN*DCOS(AN*T)
       ELSEIF(IND.EQ.0) THEN
        CCORD(1,I)=X*(1.D0+AK*DSIN(AN*T))
        VX=X*AK*AN*DCOS(AN*T)
       ENDIF

           CCORD(2,I)=Y*(1.D0+AK*DSIN(AN*T))
           VY=Y*AK*AN*DCOS(AN*T)

        VMESH(1,I)=(CCORD(1,I)-X0)/TIME
        VMESH(2,I)=(CCORD(2,I)-Y0)/TIME
       IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0) THEN
         GNODE(2,1,I)=VX
         GNODE(2,2,I)=VY
       ENDIF
  10   CONTINUE

      END
C==========================================================================
C==========================================================================
      FUNCTION DISTAN(X1,Y1,X2,Y2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

       DISTAN=DSQRT((X2-X1)**2+(Y2-Y1)**2)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE EXPAN1(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA


      IF (KKORAK.EQ.1) THEN
      DO I=MESH+1,NPT,MESH+1
        X=CORD(1,I)
        VX=RE*AMI/X
        GNODE(2,1,I)=VX
        GNODE(1,1,I)=VX
      ENDDO
      NSTAC=1
      RETURN 
      ENDIF

      IF (KKORAK.EQ.2) THEN
       NSTAC=0
       T=0.D0
      ENDIF
 
      Q=0.99D0

      DO 10 I=MESH+1,NPT,MESH+1
        VX=0.D0
        VY=0.D0
        X0=CCORD(1,I)
        Y0=CORD(2,I)
        R0=CORD(1,I)
        IF (DABS(R0).LT.1.D-10) GOTO 10
        SS=(1.D0+2.D0*(RE*AMI*T/R0**2))
        IF (SS.LT.0.D0) GOTO 10

        X=R0*DSQRT(SS)
C        X=0.5D0*(R0*DSQRT(SS)+X0)


        CCORD(1,I)=X
        VX=RE*AMI/X
C        VX=RE*AMI/(0.5D0*(X+X0))
        VMESH(1,I)=(CCORD(1,I)-X0)/TIME
        VMESH(2,I)=0.D0
C        IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0) THEN
C        IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0.AND.I.NE.NPT) THEN
C         VX0=GNODE(2,1,I)
         GNODE(2,1,I)=VX
C         GNODE(1,1,I)=0.5*(VX+VX0)
C        ENDIF

C I NACIN POMERANJA MREZE - NA JEDNAKIM RASTOJANJIMA
        DD=X/MESH
        DO J=I-1,I-MESH+1,-1
         X0=CCORD(1,J)
         CCORD(1,J)=DD*(MOD(J,MESH+1)-1)
         VMESH(1,J)=(CCORD(1,J)-X0)/TIME
        ENDDO

C II NACIN POMERANJA MREZE - NA RAZLICITIM RASTOJANJIMA
C        A1=(X*(1.D0-Q))/(1.D0-Q**MESH)
C        AN=A1
C        C=0.D0
C        DO J=I-MESH+1,I-1
C         X0=CCORD(1,J)
C         C=C+AN
C         CCORD(1,J)=C
C         AN=AN*Q 
C         VMESH(1,J)=(CCORD(1,J)-X0)/TIME
C        ENDDO

  10   CONTINUE
 
      END
C==========================================================================
C==========================================================================
      SUBROUTINE BINARY(GNODE,NPT,CCORD,CORD,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,7,*),CCORD(3,*),CORD(3,*)

      IF (KKORAK.EQ.1) THEN
       OPEN(51,FILE='BDIR', RECL=48,FORM='BINARY', ACCESS='DIRECT')
      ENDIF
      DO I=1,NPT
       WRITE(51,REC=I+NPT*(KKORAK-1)) 
C       WRITE(51) 
     &CCORD(1,I)-CORD(1,I),CCORD(2,I)-CORD(2,I),GNODE(2,1,I),
     &GNODE(2,2,I),GNODE(2,4,I),GNODE(2,6,I)
      ENDDO
      
      END
C==========================================================================
C==========================================================================
      SUBROUTINE ASCIIF(GNODE,NPT,CCORD,CORD,KKORAK,II,NEL,NET,NDIM,
     &VVREME,PRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /STORER/ ISINCR,IITIME,IFTIME
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA

      DIMENSION GNODE(2,7,*),CCORD(3,*),CORD(3,*),PRES(3,*)
      DIMENSION NEL(NDIM+1,*)

      IF (MOD((KKORAK-IITIME),ISINCR).NE.0.AND.ISINCR.NE.1) RETURN
      IF (KKORAK.LT.IITIME.OR.KKORAK.GT.IFTIME) RETURN


C      WRITE(II,*)'STEP= ', KKORAK

      WRITE(II,300) -1,KKORAK
      DO I=1,NPT
       IF(IAKIRA.EQ.4) THEN
       PP=DSQRT(PRES(1,I)**2+PRES(2,I)**2)
        WRITE(II,100)
     &CCORD(1,I)-CORD(1,I),CCORD(2,I)-CORD(2,I),GNODE(2,1,I),
     &GNODE(2,2,I),GNODE(2,4,I),GNODE(2,6,I),PP
       ELSE
        WRITE(II,100)
     &CCORD(1,I)-CORD(1,I),CCORD(2,I)-CORD(2,I),GNODE(2,1,I),
     &GNODE(2,2,I),GNODE(2,4,I),GNODE(2,6,I)
       ENDIF
      ENDDO
      WRITE(II,300) -1

C 100  FORMAT(5(1PE13.5),(1PE18.10))
 100  FORMAT(8(1PE13.5))
 101  FORMAT(4(1PE13.5),(1PE20.12),(1PE13.5))
 200  FORMAT(4(I10))
 300  FORMAT(2I6)

C 2000  FORMAT(4F10.7,2I5,F10.7)
C 2001  FORMAT(3F10.7,I5,3F10.7)
      
      END
C==========================================================================
C==========================================================================
      SUBROUTINE BIFUR1(GNODE,CORD,T,NPT,TIME,NSTAC,AMI,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*)

      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA
      COMMON /BIFURC/ DIAM,ALFA,BETA,VEL1,VEL2,PHASE1,PHASE2,IDIV1,
     &IDIV2,PENAL

      PI=4.D0*DATAN(1.D0)

      AALFA=ALFA*PI/180.D0
      IU=IDIV1

	PH1=PHASE1*PI/180.D0
	PH2=PHASE2*PI/180.D0

	  AN=2.D0*PI/PERIO1

      X1=CORD(1,NPT-2*(IU+1)+1)
      Y1=CORD(2,NPT-2*(IU+1)+1)
      X2=CORD(1,NPT-(IU+1))
      Y2=CORD(2,NPT-(IU+1))
      D=DSQRT((X2-X1)**2+(Y2-Y1)**2)
      RR=D/2.D0
      CX=0.5D0*(X1+X2)
      CY=0.5D0*(Y1+Y2)
      
        DO 10 I=1,IU-1
          II=NPT-2*(IU+1)+1+I
          X1=CORD(1,II)
          Y1=CORD(2,II)
          R=DSQRT((X1-CX)**2+(Y1-CY)**2)
          V1=VEL1*(1.D0-(R/RR)**2)
          V2=VEL2*(1.D0-(R/RR)**2)
          V1=V1*DSIN((AN*T)+PH1)
          V2=V2*DSIN((AN*T)+PH2)
          GNODE(2,1,II)=-V1*DSIN(AALFA)
          GNODE(2,2,II)=-V1*DCOS(AALFA)
          GNODE(2,1,II+IU+1)=V2*DSIN(AALFA)
          GNODE(2,2,II+IU+1)=-V2*DCOS(AALFA)
  10   CONTINUE  


      END
C==========================================================================
C==========================================================================

