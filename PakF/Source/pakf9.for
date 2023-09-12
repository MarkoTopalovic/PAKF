C=========================================================================
C=========================================================================
C    SUBROUTINE  ALE2DN
C                ALE2D1
C                ALE3DN
C                ALEF
C                FSI2D
C                GREZON
C                TCORDC
C=========================================================================
C=========================================================================
      SUBROUTINE ALE2D1(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,CCORD,CCORD0,VMESH0,GNOD0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine ALE2D1 is used for 2D ALE analysis
CE It is used global loop per elements
C

C SAMO PRIVREMENO UBACEN COMMON ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)
C      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,DPRES
C      COMMON /NASL/ NASLOV
C      COMMON /KONST/ GUSM,CC,AKT
C      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
C     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU
C      COMMON /TACNOS/ EPSTR,MAXIT
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /VISKOZ/ AMI,INDAMI
C      COMMON /TEMPER/ BETA,TETAO,FB2,FB3
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /VREPER/ NPER,NTABFT
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT
C      COMMON /NESTAC/ NSTAC,NKOR,DT

      CHARACTER*250NASLOV
      DIMENSION GNODE(2,6,*),ALEVO(*),DESNO(*),SILE(*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION SKEF(NDES,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
      DIMENSION SPSIL(2,*),PRES(3,*)
      DIMENSION VMESH(3,*),VMESH0(3,*),CCORD(3,*),GNOD0(3,*),CCORD0(3,*)

      DIMENSION X(9),Y(9),LM2(4*9),TT21(4*9),TT210(4*9),TT2100(4*9)
      DIMENSION X0(9),Y0(9)
      DIMENSION TT21A(18),TT21A0(18)
      DIMENSION R(3,3),S(3,3),W(3,3)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION AKVV2(9,9),AKMIV2(9,9),AKMIV3(9,9)
      DIMENSION A12(9,9),A21(9,9)
      DIMENSION A120(9,9),A210(9,9)
      DIMENSION AKK(9,9),AMV2(9,9),C(9,9)
      DIMENSION AJV2V2(9,9),AJV3V3(9,9),AJV2V3(9,9),AJV3V2(9,9)
      DIMENSION AKTV2(9,9),AKTV3(9,9)
      DIMENSION AKV2P(9,4),AKV3P(9,4),AKV2P1(9,4)
      DIMENSION AJK(9,9)
      DIMENSION RS2(9),RS3(9)
      DIMENSION RB2(9),RB3(9)
      DIMENSION F31(31),FALE(31)
      DIMENSION AKUAX(9,9),AKVAX(9,9),APAXIX(9,9),APAXIY(9,9)

      DIMENSION AKMV20(9,9),AKVV20(9,9),AKVP10(9,4),AMV20(9,9)
      DIMENSION AKMV30(9,9),AKV2P0(9,4),AKV3P0(9,4),C0(9,9)
      DIMENSION AKTV20(9,9),AKTV30(9,9),AKK0(9,9)
      DIMENSION AKUAX0(9,9),AKVAX0(9,9),APXIX0(9,9),APXIY0(9,9)


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


C GLAVNA PETLJA PO ELEMENTIMA
C GLOBAL LOOP PER ELEMENTS
C NBREL is counter of elements

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,4*NDIM
      TT210(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0
C=========================================================================
      DO 130 KLM=1,NDIM
      X(KLM)=CCORD(1,NEL(KLM,NBREL))
      Y(KLM)=CCORD(2,NEL(KLM,NBREL))
      X0(KLM)=CCORD0(1,NEL(KLM,NBREL))
      Y0(KLM)=CCORD0(2,NEL(KLM,NBREL))
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
       TT21A0(I)=VMESH0(1,NODE)
       TT21A0(I+NDIM)=VMESH0(2,NODE)
      ENDDO
C=======================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,2
C      IF (ID(NR,NEL(KLM,NBREL)) .NE. 0) THEN
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
      TT2100(KLM+(NR-1)*NDIM)=GNOD0(NR,NEL(KLM,NBREL))
C      ENDIF
 135  CONTINUE
      IF (KLM.LE.4) THEN
       TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
       TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
      TT21(KLM+2*NDIM+4)=GNODE(2,5,NEL(KLM,NBREL))
      TT210(KLM+2*NDIM+4)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE
C
C=======================================================================
C      IF (NDIM.GT.4) THEN
C      DO 145 NR=28,36
C      TT210(NR-5)=TT210(NR)
C 145  TT21(NR-5)=TT21(NR)
C      ENDIF
C=======================================================================
      DO I=1,3*NDIM+4
       IF(KKORAK.EQ.1.AND.ITER.EQ.1) TT210(I)=TT21(I)
      ENDDO
C=======================================================================
C NUMZAD is number of prescribed values
C
C      DO 160 NR=1,NDIM
C      DO 155 NZDT=1,NUMZAD
C
C TIMFUN is subroutine which determined values of function at currently 
C time step (at the end of step)
C
C      IF(NEL(NR,NBREL).EQ.NZAD(1,NZDT)) THEN
C       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
C     &NTABFT,IIZLAZ)
C========================================
C SAMO PRIVREMENO UBACENO ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C       FK1=PBALL
C=========================================
C      MESTO=NZAD(2,NZDT)
C        IF (MESTO.EQ.4) THEN
C        TT21(2*NDIM+4+NR)=ZADVRE(NZDT)*FK1
C        ELSE
C        TT21((MESTO-1)*NDIM+NR)=ZADVRE(NZDT)*FK1
C        ENDIF
C      ENDIF  
C 155  CONTINUE
C 160  CONTINUE
C=======================================================================

C=======================================================================
C      CALL WRR(TT21,3*NDIM+4,'T211')
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
      IF (N.LE.4) THEN
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
       AKV2P1(K,N)=0.D0
       AKV2P0(K,N)=0.D0
       AKV3P0(K,N)=0.D0
       AKVP10(K,N)=0.D0
      ENDIF
      AKVV2(K,N)=0.D0
      AKVV20(K,N)=0.D0
      APAXIX(K,N)=0.D0
      APAXIY(K,N)=0.D0
      AKUAX(K,N)=0.D0
      AKVAX(K,N)=0.D0
      APXIX0(K,N)=0.D0
      APXIY0(K,N)=0.D0
      AKUAX0(K,N)=0.D0
      AKVAX0(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKTV20(K,N)=0.D0
      AKTV30(K,N)=0.D0
      AKMIV2(K,N)=0.D0
      AKMIV3(K,N)=0.D0
      AKMV20(K,N)=0.D0
      AKMV30(K,N)=0.D0
      A12(K,N)=0.D0
      A21(K,N)=0.D0
      A120(K,N)=0.D0
      A210(K,N)=0.D0
      AMV2(K,N)=0.D0
      AMV20(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V2(K,N)=0.D0
      AKK(K,N)=0.D0
      AKK0(K,N)=0.D0
      AJK(K,N)=0.D0
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
       RS2(K)=0.D0
       RS3(K)=0.D0
       RB2(K)=0.D0
       RB3(K)=0.D0
  163 CONTINUE

       CALL CLEAR (FALE,31)
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
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ
       
       IF (NDIM.EQ.4) IPN=1
       IF (NDIM.EQ.9) IPN=4

       CALL FSI2D(FALE,ZVHX,ZVHY,NEL,NDIM,GNODE,VMESH,VMESH0,TT21,H,IPN,
     &2,GUSM,WDT,HP,NBREL,TIME,TT2100,TT210)

C       CALL ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
C       CALL AXISYF(INDAX,DEBLJP,X,HP,4)

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
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
  164 CONTINUE
C       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
C       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(0.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(0.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE

C===========================================================================
      DO 181 I=1,IBRGT
      DO 171 J=1,IBRGT

C subroutine FNTERP return to us interpolation functions...

      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0,TT210,H,HP,ZVHX,ZVHY,HV2,HV3,
     &ZVXT,ZVYT,DETJ,DETJS,X0,Y0,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ
       CALL ALEHV(TT210,TT21A0,HV2,HV3,NDIM,H)
       CALL AXISYF(INDAX,DEBLJ,X0,H,NDIM)


      IF (INDAMI.EQ.1) CALL NENJUT(ZVHX,ZVHY,TT210)
      DO  K=1,NDIM
      DO  N=1,NDIM
      IF (INDAX.EQ.1) THEN
      AKUAX0(K,N)=AKUAX0(K,N)
     &-WDT*H(K)*(ZVHX(N)/DEBLJ-H(N)/(DEBLJ**2))*AMI
      AKVAX0(K,N)=AKVAX0(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AMI
      ENDIF
      AKVV20(K,N)=AKVV20(K,N)+WDT*((H(K)*HV2*ZVHX(N)+
     1H(K)*HV3*ZVHY(N))*GUSM)
      AKTV20(K,N)=AKTV20(K,N)+WDT*(H(K)*H(N)*ZVXT
     1*GUSM*CC)
      AKTV30(K,N)=AKTV30(K,N)+WDT*(H(K)*H(N)*ZVYT
     1*GUSM*CC)
      AKMV20(K,N)=AKMV20(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI)
      AKK0(K,N)=AKK0(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AKT)
      IF (INDAX.EQ.1) THEN
       AKK0(K,N)=AKK0(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AKT
      ENDIF

      AKMV30(K,N)=AKMV20(K,N)
      AMV20(K,N)=AMV20(K,N)+WDT*H(K)*H(N)*GUSM
C      AJV2V2(K,N)=AJV2V2(K,N)+WDT*(H(K)*ZVXV2*H(N)*GUSM)
C      AJV3V3(K,N)=AJV3V3(K,N)+WDT*(H(K)*ZVYV3*H(N)*GUSM)
C      AJV2V3(K,N)=AJV2V3(K,N)+WDT*(H(K)*ZVYV2*H(N)*GUSM)
C      AJV3V2(K,N)=AJV3V2(K,N)+WDT*(H(K)*ZVXV3*H(N)*GUSM)
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P0(K,N)=AKV2P0(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P0(K,N)=AKV3P0(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKVP10(K,N)=AKVP10(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
       ENDDO  
       ENDDO  

  171 CONTINUE
  181 CONTINUE
C==========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PRITISKOM
C reduced integration for pressure
C==========================================================================
C      DO  I=1,IBRGT-1
C      DO  J=1,IBRGT-1
C      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
C     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
C     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX)
C       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
C      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
C
C      DO  K=1,NDIM
C       DO  N=1,NDIM
C        IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
C         AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
C         AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
C          IF (INDAX.EQ.1) THEN
C           AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
C          ENDIF
C        ENDIF
C       ENDDO
C      ENDDO
C
C      ENDDO
C      ENDDO
C==========================================================================
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
C      WDT1=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ*DEBLJ
      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
      DO 190 K=1,NDIM
      DO 185 N=1,NDIM
      AKMIV2(K,N)=AKMIV2(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT
      AKMIV3(K,N)=AKMIV3(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT
      A12(K,N)=A12(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
C      AJV2V3(K,N)=AJV2V3(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      A21(K,N)=A21(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
C      AJV3V2(K,N)=AJV3V2(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
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


      DO 201 I=1,IBRGT-1
      DO 196 J=1,IBRGT-1
      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0,TT210,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X0,Y0,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,
     &TAU,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,
     &INDAX,AKT,GUSM,IUPWIN,NBREL,IIZLAZ)

      CALL AXISYF(INDAX,DEBLJ,X0,H,NDIM)
      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
      DO 191 K=1,NDIM
      DO 186 N=1,NDIM
      AKMV20(K,N)=AKMV20(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT
      AKMV30(K,N)=AKMV30(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT
      A120(K,N)=A120(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      A210(K,N)=A210(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
       IF (INDAX.EQ.1) THEN
       APXIX0(K,N)=APXIX0(K,N)+PENALT*WDT*(ZVHX(K)*H(N)/DEBLJ)
       APXIY0(K,N)=APXIY0(K,N)+PENALT*WDT*(ZVHY(K)*H(N)/DEBLJ)
       ENDIF
      AKV2P0(K,N)=AKV2P0(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P0(K,N)=AKV3P0(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKVP10(K,N)=AKVP10(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
  186 CONTINUE
  191 CONTINUE
  196 CONTINUE
  201 CONTINUE
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
      DO 206 K=1,NDIM
       RS2(K)=RS2(K)+WDT*(H(K)*FS2)*NGPSIL(4,JBRPS)
       RS3(K)=RS3(K)+WDT*(H(K)*FS3)*NGPSIL(5,JBRPS)
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
         C0(I,J)=AMV20(I,J)*CC              
         AJK(I,J)=AKVV2(I,J)*CC
 254   CONTINUE
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
C=========================================================================
C PAKOVANJE MATRICE C,Kk,Jk U MATRICU SKEF
C SKEF is local stiffness matrix

      DO 263 I=1,NDIM
      DO 262 J=1,NDIM
      SKEF(I,J)=AKMIV2(I,J)+AKVV2(I,J)+AJV2V2(I,J)
      SKEF(I+NDIM,J)=AJV3V2(I,J)
      SKEF(I+NDIM,J+NDIM)=AKMIV3(I,J)+AKVV2(I,J)+AJV3V3(I,J)
      SKEF(I,J+NDIM)=AJV2V3(I,J)
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=AKK(I,J)+AJK(I,J)
      IF (NSTAC.EQ.0) THEN
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=SKEF(I+2*NDIM+4,J+2*NDIM+4)+
     &C(I,J)/TIME
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       ENDIF
      IF (J.LE.4) THEN
       SKEF(J+2*NDIM,I)=AKV2P(I,J)+AKV2P1(I,J)
       SKEF(J+2*NDIM,I+NDIM)=AKV3P(I,J)
       SKEF(I,J+2*NDIM)=AKV2P(I,J)
       SKEF(I+NDIM,J+2*NDIM)=AKV3P(I,J)
      ENDIF
       SKEF(I+2*NDIM+4,J)=AKTV2(I,J)
       SKEF(I+2*NDIM+4,J+NDIM)=AKTV3(I,J)
C LEVA STRANA USLED KONVEKCIJE
C left side which relate to convective term
      SKEF(I,J+2*NDIM+4)=SKEF(I,J+2*NDIM+4)+AMV2(I,J)*FB2*BETA
      SKEF(I+NDIM,J+2*NDIM+4)=SKEF(I+NDIM,J+2*NDIM+4)+AMV2(I,J)*FB3*BETA
      IF (INDAX.EQ.1) THEN
       SKEF(I,J)=SKEF(I,J)+AKUAX(I,J)
        IF (PENALT.GT.0.D0) THEN
         SKEF(I,J)=SKEF(I,J)+APAXIX(I,J)
         SKEF(I+NDIM,J)=SKEF(I+NDIM,J)+APAXIY(I,J)
        ENDIF
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AKVAX(I,J)
      ENDIF
  262 CONTINUE
  263 CONTINUE
C=========================================================================



      DO 290 I=1,NDIM
C ZAPREMINSKE SILE
C body forces
       F31(I)=F31(I)+RB2(I)
       F31(I+NDIM)=F31(I+NDIM)+RB3(I)

C POVRSINSKE SILE
C surface forces
C      IF (INDX.EQ.1)
      F31(I)=F31(I)+RS2(I)
C      IF (INDY.EQ.1)
      F31(I+NDIM)=F31(I+NDIM)+RS3(I)

      DO 285 J=1,NDIM
      IF (INDAX.EQ.1) THEN
       F31(I)=F31(I)-AKUAX(I,J)*TT21(J)
       F31(I+NDIM)=F31(I+NDIM)-AKVAX(I,J)*TT21(J+NDIM)
       IF (PENALT.GT.0.D0) THEN
        F31(I)=F31(I)-APAXIX(I,J)*TT21(J)
        F31(I+NDIM)=F31(I+NDIM)-APAXIY(I,J)*TT21(J)
       ENDIF
      ENDIF
      F31(I)=F31(I)-(AKMIV2(I,J)+AKVV2(I,J))*TT21(J)
      IF (J.LE.4) THEN
       F31(I)=F31(I)-AKV2P(I,J)*TT21(J+2*NDIM)
       F31(I+NDIM)=F31(I+NDIM)-AKV3P(I,J)*TT21(J+2*NDIM)
      ENDIF
      IF (I.LE.4) THEN
      F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2P(J,I)*TT21(J)-
     1AKV3P(J,I)*TT21(J+NDIM)-AKV2P1(J,I)*TT21(J)
      ENDIF


      IF (NSTAC.EQ.0) THEN
       F31(I)=F31(I)-AMV2(I,J)*(TT21(J)-TT210(J))/TIME
       F31(I+NDIM)=F31(I+NDIM)-
     &AMV2(I,J)*(TT21(J+NDIM)-TT210(J+NDIM))/TIME
       II=I+2*NDIM+4
       JJ=J+2*NDIM+4
       F31(II)=F31(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF

      F31(I)=F31(I)-A12(I,J)*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-(AKMIV3(I,J)+AKVV2(I,J))*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-A21(I,J)*TT21(J)
      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)-AKTV2(I,J)*TT21(J)
     1-AKTV3(I,J)*TT21(J+NDIM)-AKK(I,J)*TT21(J+2*NDIM+4)
C DESNA STRANA USLED UTICAJA KONVEKCIJE
c right-hand side created by natural convection
      F31(I)=F31(I)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB2
      F31(I+NDIM)=F31(I+NDIM)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB3
 285  CONTINUE							      
 290  CONTINUE

C=========================================================================
C==========================================================================



C OVDE PRIVREMENO UBACENO DA SE PRESKACE OVAJ DEO
C 23 Febr. 2004
      goto 292


      DO 291 I=1,NDIM
      DO 286 J=1,NDIM
      IF (INDAX.EQ.1) THEN
       F31(I)=F31(I)-AKUAX0(I,J)*TT210(J)
       F31(I+NDIM)=F31(I+NDIM)-AKVAX0(I,J)*TT210(J+NDIM)
       IF (PENALT.GT.0.D0) THEN
        F31(I)=F31(I)-APXIX0(I,J)*TT210(J)
        F31(I+NDIM)=F31(I+NDIM)-APXIY0(I,J)*TT210(J)
       ENDIF
      ENDIF
      F31(I)=F31(I)-(AKMV20(I,J)+AKVV20(I,J))*TT210(J)
      IF (J.LE.4) THEN
       F31(I)=F31(I)-AKV2P0(I,J)*TT210(J+2*NDIM)
       F31(I+NDIM)=F31(I+NDIM)-AKV3P0(I,J)*TT210(J+2*NDIM)
      ENDIF
      IF (I.LE.4) THEN
      F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2P0(J,I)*TT210(J)-
     1AKV3P0(J,I)*TT210(J+NDIM)-AKVP10(J,I)*TT210(J)
      ENDIF


      IF (NSTAC.EQ.0) THEN
       F31(I)=F31(I)-AMV20(I,J)*(TT210(J)-TT2100(J))/TIME
       F31(I+NDIM)=F31(I+NDIM)-
     &AMV20(I,J)*(TT210(J+NDIM)-TT2100(J+NDIM))/TIME
       II=I+2*NDIM+4
       JJ=J+2*NDIM+4
       F31(II)=F31(II)-C0(I,J)*(TT210(JJ)-TT2100(JJ))/TIME
      ENDIF

      F31(I)=F31(I)-A120(I,J)*TT210(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-(AKMV30(I,J)+AKVV20(I,J))*TT210(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-A210(I,J)*TT210(J)
      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)-AKTV20(I,J)*TT210(J)
     1-AKTV30(I,J)*TT210(J+NDIM)-AKK0(I,J)*TT210(J+2*NDIM+4)
C DESNA STRANA USLED UTICAJA KONVEKCIJE
c right-hand side created by natural convection
      F31(I)=F31(I)-AMV20(I,J)*TT210(J+2*NDIM+4)*BETA*FB2
      F31(I+NDIM)=F31(I+NDIM)-AMV20(I,J)*TT210(J+2*NDIM+4)*BETA*FB3
 286  CONTINUE							      
 291  CONTINUE

 292  CONTINUE
C=========================================================================
C==========================================================================

C==========================================================================

C==========================================================================
C==========================================================================



C============================================================================
C  FOR FLUX CALCULATION
C      IF(IPROL.EQ.1) THEN
       DO I=1,NDIM
        N=NEL(I,NBREL)
C        JJ=ID(5,N)
C        IF(JJ.NE.0) THEN
         POT=0.D0
         DO J=1,NDES
C         DO J=1+2*NDIM+4,NDES
          POT=POT+SKEF(I+2*NDIM+4,J)*TT21(J)
         ENDDO
         PRES(1,N)=PRES(1,N)+POT
C        ENDIF
       ENDDO
C      ENDIF
C============================================================================


C      CALL WRRF(SKEF,NDES*NDES,'SKEF ',IIZLAZ)      
C      CALL WRRF(F31,NDES,'F31= ',IIZLAZ)      
C       CALL SPAKDE (SILE,F31,LM2,NDES)

      DO I=1,NDES
       F31(I)=F31(I)-FALE(I)     
      ENDDO 

      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F31,MAXA,LM2,NDES,1)

C RACUNANJE SILA INTERAKCIJE KOJIMA FLUID DELUJE NA ZIDOVE
c calculation of interaction forces between fluid and solid
	 IF (NSTAC.EQ.0) THEN
        DO I=1,NDIM
          II=I+NDIM
         DO J=1,NDIM
          JJ=J+NDIM
          SKEF(I,J)=SKEF(I,J)-AMV2(I,J)*TT210(J)/TIME
          SKEF(II,JJ)=SKEF(II,JJ)-AMV2(I,J)*TT210(JJ)/TIME
         ENDDO
        ENDDO
       ENDIF

       FACAXY=1.D0
       IF(INDAX.EQ.1) FACAXY=8.D0*DATAN(1.D0)

       DO I=1,NDIM
       N=NEL(I,NBREL)
        DO J=1,NDES
         P1=SKEF(I,J)*TT21(J)
         P2=SKEF(I+NDIM,J)*TT21(J)
         SPSIL(1,N)=SPSIL(1,N)-P1/FACAXY
         SPSIL(2,N)=SPSIL(2,N)-P2/FACAXY
       ENDDO
      ENDDO



C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C end of loop per elements
C=======================================================================
 400  CONTINUE

C End of subroutine
      END
C=======================================================================
C=========================================================================
C=========================================================================
      SUBROUTINE ALE2DN(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,CCORD,CCORD0,VMESH0,GNOD0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine ALE2DN is used for 2D ALE analysis
CE It is used global loop per elements
C

C SAMO PRIVREMENO UBACEN COMMON ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)
C      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,DPRES
C      COMMON /NASL/ NASLOV
C      COMMON /KONST/ GUSM,CC,AKT
C      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
C     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU
C      COMMON /TACNOS/ EPSTR,MAXIT
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /VISKOZ/ AMI,INDAMI
C      COMMON /TEMPER/ BETA,TETAO,FB2,FB3
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /VREPER/ NPER,NTABFT
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT
C      COMMON /NESTAC/ NSTAC,NKOR,DT

      CHARACTER*250NASLOV
      DIMENSION GNODE(2,6,*),ALEVO(*),DESNO(*),SILE(*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION SKEF(NDES,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
      DIMENSION SPSIL(2,*),PRES(3,*)
      DIMENSION VMESH(3,*),VMESH0(3,*),CCORD(3,*),GNOD0(3,*),CCORD0(3,*)

      DIMENSION X(9),Y(9),LM2(4*9),TT21(4*9),TT210(4*9),TT2100(4*9)
      DIMENSION X0(9),Y0(9)
      DIMENSION TT21A(18),TT21A0(18)
      DIMENSION R(3,3),S(3,3),W(3,3)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION AKVV2(9,9),AKMIV2(9,9),AKMIV3(9,9)
      DIMENSION A12(9,9),A21(9,9)
      DIMENSION A120(9,9),A210(9,9)
      DIMENSION AKK(9,9),AMV2(9,9),C(9,9)
      DIMENSION AJV2V2(9,9),AJV3V3(9,9),AJV2V3(9,9),AJV3V2(9,9)
      DIMENSION AKTV2(9,9),AKTV3(9,9)
      DIMENSION AKV2P(9,4),AKV3P(9,4),AKV2P1(9,4)
      DIMENSION AJK(9,9)
      DIMENSION RS2(9),RS3(9)
      DIMENSION RB2(9),RB3(9)
      DIMENSION F31(31),FALE(31)
      DIMENSION AKUAX(9,9),AKVAX(9,9),APAXIX(9,9),APAXIY(9,9)

      DIMENSION AKMV20(9,9),AKVV20(9,9),AKVP10(9,4),AMV20(9,9)
      DIMENSION AKMV30(9,9),AKV2P0(9,4),AKV3P0(9,4),C0(9,9)
      DIMENSION AKTV20(9,9),AKTV30(9,9),AKK0(9,9)
      DIMENSION AKUAX0(9,9),AKVAX0(9,9),APXIX0(9,9),APXIY0(9,9)


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


C GLAVNA PETLJA PO ELEMENTIMA
C GLOBAL LOOP PER ELEMENTS
C NBREL is counter of elements

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,4*NDIM
      TT210(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0
C=========================================================================
      DO 130 KLM=1,NDIM
      X(KLM)=CCORD(1,NEL(KLM,NBREL))
      Y(KLM)=CCORD(2,NEL(KLM,NBREL))
      X0(KLM)=CCORD0(1,NEL(KLM,NBREL))
      Y0(KLM)=CCORD0(2,NEL(KLM,NBREL))
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
       TT21A0(I)=VMESH0(1,NODE)
       TT21A0(I+NDIM)=VMESH0(2,NODE)
      ENDDO
C=======================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,2
C      IF (ID(NR,NEL(KLM,NBREL)) .NE. 0) THEN
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
      TT2100(KLM+(NR-1)*NDIM)=GNOD0(NR,NEL(KLM,NBREL))
C      ENDIF
 135  CONTINUE
      IF (KLM.LE.4) THEN
       TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
       TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
      TT21(KLM+2*NDIM+4)=GNODE(2,5,NEL(KLM,NBREL))
      TT210(KLM+2*NDIM+4)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE
C
C=======================================================================
C      IF (NDIM.GT.4) THEN
C      DO 145 NR=28,36
C      TT210(NR-5)=TT210(NR)
C 145  TT21(NR-5)=TT21(NR)
C      ENDIF
C=======================================================================
      DO I=1,3*NDIM+4
       IF(KKORAK.EQ.1.AND.ITER.EQ.1) TT210(I)=TT21(I)
      ENDDO
C=======================================================================
C NUMZAD is number of prescribed values
C
C      DO 160 NR=1,NDIM
C      DO 155 NZDT=1,NUMZAD
C
C TIMFUN is subroutine which determined values of function at currently 
C time step (at the end of step)
C
C      IF(NEL(NR,NBREL).EQ.NZAD(1,NZDT)) THEN
C       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
C     &NTABFT,IIZLAZ)
C========================================
C SAMO PRIVREMENO UBACENO ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C       FK1=PBALL
C=========================================
C      MESTO=NZAD(2,NZDT)
C        IF (MESTO.EQ.4) THEN
C        TT21(2*NDIM+4+NR)=ZADVRE(NZDT)*FK1
C        ELSE
C        TT21((MESTO-1)*NDIM+NR)=ZADVRE(NZDT)*FK1
C        ENDIF
C      ENDIF  
C 155  CONTINUE
C 160  CONTINUE
C=======================================================================

C=======================================================================
C      CALL WRR(TT21,3*NDIM+4,'T211')
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
      IF (N.LE.4) THEN
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
       AKV2P1(K,N)=0.D0
       AKV2P0(K,N)=0.D0
       AKV3P0(K,N)=0.D0
       AKVP10(K,N)=0.D0
      ENDIF
      AKVV2(K,N)=0.D0
      AKVV20(K,N)=0.D0
      APAXIX(K,N)=0.D0
      APAXIY(K,N)=0.D0
      AKUAX(K,N)=0.D0
      AKVAX(K,N)=0.D0
      APXIX0(K,N)=0.D0
      APXIY0(K,N)=0.D0
      AKUAX0(K,N)=0.D0
      AKVAX0(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKTV20(K,N)=0.D0
      AKTV30(K,N)=0.D0
      AKMIV2(K,N)=0.D0
      AKMIV3(K,N)=0.D0
      AKMV20(K,N)=0.D0
      AKMV30(K,N)=0.D0
      A12(K,N)=0.D0
      A21(K,N)=0.D0
      A120(K,N)=0.D0
      A210(K,N)=0.D0
      AMV2(K,N)=0.D0
      AMV20(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V2(K,N)=0.D0
      AKK(K,N)=0.D0
      AKK0(K,N)=0.D0
      AJK(K,N)=0.D0
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
       RS2(K)=0.D0
       RS3(K)=0.D0
       RB2(K)=0.D0
       RB3(K)=0.D0
  163 CONTINUE

       CALL CLEAR (FALE,31)
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
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ
       
       IF (NDIM.EQ.4) IPN=1
       IF (NDIM.EQ.9) IPN=4

       CALL FSI2D(FALE,ZVHX,ZVHY,NEL,NDIM,GNODE,VMESH,VMESH0,TT21,H,IPN,
     &2,GUSM,WDT,HP,NBREL,TIME,TT2100,TT210)

       CALL ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
C       CALL AXISYF(INDAX,DEBLJP,X,HP,4)

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
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
  164 CONTINUE
C       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
C       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(0.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(0.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE

C===========================================================================
      DO 181 I=1,IBRGT
      DO 171 J=1,IBRGT

C subroutine FNTERP return to us interpolation functions...

      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0,TT210,H,HP,ZVHX,ZVHY,HV2,HV3,
     &ZVXT,ZVYT,DETJ,DETJS,X0,Y0,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ
       CALL ALEHV(TT210,TT21A0,HV2,HV3,NDIM,H)
       CALL AXISYF(INDAX,DEBLJ,X0,H,NDIM)


      IF (INDAMI.EQ.1) CALL NENJUT(ZVHX,ZVHY,TT210)
      DO  K=1,NDIM
      DO  N=1,NDIM
      IF (INDAX.EQ.1) THEN
      AKUAX0(K,N)=AKUAX0(K,N)
     &-WDT*H(K)*(ZVHX(N)/DEBLJ-H(N)/(DEBLJ**2))*AMI
      AKVAX0(K,N)=AKVAX0(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AMI
      ENDIF
      AKVV20(K,N)=AKVV20(K,N)+WDT*((H(K)*HV2*ZVHX(N)+
     1H(K)*HV3*ZVHY(N))*GUSM)
      AKTV20(K,N)=AKTV20(K,N)+WDT*(H(K)*H(N)*ZVXT
     1*GUSM*CC)
      AKTV30(K,N)=AKTV30(K,N)+WDT*(H(K)*H(N)*ZVYT
     1*GUSM*CC)
      AKMV20(K,N)=AKMV20(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI)
      AKK0(K,N)=AKK0(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AKT)
      IF (INDAX.EQ.1) THEN
       AKK0(K,N)=AKK0(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AKT
      ENDIF

      AKMV30(K,N)=AKMV20(K,N)
      AMV20(K,N)=AMV20(K,N)+WDT*H(K)*H(N)*GUSM
C      AJV2V2(K,N)=AJV2V2(K,N)+WDT*(H(K)*ZVXV2*H(N)*GUSM)
C      AJV3V3(K,N)=AJV3V3(K,N)+WDT*(H(K)*ZVYV3*H(N)*GUSM)
C      AJV2V3(K,N)=AJV2V3(K,N)+WDT*(H(K)*ZVYV2*H(N)*GUSM)
C      AJV3V2(K,N)=AJV3V2(K,N)+WDT*(H(K)*ZVXV3*H(N)*GUSM)
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P0(K,N)=AKV2P0(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P0(K,N)=AKV3P0(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKVP10(K,N)=AKVP10(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
       ENDDO  
       ENDDO  

  171 CONTINUE
  181 CONTINUE
C==========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PRITISKOM
C reduced integration for pressure
C==========================================================================
C      DO  I=1,IBRGT-1
C      DO  J=1,IBRGT-1
C      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
C     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
C     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX)
C       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
C      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
C
C      DO  K=1,NDIM
C       DO  N=1,NDIM
C        IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
C         AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
C         AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
C          IF (INDAX.EQ.1) THEN
C           AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
C          ENDIF
C        ENDIF
C       ENDDO
C      ENDDO
C
C      ENDDO
C      ENDDO
C==========================================================================
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
C      WDT1=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ*DEBLJ
      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
      DO 190 K=1,NDIM
      DO 185 N=1,NDIM
      AKMIV2(K,N)=AKMIV2(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT
      AKMIV3(K,N)=AKMIV3(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT
      A12(K,N)=A12(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
C      AJV2V3(K,N)=AJV2V3(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      A21(K,N)=A21(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
C      AJV3V2(K,N)=AJV3V2(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
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


      DO 201 I=1,IBRGT-1
      DO 196 J=1,IBRGT-1
      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0,TT210,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X0,Y0,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,
     &TAU,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,
     &INDAX,AKT,GUSM,IUPWIN,NBREL,IIZLAZ)

      CALL AXISYF(INDAX,DEBLJ,X0,H,NDIM)
      WDT=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
      DO 191 K=1,NDIM
      DO 186 N=1,NDIM
      AKMV20(K,N)=AKMV20(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT
      AKMV30(K,N)=AKMV30(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT
      A120(K,N)=A120(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      A210(K,N)=A210(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
       IF (INDAX.EQ.1) THEN
       APXIX0(K,N)=APXIX0(K,N)+PENALT*WDT*(ZVHX(K)*H(N)/DEBLJ)
       APXIY0(K,N)=APXIY0(K,N)+PENALT*WDT*(ZVHY(K)*H(N)/DEBLJ)
       ENDIF
      AKV2P0(K,N)=AKV2P0(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P0(K,N)=AKV3P0(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKVP10(K,N)=AKVP10(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
  186 CONTINUE
  191 CONTINUE
  196 CONTINUE
  201 CONTINUE
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
      DO 206 K=1,NDIM
       RS2(K)=RS2(K)+WDT*(H(K)*FS2)*NGPSIL(4,JBRPS)
       RS3(K)=RS3(K)+WDT*(H(K)*FS3)*NGPSIL(5,JBRPS)
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
         C0(I,J)=AMV20(I,J)*CC              
         AJK(I,J)=AKVV2(I,J)*CC
 254   CONTINUE
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
C=========================================================================
C PAKOVANJE MATRICE C,Kk,Jk U MATRICU SKEF
C SKEF is local stiffness matrix

      DO 263 I=1,NDIM
      DO 262 J=1,NDIM
      SKEF(I,J)=AKMIV2(I,J)+AKVV2(I,J)+AJV2V2(I,J)
      SKEF(I+NDIM,J)=AJV3V2(I,J)
      SKEF(I+NDIM,J+NDIM)=AKMIV3(I,J)+AKVV2(I,J)+AJV3V3(I,J)
      SKEF(I,J+NDIM)=AJV2V3(I,J)
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=AKK(I,J)+AJK(I,J)
      IF (NSTAC.EQ.0) THEN
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=SKEF(I+2*NDIM+4,J+2*NDIM+4)+
     &C(I,J)/TIME
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       ENDIF
      IF (J.LE.4) THEN
       SKEF(J+2*NDIM,I)=AKV2P(I,J)+AKV2P1(I,J)
       SKEF(J+2*NDIM,I+NDIM)=AKV3P(I,J)
       SKEF(I,J+2*NDIM)=AKV2P(I,J)
       SKEF(I+NDIM,J+2*NDIM)=AKV3P(I,J)
      ENDIF
       SKEF(I+2*NDIM+4,J)=AKTV2(I,J)
       SKEF(I+2*NDIM+4,J+NDIM)=AKTV3(I,J)
C LEVA STRANA USLED KONVEKCIJE
C left side which relate to convective term
      SKEF(I,J+2*NDIM+4)=SKEF(I,J+2*NDIM+4)+AMV2(I,J)*FB2*BETA
      SKEF(I+NDIM,J+2*NDIM+4)=SKEF(I+NDIM,J+2*NDIM+4)+AMV2(I,J)*FB3*BETA
      IF (INDAX.EQ.1) THEN
       SKEF(I,J)=SKEF(I,J)+AKUAX(I,J)
        IF (PENALT.GT.0.D0) THEN
         SKEF(I,J)=SKEF(I,J)+APAXIX(I,J)
         SKEF(I+NDIM,J)=SKEF(I+NDIM,J)+APAXIY(I,J)
        ENDIF
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AKVAX(I,J)
      ENDIF
  262 CONTINUE
  263 CONTINUE
C=========================================================================



      DO 290 I=1,NDIM
C ZAPREMINSKE SILE
C body forces
       F31(I)=F31(I)+RB2(I)
       F31(I+NDIM)=F31(I+NDIM)+RB3(I)

C POVRSINSKE SILE
C surface forces
C      IF (INDX.EQ.1)
      F31(I)=F31(I)+RS2(I)
C      IF (INDY.EQ.1)
      F31(I+NDIM)=F31(I+NDIM)+RS3(I)

      DO 285 J=1,NDIM
      IF (INDAX.EQ.1) THEN
       F31(I)=F31(I)-AKUAX(I,J)*TT21(J)
       F31(I+NDIM)=F31(I+NDIM)-AKVAX(I,J)*TT21(J+NDIM)
       IF (PENALT.GT.0.D0) THEN
        F31(I)=F31(I)-APAXIX(I,J)*TT21(J)
        F31(I+NDIM)=F31(I+NDIM)-APAXIY(I,J)*TT21(J)
       ENDIF
      ENDIF
      F31(I)=F31(I)-(AKMIV2(I,J)+AKVV2(I,J))*TT21(J)
      IF (J.LE.4) THEN
       F31(I)=F31(I)-AKV2P(I,J)*TT21(J+2*NDIM)
       F31(I+NDIM)=F31(I+NDIM)-AKV3P(I,J)*TT21(J+2*NDIM)
      ENDIF
      IF (I.LE.4) THEN
      F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2P(J,I)*TT21(J)-
     1AKV3P(J,I)*TT21(J+NDIM)-AKV2P1(J,I)*TT21(J)
      ENDIF


      IF (NSTAC.EQ.0) THEN
       F31(I)=F31(I)-AMV2(I,J)*(TT21(J)-TT210(J))/TIME
       F31(I+NDIM)=F31(I+NDIM)-
     &AMV2(I,J)*(TT21(J+NDIM)-TT210(J+NDIM))/TIME
       II=I+2*NDIM+4
       JJ=J+2*NDIM+4
       F31(II)=F31(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF

      F31(I)=F31(I)-A12(I,J)*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-(AKMIV3(I,J)+AKVV2(I,J))*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-A21(I,J)*TT21(J)
      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)-AKTV2(I,J)*TT21(J)
     1-AKTV3(I,J)*TT21(J+NDIM)-AKK(I,J)*TT21(J+2*NDIM+4)
C DESNA STRANA USLED UTICAJA KONVEKCIJE
c right-hand side created by natural convection
      F31(I)=F31(I)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB2
      F31(I+NDIM)=F31(I+NDIM)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB3
 285  CONTINUE							      
 290  CONTINUE

C=========================================================================
C==========================================================================






      DO 291 I=1,NDIM
      DO 286 J=1,NDIM
      IF (INDAX.EQ.1) THEN
       F31(I)=F31(I)-AKUAX0(I,J)*TT210(J)
       F31(I+NDIM)=F31(I+NDIM)-AKVAX0(I,J)*TT210(J+NDIM)
       IF (PENALT.GT.0.D0) THEN
        F31(I)=F31(I)-APXIX0(I,J)*TT210(J)
        F31(I+NDIM)=F31(I+NDIM)-APXIY0(I,J)*TT210(J)
       ENDIF
      ENDIF
      F31(I)=F31(I)-(AKMV20(I,J)+AKVV20(I,J))*TT210(J)
      IF (J.LE.4) THEN
       F31(I)=F31(I)-AKV2P0(I,J)*TT210(J+2*NDIM)
       F31(I+NDIM)=F31(I+NDIM)-AKV3P0(I,J)*TT210(J+2*NDIM)
      ENDIF
      IF (I.LE.4) THEN
      F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2P0(J,I)*TT210(J)-
     1AKV3P0(J,I)*TT210(J+NDIM)-AKVP10(J,I)*TT210(J)
      ENDIF


      IF (NSTAC.EQ.0) THEN
       F31(I)=F31(I)-AMV20(I,J)*(TT210(J)-TT2100(J))/TIME
       F31(I+NDIM)=F31(I+NDIM)-
     &AMV20(I,J)*(TT210(J+NDIM)-TT2100(J+NDIM))/TIME
       II=I+2*NDIM+4
       JJ=J+2*NDIM+4
       F31(II)=F31(II)-C0(I,J)*(TT210(JJ)-TT2100(JJ))/TIME
      ENDIF

      F31(I)=F31(I)-A120(I,J)*TT210(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-(AKMV30(I,J)+AKVV20(I,J))*TT210(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-A210(I,J)*TT210(J)
      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)-AKTV20(I,J)*TT210(J)
     1-AKTV30(I,J)*TT210(J+NDIM)-AKK0(I,J)*TT210(J+2*NDIM+4)
C DESNA STRANA USLED UTICAJA KONVEKCIJE
c right-hand side created by natural convection
      F31(I)=F31(I)-AMV20(I,J)*TT210(J+2*NDIM+4)*BETA*FB2
      F31(I+NDIM)=F31(I+NDIM)-AMV20(I,J)*TT210(J+2*NDIM+4)*BETA*FB3
 286  CONTINUE							      
 291  CONTINUE

C=========================================================================
C==========================================================================

C==========================================================================

C==========================================================================
C==========================================================================



C============================================================================
C  FOR FLUX CALCULATION
C      IF(IPROL.EQ.1) THEN
       DO I=1,NDIM
        N=NEL(I,NBREL)
C        JJ=ID(5,N)
C        IF(JJ.NE.0) THEN
         POT=0.D0
         DO J=1,NDES
C         DO J=1+2*NDIM+4,NDES
          POT=POT+SKEF(I+2*NDIM+4,J)*TT21(J)
         ENDDO
         PRES(1,N)=PRES(1,N)+POT
C        ENDIF
       ENDDO
C      ENDIF
C============================================================================


C      CALL WRRF(SKEF,NDES*NDES,'SKEF ',IIZLAZ)      
C      CALL WRRF(F31,NDES,'F31= ',IIZLAZ)      
C       CALL SPAKDE (SILE,F31,LM2,NDES)

      DO I=1,NDES
       F31(I)=F31(I)-FALE(I)     
      ENDDO 

      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F31,MAXA,LM2,NDES,1)

C RACUNANJE SILA INTERAKCIJE KOJIMA FLUID DELUJE NA ZIDOVE
c calculation of interaction forces between fluid and solid
	 IF (NSTAC.EQ.0) THEN
        DO I=1,NDIM
          II=I+NDIM
         DO J=1,NDIM
          JJ=J+NDIM
          SKEF(I,J)=SKEF(I,J)-AMV2(I,J)*TT210(J)/TIME
          SKEF(II,JJ)=SKEF(II,JJ)-AMV2(I,J)*TT210(JJ)/TIME
         ENDDO
        ENDDO
       ENDIF

       FACAXY=1.D0
       IF(INDAX.EQ.1) FACAXY=8.D0*DATAN(1.D0)

       DO I=1,NDIM
       N=NEL(I,NBREL)
        DO J=1,NDES
         P1=SKEF(I,J)*TT21(J)
         P2=SKEF(I+NDIM,J)*TT21(J)
         SPSIL(1,N)=SPSIL(1,N)-P1/FACAXY
         SPSIL(2,N)=SPSIL(2,N)-P2/FACAXY
       ENDDO
      ENDDO



C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C end of loop per elements
C=======================================================================
 400  CONTINUE

C End of subroutine
      END
C=======================================================================
C=========================================================================
      SUBROUTINE ALE3DN(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine ALE3DN is used for 3D ALE analysis
CE It is used global loop per elements
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21),TT21A(3*21)
      DIMENSION H(21),HP(8)
      DIMENSION AKVV1(21,21),AKMIV1(21,21),AMV2(21,21),AJV1V1(21,21)
      DIMENSION AJV1V2(21,21),AJV1V3(21,21),AJV2V1(21,21),AJV2V2(21,21)
      DIMENSION AJV2V3(21,21),AJV3V1(21,21),AJV3V2(21,21),AJV3V3(21,21)
      DIMENSION AKTV1(21,21),AKTV2(21,21),AKTV3(21,21),AKK(21,21)
      DIMENSION AKV1P(21,8),AKV2P(21,8),AKV3P(21,8),RS1(21),RS2(21)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(8,8),PENXY(8,8),PENXZ(8,8)
      DIMENSION PENYX(8,8),PENYY(8,8),PENYZ(8,8)
      DIMENSION PENZX(8,8),PENZY(8,8),PENZZ(8,8)
      DIMENSION C(21,21),VSTAR(3)
      DATA NREF/0,1,3,6,10,15/



 
      DATA WGT/2.,1.00000000000000,1.00000000000000,
     10.55555555555556,0.88888888888889,0.55555555555556,
     20.34785484513745,0.65214515486255,0.65214515486255,
     30.34785484513745,0.23692688505619,0.47862867049937,
     40.56888888888889,0.47862867049937,0.23692688505619/
 
      DATA XG/0.,-0.5773502691896,0.57735026918963,
     1-0.7745966692415,0.,0.77459666924148,-0.8611363115941,
     2-0.3399810435849,0.33998104358486,0.86113631159405,
     3-0.9061798459387,-0.5384693101057,0.,0.53846931010568,
     40.90617984593866/


C      R(3,1)=-0.7745966692415
C      R(3,2)=0.0
C      R(3,3)=0.77459666924148
C      S(3,1)=-0.7745966692415
C      S(3,2)=0.0
C      S(3,3)=0.77459666924148
C      W(3,1)=0.55555555555556
C      W(3,2)=0.88888888888889
C      W(3,3)=0.55555555555556


C      R(2,1)=0.57735026918963
C      R(2,2)=-0.5773502691896
C      S(2,1)=0.57735026918963
C      S(2,2)=-0.5773502691896
C      W(2,1)=1.000000000000000
C      W(2,2)=1.000000000000000

C      R(1,1)=0.0
C      S(1,1)=0.0
C      W(1,1)=2.0
         
C       PENALT=1.D9








C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET

      DO 125 I=1,92
      TT210(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0


C=========================================================================
      DO 130 KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      LM2(KLM+2*NDIM)=ID(3,NEL(KLM,NBREL))
      IF (KLM.LE.8) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
      LM2(KLM+3*NDIM+8)=ID(5,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
C=========================================================================
      IF (IALE.EQ.1) THEN
       DO I=1,NDIM
        NODE=NEL(I,NBREL)
        TT21A(I)=VMESH(1,NODE)
        TT21A(I+NDIM)=VMESH(2,NODE)
        TT21A(I+2*NDIM)=VMESH(3,NODE)
       ENDDO
      ENDIF
C=========================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
C      IF (KLM.LE.8.AND.ID(4,NEL(KLM,NBREL)).NE.0) THEN
      IF (KLM.LE.8) THEN
        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
        TT21(KLM+3*NDIM+8)=GNODE(2,5,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM+8)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE
C=======================================================================
C      DO 160 NR=1,NDIM
C      DO 155 NZDT=1,NUMZAD
C      IF(NEL(NR,NBREL).EQ.NZAD(1,NZDT)) THEN
C       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
C     &NTABFT,IIZLAZ)
C      MESTO=NZAD(2,NZDT)
C       IF ((MESTO.EQ.4).AND.(NR.GT.8)) GOTO 160
C       IF (MESTO.EQ.5) THEN
C        TT21(3*NDIM+8+NR)=ZADVRE(NZDT)*FK1
C       ELSE
C        TT21((MESTO-1)*NDIM+NR)=ZADVRE(NZDT)*FK1
C       ENDIF
C      ENDIF  
C 155  CONTINUE
C 160  CONTINUE
C=======================================================================
C======================================================================= 
C RACUNANJE PRITISAKA 
C      IF (INDPR.EQ.1) THEN
C      K=0
C      DO 450 I=1,IBRGT-1
C      DO 450 J=1,IBRGT-1
C      K=K+1
C      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),1)
C      PRIT(K,NBREL)=PRESS
C450   CONTINUE
C      GOTO 400
C      ENDIF
C
C=======================================================================
C      CALL WRR(TT21,3*NDIM+4,'T211')
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
  
      IF ((K.LE.8).AND.(N.LE.8)) THEN
       PENXX(K,N)=0.D0
       PENXY(K,N)=0.D0
       PENXZ(K,N)=0.D0
       PENYX(K,N)=0.D0
       PENYY(K,N)=0.D0
       PENYZ(K,N)=0.D0
       PENZX(K,N)=0.D0
       PENZY(K,N)=0.D0
       PENZZ(K,N)=0.D0
      ENDIF
      IF (N.LE.8) THEN
       AKV1P(K,N)=0.D0
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
      ENDIF
      AKVV1(K,N)=0.D0
      AKTV1(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKMIV1(K,N)=0.D0
C      AKMIV3(K,N)=0.D0
C      A12(K,N)=0.
C      A21(K,N)=0.
      AMV2(K,N)=0.D0
      C(K,N)=0.D0
      AJV1V1(K,N)=0.D0
      AJV1V2(K,N)=0.D0
      AJV1V3(K,N)=0.D0
      AJV2V1(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V1(K,N)=0.D0
      AJV3V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AKK(K,N)=0.D0
C      AJK(K,N)=0.D0
  162 CONTINUE
  163 CONTINUE


C POVRSINSKE SILE I ZAPREMINSKE SILE:
      DO 199 I=1,NDIM
      RS1(I)=0.D0
      RS2(I)=0.D0
      RS3(I)=0.D0
      RB1(I)=0.D0
      RB2(I)=0.D0
      RB3(I)=0.D0
 199  CONTINUE

C===========================================================================

      NGAUSX=IBRGT
      NGAUSY=IBRGT
      NGAUSZ=IBRGT
C INTEGRACIJA U GAUSOVIM TACKAMA
C      DO 180 I=1,IBRGT
C      DO 170 J=1,IBRGT
C      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0)
C      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

      ZAPRE=0.D0

      DO 180 NGX=1,NGAUSX
      JR=NREF(NGAUSX) + NGX
      R = XG(JR)
 
      DO 175 NGY=1,NGAUSY
      JS=NREF(NGAUSY) + NGY
      S = XG(JS)
 
      DO 170 NGZ=1,NGAUSZ
      JT=NREF(NGAUSZ) + NGZ
      T = XG(JT)
 
      WT=WGT(JR)*WGT(JS)*WGT(JT)
 
       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)

      CALL ALEF(FALE,PJ,NEL,NDIM,GNODE,VMESH,VMESH0,TT21,H,IPN,NETIP,
     &GUSM,WDT,HP,TIME,NBREL,VSTAR,TT210)

       IF (IALE.EQ.1) CALL ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
 
      WDT=WT*DET1
      ZAPRE=ZAPRE+WDT
  
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM
      AKVV1(K,N)=AKVV1(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM)
      AKTV1(K,N)=AKTV1(K,N)+WDT*(H(K)*H(N)*ZVXT*GUSM*CC)
      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVYT*GUSM*CC)
      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVZT*GUSM*CC)
      AKMIV1(K,N)=AKMIV1(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
      AKK(K,N)=AKK(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AKT)
C      AKMIV3(K,N)=AKMIV1(K,N)
      AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM
      AJV1V1(K,N)=AJV1V1(K,N)+H(K)*HXU*H(N)*GUSM*WDT
      AJV1V2(K,N)=AJV1V2(K,N)+H(K)*HYU*H(N)*GUSM*WDT
      AJV1V3(K,N)=AJV1V3(K,N)+H(K)*HZU*H(N)*GUSM*WDT
      AJV2V1(K,N)=AJV2V1(K,N)+H(K)*HXV*H(N)*GUSM*WDT
      AJV2V2(K,N)=AJV2V2(K,N)+H(K)*HYV*H(N)*GUSM*WDT
      AJV2V3(K,N)=AJV2V3(K,N)+H(K)*HZV*H(N)*GUSM*WDT
      AJV3V1(K,N)=AJV3V1(K,N)+H(K)*HXW*H(N)*GUSM*WDT
      AJV3V2(K,N)=AJV3V2(K,N)+H(K)*HYW*H(N)*GUSM*WDT
      AJV3V3(K,N)=AJV3V3(K,N)+H(K)*HZW*H(N)*GUSM*WDT
       IF (N.LE.8.AND.PENALT.LT.1.D0) THEN
       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
      ENDIF
 164  CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(0.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(0.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  175 CONTINUE
  180 CONTINUE



      IF (PENALT.GT.0.D0) THEN
       NGAUSX=IBRGT-1
       NGAUSY=IBRGT-1
       NGAUSZ=IBRGT-1

       DO  NGX=1,NGAUSX
       JR=NREF(NGAUSX) + NGX
       R = XG(JR)
 
       DO  NGY=1,NGAUSY
       JS=NREF(NGAUSY) + NGY
       S = XG(JS)
 
       DO  NGZ=1,NGAUSZ
       JT=NREF(NGAUSZ) + NGZ
       T = XG(JT)
 
       WT=WGT(JR)*WGT(JS)*WGT(JT)
 
       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)

 
        WDT=WT*DET1
  
        DO  K=1,NDIM
        DO  N=1,NDIM
C      IF (N.LE.8.AND.PENALT.LT.1.D0) THEN
C       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
C      ENDIF
C         IF (PENALT.GT.0.D0) THEN
          PENXX(K,N)=PENXX(K,N)+WDT*PJ(1,K)*PJ(1,N)*PENALT
          PENXY(K,N)=PENXY(K,N)+WDT*PJ(1,K)*PJ(2,N)*PENALT
          PENXZ(K,N)=PENXZ(K,N)+WDT*PJ(1,K)*PJ(3,N)*PENALT
          PENYX(K,N)=PENYX(K,N)+WDT*PJ(2,K)*PJ(1,N)*PENALT
          PENYY(K,N)=PENYY(K,N)+WDT*PJ(2,K)*PJ(2,N)*PENALT
          PENYZ(K,N)=PENYZ(K,N)+WDT*PJ(2,K)*PJ(3,N)*PENALT
          PENZX(K,N)=PENZX(K,N)+WDT*PJ(3,K)*PJ(1,N)*PENALT
          PENZY(K,N)=PENZY(K,N)+WDT*PJ(3,K)*PJ(2,N)*PENALT
          PENZZ(K,N)=PENZZ(K,N)+WDT*PJ(3,K)*PJ(3,N)*PENALT
C         ENDIF
        ENDDO
        ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDIF



C       WRITE (IIZLAZ,*)'ZAPREMINA=',ZAPRE
C===========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PENALTY FAKTOROM
C     IF (PENALT.GT.0.D0) THEN
C     DO 200 I=1,IBRGT-1
C     DO 195 J=1,IBRGT-1
C     CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0)
C     WDT1=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
C     DO 190 K=1,NDIM
C     DO 185 N=1,NDIM
C     AKMIV1(K,N)=AKMIV1(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT1
C     AKMIV3(K,N)=AKMIV3(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT1
C     A12(K,N)=A12(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT1
C     AJV2V3(K,N)=AJV2V3(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT1
C     A21(K,N)=A21(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT1
C     AJV3V2(K,N)=AJV3V2(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT1
C 185 CONTINUE
C 190 CONTINUE
C 195 CONTINUE
C 200 CONTINUE
C     ENDIF
C===========================================================================

C======================================================================= 
C POVRSINSKE SILE
      

      DO 250 JBRPS=1,MAXSIL
      IF (NGPSIL(1,JBRPS).EQ.NBREL) THEN 
      CALL STRANA(NEL,NDIM,NGPSIL(1,JBRPS),NPOV)
C      IF(JBRPS.LE.4) NPOV=1
C      IF(JBRPS.GT.4) NPOV=2
C       WRITE(IIZLAZ,*)'STRANA=',NPOV
       POVRS=0.D0
C
CS  PETLJA PO GAUSOVIM TACKAMA
CE  GAUSS POINTS LOOP
C
C   40 IF(NPOV.GT.2) GO TO 45
      IF(NPOV.GT.2) GO TO 45
      KFIX=1
      NGXP=NGAUSY
      NGYP=NGAUSZ
      R=1.
      IF(NPOV.EQ.2) R=-1.
      GO TO 60
C
   45 IF(NPOV.GT.4) GO TO 50
      KFIX=2
      NGXP=NGAUSX
      NGYP=NGAUSZ
      S=1.
      IF(NPOV.EQ.4) S=-1.
      GO TO 60
C
   50 NGXP=NGAUSX
      NGYP=NGAUSY
      KFIX=3
      T=1.
      IF(NPOV.EQ.6) T=-1.
C
   60 DO 210 NGX=1,NGXP
      JR=NREF(NGXP) + NGX
      XX = XG(JR)
      WX=WGT(JR)
      DO 210 NGY=1,NGYP
      JS=NREF(NGYP) + NGY
      YY = XG(JS)
      WY=WGT(JS)
      GO TO (71,72,73),KFIX
   71 S=XX
      T=YY
      GO TO 75
   72 R=XX
      T=YY
      GO TO 75
   73 R=XX
      S=YY
   75 WT=WX*WY
C
C      WRITE(IIZLAZ,*)'POVRS KFIX=',KFIX
C       KFIX=1
       CALL JACT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)
C
      WDT=WT*DET
      POVRS=POVRS+WDT

C      WRITE(IIZLAZ,*)'SF1=',SF1
      DO K=1,NDIM
       RS1(K)=RS1(K)+(H(K)*SF1)*WDT*NGPSIL(6,JBRPS)
       RS2(K)=RS2(K)+(H(K)*SF2)*WDT*NGPSIL(7,JBRPS)
       RS3(K)=RS3(K)+(H(K)*SF3)*WDT*NGPSIL(8,JBRPS)
      ENDDO
C       WRITE(IIZLAZ,*)'NPOV',NPOV,'NBREL=',NBREL
C       CALL IWRRF(NGPSIL(1,JBRPS),5,'NGPSI',IIZLAZ)
C       CALL WRRF(RS3,21,'RS3=',IIZLAZ)


  210 CONTINUE    
C       WRITE(IIZLAZ,*)'POVRSINA',JBRPS,'=',POVRS
      ENDIF
  250 CONTINUE


C INCIJALIZACIJA MATRICE SKEF I F92
      DO 260 I=1,NDES
      DO 258 J=1,NDES
      SKEF(I,J)=0.D0
 258  CONTINUE
       F92(I)=0.D0
 260  CONTINUE
C=========================================================================
C=========================================================================
C PAKOVANJE MATRICE C,Kk,Jk U MATRICU SKEF
      DO 263 I=1,NDIM
      DO 262 J=1,NDIM
C       AJK(I,J)=AKVV1(I,J)*CC
       SKEF(I+3*NDIM+8,J+3*NDIM+8)=AKK(I,J)+AKVV1(I,J)*CC
      SKEF(I+3*NDIM+8,J)=AKTV1(I,J)
      SKEF(I+3*NDIM+8,J+NDIM)=AKTV2(I,J)
      SKEF(I+3*NDIM+8,J+2*NDIM)=AKTV3(I,J)
      IF (NSTAC.EQ.0) THEN
      SKEF(I+3*NDIM+8,J+3*NDIM+8)=SKEF(I+3*NDIM+8,J+3*NDIM+8)+
     &AMV2(I,J)*CC/TIME
      ENDIF
 262  CONTINUE
 263  CONTINUE
C=========================================================================
C     CALL NUL(FPOM,21)
C
C     DO I=1,NDIM
C     DO J=1,NDIM
C       FPOM(I)=FPOM(I)+AKMIV1(I,J)*TT21(J)
C      IF (J.LE.8) FPOM(I)=FPOM(I)+AKV1P(I,J)*TT21(J+3*NDIM)
C     ENDDO
C       FPOM(I)=FPOM(I)-RS1(I)
C     ENDDO
C
C     CALL WRR(FPOM,21,'FPOM=')

C LEVA STRANA:
      DO 270 I=1,NDIM
      DO 265 J=1,NDIM
       SKEF(I,J)=AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J)
       SKEF(I,J+NDIM)=AJV1V2(I,J)
       SKEF(I,J+2*NDIM)=AJV1V3(I,J)
      IF (J.LE.8) THEN
       SKEF(I,J+3*NDIM)=AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AKV3P(I,J)
      ENDIF
      IF (I.LE.8) THEN
       SKEF(I+3*NDIM,J)=AKV1P(J,I)
       SKEF(I+3*NDIM,J+NDIM)=AKV2P(J,I)
       SKEF(I+3*NDIM,J+2*NDIM)=AKV3P(J,I)
      ENDIF

       SKEF(I+NDIM,J)=AJV2V1(I,J)
       SKEF(I+NDIM,J+NDIM)=AKMIV1(I,J)+AKVV1(I,J)+AJV2V2(I,J)
       SKEF(I+NDIM,J+2*NDIM)=AJV2V3(I,J)

      SKEF(I+2*NDIM,J)=AJV3V1(I,J)
      SKEF(I+2*NDIM,J+NDIM)=AJV3V2(I,J)
      SKEF(I+2*NDIM,J+2*NDIM)=AKMIV1(I,J)+AKVV1(I,J)+AJV3V3(I,J)
      IF (NSTAC.EQ.0) THEN
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       SKEF(I+2*NDIM,J+2*NDIM)=SKEF(I+2*NDIM,J+2*NDIM)+AMV2(I,J)/TIME
      ENDIF
 265  CONTINUE
 270  CONTINUE
       
      IF (PENALT.GT.0.D0) THEN
      DO I=1,NDIM
        I1=I
        I2=I+NDIM
        I3=I+2*NDIM
       DO J=1,NDIM
        J1=J
        J2=J+NDIM
        J3=J+2*NDIM
        SKEF(I1,J1)=SKEF(I1,J1)+PENXX(I,J)
        SKEF(I1,J2)=SKEF(I1,J2)+PENXY(I,J)
        SKEF(I1,J3)=SKEF(I1,J3)+PENXZ(I,J)
        SKEF(I2,J1)=SKEF(I2,J1)+PENYX(I,J)
        SKEF(I2,J2)=SKEF(I2,J2)+PENYY(I,J)
        SKEF(I2,J3)=SKEF(I2,J3)+PENYZ(I,J)
        SKEF(I3,J1)=SKEF(I3,J1)+PENZX(I,J)
        SKEF(I3,J2)=SKEF(I3,J2)+PENZY(I,J)
        SKEF(I3,J3)=SKEF(I3,J3)+PENZZ(I,J)
        F92(I1)=F92(I1)-PENXX(I,J)*TT21(J1)
        F92(I1)=F92(I1)-PENXY(I,J)*TT21(J2)
        F92(I1)=F92(I1)-PENXZ(I,J)*TT21(J3)
        F92(I2)=F92(I2)-PENYX(I,J)*TT21(J1)
        F92(I2)=F92(I2)-PENYY(I,J)*TT21(J2)
        F92(I2)=F92(I2)-PENYZ(I,J)*TT21(J3)
        F92(I3)=F92(I3)-PENZX(I,J)*TT21(J1)
        F92(I3)=F92(I3)-PENZY(I,J)*TT21(J2)
        F92(I3)=F92(I3)-PENZZ(I,J)*TT21(J3)
       ENDDO
      ENDDO
      ENDIF

C     DO 272 I=1,NDIM
C     DO 271 J=1,NDIM
C     SKEF(I+2*NDIM+4,J)=AKTV2(I,J)
C     SKEF(I+2*NDIM+4,J+NDIM)=AKTV3(I,J)
C271  CONTINUE
C272  CONTINUE

C LEVA STRANA USLED KONVEKCIJE
C     DO 282 I=1,NDIM
C     DO 281 J=1,NDIM
C     SKEF(I,J+2*NDIM+4)=SKEF(I,J+2*NDIM+4)+AMV2(I,J)*FB2*BETA
C     SKEF(I+NDIM,J+2*NDIM+4)=SKEF(I+NDIM,J+2*NDIM+4)+AMV2(I,J)*FB3*BETA
C281  CONTINUE
C282  CONTINUE

C======================================================================
C      IF (INDSIM(AKMIV1,NDIM).EQ.0) THEN
C        WRITE(*,*)'MATRICA AKMIV1 NIJE SIMETRICNA'
C        STOP
C      ENDIF
C      IF (INDSIM(AKVV1,NDIM).EQ.0) THEN
C        WRITE(*,*)'MATRICA AKVV1 NIJE SIMETRICNA'
C        STOP
C      ENDIF

C DESNA STRANA:
      DO K=0,2
      DO 290 I=1,NDIM
        II=I+K*NDIM
      DO 285 J=1,NDIM
        JJ=J+K*NDIM
        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
      IF (NSTAC.EQ.0) THEN
        F92(II)=F92(II)-AMV2(I,J)*(TT21(JJ)-TT210(JJ))/TIME
C       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
C     F92(I)=F92(I)-A12(I,J)*TT21(J+NDIM)
 285  CONTINUE							      
 290  CONTINUE
      ENDDO

C=========================================================================
      DO 294 I=1,NDIM
      DO 292 J=1,NDIM
      II=I+3*NDIM+8
      JJ=J+3*NDIM+8
       F92(II)=F92(II)-AKTV1(I,J)*TT21(J)-AKTV2(I,J)*TT21(J+NDIM)
     &-AKTV3(I,J)*TT21(J+2*NDIM)-AKK(I,J)*TT21(JJ)
      IF (NSTAC.EQ.0) THEN
       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
 292  CONTINUE
 294  CONTINUE

C==========================================================================
C       CALL NUL(FPOM,21)
C DESNA STRANA USLED PRITISKA
      DO 310 I=1,NDIM
      DO 300 J=1,NDIM
      IF (J.LE.8) THEN
       F92(I)=F92(I)-AKV1P(I,J)*TT21(J+3*NDIM)
C       FPOM(I)=FPOM(I)+AKV1P(I,J)*TT21(J+3*NDIM)
       F92(I+NDIM)=F92(I+NDIM)-AKV2P(I,J)*TT21(J+3*NDIM)
       F92(I+2*NDIM)=F92(I+2*NDIM)-AKV3P(I,J)*TT21(J+3*NDIM)
      ENDIF
      IF (I.LE.8) THEN
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV1P(J,I)*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV2P(J,I)*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV3P(J,I)*TT21(J+2*NDIM)
      ENDIF
 300  CONTINUE
 310  CONTINUE

C==========================================================================
C ZAPREMINSKE SILE SA DESNE STRANE
C      DO 315 I=1,NDIM
C       F92(I)=F92(I)+RB1(I)
C       F92(I+NDIM)=F92(I+NDIM)+RB2(I)
C       F92(I+2*NDIM)=F92(I+2*NDIM)+RB3(I)
C 315  CONTINUE
C==========================================================================
C       CALL WRR(RS1,21,'RS1= ')
C       CALL WRR(FPOM,21,'FPOM=')
C==========================================================================
C POVRSINSKE SILE SA DESNE STRANE
      DO 320 I=1,NDIM
       F92(I)=F92(I)+RS1(I)
       F92(I+NDIM)=F92(I+NDIM)+RS2(I)
       F92(I+2*NDIM)=F92(I+2*NDIM)+RS3(I)
320   CONTINUE
C==========================================================================





      CALL SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)

C      CALL WRRF(F92,NDES,'F92= ',IIZLAZ)      
C      CALL WRRF(SKEF,NDES*NDES,'SKEF=',IIZLAZ)      
      
      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F92,MAXA,LM2,NDES,1)



C RACUNANJE SILA INTERAKCIJE KOJIMA FLUID DELUJE NA ZIDOVE
c calculation of interaction forces between fluid and solid


	 IF (NSTAC.EQ.0) THEN
       DO K=0,2
        DO I=1,NDIM
          II=I+K*NDIM
         DO J=1,NDIM
          JJ=J+K*NDIM
          SKEF(II,JJ)=SKEF(II,JJ)-AMV2(I,J)*TT210(JJ)/TIME
         ENDDO
        ENDDO
       ENDDO
       ENDIF


       DO I=1,NDIM
       N=NEL(I,NBREL)
        DO J=1,NDES
         P1=SKEF(I,J)*TT21(J)
         P2=SKEF(I+NDIM,J)*TT21(J)
         P3=SKEF(I+2*NDIM,J)*TT21(J)
         SPSIL(1,N)=SPSIL(1,N)-P1
         SPSIL(2,N)=SPSIL(2,N)-P2
         SPSIL(3,N)=SPSIL(3,N)-P3
       ENDDO
      ENDDO



C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE

C End of subroutine

      END
C=======================================================================
C==========================================================================
      SUBROUTINE ALEF(FALE,PJ,NEL,NDIM,GNODE,VMESH,VMESH0,TT21,H,IPN,IN,
     &RO,WDT,HP,TIME,NBREL,VSTAR,TT210)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION NEL(NDIM+1,*)

      DIMENSION FALE(*),PJ(3,*),GNODE(2,6,*),VMESH(3,*),VMESH0(3,*)
      DIMENSION TT21(*),H(*),HP(*),VSTAR(3),TT210(*)

     
      DIMENSION DVM(3),DUM(3,3),VD(3,3),V(3,21),C(3),PP(8),V0(3,21)
      DIMENSION DVD(3,3),VM(3)


C=========================================================================
C DODATAK ZA ALE FORMULACIJU
C
C INITIALISATION VALUES
       DO I=1,IN
        DVM(I)=0.D0
         DO J=1,IN
          DUM(I,J)=0.D0
          VD(I,J)=0.D0
          DVD(I,J)=0.D0
         ENDDO
       ENDDO
       P=0.D0
       CALL CLEAR(C,3)
       CALL CLEAR(VM,3)

       DO JJ=1,NDIM
       NODE=NEL(JJ,NBREL)
        DO I=1,IN
          V(I,JJ)=TT21(JJ+(I-1)*NDIM)
          V0(I,JJ)=TT210(JJ+(I-1)*NDIM)
          C(I)=C(I)+H(JJ)*(GNODE(2,I,NODE)-VMESH(I,NODE))
	    VM(I)=VM(I)+H(JJ)*VMESH(I,NODE)
          IF (JJ.LE.IPN) PP(JJ)=TT21(JJ+IN*NDIM)
          DVM(I)=DVM(I)+H(JJ)*(VMESH(I,NODE)-VMESH0(I,NODE))
         DO J=1,IN
          DUM(I,J)=DUM(I,J)+PJ(J,JJ)*VMESH(I,NODE)*TIME
          VD(I,J)=VD(I,J)+PJ(J,JJ)*V(I,JJ)
          DVD(I,J)=DVD(I,J)+PJ(J,JJ)*(V(I,JJ)-V0(I,JJ))
         ENDDO
        IF (JJ.LE.IPN) P=P+H(JJ)*PP(JJ)
        ENDDO
       ENDDO
       
        DIVUM=0.D0
        DO K=1,IN
         DIVUM=DIVUM+DUM(K,K)
        ENDDO



C ADDITION FOR TERM FALE1
    
       DO JJ=1,NDIM
        DO I=1,IN
         IDIM=JJ+NDIM*(I-1)
         FALE(IDIM)=FALE(IDIM)+RO*H(JJ)*VSTAR(I)*DIVUM*WDT
        ENDDO
       ENDDO

C ADDITION FOR TERM FALE2

       DO JJ=1,NDIM
        DO I=1,IN
         IDIM=JJ+NDIM*(I-1)
        ALE2=0.D0
        DO J=1,IN
         DO K=1,IN
       ALE2=ALE2-DVM(J)*VD(I,J)+C(J)*(VD(I,J)*DUM(K,K)-VD(I,K)*DUM(K,J))
         ENDDO
       ALE2=ALE2-DVD(I,J)*VM(J)
        ENDDO
         FALE(IDIM)=FALE(IDIM)+RO*H(JJ)*ALE2*WDT
        ENDDO
       ENDDO

C ADDITION FOR TERM FALE3

       DO JJ=1,NDIM
        DO I=1,IN
         IDIM=JJ+NDIM*(I-1)
        ALE3=0.D0
         DO K=1,IN
         ALE3=ALE3+PJ(K,JJ)*DUM(K,I)*P-PJ(I,JJ)*P*DUM(K,K)
C         ALE3=ALE3-PJ(I,JJ)*P*DUM(K,K)
         ENDDO
         FALE(IDIM)=FALE(IDIM)+ALE3*WDT
        ENDDO
       ENDDO

C ADDITION FOR TERM FALE4

       DO JJ=1,NDIM
        DO I=1,IN
         IDIM=JJ+NDIM*(I-1)
        ALE4=0.D0
         DO J=1,IN
          DO K=1,IN
           ALE4=ALE4-PJ(K,JJ)*DUM(K,J)*VD(I,J)+
     &PJ(J,JJ)*(VD(I,J)*DUM(K,K)-VD(I,K)*DUM(K,J))
C     &PJ(J,JJ)*(VD(I,J)*DUM(K,K))
          ENDDO
         ENDDO
         FALE(IDIM)=FALE(IDIM)+ALE4*WDT
        ENDDO
       ENDDO



C ADDITION FOR TERM FALE5

       DO JJ=1,IPN
        DO I=1,IN
         IDIM=JJ+NDIM*(I-1)
        ALE5=0.D0
         DO K=1,IN
          ALE5=ALE5+VD(I,I)*DUM(K,K)-VD(I,K)*DUM(K,I)
         ENDDO
         FALE(IDIM)=FALE(IDIM)+HP(JJ)*ALE5*WDT
        ENDDO
       ENDDO


C=========================================================================


      END
C==========================================================================
       SUBROUTINE FSI2D(FALE,ZVHX,ZVHY,NEL,NDIM,GNODE,VMESH,VMESH0,TT21,
     &H,IPN,NETIP,RO,WDT,HP,NBREL,TIME,TT2100,TT210)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION NEL(NDIM+1,*)

      DIMENSION FALE(*),ZVHX(*),ZVHY(*),GNODE(2,6,*),VMESH(3,*)
      DIMENSION TT21(*),H(*),HP(*),VMESH0(3,*),TT2100(*),TT210(*)
      DIMENSION PJ(3,21),VSTAR(3)

             
       DO I=1,3
         VSTAR(I)=0.D0
        DO J=1,21
         PJ(I,J)=0.D0
        ENDDO
       ENDDO

       DO I=1,NDIM
        PJ(1,I)=ZVHX(I)
        PJ(2,I)=ZVHY(I)
	  VSTAR(1)=VSTAR(1)+H(I)*(TT210(I)-TT2100(I))/TIME
	  VSTAR(2)=VSTAR(2)+H(I)*(TT210(I+NDIM)-TT2100(I+NDIM))/TIME
	  VSTAR(3)=VSTAR(3)+H(I)*(TT210(I+2*NDIM)-TT2100(I+2*NDIM))/TIME
       ENDDO



      CALL ALEF(FALE,PJ,NEL,NDIM,GNODE,VMESH,VMESH0,TT21,H,IPN,NETIP,
     &RO,WDT,HP,TIME,NBREL,VSTAR,TT210)


      END
C==========================================================================
C=========================================================================
      SUBROUTINE GREZON(VMESH,CCORD,CCORD0,NPT,NETIP,TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VMESH(3,*),CCORD(3,*),CCORD0(3,*)

       DO I=1,NPT
        DO J=1,NETIP
         VMESH(J,I)=(CCORD(J,I)-CCORD0(J,I))/TIME
        ENDDO
       ENDDO
      END
C==========================================================================
C=========================================================================
      SUBROUTINE TCORDC(CCORD,CCORD0,NPT,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CCORD(3,*),CCORD0(3,*)

       DO I=1,NPT
        DO J=1,NETIP
         CCORD0(J,I)=CCORD(J,I)
        ENDDO
       ENDDO
      END
C==========================================================================

