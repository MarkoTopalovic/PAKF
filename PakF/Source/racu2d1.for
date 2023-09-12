C=========================================================================
      SUBROUTINE RACU2D1(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,AF,INDEL,NJUTRA,ISYMMS,NZAD,ZADVRE,
     &IPASS,NUMPASS,METOD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine RACU2D is used for 2D analysis
CE It is used global loop per elements
C

C SAMO PRIVREMENO UBACEN COMMON ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,DPRES
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
	COMMON /ALPHA_SEG/ ALPHAU,ALPHAV,ALPHAW,ALPHAP

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
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*),NZAD(3,*)
      DIMENSION SKEF(NDES,*),zadvre(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
      DIMENSION SPSIL(2,*),PRES(3,*),VMESH(3,*)

      DIMENSION X(9),Y(9),LM2(4*9),TT21(4*9),TT210(4*9),TT21A(2*9)
      DIMENSION R(3,3),S(3,3),W(3,3)
C      DIMENSION R1(2),S1(2),W1(2)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4),ZVHXP(9),ZVHYP(9)
      DIMENSION AKVV2(9,9),TT21P(36)
      DIMENSION AKMIV2(9,9)
      DIMENSION AKMIV2P(9,9),AKMIV3P(9,9),AKVV2P(9,9),AKVV3P(9,9)
	DIMENSION AKIIX(9,9),AKIIY(9,9)
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
      DIMENSION AKV2P(9,4),AKV2PT(9,4),AKV3PT(9,4)
      DIMENSION AKV3P(9,4),AKV2P1(9,4)
      DIMENSION AJK(9,9)
      DIMENSION RS2(9)
      DIMENSION RS3(9)
      DIMENSION RB2(9)
      DIMENSION RB3(9)
      DIMENSION F31(31)
      DIMENSION AKUAX(9,9),AKVAX(9,9),APAXIX(9,9),APAXIY(9,9)

C	DIMENSION SKE(31*16)
	DIMENSION SKE(31*16)

C Coorection for UPWIND, stability, Taylor, Oct, 21, 2004 Implemented
	DIMENSION ATAYL(9,9),AKPP(9,9)







      NDIMP=4



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
        

C=========================================================================
C      DO NODE=1,NPT
C       DO NZDT=1,NUMZAD
C        IF(NODE.EQ.NZAD(1,NZDT)) THEN
C         CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
C     &NTABFT,IIZLAZ)
C         MESTO=NZAD(2,NZDT)
C	   GNODE(2,MESTO,NODE)=ZADVRE(NZDT)*FK1
Cc	   GNODE(1,MESTO,NODE)=ZADVRE(NZDT)*FK1
C        ENDIF  
C	 ENDDO
C	ENDDO
c=======================================================================

c      if (kkorak.eq.1.and. iter.eq.0) 
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,2)

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,4*NDIM
      TT210(I)=0.D0
	TT21P(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0
C=========================================================================
      IF (IALE.EQ.1) THEN
       DO I=1,NDIM
        NODE=NEL(I,NBREL)
        TT21A(I)=VMESH(1,NODE)
        TT21A(I+NDIM)=VMESH(2,NODE)
       ENDDO
      ENDIF
C=========================================================================
      DO 130 KLM=1,NDIM
      X(KLM)=CORD(1,NEL(KLM,NBREL))
      Y(KLM)=CORD(2,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      IF (KLM.LE.4) LM2(KLM+2*NDIM)=ID(4,NEL(KLM,NBREL))
      LM2(KLM+2*NDIM+4)=ID(5,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,2
C      IF (ID(NR,NEL(KLM,NBREL)) .NE. 0) THEN
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
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
      DO 160 NR=1,NDIM
      DO 155 NZDT=1,NUMZAD
C
C TIMFUN is subroutine which determined values of function at currently 
C time step (at the end of step)
C
      IF(NEL(NR,NBREL).EQ.NZAD(1,NZDT)) THEN
       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
     &NTABFT,IIZLAZ)
C========================================
C SAMO PRIVREMENO UBACENO ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C       FK1=PBALL
C=========================================
      MESTO=NZAD(2,NZDT)
        IF (MESTO.EQ.4) THEN
        TT21P(2*NDIM+NR)=ZADVRE(NZDT)*FK1
        ELSE
        TT21P((MESTO-1)*NDIM+NR)=ZADVRE(NZDT)*FK1
        ENDIF
      ENDIF  
 155  CONTINUE
 160  CONTINUE
C=======================================================================

C=======================================================================
C      CALL WRR(TT21,3*NDIM+4,'T211')
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
      IF (N.LE.4) THEN
      AKV2P(K,N)=0.D0
      AKV3P(K,N)=0.D0
      AKV2PT(K,N)=0.D0
      AKV3PT(K,N)=0.D0
      AKV2P1(K,N)=0.D0
      ENDIF
      AKVV2(K,N)=0.D0
	AKPP(K,N)=0.D0
      APAXIX(K,N)=0.D0
      APAXIY(K,N)=0.D0
      AKUAX(K,N)=0.D0
      AKVAX(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKMIV2(K,N)=0.D0
      AKMIV3(K,N)=0.D0
      AKMIV2P(K,N)=0.D0
      AKMIV3P(K,N)=0.D0
	AKVV2P(K,N)=0.D0
	AKVV3P(K,N)=0.D0
      A12(K,N)=0.D0
      A21(K,N)=0.D0
      AMV2(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V2(K,N)=0.D0
	AKIIX(K,N)=0.D0
	AKIIY(K,N)=0.D0
      AKK(K,N)=0.D0
      AJK(K,N)=0.D0
c	ATAYL(K,N)=0.D0
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
C SURFACE FORCES AND BODY FORCES
       RS2(K)=0.D0
       RS3(K)=0.D0
       RB2(K)=0.D0
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
       IF (IALE.EQ.1) CALL ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
C       IF (IALE.EQ.1) CALL ALEHV4(TT21,TT21A,HV2,HV3,NDIM,H,ZVXV2,
C     &ZVYV3,ZVYV2,ZVXV3,ZVHX,ZVHY)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)

      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

      IF (INDAMI.EQ.1) CALL NENJUT(ZVHX,ZVHY,TT210)
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM
      IF (INDAX.EQ.1) THEN
      AKUAX(K,N)=AKUAX(K,N)-WDT*H(K)*(ZVHX(N)/DEBLJ-H(N)/(DEBLJ**2))*AMI
      AKVAX(K,N)=AKVAX(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AMI
      ENDIF


      AKVV2(K,N)=AKVV2(K,N)+WDT*((H(K)*HV2*ZVHX(N)+
     1H(K)*HV3*ZVHY(N))*GUSM)

C CONVECTIVE-DIFFUSION EXAMPLE WHERE V=1.0
c      AKVV2(K,N)=AKVV2(K,N)+WDT*((H(K)*ZVHX(N)+
c     1H(K)*ZVHY(N))*GUSM)

	AKVV2P(K,N)=AKVV2P(K,N)+WDT*ZVHX(K)*
     1(HV2*ZVHX(N)+HV3*ZVHY(N))*GUSM
	AKVV3P(K,N)=AKVV3P(K,N)+WDT*ZVHY(K)*
     1(HV2*ZVHX(N)+HV3*ZVHY(N))*GUSM
     
c=================================================================

      AKMIV2P(K,N)=AKMIV2P(K,N)+WDT*ZVHX(K)*(ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI
      AKMIV3P(K,N)=AKMIV3P(K,N)+WDT*ZVHY(K)*(ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI
c=================================================================
c=================================================================





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
      AJV2V2(K,N)=AJV2V2(K,N)+WDT*(H(K)*ZVXV2*H(N)*GUSM)
      AJV3V3(K,N)=AJV3V3(K,N)+WDT*(H(K)*ZVYV3*H(N)*GUSM)
      AJV2V3(K,N)=AJV2V3(K,N)+WDT*(H(K)*ZVYV2*H(N)*GUSM)
      AJV3V2(K,N)=AJV3V2(K,N)+WDT*(H(K)*ZVXV3*H(N)*GUSM)
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
C     AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
C      AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
C       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
c      ATAYL(K,N)=ATAYL(K,N)+WDT*(ZVHX(K)*ZVHX(N)*HV2**2+
c     1ZVHY(K)*ZVHY(N)*HV3**2)*GUSM*TIME/2.D0
c      ATAYL(K,N)=ATAYL(K,N)+WDT*(ZVHX(K)*ZVHY(N)*HV2*HV3+
c     1ZVHY(K)*ZVHX(N)*HV3*HV2)*GUSM*TIME/2.D0
  164 CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
C       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(0.D0+BETA*TETAO)*WDT
C       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(0.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE

C===========================================================================
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
c=================================================================
c=================================================================
C SEGREGATED PROCEDURE:
        DO II=1,NDIM
	   DO JJ=1,NDIM
	    AKIIX(II,II)=AKIIX(II,II)+
     &DABS(AKMIV2(II,JJ))+DABS(AKVV2(II,JJ))+DABS(AJV2V2(II,JJ))
	    AKIIY(II,II)=AKIIY(II,II)+
     &DABS(AKMIV3(II,JJ))+DABS(AKVV2(II,JJ))+DABS(AJV3V3(II,JJ))
C	    AKIIX(II,II)=AKIIX(II,II)+
C     &(AKMIV2(II,JJ)+AKVV2(II,JJ)+AJV2V2(II,JJ))
C	    AKIIY(II,II)=AKIIY(II,II)+
C     &(AKMIV3(II,JJ)+AKVV2(II,JJ)+AJV3V3(II,JJ))
	   ENDDO
C	    AKIIX(II,II)=1.D0
C	    AKIIY(II,II)=1.D0
C	    AKIIX(II,II)=AKIIX(II,II)/(1.D0/NDIM)
C	    AKIIY(II,II)=AKIIY(II,II)/(1.D0/NDIM)
C	    AKMIV2(II,II)=AKMIV2(II,II)+AKIIX(II,II)
C	    AKMIV3(II,II)=AKMIV3(II,II)+AKIIY(II,II)
c	    AKIIX(II,II)=1.D0/AKIIX(II,II)
c	    AKIIY(II,II)=1.D0/AKIIY(II,II)
	  ENDDO
c=================================================================
c=================================================================

C===========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PENALTY FAKTOROM
C reduced integration if penalty function is defined

      IF (PENALT.LT.1.D0) THEN
	IIBRGT=IBRGT-1
      DO I=1,IIBRGT
      DO J=1,IIBRGT
      CALL FNTERP(R(IIBRGT,I),S(IIBRGT,J),0,TT21,H,HP,ZVHXP,ZVHYP,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,4,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
      CALL FNTERP(R(IIBRGT,I),S(IIBRGT,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)
      WDT=W(IIBRGT,I)*W(IIBRGT,J)*DETJ
      DO K=1,NDIM
      DO N=1,NDIMP
c        AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
c        AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
        AKV2P(K,N)=AKV2P(K,N)-WDT*(H(K)*ZVHXP(N))
        AKV3P(K,N)=AKV3P(K,N)-WDT*(H(K)*ZVHYP(N))
        AKV2PT(K,N)=AKV2PT(K,N)+WDT*(-ZVHX(K)*HP(N))
        AKV3PT(K,N)=AKV3PT(K,N)+WDT*(-ZVHY(K)*HP(N))
       IF (INDAX.EQ.1) THEN
C        AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
        AKV2P1(K,N)=AKV2P1(K,N)+WDT*(H(K)*HP(N)/DEBLJ)
	 ENDIF
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDIF


      IF (PENALT.LT.1.D0.AND.METOD.EQ.4) THEN
      DO K=1,NDIMP
       DO N=1,NDIMP
	  DO K1=1,NDIM
	  AKPP(K,N)=AKPP(K,N)+
     &AKV2P(K1,K)*(1.D0/AKIIX(K1,K1))*AKV2P(K1,N)+
     &AKV3P(K1,K)*(1.D0/AKIIY(K1,K1))*AKV3P(K1,N)
       ENDDO  
	 ENDDO
	ENDDO


	IIBRGT=IBRGT-1
      DO I=1,IIBRGT
      DO J=1,IIBRGT
      CALL FNTERP(R(IIBRGT,I),S(IIBRGT,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,4,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
      WDT=W(IIBRGT,I)*W(IIBRGT,J)*DETJ
      DO K=1,NDIMP
       DO N=1,NDIMP
	  DO K1=1,NDIM
c	  AKPP(K,N)=AKPP(K,N)+WDT*(
c     &HP(K)*ZVHX(K1)*(1.D0/AKIIX(K1,K1))*ZVHX(K1)*HP(N)+
c     &HP(K)*ZVHY(K1)*(1.D0/AKIIY(K1,K1))*ZVHY(K1)*HP(N))
c	  AKPP(K,N)=AKPP(K,N)+WDT*(ZVHX(K)*(1.D0/AKIIX(K1,K1))*ZVHX(N)
c     &+ZVHY(K)*(1.D0/AKIIY(K1,K1))*ZVHY(N))
c	  AKPP(K,N)=AKPP(K,N)+WDT*(ZVHX(K)*ZVHX(N)
c     &+ZVHY(K)*ZVHY(N))
       ENDDO  
	 ENDDO
	ENDDO
	ENDDO
	ENDDO


	ENDIF


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
      AJV2V3(K,N)=AJV2V3(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT
      A21(K,N)=A21(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
      AJV3V2(K,N)=AJV3V2(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT
       IF (INDAX.EQ.1) THEN
       APAXIX(K,N)=APAXIX(K,N)+PENALT*WDT*(ZVHX(K)*H(N)/DEBLJ)
       APAXIY(K,N)=APAXIY(K,N)+PENALT*WDT*(ZVHY(K)*H(N)/DEBLJ)
       ENDIF
C      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
C      ENDIF
C      IF (INDAX.EQ.1) THEN
C       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
C      ENDIF
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
C        NDIM=4


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
      WDT=W(IBRGT-1,I)*DETJS*DEBLJ
c      WDT=W(IBRGT-1,I)*DETJS
C      IF(NGPSIL(4,JBRPS).EQ.3) GOTO 207
c===============================================================================
C zadavanje pritiska za desne strane, Fica, July 16, 2006
C      DO K=1,NDIM
c 	ivr=NGPSIL(4,JBRPS)
c       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(ivr),ivr,
c     &NTABFT,IIZLAZ)
cc       write(iizlaz,*)'eleme=',nbrel,'  sf3=  ',sf3,sf2,sf1
cc       write(iizlaz,*)'eleme=',nbrel,'  wdt=  ',wdt
c       RS2(K)=RS2(K)+(H(K)*FS2)*WDT*fk1
c       RS3(K)=RS3(K)+(H(K)*FS3)*WDT*fk1
C      ENDDO
c===============================================================================
      DO 206 K=1,NDIM 
Cc	 IF (NBREL.LT.900) GOTO 206
Cc       RS2(K)=RS2(K)+WDT*(H(K)*FS2)*NGPSIL(4,JBRPS)*AKT*(1.0D0)
C=====================================================================================
C Wada and Karino example, LDL transfer from blood to arterial wall, April 14, 2005
C       RS2(K)=RS2(K)+WDT*H(K)*20.0*(4.d-6-2.d-8)*TT21(2*NDIM+4+K)
C=====================================================================================
C=====================================================================================
C Tada and Tarbell example, IEL pores in medial mass transport, May 27, 2005
c       RS2(K)=RS2(K)+WDT*H(K)*FS2*1.D6*2.2D-6*TT21(2*NDIM+4+K)
c       RS2(K)=RS2(K)+WDT*H(K)*2.2D-6*TT21(2*NDIM+4+K)
C=====================================================================================
cC       RS3(K)=RS3(K)+WDT*(H(K)*FS3)*NGPSIL(5,JBRPS)*AKT
c       fs2=1.d0
       RS2(K)=RS2(K)+WDT*(H(K)*FS2)*NGPSIL(4,JBRPS)
       RS3(K)=RS3(K)+WDT*(H(K)*FS3)*NGPSIL(5,JBRPS)
 206  CONTINUE
 207   TTAU=TTAU+WDT*TAU
C 207   WRITE(IIZLAZ,*)'NBREL= ',NBREL,'TAU= ',TAU
 225  CONTINUE

C       IF(NGPSIL(4,JBRPS).EQ.3) THEN
C        WRITE(IIZLAZ,*)'ELEM',NBREL,'TTAU',TTAU
C       ENDIF


      IF (PENALT.LE.1.D0) IBRGT=IBRGT-1
C	NDIM=9
	nparam=0
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
C=========================================================================
C PAKOVANJE MATRICE C,Kk,Jk U MATRICU SKEF
C SKEF is local stiffness matrix

      DO 263 I=1,NDIM
      DO 262 J=1,NDIM
      SKEF(I,J)=AF*(AKMIV2(I,J)+AKVV2(I,J)+AJV2V2(I,J))
C      SKEF(I,J)=AF*(AKMIV2(I,J)+AKVV2(I,J))
      SKEF(I+NDIM,J)=AF*AJV3V2(I,J)
      SKEF(I+NDIM,J+NDIM)=AF*(AKMIV3(I,J)+AKVV2(I,J)+AJV3V3(I,J))
C      SKEF(I+NDIM,J+NDIM)=AF*(AKMIV3(I,J)+AKVV2(I,J))
      SKEF(I,J+NDIM)=AF*AJV2V3(I,J)
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=AF*(AKK(I,J)+AJK(I,J))
      IF (NSTAC.EQ.0) THEN
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=SKEF(I+2*NDIM+4,J+2*NDIM+4)+
     &C(I,J)/TIME
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       ENDIF
      IF (J.LE.4) THEN
       SKEF(J+2*NDIM,I)=AF*(AKV2P(I,J)+AKV2P1(I,J))
       SKEF(J+2*NDIM,I+NDIM)=AF*(AKV3P(I,J))
       SKEF(I,J+2*NDIM)=AF*(AKV2P(I,J))
       SKEF(I+NDIM,J+2*NDIM)=AF*(AKV3P(I,J))
      ENDIF
       SKEF(I+2*NDIM+4,J)=AF*(AKTV2(I,J))
       SKEF(I+2*NDIM+4,J+NDIM)=AF*(AKTV3(I,J))
C LEVA STRANA USLED KONVEKCIJE
C left side which relate to convective term
      SKEF(I,J+2*NDIM+4)=SKEF(I,J+2*NDIM+4)+AMV2(I,J)*FB2*BETA
      SKEF(I+NDIM,J+2*NDIM+4)=SKEF(I+NDIM,J+2*NDIM+4)+AMV2(I,J)*FB3*BETA
      IF (INDAX.EQ.1) THEN
       SKEF(I,J)=SKEF(I,J)+AF*AKUAX(I,J)
        IF (PENALT.GT.0.D0) THEN
         SKEF(I,J)=SKEF(I,J)+AF*APAXIX(I,J)
         SKEF(I+NDIM,J)=SKEF(I+NDIM,J)+AF*APAXIY(I,J)
        ENDIF
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AF*AKVAX(I,J)
      ENDIF
  262 CONTINUE
  263 CONTINUE
C=========================================================================
 
      IF (IPASS.NE.0.AND.METOD.EQ.4) THEN
      DO I=1,NDIM
	  II=I+NDIM
       DO J=1,NDIM
	   JJ=J+NDIM
        SKEF(I,J)=SKEF(I,J)+AKIIX(I,J)*ALPHAU/(1.D0-ALPHAU)
        SKEF(II,JJ)=SKEF(II,JJ)+AKIIY(I,J)*ALPHAV/(1.D0-ALPHAV)
	  F31(I)=F31(I)+AKIIX(I,J)*TT21(J)*ALPHAU/(1.D0-ALPHAU)
	  F31(II)=F31(II)+AKIIY(I,J)*TT21(JJ)*ALPHAV/(1.D0-ALPHAV)

	 ENDDO
	ENDDO
	ENDIF


      IF ((IPASS.EQ.0.OR.IPASS.EQ.4).AND.METOD.EQ.4) THEN
      DO I=1,NDIMP 
       DO J=1,NDIMP 
        SKEF(I+2*NDIM,J+2*NDIM)=AKPP(I,J)
	 ENDDO
	ENDDO
	ENDIF


C       QQ=0.D0
C       IF (NBREL.EQ.1240) QQ=1.D2
C       IF (NBREL.EQ.6760) QQ=-1.D2

      DO 290 I=1,NDIM
      

C ZAPREMINSKE SILE
C body forces
C       F31(I)=F31(I)+RB2(I)
C       F31(I+NDIM)=F31(I+NDIM)+RB3(I)

C POVRSINSKE SILE
C surface forces
C      IF (INDX.EQ.1)

      F31(I)=F31(I)+RS2(I)

c==========================================================
c samo za zadati flux od temperature......
C      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)+RS2(I)
c==========================================================


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


c      F31(I)=F31(I)-(AKMIV2(I,J)+AKVV2(I,J))*TT21(J)
      F31(I)=F31(I)-SKEF(I,J)*TT21P(J)



      IF (J.LE.4) THEN
c       F31(I)=F31(I)-AKV2P(I,J)*TT21(J+2*NDIM)
c       F31(I+NDIM)=F31(I+NDIM)-AKV3P(I,J)*TT21(J+2*NDIM)
       F31(I)=F31(I)+AKV2P(I,J)*TT21(J+2*NDIM)
       F31(I+NDIM)=F31(I+NDIM)+AKV3P(I,J)*TT21(J+2*NDIM)
      ENDIF
      IF (I.LE.4) THEN
	IF (METOD.EQ.4.AND.IPASS.NE.0) THEN
c      F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2P(J,I)*TT21(J)-
c     1AKV3P(J,I)*TT21(J+NDIM)-AKV2P1(J,I)*TT21(J)
	ENDIF
	IF (METOD.EQ.1) THEN
      F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2P(J,I)*TT21(J)-
     1AKV3P(J,I)*TT21(J+NDIM)-AKV2P1(J,I)*TT21(J)
	ENDIF
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
 
 
c      F31(I+NDIM)=F31(I+NDIM)-(AKMIV3(I,J)+AKVV2(I,J))*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-SKEF(I+NDIM,J+NDIM)*TT21P(J+NDIM)
 
 
      F31(I+NDIM)=F31(I+NDIM)-A21(I,J)*TT21(J)
      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)-AKTV2(I,J)*TT21(J)
     1-AKTV3(I,J)*TT21(J+NDIM)-AKK(I,J)*TT21(J+2*NDIM+4)
C DESNA STRANA USLED UTICAJA KONVEKCIJE
c right-hand side created by natural convection
      F31(I)=F31(I)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB2
      F31(I+NDIM)=F31(I+NDIM)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB3
c Correction Oct 21, 2004, N. Filipovic, only for Burger's non-viscous example
c         F31(I)=F31(I)-ATAYL(I,J)*TT210(J)
c         F31(I+NDIM)=F31(I+NDIM)-ATAYL(I,J)*TT210(J+NDIM)
 285  CONTINUE							      
 290  CONTINUE

C=========================================================================
C==========================================================================

C==========================================================================

C==========================================================================
C==========================================================================


C==========================================================================
C Poisson equation
C==========================================================================
	IF (IPASS.EQ.0) THEN
  
	DO J=1,NDIM
	 DO K=1,NDIM
	  RB2(J)=RB2(J)-(AKVV2(J,K)+AKMIV2(J,K))*TT21(K)
 	  RB3(J)=RB3(J)-(AKVV2(J,K)+AKMIV2(J,K))*TT21(K+NDIM)
C	  RB2(J)=RB2(J)-(AJV2V2(J,K)*TT21(K))-AJV2V3(J,K)*TT21(K+NDIM)
C	  RB3(J)=RB3(J)-(AJV3V2(J,K)*TT21(K))-AJV3V3(J,K)*TT21(K+NDIM)
	 ENDDO
	ENDDO

      DO I=1,NDIMP
	II=I+2*NDIM
	DO J=1,NDIM
        F31(II)=F31(II)-AKV2P(J,I)*(1.D0/AKIIX(J,J))*RB2(J)
        F31(II)=F31(II)-AKV3P(J,I)*(1.D0/AKIIY(J,J))*RB3(J)
c        F31(II)=F31(II)+AKV2P(J,I)*(1.D0/AKIIX(J,J))*RB2(J)
c        F31(II)=F31(II)+AKV3P(J,I)*(1.D0/AKIIY(J,J))*RB3(J)
	ENDDO
	ENDDO

      



      DO 311 I=1,NDIMP
      DO 301 J=1,NDIM

C       F31(I+2*NDIM)=F31(I+2*NDIM)-(AKVV2P(J,I)+AKMIV2P(J,I))
C     &*TT21(J)
C       F31(I+2*NDIM)=F31(I+2*NDIM)-(AKVV3P(J,I)+AKMIV3P(J,I))
C     &*TT21(J+NDIM)

       ALAMBDA=0.D0
C	 IF (NUMPASS.EQ.3) ALAMBDA=0.1D0
c       F31(I+2*NDIM)=F31(I+2*NDIM)+AKV2P(J,I)*TT21(J)*ALAMBDA
c       F31(I+2*NDIM)=F31(I+2*NDIM)+AKV3P(J,I)*TT21(J+NDIM)*ALAMBDA

C       F31(I+2*NDIM)=F31(I+2*NDIM)+1.D0
       IF (I.LE.NDIMP.AND.J.LE.NDIMP) THEN
        F31(I+2*NDIM)=F31(I+2*NDIM)-AKPP(I,J)*TT21P(J+2*NDIM)
	 ENDIF
C==========================================================================
 301  CONTINUE
C       F92(I+3*NDIM)=F92(I+3*NDIM)+RSVP(I)
 311  CONTINUE
      ENDIF
C==========================================================================

	IF (IPASS.EQ.4) THEN
  
      DO I=1,NDIMP
      DO J=1,NDIM
c       F31(I+2*NDIM)=F31(I+2*NDIM)-AKV2PT(J,I)*TT21(J)
c       F31(I+2*NDIM)=F31(I+2*NDIM)-AKV3PT(J,I)*TT21(J+NDIM)
       F31(I+2*NDIM)=F31(I+2*NDIM)+AKV2PT(J,I)*TT21(J)
       F31(I+2*NDIM)=F31(I+2*NDIM)+AKV3PT(J,I)*TT21(J+NDIM)
C        F31(I+2*NDIM)=F31(I+2*NDIM)-AKPP(I,J)*(TT21P(J+2*NDIM))

C==========================================================================
      ENDDO
c       F92(I+3*NDIM)=F92(I+3*NDIM)+RSVP(I)
      ENDDO
C	 PRIT(1,NBREL)=PRIT(1,NBREL)+F92(1+3*NDIM)/AKPP(1,1)/1.D-15
      ENDIF


	IF (IPASS.EQ.5) THEN
  
      DO I=1,NDIM
	NODE=NEL(I,NBREL)
      DO J=1,NDIMP
C	DO K=1,NDIM
	IF (LM2(I).NE.0) THEN
       GNODE(2,1,NODE)=GNODE(2,1,NODE)+
c     &(AKVV1P(I,J)+AKMIV1P(I,J))*TT21(J+4*NDIM)
     &(AKV2P(I,J))*(1.D0/AKIIX(I,I))*TT21(J+2*NDIM+4)
	ENDIF
	IF (LM2(I+NDIM).NE.0) THEN
       GNODE(2,2,NODE)=GNODE(2,2,NODE)+
c     &(AKVV2P(I,J)+AKMIV2P(I,J))*TT21(J+4*NDIM)
     &(AKV3P(I,J))*(1.D0/AKIIY(I,I))*TT21(J+2*NDIM+4)
	ENDIF
     
C==========================================================================
C      ENDDO
      ENDDO
      ENDDO
      ENDIF
C==========================================================================



C============================================================================
C  FOR FLUX CALCULATION
C      IF(IPROL.EQ.1) THEN
C       DO I=1,NDIM
C        N=NEL(I,NBREL)
C         POT=0.D0
C         DO J=1,NDES
C          POT=POT+SKEF(I+2*NDIM+4,J)*TT21(J)
C         ENDDO
C=========================================
C OVDE PRIVREMENO UKINUTO RACUNANJE FLUXA
C         PRES(1,N)=PRES(1,N)+POT
C=========================================
C       ENDDO
C============================================================================
C SHEAR STRESS CALCULATION FOR 2D PROBLEM
c      CALL SSTRE2(NEL,ID,NPARAM,TT21,H,HP,X,Y,R,S,W,ZVHX,ZVHY,
c     &PRES,INDEL,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,
c     &INDAX,AKT,GUSM,IUPWIN,NBREL,IIZLAZ,IBRGT)
C============================================================================
C============================================================================
C============================================================================
C adding of pressure boundary condition at outlet of the tube, non-zero traction boundary condition
c      if (nbrel.ge.121) then

c       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(1),1,
c     &NTABFT,IIZLAZ)

c	 do i=1,ndim
c	  node=nel(i,nbrel)
c	  if (cord(1,node).eq.20.d0) then
C         f31(i)=f31(i)-1.d6*fk1
c         f31(i)=f31(i)-0.d0*fk1
c	  endif
c	 enddo
c	endif
C============================================================================
C============================================================================
C============================================================================
C    
C
c      CALL WRRF(SKEF,NDES*NDES,'SKEF ',IIZLAZ)      
c	stop
C      CALL WRRF(F31,NDES,'F31= ',IIZLAZ)      


c=================================================================================================
C Approximation, making global symmetry matrix from unsymmetrix by averaged out of diagonal terms
	IF (ISYMMS.EQ.1)THEN
	   DO I=1,NDES
	     DO J=1,NDES
	      SKEF(I,J)=0.5D0*(SKEF(I,J)+SKEF(J,I))  
	   ENDDO
	  ENDDO

C       DO I=1,NDES*(NDES+1)*0.5
C        SKE(I)=0.D0
C       ENDDO
       CALL PSKEFN(SKEF,SKE,NDES)	  
	ENDIF
c=================================================================================================

c      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
c	  CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,2)
c	  GOTO 400
c      ENDIF


	IF(NJUTRA.EQ.1.AND.ITER.GT.0) THEN
       CALL SPAKDE (SILE,F31,LM2,NDES)
	ELSE
       IF (ISYMMS.EQ.1) THEN
	    CALL MATSTE (ALEVO,MAXA,SILE,SKE,F31,LM2,NDES,1)
	 ELSE
          CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F31,MAXA,LM2,NDES,1)
	 ENDIF
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
C
C       FACAXY=1.D0
C       IF(INDAX.EQ.1) FACAXY=8.D0*DATAN(1.D0)
C
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
c      write(iizlaz,*) jedn
c     do i=1,jedn
c	 write(iizlaz,*)'sile= ',sile(i)
c	enddo 
c	stop

c      CALL MUMPSRIGHT(SILE,JEDN)
	


C End of subroutine
      END
C=======================================================================
