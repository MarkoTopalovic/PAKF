C=========================================================================
      SUBROUTINE AtheroSun2D(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,AF,INDEL,NJUTRA,ISYMMS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine RACU2D is used for 2D analysis
CE It is used global loop per elements
C

C SAMO PRIVREMENO UBACEN COMMON ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,DPRES
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)
	COMMON /MUMPS/ IMUMPS, MUFILE,MUFILE2
      COMMON /ATHERO/ HM,ZVXM,ZVYM,HOX
      COMMON /ATHERO1/ IATHERO,AK1,GAMA,O_THRES,ALAMBDA,DD,AKleg,Rw,ALp,
     &Dpw,Pp,SigmaF,PlaqMove
     

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

      CHARACTER*80 NASLOV
      DIMENSION GNODE(2,11,*),ALEVO(*),DESNO(*),SILE(*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(11,*),NGPSIL(8,*),MAXA(*)
      DIMENSION SKEF(NDES,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
      DIMENSION SPSIL(2,*),PRES(3,*),VMESH(3,*)

      DIMENSION X(9),Y(9),LM2(9*9),TT21(9*9),TT210(9*9),TT21A(9*9)
      DIMENSION R(3,3),S(3,3),W(3,3)
C      DIMENSION R1(2),S1(2),W1(2)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION AKVV2(9,9)
      DIMENSION AKMIV2(9,9)
      DIMENSION AKMIV3(9,9)
      DIMENSION A12(9,9)
      DIMENSION A21(9,9)
      DIMENSION AKK(9,9),AKKD2(9,9)
      DIMENSION AMV2(9,9)
      DIMENSION AK1Q(9,9)
      DIMENSION C(9,9)
      DIMENSION AJV2V2(9,9)
      DIMENSION AJV3V3(9,9)
      DIMENSION AJV2V3(9,9)
      DIMENSION AJV3V2(9,9)
      DIMENSION AKTV2(9,9)
      DIMENSION AKTV3(9,9)
      DIMENSION AKOXV2(9,9)
      DIMENSION AKOXV3(9,9)

      DIMENSION AKDIVV(9,9)
      DIMENSION AKMV2(9,9)
      DIMENSION AKMV3(9,9)
      DIMENSION AKTMV2(9,9)
      DIMENSION AKTMV3(9,9)
      DIMENSION AFIFI(9,9)

      DIMENSION AK1OX(9,9)
      DIMENSION AK1M(9,9)


      DIMENSION AKV2P(9,4)
      DIMENSION AKV3P(9,4),AKV2P1(9,4)
      DIMENSION AJK(9,9),AJKLOG(9,9)
      DIMENSION RS2(9)
      DIMENSION RS3(9)
      DIMENSION RB2(9)
      DIMENSION FSS(9)
      DIMENSION FTSS(9),FTLDL(9),FFIFI(9)
      DIMENSION FK1OXM(9),FRWCW(9),FUWX(9),FUWY(9)
      DIMENSION FGAMA(9)
      DIMENSION RB3(9)
      DIMENSION F31(81)
      DIMENSION AKUAX(9,9),AKVAX(9,9),APAXIX(9,9),APAXIY(9,9)

C	DIMENSION SKE(31*16)
	DIMENSION SKE(31*16)

C Coorection for UPWIND, stability, Taylor, Oct, 21, 2004 Implemented
	DIMENSION ATAYL(9,9)
	DIMENSION DD(5)




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
        

C===========================================================================
C Plaque growing
C===========================================================================
      AKONST=1.D-18/AMI
C      AK1=1.D0
C     GAMA=1.D0
C	O_THRES=0.0D0
C	ALAMBDA=1.D-5
	AKT=DD(1)
C	DD(2)=1.D-3
C	DD(3)=1.D-5
C	DD(4)=1.D-3
C	DD(5)=AMI
C	DD(6)=AMI
c	DD(2)=1.D0
c	DD(3)=1.D0
c	DD(4)=1.D0

C=======================================================================
C loop for shear stress calculation
      goto 159
C=======================================================================
      do nbrel=1,net

      DO 126 I=1,4*NDIM
      TT210(I)=0.D0
      LM2(I)=0
 126  TT21(I)=0.D0
C=========================================================================
C=========================================================================
      DO 131 KLM=1,NDIM
      X(KLM)=CORD(1,NEL(KLM,NBREL))
      Y(KLM)=CORD(2,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      IF (KLM.LE.4) LM2(KLM+2*NDIM)=ID(4,NEL(KLM,NBREL))
	DO J=2,6
        LM2(KLM+J*NDIM+4)=ID(J+3,NEL(KLM,NBREL))
	ENDDO
 131  CONTINUE
C=======================================================================
      DO 141 KLM=1,NDIM
      DO 136 NR=1,2
C      IF (ID(NR,NEL(KLM,NBREL)) .NE. 0) THEN
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
C      ENDIF
 136  CONTINUE
      IF (KLM.LE.4) THEN
C       TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
C       TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
	DO J=2,6
       TT21(KLM+J*NDIM+4)=GNODE(2,J+3,NEL(KLM,NBREL))
       TT210(KLM+J*NDIM+4)=GNODE(1,J+3,NEL(KLM,NBREL))
	ENDDO
 141  CONTINUE
C
      DO I=1,7*NDIM+4
       IF(KKORAK.EQ.1.AND.ITER.EQ.1) TT210(I)=TT21(I)
      ENDDO

C============================================================================
C SHEAR STRESS CALCULATION FOR 2D PROBLEM
      CALL SSTRE2(NEL,ID,NPARAM,TT21,H,HP,X,Y,R,S,W,ZVHX,ZVHY,
     &PRES,INDEL,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,
     &INDAX,AKT,GUSM,IUPWIN,NBREL,IIZLAZ,IBRGT)
C============================================================================
      enddo



c mesh moving because plaque developing
c      CALL PLAQGROW(CORD,GNODE,NPT,TIME,AK1,NEL,NDIM) 

c      if (kkorak.eq.1.and. iter.eq.0) 
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,2)

159   continue
      DO 400 NBREL=1,NET
     

	MAT=NEL(NDIM+1,NBREL)
C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,4*NDIM
      TT210(I)=0.D0
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
	DO J=2,6
        LM2(KLM+J*NDIM+4)=ID(J+3,NEL(KLM,NBREL))
	ENDDO
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
C       TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
C       TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
	DO J=2,6
       TT21(KLM+J*NDIM+4)=GNODE(2,J+3,NEL(KLM,NBREL))
       TT210(KLM+J*NDIM+4)=GNODE(1,J+3,NEL(KLM,NBREL))
	ENDDO
 140  CONTINUE
C
C=======================================================================
C      IF (NDIM.GT.4) THEN
C      DO 145 NR=28,36
C      TT210(NR-5)=TT210(NR)
C 145  TT21(NR-5)=TT21(NR)
C      ENDIF
C=======================================================================
      DO I=1,7*NDIM+4
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
      ENDIF
      AKVV2(K,N)=0.D0
      APAXIX(K,N)=0.D0
      APAXIY(K,N)=0.D0
      AKUAX(K,N)=0.D0
      AKVAX(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0

      AKDIVV(K,N)=0.D0
      AKMV2(K,N)=0.D0
      AKMV3(K,N)=0.D0
      AKTMV2(K,N)=0.D0
      AKTMV3(K,N)=0.D0

      AK1OX(K,N)=0.D0
      AK1M(K,N)=0.D0
      AKOXV2(K,N)=0.D0
      AKOXV3(K,N)=0.D0

      AKMIV2(K,N)=0.D0
      AKMIV3(K,N)=0.D0
      A12(K,N)=0.D0
      A21(K,N)=0.D0
      AMV2(K,N)=0.D0
      AK1Q(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V2(K,N)=0.D0
      AKK(K,N)=0.D0
      AFIFI(K,N)=0.D0
      AKKD2(K,N)=0.D0
      AJK(K,N)=0.D0
      AJKLOG(K,N)=0.D0
c	ATAYL(K,N)=0.D0
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
C SURFACE FORCES AND BODY FORCES
       RS2(K)=0.D0
       RS3(K)=0.D0
       RB2(K)=0.D0
       FSS(K)=0.D0
       FTSS(K)=0.D0
       FTLDL(K)=0.D0
       FFIFI(K)=0.D0
       FK1OXM(K)=0.D0
       FRWCW(K)=0.D0
       FUWX(K)=0.D0
       FUWY(K)=0.D0
       FGAMA(K)=0.D0
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

      IF (MAT.EQ.1) THEN
       AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVXT)
       AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVYT)
	ENDIF

      ZVXOX=DOT(ZVHX,TT21(1+3*NDIM+4),NDIM)
      ZVYOX=DOT(ZVHY,TT21(1+3*NDIM+4),NDIM)
	AKLOG=AK1

      IF (MAT.EQ.2) THEN
       AKOXV2(K,N)=AKOXV2(K,N)+WDT*H(K)*H(N)*ZVXOX*AKLOG
       AKOXV3(K,N)=AKOXV3(K,N)+WDT*H(K)*H(N)*ZVYOX*AKLOG
      AFIFI(K,N)=AFIFI(K,N)+AKONST*(ZVHX(K)*ZVHX(N)+ZVHY(K)*ZVHY(N))*WDT
	ENDIF


      DIVV=ZVXV2+ZVYV3
C      AKDIVV(K,N)=AKDIVV(K,N)+WDT*(H(K)*DIVV*H(N))
C      AKTMV2(K,N)=AKTMV2(K,N)+WDT*(H(K)*HM*ZVHX(N))
C      AKTMV3(K,N)=AKTMV3(K,N)+WDT*(H(K)*HM*ZVHY(N))
C      AKMV2(K,N)=AKMV2(K,N)+WDT*(H(K)*ZVXM*H(N))
C      AKMV3(K,N)=AKMV3(K,N)+WDT*(H(K)*ZVYM*H(N))
     
c      AK1OX(K,N)=AK1OX(K,N)+AK1*(H(K)*HOX*H(N))*WDT
c      AK1M(K,N)=AK1M(K,N)+AK1*(H(K)*HM*H(N))*WDT



      AKMIV2(K,N)=AKMIV2(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI)

	IF (MAT.EQ.1) THEN
      AKK(K,N)=AKK(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N)))*DD(1)
      IF (INDAX.EQ.1) THEN
       AKK(K,N)=AKK(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*DD(1)
      ENDIF
	ENDIF

      IF (MAT.EQ.2) THEN
      AKKD2(K,N)=AKKD2(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N)))*DD(2)
      IF (INDAX.EQ.1) THEN
       AKKD2(K,N)=AKKD2(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*DD(2)
      ENDIF
      ENDIF

      AKMIV3(K,N)=AKMIV2(K,N)
      AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM
  
C===========================================================================
C Plaque growing
C===========================================================================
C      HQOX=DOT(H,TT21(2*NDIM+4+1),NDIM)
C      AK1Q(K,N)=AK1Q(K,N)+WDT*H(K)*H(N)*HQOX*AK1
C===========================================================================
  
      AJV2V2(K,N)=AJV2V2(K,N)+WDT*(H(K)*ZVXV2*H(N)*GUSM)
      AJV3V3(K,N)=AJV3V3(K,N)+WDT*(H(K)*ZVYV3*H(N)*GUSM)
      AJV2V3(K,N)=AJV2V3(K,N)+WDT*(H(K)*ZVYV2*H(N)*GUSM)
      AJV3V2(K,N)=AJV3V2(K,N)+WDT*(H(K)*ZVXV3*H(N)*GUSM)
      IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
c      ATAYL(K,N)=ATAYL(K,N)+WDT*(ZVHX(K)*ZVHX(N)*HV2**2+
c     1ZVHY(K)*ZVHY(N)*HV3**2)*GUSM*TIME/2.D0
c      ATAYL(K,N)=ATAYL(K,N)+WDT*(ZVHX(K)*ZVHY(N)*HV2*HV3+
c     1ZVHY(K)*ZVHX(N)*HV3*HV2)*GUSM*TIME/2.D0
  164 CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
C===========================================================================
C Plaque growing
C===========================================================================
	 OOX=TT21(K+3*NDIM+4)
C	 AMM=TT21(K+4*NDIM+4)
C	 SS=TT21(K+5*NDIM+4)
	 NN=NEL(K,NBREL)
	 TSS=DSQRT(PRES(1,NN)**2+PRES(2,NN)**2+PRES(3,NN)**2)
C       FSS(K)=FSS(K)+H(K)*(SS/(1.D0+SS))*WDT
Cc       FTSS(K)=FTSS(K)+H(K)*(TSS*OOX)*WDT
Cc       OD=1.39*TSS**(-0.118D0)
Cc	IF (DABS(TSS).LT.1.D-10) OD=1.D0
       tauOx=1.d-4
Cc	 if (kkorak.gt.1) tauOx=OOX
C       FTSS(K)=FTSS(K)+H(K)*(TSS*tauOx)*WDT
Cc       FTSS(K)=FTSS(K)+H(K)*(OD)*WDT
C       FK1OXM(K)=FK1OXM(K)+H(K)*(AK1*OOX*AMM)*WDT
	 RW=GAMA
c	 RW=1.D-6
       
	 FRWCW(K)=FRWCW(K)+H(K)*(RW*OOX)*WDT
	 
	 PRESSURE=GNODE(2,4,NEL(K,NBREL))-9.331D3
C	 PRESSURE=PRESS-9.331D3
	 AKAPA=1.D-15
C      FUWX(K)=FUWX(K)+ZVHX(K)*(AKAPA*PRESSURE)*WDT
C       FUWY(K)=FUWY(K)+ZVHY(K)*(AKAPA*PRESSURE)*WDT
C       FGAMA(K)=FGAMA(K)+H(K)*(GAMA*(OOX-O_THRES))*WDT
C===========================================================================
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
C===========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PENALTY FAKTOROM
C reduced integration if penalty function is defined

C      IF (PENALT.GT.0.D0) THEN
      IF (PENALT.GT.0.D0.AND.MAT.EQ.1) THEN
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
       AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
       AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
C      ENDIF
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
c	goto 251

      DO 250 JBRPS=1,MAXSIL
      IF (NBREL.EQ.NGPSIL(1,JBRPS)) THEN
C       WRITE(IIZLAZ,*)'IND,ELEM',NGPSIL(4,JBRPS),NBREL
        IF (PENALT.LE.1.D0) IBRGT=IBRGT+1
C        INDX=NGPSIL(6,JBRPS)
C        INDY=NGPSIL(7,JBRPS)

C        IF (NGPSIL(4,JBRPS).EQ.0) THEN
        INDX=1
	  
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
c       IF(NGPSIL(4,JBRPS).EQ.3) NPARAM=2
       
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
c===============================================================================
      DO 206 K=1,NDIM 
 	ivr=NGPSIL(4,JBRPS)
	if (ivr.eq.4) then
       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(ivr),ivr,
     &NTABFT,IIZLAZ)
cc       write(iizlaz,*)'eleme=',nbrel,'  sf3=  ',sf3,sf2,sf1
cc       write(iizlaz,*)'eleme=',nbrel,'  wdt=  ',wdt
c       fk1=9331.d1
       RS2(K)=RS2(K)+(H(K)*FS2)*WDT*fk1
       RS3(K)=RS3(K)+(H(K)*FS3)*WDT*fk1
	 else
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
C===========================================================================
C===========================================================================
C===========================================================================
C Plaque growing
C===========================================================================
	 ALDL=TT21(K+2*NDIM+4)
	 OOX=TT21(K+3*NDIM+4)
C	 AMM=TT21(K+4*NDIM+4)
	 SS=TT21(K+5*NDIM+4)
	 NN=NEL(K,NBREL)
	 TSS=DSQRT(PRES(1,NN)**2+PRES(2,NN)**2+PRES(3,NN)**2)
       FSS(K)=FSS(K)+H(K)*(SS/(1.D0+SS))*WDT
c       FTSS(K)=FTSS(K)+H(K)*(TSS*OOX)*WDT
       OD=1.39*TSS**(-0.118D0)
	IF (DABS(TSS).LT.1.D-10) OD=1.D0
       tauOx=1.d-6
	AJs=tauOx*1.d-6*TSS
	PP=DD(3)
      SIGMAF=DD(4)
	ALP=3.D-12
      ALP=O_THRES


      PRESSURE=GNODE(2,4,NEL(K,NBREL))
      Pw=GNODE(2,7,NEL(K,NBREL))

C	PRESSURE=9531.D0
C      DPW=PRESSURE-ALAMBDA
C      DPW=-DPW
      DPW=PRESSURE-Pw
      AJv=ALP*DPW
c      write(iizlaz,*)'AJv= ',AJv
C      AJv=0.D0

      AJs=PP*(ALDL-OOX)+(1.D0-SIGMAF)*AJv*(ALDL+OOX)*0.50D0
C      AJs=PP*(ALDL-OOX)
C      AJs=0.D0
C	AJs=AJS
c       tauOx=1.d5
c	 TSS=4.8D-2-TSS
c	 if (kkorak.gt.1) tauOx=OOX


       FTSS(K)=FTSS(K)+H(K)*(AJs)*WDT
       FTLDL(K)=FTLDL(K)+H(K)*(AJs)*WDT
       FFIFI(K)=FFIFI(K)+H(K)*(AJv)*WDT
c===============================================================
c for nikola only
C	 TTEMP=TT21(K+2*NDIM+4)
C       FTSS(K)=FTSS(K)+H(K)*(ak1)*WDT
c===============================================================



c       FTSS(K)=FTSS(K)+H(K)*(TSS)*WDT
C       FTSS(K)=FTSS(K)+H(K)*(OD*tauOx)*WDT
c       FK1OXM(K)=FK1OXM(K)+H(K)*(AK1*OOX*AMM)*WDT
c       FGAMA(K)=FGAMA(K)+H(K)*(GAMA*(OOX-O_THRES))*WDT
C===========================================================================
C===========================================================================
C===========================================================================

c       RS2(K)=RS2(K)+WDT*(H(K)*FS2)*NGPSIL(4,JBRPS)
c       RS3(K)=RS3(K)+WDT*(H(K)*FS3)*NGPSIL(5,JBRPS)
        endif
 206  CONTINUE
 207   TTAU=TTAU+WDT*TAU
C 207   WRITE(IIZLAZ,*)'NBREL= ',NBREL,'TAU= ',TAU
 225  CONTINUE

C       IF(NGPSIL(4,JBRPS).EQ.3) THEN
C        WRITE(IIZLAZ,*)'ELEM',NBREL,'TTAU',TTAU
C       ENDIF


      IF (PENALT.LE.1.D0) IBRGT=IBRGT-1
	nparam=0
      ENDIF
 250  CONTINUE

c 251  CONTINUE
C C matrix include heat conduction


      DO 255 I=1,NDIM
      DO 254 J=1,NDIM

      C(I,J)=AMV2(I,J)*CC/gusm              
      if (MAT.EQ.1) then
	  AJK(I,J)=AKVV2(I,J)*CC/gusm
	endif

      if (MAT.EQ.2) then
         AJKLOG(I,J)=AKVV2(I,J)*AKLOG*CC/gusm
	ENDIF

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
      SKEF(I+NDIM,J)=AF*AJV3V2(I,J)
      SKEF(I+NDIM,J+NDIM)=AF*(AKMIV3(I,J)+AKVV2(I,J)+AJV3V3(I,J))
      SKEF(I,J+NDIM)=AF*AJV2V3(I,J)
C=========================================================================      
C Plaque growing   
C=========================================================================   
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=AF*(AKK(I,J)+AJK(I,J))
	SKEF(I+3*NDIM+4,J+3*NDIM+4)=AF*(AKKD2(I,J)+AJKLOG(I,J))
cc	SKEF(I+4*NDIM+4,J+4*NDIM+4)=AF*(AKK(I,J)*DD(3))
cc	SKEF(I+5*NDIM+4,J+5*NDIM+4)=AF*(AKK(I,J)*DD(4))
c     SKEF(I+6*NDIM+4,J+6*NDIM+4)=AF*(AKK(I,J)*DD(5))
C=========================================================================      
	SKEF(I+3*NDIM+4,J+3*NDIM+4)=SKEF(I+3*NDIM+4,J+3*NDIM+4)+
     &AF*(RW*AMV2(I,J)/GUSM)
C=========================================================================      
C Plaque growing   
C=========================================================================   
C	SKEF(I+2*NDIM+4,J+3*NDIM+4)=AK1Q(I,J)
C	SKEF(I+3*NDIM+4,J+3*NDIM+4)=AK1Q(I,J)
C	SKEF(I+4*NDIM+4,J+3*NDIM+4)=-AK1Q(I,J)
C=========================================================================      
C=========================================================================      
C Plaque growing   
C=========================================================================   
c	IF (MAT.EQ.2) THEN
c	SKEF(I+5*NDIM+4,J+5*NDIM+4)=SKEF(I+5*NDIM+4,J+5*NDIM+4)+
c     &ALAMBDA*AMV2(I,J)
c	ENDIF
C=========================================================================      


      IF (NSTAC.EQ.0) THEN
C=========================================================================   
C Plaque growing   
C=========================================================================   
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=SKEF(I+2*NDIM+4,J+2*NDIM+4)+
     &C(I,J)/TIME
c	 IF(MAT.EQ.2) THEN
      SKEF(I+3*NDIM+4,J+3*NDIM+4)=SKEF(I+3*NDIM+4,J+3*NDIM+4)+
     &C(I,J)/TIME
c      SKEF(I+4*NDIM+4,J+4*NDIM+4)=SKEF(I+4*NDIM+4,J+4*NDIM+4)+
c     &C(I,J)/TIME
c      SKEF(I+5*NDIM+4,J+5*NDIM+4)=SKEF(I+5*NDIM+4,J+5*NDIM+4)+
c     &C(I,J)/TIME
c      SKEF(I+6*NDIM+4,J+6*NDIM+4)=SKEF(I+6*NDIM+4,J+6*NDIM+4)+
c     &C(I,J)/TIME
c      ENDIF
C=========================================================================      
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       ENDIF
      IF (J.LE.4) THEN
       SKEF(J+2*NDIM,I)=AF*(AKV2P(I,J)+AKV2P1(I,J))
       SKEF(J+2*NDIM,I+NDIM)=AF*(AKV3P(I,J))
       SKEF(I,J+2*NDIM)=AF*(AKV2P(I,J))
       SKEF(I+NDIM,J+2*NDIM)=AF*(AKV3P(I,J))
      ENDIF

C=========================================================================      
       SKEF(I+2*NDIM+4,J)=AF*(AKTV2(I,J))
       SKEF(I+2*NDIM+4,J+NDIM)=AF*(AKTV3(I,J))
C=========================================================================      

C=========================================================================      
C ATHEROSCLEROSIS div(v*M)
C=========================================================================  
C=========================================================================      
C       SKEF(I+3*NDIM+4,J+4*NDIM+4)=AF*(AKOXV2(I,J))
C       SKEF(I+3*NDIM+4,J+5*NDIM+4)=AF*(AKOXV3(I,J))
C=========================================================================      
      SKEF(I+4*NDIM+4,J+4*NDIM+4)=SKEF(I+4*NDIM+4,J+4*NDIM+4)+
     &AF*AFIFI(I,J)

C       SKEF(I+4*NDIM+4,J)=SKEF(I+4*NDIM+4,J)+AF*(AMV2(I,J)/gusm)
C       SKEF(I+5*NDIM+4,J)=SKEF(I+5*NDIM+4,J)+AF*(AMV2(I,J)/gusm)

 
c       SKEF(I+4*NDIM+4,J)=SKEF(I+4*NDIM+4,J)+AF*(AKMV2(I,J))
c       SKEF(I+4*NDIM+4,J+NDIM)=SKEF(I+4*NDIM+4,J+NDIM)+AF*(AKMV3(I,J))
c       SKEF(I+4*NDIM+4,J+4*NDIM+4)=SKEF(I+4*NDIM+4,J+4*NDIM+4)+
c     &AF*(AJK(I,J))
c       SKEF(I+4*NDIM+4,J)=SKEF(I+4*NDIM+4,J)+AF*(AKTMV2(I,J))
c       SKEF(I+4*NDIM+4,J+NDIM)=SKEF(I+4*NDIM+4,J+NDIM)+AF*(AKTMV3(I,J))
c       SKEF(I+4*NDIM+4,J+4*NDIM+4)=SKEF(I+4*NDIM+4,J+4*NDIM+4)+
c     &AF*(AKDIVV(I,J))
C=========================================================================      


C=========================================================================      
C ATHEROSCLEROSIS k1*Ox*M
C=========================================================================      
c       SKEF(I+3*NDIM+4,J+4*NDIM+4)=SKEF(I+3*NDIM+4,J+4*NDIM+4)+
c     &AF*(AK1OX(I,J))
c       SKEF(I+3*NDIM+4,J+3*NDIM+4)=SKEF(I+3*NDIM+4,J+3*NDIM+4)+
c     &AF*(AK1M(I,J))
       
c       SKEF(I+4*NDIM+4,J+4*NDIM+4)=SKEF(I+4*NDIM+4,J+4*NDIM+4)+
c     &AF*(AK1OX(I,J))
c       SKEF(I+4*NDIM+4,J+3*NDIM+4)= SKEF(I+4*NDIM+4,J+3*NDIM+4)+
c     &AF*(AK1M(I,J))

c       SKEF(I+5*NDIM+4,J+4*NDIM+4)=SKEF(I+5*NDIM+4,J+4*NDIM+4)-
c     &AF*(AK1OX(I,J))
c       SKEF(I+5*NDIM+4,J+3*NDIM+4)=SKEF(I+5*NDIM+4,J+3*NDIM+4)-
c     &AF*(AK1M(I,J))
c	ENDIF
C=========================================================================      



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


C       QQ=0.D0
C       IF (NBREL.EQ.1240) QQ=1.D2
C       IF (NBREL.EQ.6760) QQ=-1.D2

C      DO I=1,NDIM
C       write(iizlaz,*)'element fk1oxm= ', nbrel,fk1oxm(i)
C	enddo

      DO 290 I=1,NDIM
C=========================================================================   
C Plaque growing   
C=========================================================================   
c      if (mod(nbrel-10,20).eq.0.and.(i.eq.3.or.i.eq.4)) then
c      if (indx.eq.1.and.(i.eq.3.or.i.eq.4)) then
C      if (indx.eq.1) then
c	 write(iizlaz,*) ftss(i),FSS(I),FGAMA(I)
c       if (mod(kkorak,3).eq.1) then
        F31(I+2*NDIM+4)=F31(I+2*NDIM+4)+FTLDL(I)
        F31(I+3*NDIM+4)=F31(I+3*NDIM+4)+FTSS(I)
        F31(I+4*NDIM+4)=F31(I+4*NDIM+4)+FFIFI(I)
c	 endif
c       if (mod(kkorak,3).eq.0) then
c        F31(I+4*NDIM+4)=F31(I+4*NDIM+4)+FSS(I)
c	 endif
c       if (mod(kkorak,3).eq.2) then
c        F31(I+5*NDIM+4)=F31(I+5*NDIM+4)+FGAMA(I)
c	 endif
C	endif

c       IF (MAT.EQ.2) THEN
c       F31(I+3*NDIM+4)=F31(I+3*NDIM+4)-FK1OXM(I)
c       F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-FK1OXM(I)
c       F31(I+5*NDIM+4)=F31(I+5*NDIM+4)+FK1OXM(I)
c       F31(I+5*NDIM+4)=F31(I+5*NDIM+4)-ALAMBDA*TT21(I+5*NDIM+4)
c	ENDIF



       F31(I+3*NDIM+4)=F31(I+3*NDIM+4)-FRWCW(I)
C       F31(I+3*NDIM+4)=F31(I+3*NDIM+4)-FRWCW(I)

      
C ZAPREMINSKE SILE
C body forces
       F31(I)=F31(I)+RB2(I)
       F31(I+NDIM)=F31(I+NDIM)+RB3(I)

c==========================================================
c==========================================================
C POVRSINSKE SILE
C surface forces
c==========================================================
C      IF (INDX.EQ.1)
      F31(I)=F31(I)+RS2(I)
c==========================================================
c samo za zadati flux od temperature......
C      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)+RS2(I)
c==========================================================
C      IF (INDY.EQ.1)
      F31(I+NDIM)=F31(I+NDIM)+RS3(I)
c==========================================================
c==========================================================

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
C=========================================================================   
C Plaque growing   
C=========================================================================   
c	 DO KK=2,5
 	 DO KK=2,3
        II=I+KK*NDIM+4
        JJ=J+KK*NDIM+4
        F31(II)=F31(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
	 ENDDO
      ENDIF
C=========================================================================   


C=========================================================================   
C Plaque growing   
C=========================================================================   
	 DO KK=3,3
c	 DO KK=2,5
        II=I+KK*NDIM+4
        JJ=J+KK*NDIM+4
        F31(II)=F31(II)-AKKD2(I,J)*TT21(JJ)
	 ENDDO
	 DO KK=2,2
c	 DO KK=2,5
        II=I+KK*NDIM+4
        JJ=J+KK*NDIM+4
        F31(II)=F31(II)-AKK(I,J)*TT21(JJ)
	 ENDDO
C=========================================================================   


      F31(I)=F31(I)-A12(I,J)*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-(AKMIV3(I,J)+AKVV2(I,J))*TT21(J+NDIM)
      F31(I+NDIM)=F31(I+NDIM)-A21(I,J)*TT21(J)

      F31(I+2*NDIM+4)=F31(I+2*NDIM+4)-AKTV2(I,J)*TT21(J)
     1-AKTV3(I,J)*TT21(J+NDIM)
	


C=========================================================================   
C Plaque growing   
C========================================================================= 
c      IF (MAT.EQ.2) THEN  
C      F31(I+3*NDIM+4)=F31(I+3*NDIM+4)-AKOXV2(I,J)*TT21(J+4*NDIM+4)
C     1-AKOXV3(I,J)*TT21(J+5*NDIM+4)


C      F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-AMV2(I,J)*TT21(J+4*NDIM+4)
C      F31(I+5*NDIM+4)=F31(I+5*NDIM+4)-AMV2(I,J)*TT21(J+5*NDIM+4)
     
      F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-AFIFI(I,J)*TT21(J+4*NDIM+4)

C      F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-FUWX(I)
C      F31(I+5*NDIM+4)=F31(I+5*NDIM+4)-FUWY(I)



c      F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-AKMV2(I,J)*TT21(J)
c     1-AKMV3(I,J)*TT21(J+NDIM)
c      F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-AKTMV2(I,J)*TT21(J)
c     1-AKTMV3(I,J)*TT21(J+NDIM)
c      F31(I+4*NDIM+4)=F31(I+4*NDIM+4)-AKDIVV(I,J)*TT21(J+4*NDIM+4)
c	ENDIF
     



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
      CALL SSTRE2(NEL,ID,NPARAM,TT21,H,HP,X,Y,R,S,W,ZVHX,ZVHY,
     &PRES,INDEL,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,
     &INDAX,AKT,GUSM,IUPWIN,NBREL,IIZLAZ,IBRGT)
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
         
c      CALL MUMPSRIGHT(SILE,JEDN)
	

c       if (iter.eq.2) stop
C End of subroutine
      END
C=======================================================================
C=========================================================================
