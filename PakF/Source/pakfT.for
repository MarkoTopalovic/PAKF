C=========================================================================
C MODIFIKOVAO I KREIRAO ALEKSANDAR NIKOLIC, AVGUST 2016-FEBRUAR 2017.
C PROGRAM ZA TURBULENTNO STRUJANJE FLUIDA, K-OMEGA TURBULENTNI MODEL, 2D
C UBACENA SUBROUTINA U FICIN SOLVER
C=========================================================================
      SUBROUTINE RACU2DT(GNODE,NEL,ID,CORD,SKEF,TABF,ITFMAX,
     &NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,TETAO,FB2,FB3,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,INDAMI,BETA,NDES,IDPRIT,
     &IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,NPER,NTABFT,PRES,VMESH,
     &IALE,AF,NJUTRA,ISYMMS,DELTAL,PRES1,TAU,VOSI,DT,IBKOR)
     
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)

CE Subroutine RACU2DT is used for 2D analysis in turbulent mode
CE It is used global loop per elements

	COMMON /MUMPS_PAK/ IMUMPS, MUFILE, MUFILE2
	COMMON /TRANSP/ INDFL,LID1
	COMMON /TURB/ ITURB
      COMMON /POCK/ POCK
      COMMON /POCO/ POCO
      COMMON /IZID/ IZID
	
      CHARACTER*250 NASLOVF
      DIMENSION GNODE(2,7,*),ALEVO(*),DESNO(*),SILE(*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(7,*),NGPSIL(8,*),MAXA(*)
      DIMENSION SKEF(NDES,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
      DIMENSION SPSIL(2,*),PRES(3,*),VMESH(3,*)
      DIMENSION X(9),Y(9),LM2(7*9),TT21(7*9),TT210(7*9),TT21A(2*9)
      DIMENSION R(3,3),S(3,3),W(3,3)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION AKVV2(9,9)
      DIMENSION AKMIV2(9,9)
      DIMENSION AKMIV3(9,9)
C-----------------------------------------
C Dodate matrice za turbulentno strujanje
      DIMENSION AKMIV2T(9,9)
      DIMENSION AKMIV3T(9,9)
      DIMENSION AMK(9,9)
      DIMENSION AKVK(9,9)
      DIMENSION AKMK(9,9)
      DIMENSION AKBETAK(9,9)
      DIMENSION AKMOMEGA(9,9)
      DIMENSION AKBETAOM(9,9)

      DIMENSION AKKV(9,9)
      DIMENSION AKVVK1(9,9)
      DIMENSION AKVVK2(9,9)
      DIMENSION AKVVK3(9,9)
      DIMENSION AKVVK4(9,9)
      DIMENSION AKOMV(9,9)
      DIMENSION AKVVOM1(9,9)
      DIMENSION AKVVOM2(9,9)
      DIMENSION AKVVOM3(9,9)
      DIMENSION AKVVOM4(9,9)
      
      DIMENSION FSK(9)
      DIMENSION FSOMEGA(9)
      DIMENSION AKVV11(9,9)
      DIMENSION AKVV12(9,9)
      DIMENSION AKVV13(9,9)
      DIMENSION AKVV14(9,9)
      DIMENSION AKVV21(9,9)
      DIMENSION AKVV22(9,9)
      DIMENSION AKVV23(9,9)
      DIMENSION AKVV24(9,9)
C-----------------------------------------
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
      DIMENSION RB2(9)
      DIMENSION RB3(9)
C------------------------------
C PROMENJENO F31(31) U F31(60)
C------------------------------
      DIMENSION F31(60)
      DIMENSION AKUAX(9,9),AKVAX(9,9),APAXIX(9,9),APAXIY(9,9)
	DIMENSION SKE(31*16)

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
C------------------------------------------- 
C     KONSTANTE TURBULENTNOG MODELA K-OMEGA
C-------------------------------------------  
      ALFAZV=1.0
      ALFAK=0.55555555555556
      BETAOM=0.075
      BETAK=0.09
      SIGMAT=0.5
      SIGMAZV=0.5
      AMIT0=1.0
      TETA210=1.0
C----------------------------------------- 
      IF (IMUMPS.EQ.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN 
#if(MUMPS_CLUSTER)
      if (kkorak.eq.0.and. iter.eq.0) 
     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,6)
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,5)
C	REWIND(MUFILE2)
#else
        CALL sparseassembler_init(0)
#endif
	ENDIF
	
C GLAVNA PETLJA PO ELEMENTIMA
C GLOBAL LOOP PER ELEMENTS
C NBREL is counter of elements


      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,6*NDIM
      TT210(I)=0.D0
      LM2(I)=0
      TT21(I)=0.D0
 125  CONTINUE
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
C=======================================================================
C     ZA TURBULENTNO
      LM2(KLM+4*NDIM)=ID(6,NEL(KLM,NBREL))
      LM2(KLM+5*NDIM)=ID(7,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
C     ZA BRZINE, DVE KOMPONENTE
      DO 140 KLM=1,NDIM
      DO 135 NR=1,2
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
C=======================================================================
C     ZA PRITISKE 
      IF (KLM.LE.4) THEN
       TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
       TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
C=======================================================================
C     ZA TEMPERATURE
      TT21(KLM+2*NDIM+4)=GNODE(2,5,NEL(KLM,NBREL))
      TT210(KLM+2*NDIM+4)=GNODE(1,5,NEL(KLM,NBREL))
C=======================================================================
C     ZA KT TURBULENTNO
      TT21(KLM+4*NDIM)=GNODE(2,6,NEL(KLM,NBREL))
      TT210(KLM+4*NDIM)=GNODE(1,6,NEL(KLM,NBREL))
C=======================================================================
C     ZA OMT TURBULENTNO
      TT21(KLM+5*NDIM)=GNODE(2,7,NEL(KLM,NBREL))
      TT210(KLM+5*NDIM)=GNODE(1,7,NEL(KLM,NBREL))
C=======================================================================
 140  CONTINUE

C--------------------------------------------------------------------
      DO I=1,6*NDIM
       IF(KKORAK.EQ.1.AND.ITER.EQ.1) TT210(I)=TT21(I)
      ENDDO

C==============================================================================
C TURBULENTNA DINAMICKA VISKOZNOST IZ PRETHODNE ITERACIJE (AMIT) I ODNOS TETA21
C------------------------------------------------------------------------------
C AMIT I TETA21
C==============================================================================
      DO 150 I=1,NDIM
      AMIT=ALFAZV*GUSM*(GNODE(1,6,NEL(I,NBREL))/GNODE(1,7,NEL(I,NBREL)))
      TETA21=GNODE(1,7,NEL(I,NBREL))/GNODE(1,6,NEL(I,NBREL))
  150   CONTINUE
C==============================================================================
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
      AKMIV2(K,N)=0.D0
      AKMIV3(K,N)=0.D0
C-------------------------------------------
C POCETNE VREDNOSTI TURBULENTNIH MATRICA
      AKMIV2T(K,N)=0.D0
      AKMIV3T(K,N)=0.D0
      AMK(K,N)=0.D0
      AKVK(K,N)=0.D0
      AKMK(K,N)=0.D0
      AKBETAK(K,N)=0.D0
      AKMOMEGA(K,N)=0.D0
      AKBETAOM(K,N)=0.D0

      AKKV(K,N)=0.D0
      AKVVK1(K,N)=0.D0
      AKVVK2(K,N)=0.D0
      AKVVK3(K,N)=0.D0
      AKVVK4(K,N)=0.D0
      AKOMV(K,N)=0.D0
      AKVVOM1(K,N)=0.D0
      AKVVOM2(K,N)=0.D0
      AKVVOM3(K,N)=0.D0
      AKVVOM4(K,N)=0.D0
      
      FSK(K)=0.D0
      FSOMEGA(K)=0.D0
      AKVV11(K,N)=0.D0
      AKVV12(K,N)=0.D0
      AKVV13(K,N)=0.D0
      AKVV14(K,N)=0.D0
      AKVV21(K,N)=0.D0
      AKVV22(K,N)=0.D0
      AKVV23(K,N)=0.D0
      AKVV24(K,N)=0.D0
C-------------------------------------------
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

C subroutine FNTERPT return to us interpolation functions...

      CALL FNTERPT(R(IBRGT,I),S(IBRGT,J),1,TT21,H,HP,ZVHX,ZVHY,HV2,HV3,
     &ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVXOM,ZVYOM,NBREL,
     &IIZLAZ,GNODE,NEL,CORD,ID,TT210)
     
       IF (IALE.EQ.1) CALL ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
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
      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVXT
     1*GUSM*CC)
      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVYT
     1*GUSM*CC)
C---------------------------------------------------------------
C     DEAKTIVIRANA STARA MATRICA 

      AKMIV2(K,N)=AKMIV2(K,N)+WDT*((ZVHX(K)*ZVHX(N)+ZVHY(K)*ZVHY(N))
     &*(AMI+AMIT))

C     KREIRANA NOVA MATRICA AKMIV2T 
C---------------------------------------------------------------
C      AKMIV2T(K,N)=AKMIV2T(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
C     1ZVHY(K)*ZVHY(N))*(AMI+AMIT))
C---------------------------------------------------------------
      AKK(K,N)=AKK(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AKT)
      IF (INDAX.EQ.1) THEN
       AKK(K,N)=AKK(K,N)-WDT*(H(K)*ZVHX(N)/DEBLJ)*AKT
      ENDIF
C---------------------------------------------------------------
C     ZAMENJENO SA NOVIM IZRAZOM      
      AKMIV3(K,N)=AKMIV2(K,N)
C---------------------------------------------------------------
C      AKMIV3T(K,N)=AKMIV2T(K,N)
      AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM
      AJV2V2(K,N)=AJV2V2(K,N)+WDT*(H(K)*ZVXV2*H(N)*GUSM)
      AJV3V3(K,N)=AJV3V3(K,N)+WDT*(H(K)*ZVYV3*H(N)*GUSM)
      AJV2V3(K,N)=AJV2V3(K,N)+WDT*(H(K)*ZVYV2*H(N)*GUSM)
      AJV3V2(K,N)=AJV3V2(K,N)+WDT*(H(K)*ZVXV3*H(N)*GUSM)
C--------------------------------------------------------------------
C MATRICE IZ TURBULENTNOG MODELA K-OMEGA:
C--------------------------------------------------------------------
C     LEVA STRANA JEDNACINE 4.182
C     4.183 
      AMK(K,N)=AMK(K,N)+WDT*H(K)*H(N)*GUSM
C     4.184
      AKVK(K,N)=AKVK(K,N)+WDT*((H(K)*HV2*ZVHX(N)+
     1H(K)*HV3*ZVHY(N))*GUSM)  
C     4.185      
      AKMK(K,N)=AKMK(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     &ZVHY(K)*ZVHY(N))*(AMI+SIGMAZV*AMIT))
C     4.186      
      AKBETAK(K,N)=AKBETAK(K,N)+WDT*(2*GUSM*BETAK*TETA21*H(K)*H(N)*HKT)
C     4.187
      AKMOMEGA(K,N)=AKMOMEGA(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     &ZVHY(K)*ZVHY(N))*(AMI+SIGMAT*AMIT))
C     4.188
      AKBETAOM(K,N)=AKBETAOM(K,N)+WDT*(2*GUSM*BETAOM*H(K)*H(N)*HOM)
C--------------------------------------------------------------------
C     4.189
      AKKV(K,N)=AKKV(K,N)+WDT*(H(K)*ZVXK*H(N)+H(K)*ZVYK*H(N))*GUSM
C     4.190
      AKVVK1(K,N)=AKVVK1(K,N)+WDT*(-4*AMIT*((H(K)*ZVXV2*ZVHX(N))+
     &(H(K)*ZVXV2*ZVHY(N))))
      AKVVK2(K,N)=AKVVK2(K,N)+WDT*(-4*AMIT*((H(K)*ZVXV3*ZVHX(N))+
     &(H(K)*ZVXV3*ZVHY(N))))
      AKVVK3(K,N)=AKVVK3(K,N)+WDT*(-4*AMIT*((H(K)*ZVYV2*ZVHX(N))+
     &(H(K)*ZVYV2*ZVHY(N))))
      AKVVK4(K,N)=AKVVK4(K,N)+WDT*(-4*AMIT*((H(K)*ZVYV3*ZVHX(N))+
     &(H(K)*ZVYV3*ZVHY(N))))
C     4.191
      AKOMV(K,N)=AKOMV(K,N)+WDT*(H(K)*ZVXOM*H(N)+H(K)*ZVYOM*H(N))*GUSM
C     4.192
      AKVVOM1(K,N)=AKVVOM1(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV2*ZVHX(N))+(H(K)*ZVXV2*ZVHY(N))))
      AKVVOM2(K,N)=AKVVOM2(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV3*ZVHX(N))+(H(K)*ZVXV3*ZVHY(N))))
      AKVVOM3(K,N)=AKVVOM3(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV2*ZVHX(N))+(H(K)*ZVYV2*ZVHY(N))))
      AKVVOM4(K,N)=AKVVOM4(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV3*ZVHX(N))+(H(K)*ZVYV3*ZVHY(N))))
C--------------------------------------------------------------------
C DESNA STRANA JEDNACINE 4.182
C------------------------------------------------------
C     4.198
      AKVV11(K,N)=AKVV11(K,N)+WDT*(2*AMIT*((H(K)*ZVXV2*ZVHX(N))+
     &(H(K)*ZVXV2*ZVHY(N))))
      AKVV12(K,N)=AKVV12(K,N)+WDT*(2*AMIT*((H(K)*ZVXV3*ZVHX(N))+
     &(H(K)*ZVXV3*ZVHY(N))))
      AKVV13(K,N)=AKVV13(K,N)+WDT*(2*AMIT*((H(K)*ZVYV2*ZVHX(N))+
     &(H(K)*ZVYV2*ZVHY(N))))
      AKVV14(K,N)=AKVV14(K,N)+WDT*(2*AMIT*((H(K)*ZVYV3*ZVHX(N))+
     &(H(K)*ZVYV3*ZVHY(N))))
C     4.199
      AKVV21(K,N)=AKVV21(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV2*ZVHX(N))+(H(K)*ZVXV2*ZVHY(N))))
      AKVV22(K,N)=AKVV22(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV3*ZVHX(N))+(H(K)*ZVXV3*ZVHY(N))))
      AKVV23(K,N)=AKVV23(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV2*ZVHX(N))+(H(K)*ZVYV2*ZVHY(N))))
      AKVV24(K,N)=AKVV24(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV3*ZVHX(N))+(H(K)*ZVYV3*ZVHY(N))))
C--------------------------------------------------------------------
            IF (N.LE.4.AND.PENALT.LT.1.D0) THEN
      AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
      AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
      IF (INDAX.EQ.1) THEN
       AKV2P1(K,N)=AKV2P1(K,N)+WDT*(-H(K)*HP(N)/DEBLJ)
      ENDIF
      ENDIF
  164 CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(0.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(0.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE
C===========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PENALTY FAKTOROM
C reduced integration if penalty function is defined

      IF (PENALT.GT.0.D0) THEN
      DO 200 I=1,IBRGT-1
      DO 195 J=1,IBRGT-1
      CALL FNTERPT(R(IBRGT-1,I),S(IBRGT-1,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVXOM,ZVYOM,NBREL,
     &IIZLAZ,GNODE,NEL,CORD,ID,TT210)
         
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
C POVRSINSKE SILE
C surface forces

       INDX=0
       INDY=0

      DO 250 JBRPS=1,MAXSIL
      IF (NBREL.EQ.NGPSIL(1,JBRPS)) THEN
        IF (PENALT.LE.1.D0) IBRGT=IBRGT+1
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
      IF(NGPSIL(4,JBRPS).EQ.3) GOTO 207
      DO 206 K=1,NDIM
       RS2(K)=RS2(K)+(H(K)*FS2)*WDT
       RS3(K)=RS3(K)+(H(K)*FS3)*WDT
C------------------------------------------------------
C DESNA STRANA JEDNACINE 4.182
C------------------------------------------------------
C     4.196
      FSK(K+4*NDIM)=FSK(K+4*NDIM)+WDT*(AMI+SIGMAZV*AMIT)*H(K)*
     &FSK1*NGPSIL(4,JBRPS)
C     4.197
      FSOMEGA(K+5*NDIM)=FSOMEGA(K+5*NDIM)+WDT*(AMI+SIGMAT*AMIT)*H(K)*
     &FSOMEGA1*NGPSIL(4,JBRPS)
C------------------------------------------------------
C------------------------------------------------------
C------------------------------------------------------
 206  CONTINUE
 207   WRITE(IIZLAZ,*)'NBREL= ',NBREL,'TAU= ',TAU
 225  CONTINUE

      IF (PENALT.LE.1.D0) IBRGT=IBRGT-1
      ENDIF
 250  CONTINUE
C------------------------------------------------------
C C matrix include heat conduction
C------------------------------------------------------
C   FALI MATRICA JTT

      DO 255 I=1,NDIM
      DO 254 J=1,NDIM

      C(I,J)=AMV2(I,J)*CC              
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
C------------------------------------------------------
C      SKEF(I,J)=AKMIV2(I,J)+AKVV2(I,J)+AJV2V2(I,J)
C     PROMENJENO UMESTO AKMIV2 STAVLJENO AKMIV2T 
C------------------------------------------------------------------ 
C MNOZENJE SA BRZINAMA
C------------------------------------------------------------------
      SKEF(I,J)=AKMIV2(I,J)+AKVV2(I,J)+AJV2V2(I,J)
      SKEF(I+NDIM,J)=AJV3V2(I,J)
      SKEF(I+NDIM,J+NDIM)=AKMIV3(I,J)+AKVV2(I,J)+AJV3V3(I,J)
      SKEF(I,J+NDIM)=AJV2V3(I,J)
C------------------------------------------------------------------
C MNOZENJE SA TEMPERATURAMA
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=AKK(I,J)+AJK(I,J)
      IF (NSTAC.EQ.0) THEN
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=SKEF(I+2*NDIM+4,J+2*NDIM+4)+
     &C(I,J)/TIME
      SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
      SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
      ENDIF
C------------------------------------------------------------------
      IF (J.LE.4) THEN
       SKEF(J+2*NDIM,I)=AKV2P(I,J)+AKV2P1(I,J)
       SKEF(J+2*NDIM,I+NDIM)=AKV3P(I,J)
       SKEF(I,J+2*NDIM)=AKV2P(I,J)
       SKEF(I+NDIM,J+2*NDIM)=AKV3P(I,J)
      ENDIF
      
       SKEF(I+2*NDIM+4,J)=AKTV2(I,J)
       SKEF(I+2*NDIM+4,J+NDIM)=AKTV3(I,J)
C=========================================================================
C K-OMEGA TURBULENTNI MODEL
C------------------------------------------------------
C------------------------------------------------------
C LEVA STRANA K JEDNACINE
C------------------------------------------------------
C
      SKEF(I+4*NDIM,J+4*NDIM)=AKVK(I,J)+AKMK(I,J)+AKBETAK(I,J) 
      SKEF(I+4*NDIM,J)=SKEF(I+4*NDIM,J)+AKKV(I,J)+AKVVK1(I,J)+
     &AKVVK2(I,J)+AKVVK3(I,J)+AKVVK4(I,J)
             
      IF (NSTAC.EQ.0) THEN
      SKEF(I+4*NDIM,J+4*NDIM)=SKEF(I+4*NDIM,J+4*NDIM)+AMK(I,J)/TIME
      ENDIF
C------------------------------------------------------
C LEVA STRANA OMEGA JEDNACINE
C------------------------------------------------------
      SKEF(I+5*NDIM,J+5*NDIM)=AKVK(I,J)+AKMOMEGA(I,J)+AKBETAOM(I,J)
      SKEF(I+5*NDIM,J)=SKEF(I+5*NDIM,J)+AKOMV(I,J)+AKVVOM1(I,J)+
     &AKVVOM2(I,J)+AKVVOM3(I,J)+AKVVOM4(I,J)

      IF (NSTAC.EQ.0) THEN
      SKEF(I+5*NDIM,J+5*NDIM)=SKEF(I+5*NDIM,J+5*NDIM)+AMK(I,J)/TIME
      ENDIF
C=========================================================================     
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
C------------------------------------------------------
      F31(I)=F31(I)-(AKMIV2(I,J)+AKVV2(I,J))*TT21(J) 
C------------------------------------------------------
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
C-------------------------------------------------------------------- 
C K-OMEGA TURBULENTNI MODEL
C------------------------------------------------------
C DESNA STRANA K JEDNACINE
C------------------------------------------------------
C     4.193
      F31(I+4*NDIM)=F31(I+4*NDIM)-(AKVK(I,J)+AKMK(I,J)+
     &0.5*AKBETAK(I,J))*TT21(J+4*NDIM)
C------------------------------------------------------
C     JEDNACINA  4.198
C------------------------------------------------------
      F31(I+4*NDIM)=F31(I+4*NDIM)+(AKVV11(I,J)+AKVV12(I,J)+AKVV13(I,J)+
     &AKVV14(I,J))*TT21(J)
     
      IF (NSTAC.EQ.0) THEN
      II4=I+4*NDIM
      JJ4=J+4*NDIM
      F31(II4)=F31(II4)-AMK(I,J)*(TT21(JJ4)-TT210(JJ4))/TIME
      ENDIF
      
      F31(I+4*NDIM)=F31(I+4*NDIM)+FSK(I)
C------------------------------------------------------
C DESNA STRANA OMEGA JEDNACINE
C------------------------------------------------------
C     4.194
      F31(I+5*NDIM)=F31(I+5*NDIM)-(AKVK(I,J)+AKMOMEGA(I,J)+
     &0.5*AKBETAOM(I,J))*TT21(J+5*NDIM)
C------------------------------------------------------
C     JEDNACINA 4.199
C------------------------------------------------------
      F31(I+5*NDIM)=F31(I+5*NDIM)+(AKVV21(I,J)+AKVV22(I,J)+AKVV23(I,J)+
     &AKVV24(I,J))*TT21(J)
     
      IF (NSTAC.EQ.0) THEN
      II5=I+5*NDIM
      JJ5=J+5*NDIM
      F31(II5)=F31(II5)-AMK(I,J)*(TT21(JJ5)-TT210(JJ5))/TIME
      ENDIF

      F31(I+5*NDIM)=F31(I+5*NDIM)+FSOMEGA(I)
C------------------------------------------------------
     
C DESNA STRANA USLED UTICAJA KONVEKCIJE
c right-hand side created by natural convection
      F31(I)=F31(I)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB2
      F31(I+NDIM)=F31(I+NDIM)-AMV2(I,J)*TT21(J+2*NDIM+4)*BETA*FB3
 285  CONTINUE							      
 290  CONTINUE
 
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
C      CALL SSTRE2(NEL,ID,NPARAM,TT21,H,HP,X,Y,R,S,W,ZVHX,ZVHY,
C     &PRES,INDEL,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,
C     &INDAX,AKT,GUSM,IUPWIN,NBREL,IIZLAZ,IBRGT)
C============================================================================
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
C	  CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,6)
C        CALL SPAKDE (SILE,F92,LM2,NDES)
C	  GOTO 400
C      ENDIF
      
      IF (IMUMPS.EQ.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN  
C	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,7)
#if(MUMPS_CLUSTER)
	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,7)
#else
       CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
#endif
       CALL SPAKDE (SILE,F31,LM2,NDES)
	 GOTO 400
      ENDIF

	IF(NJUTRA.EQ.1.AND.ITER.GT.0) THEN
       CALL SPAKDE (SILE,F31,LM2,NDES)
	ELSE
       IF (ISYMMS.EQ.1) THEN
	    CALL MATSTE (ALEVO,MAXA,SILE,SKE,F31,LM2,NDES,1)
	 ELSE
          CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F31,MAXA,LM2,NDES,1)
	 ENDIF
	ENDIF

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C end of loop per elements
C=======================================================================
 400  CONTINUE
C       CALL MUMPSRIGHT(SILE,JEDN)

C End of subroutine
      END
C=======================================================================   

C=========================================================================
C MODIFIKOVAO I KREIRAO ALEKSANDAR NIKOLIC, FEBRUAR 2017.
C PROGRAM ZA TURBULENTNO STRUJANJE FLUIDA, K-OMEGA TURBULENTNI MODEL, 3D
C=========================================================================
C=========================================================================

      SUBROUTINE RACU3DT(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,DELTAL,PRES1,TAU,VOSI,DT,
     &IBKOR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine RACU3DT is used for 3D analysis in turbulent K-OMEGA mode
CE It is used global loop per elements
C
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE, MUFILE2
	COMMON /TRANSP/ INDFL,LID1
		
      DIMENSION GNODE(2,7,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(7,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)
      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)
      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21)
      DIMENSION H(21),HP(8),TT21A(3*21)
      DIMENSION AKVV1(21,21),AKMIV1(21,21),AMV2(21,21),AJV1V1(21,21)
      DIMENSION AJV1V2(21,21),AJV1V3(21,21),AJV2V1(21,21),AJV2V2(21,21)
      DIMENSION AJV2V3(21,21),AJV3V1(21,21),AJV3V2(21,21),AJV3V3(21,21)
      DIMENSION AKTV1(21,21),AKTV2(21,21),AKTV3(21,21),AKK(21,21)
      DIMENSION AKV1P(21,8),AKV2P(21,8),AKV3P(21,8),RS1(21),RS2(21)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3)
      DIMENSION XG(55),WGT(55),NREF(11)
      DIMENSION PENXX(8,8),PENXY(8,8),PENXZ(8,8)
      DIMENSION PENYX(8,8),PENYY(8,8),PENYZ(8,8)
      DIMENSION PENZX(8,8),PENZY(8,8),PENZZ(8,8)
      DIMENSION C(21,21)
      DIMENSION PRES1(100,3,*)
      DIMENSION TAU(100,2,*)
      DIMENSION VOSI(100,*)
C-----------------------------------------
C Dodate matrice za turbulentno strujanje
      DIMENSION AMK(21,21)
      DIMENSION AKVK(21,21)
      DIMENSION AKMK(21,21)
      DIMENSION AKBETAK(21,21)
      DIMENSION AKMOMEGA(21,21)
      DIMENSION AKBETAOM(21,21)
      DIMENSION AKKV(21,21)
      DIMENSION AKOMV(21,21)
      
      DIMENSION AKVVK1(21,21)
      DIMENSION AKVVK2(21,21)
      DIMENSION AKVVK3(21,21)
      DIMENSION AKVVK4(21,21)
      DIMENSION AKVVK5(21,21)
      DIMENSION AKVVK6(21,21)
      DIMENSION AKVVK7(21,21)
      DIMENSION AKVVK8(21,21)
      DIMENSION AKVVK9(21,21)

      DIMENSION AKVVOM1(21,21)
      DIMENSION AKVVOM2(21,21)
      DIMENSION AKVVOM3(21,21)
      DIMENSION AKVVOM4(21,21)
      DIMENSION AKVVOM5(21,21)
      DIMENSION AKVVOM6(21,21)
      DIMENSION AKVVOM7(21,21)
      DIMENSION AKVVOM8(21,21)
      DIMENSION AKVVOM9(21,21)

      DIMENSION AKVV11(21,21)
      DIMENSION AKVV12(21,21)
      DIMENSION AKVV13(21,21)
      DIMENSION AKVV14(21,21)
      DIMENSION AKVV15(21,21)
      DIMENSION AKVV16(21,21)
      DIMENSION AKVV17(21,21)
      DIMENSION AKVV18(21,21)
      DIMENSION AKVV19(21,21)
      
      DIMENSION AKVV21(21,21)
      DIMENSION AKVV22(21,21)
      DIMENSION AKVV23(21,21)
      DIMENSION AKVV24(21,21)
      DIMENSION AKVV25(21,21)
      DIMENSION AKVV26(21,21)
      DIMENSION AKVV27(21,21)
      DIMENSION AKVV28(21,21)
      DIMENSION AKVV29(21,21)
      
      DIMENSION FSK(21)
      DIMENSION FSOMEGA(21)
C-----------------------------------------
      
      DATA NREF/0,1,3,6,10,15,21,28,36,45,55/
      DATA WGT/            2.D0,               1.D0,               1.D0,
     1       .555555555555556D0, .888888888888889D0, .555555555555556D0,
     2       .347854845137454D0, .652145154862546D0, .652145154862546D0,
     3       .347854845137454D0, .236926885056189D0, .478628670499366D0,
     4       .568888888888889D0, .478628670499366D0, .236926885056189D0,
     5       .171324492379170D0, .360761573048139D0, .467913934572691D0,
     6       .467913934572691D0, .360761573048139D0, .171324492379170D0,
     7       .129484966168870D0, .279705391489277D0, .381830050505119D0,
     8       .417959183673469D0, .381830050505119D0, .279705391489277D0,
     9       .129484966168870D0, .101228536290376D0, .222381034453374D0,
     9       .313706645877887D0, .362683783378362D0, .362683783378362D0,
     1       .313706645877887D0, .222381034453374D0, .101228536290376D0,
     2       .081274388361574D0, .180648160694857D0, .260610696402935D0,
     3       .312347077040003D0, .330239355001260D0, .312347077040003D0,
     4       .260610696402935D0, .180648160694857D0, .081274388361574D0,
     5       .066671344308688D0, .149451349150581D0, .219086362515982D0,
     6       .269266719309996D0, .295524224714753D0, .295524224714753D0,
     7       .269266719309996D0, .219086362515982D0, .149451349150581D0,
     8       .066671344308688D0/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0,
     5      -.932469514203152D0,-.661209386466265D0,-.238619186083197D0,
     6       .238619186083197D0, .661209386466265D0, .932469514203152D0,
     7      -.949107912342759D0,-.741531185599394D0,-.405845151377397D0,
     8                     0.D0, .405845151377397D0, .741531185599394D0,
     9       .949107912342759D0,-.960289856497536D0,-.796666477413627D0,
     9      -.525532409916329D0,-.183434642495650D0, .183434642495650D0,
     1       .525532409916329D0, .796666477413627D0, .960289856497536D0,
     2      -.968160239507626D0,-.836031107326636D0,-.613371432700590D0,
     3      -.324253423403809D0,               0.D0, .324253423403809D0,
     4       .613371432700590D0, .836031107326636D0, .968160239507626D0,
     5      -.973906528517172D0,-.865063366688985D0,-.679409568299024D0,
     6      -.433395394129247D0,-.148874338981631D0, .148874338981631D0,
     7       .433395394129247D0, .679409568299024D0, .865063366688985D0,
     8       .973906528517172D0/

C------------------------------------------- 
C     KONSTANTE TURBULENTNOG MODELA K-OMEGA
C-------------------------------------------  
      ALFAZV=1.0
      ALFAK=0.55555555555556
      BETAOM=0.075
      BETAK=0.09
      SIGMAT=0.5
      SIGMAZV=0.5
      AMIT0=1.0
      TETA210=1.0
C-----------------------------------------
      IF (IMUMPS.EQ.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN 
#if(MUMPS_CLUSTER)
      if (kkorak.eq.0.and. iter.eq.0) 
C     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)
     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,6)
C	REWIND(MUFILE2)
#else
        CALL sparseassembler_init(0)
#endif
	ENDIF

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
C  ZA TURBULENTNO
      LM2(KLM+5*NDIM)=ID(6,NEL(KLM,NBREL))
      LM2(KLM+6*NDIM)=ID(7,NEL(KLM,NBREL))
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
C ZA BRZINE 3 KOMPONENTE
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
C=========================================================================
C ZA PRITISKE
      IF (KLM.LE.8) THEN
        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
C=========================================================================
C ZA TEMPERATURE
        TT21(KLM+3*NDIM+8)=GNODE(2,5,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM+8)=GNODE(1,5,NEL(KLM,NBREL))
C=========================================================================
C ZA KT TURBULENTNO
        TT21(KLM+5*NDIM)=GNODE(2,6,NEL(KLM,NBREL))
        TT210(KLM+5*NDIM)=GNODE(1,6,NEL(KLM,NBREL))
C=========================================================================
C ZA OMT TURBULENTNO
        TT21(KLM+6*NDIM)=GNODE(2,7,NEL(KLM,NBREL))
        TT210(KLM+6*NDIM)=GNODE(1,7,NEL(KLM,NBREL))
 140  CONTINUE
C==============================================================================
C==============================================================================
C TURBULENTNA DINAMICKA VISKOZNOST IZ PRETHODNE ITERACIJE (AMIT) I ODNOS TETA21
C------------------------------------------------------------------------------
C AMIT I TETA21
C==============================================================================
      DO 150 I=1,NDIM
      AMIT=ALFAZV*GUSM*(GNODE(1,6,NEL(I,NBREL))/GNODE(1,7,NEL(I,NBREL)))
      TETA21=GNODE(1,7,NEL(I,NBREL))/GNODE(1,6,NEL(I,NBREL))
  150   CONTINUE
C==============================================================================

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
C------------------------------------------------------------------------------
C POCETNE VREDNOSTI TURBULENTNIH MATRICA
      AMK(K,N)=0.D0
      AKVK(K,N)=0.D0
      AKMK(K,N)=0.D0
      AKBETAK(K,N)=0.D0
      AKMOMEGA(K,N)=0.D0
      AKBETAOM(K,N)=0.D0
      AKKV(K,N)=0.D0
      AKOMV(K,N)=0.D0
      
      AKVVK1(K,N)=0.D0
      AKVVK2(K,N)=0.D0
      AKVVK3(K,N)=0.D0
      AKVVK4(K,N)=0.D0
      AKVVK5(K,N)=0.D0
      AKVVK6(K,N)=0.D0
      AKVVK7(K,N)=0.D0
      AKVVK8(K,N)=0.D0
      AKVVK9(K,N)=0.D0

      AKVVOM1(K,N)=0.D0
      AKVVOM2(K,N)=0.D0
      AKVVOM3(K,N)=0.D0
      AKVVOM4(K,N)=0.D0
      AKVVOM5(K,N)=0.D0
      AKVVOM6(K,N)=0.D0
      AKVVOM7(K,N)=0.D0
      AKVVOM8(K,N)=0.D0
      AKVVOM9(K,N)=0.D0

      AKVV11(K,N)=0.D0
      AKVV12(K,N)=0.D0
      AKVV13(K,N)=0.D0
      AKVV14(K,N)=0.D0
      AKVV15(K,N)=0.D0
      AKVV16(K,N)=0.D0
      AKVV17(K,N)=0.D0
      AKVV18(K,N)=0.D0
      AKVV19(K,N)=0.D0
      
      AKVV21(K,N)=0.D0
      AKVV22(K,N)=0.D0
      AKVV23(K,N)=0.D0
      AKVV24(K,N)=0.D0
      AKVV25(K,N)=0.D0
      AKVV26(K,N)=0.D0
      AKVV27(K,N)=0.D0
      AKVV28(K,N)=0.D0
      AKVV29(K,N)=0.D0
      
      FSK(K)=0.D0
      FSOMEGA(K)=0.D0
C------------------------------------------------------------------------------
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
 
       CALL JACTT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     &,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,
     &FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVZK,ZVXOM,ZVYOM,ZVZOM,
     &ZVXV1,ZVXV2,ZVXV3,ZVYV1,ZVYV2,ZVYV3,ZVZV1,ZVZV2,ZVZV3,CORD,ID,
     &PRITISAK,PENALT,PRIT,IDPRIT,GUSM)

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
C---------------------------------------------------------------
C DODATO AMIT U AKMIV1
      AKMIV1(K,N)=AKMIV1(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*(AMI+AMIT))
C---------------------------------------------------------------
      AKK(K,N)=AKK(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AKT)
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
C--------------------------------------------------------------------
C MATRICE IZ TURBULENTNOG MODELA K-OMEGA:
C--------------------------------------------------------------------
C     LEVA STRANA JEDNACINE 4.182
C     4.183 
      AMK(K,N)=AMK(K,N)+WDT*H(K)*H(N)*GUSM
C     4.184
      AKVK(K,N)=AKVK(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM) 
C     4.185      
      AKMK(K,N)=AKMK(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*(AMI+SIGMAZV*AMIT))
C     4.186      
      AKBETAK(K,N)=AKBETAK(K,N)+WDT*(2*GUSM*BETAK*TETA21*H(K)*H(N)*HKT)
C     4.187
      AKMOMEGA(K,N)=AKMOMEGA(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*(AMI+SIGMAT*AMIT))
C     4.188
      AKBETAOM(K,N)=AKBETAOM(K,N)+WDT*(2*GUSM*BETAOM*H(K)*H(N)*HOM)
C--------------------------------------------------------------------
C     4.189
      AKKV(K,N)=AKKV(K,N)+WDT*(H(K)*ZVXK*H(N)+H(K)*ZVYK*H(N)+
     &H(K)*ZVZK*H(N))*GUSM
C     4.190
      AKVVK1(K,N)=AKVVK1(K,N)+WDT*(-4*AMIT*((H(K)*HXU*PJ(1,N))+
     &(H(K)*HXU*PJ(2,N))+(H(K)*HXU*PJ(3,N))))
      AKVVK2(K,N)=AKVVK2(K,N)+WDT*(-4*AMIT*((H(K)*HYU*PJ(1,N))+
     &(H(K)*HYU*PJ(2,N))+(H(K)*HYU*PJ(3,N))))
      AKVVK3(K,N)=AKVVK3(K,N)+WDT*(-4*AMIT*((H(K)*HZU*PJ(1,N))+
     &(H(K)*HZU*PJ(2,N))+(H(K)*HZU*PJ(3,N))))
           
      AKVVK4(K,N)=AKVVK4(K,N)+WDT*(-4*AMIT*((H(K)*HXV*PJ(1,N))+
     &(H(K)*HXV*PJ(2,N))+(H(K)*HXV*PJ(3,N))))
      AKVVK5(K,N)=AKVVK5(K,N)+WDT*(-4*AMIT*((H(K)*HYV*PJ(1,N))+
     &(H(K)*HYV*PJ(2,N))+(H(K)*HYV*PJ(3,N))))
      AKVVK6(K,N)=AKVVK6(K,N)+WDT*(-4*AMIT*((H(K)*HZV*PJ(1,N))+
     &(H(K)*HZV*PJ(2,N))+(H(K)*HZV*PJ(3,N))))
     
      AKVVK7(K,N)=AKVVK7(K,N)+WDT*(-4*AMIT*((H(K)*HXW*PJ(1,N))+
     &(H(K)*HXW*PJ(2,N))+(H(K)*HXW*PJ(3,N))))
      AKVVK8(K,N)=AKVVK8(K,N)+WDT*(-4*AMIT*((H(K)*HYW*PJ(1,N))+
     &(H(K)*HYW*PJ(2,N))+(H(K)*HYW*PJ(3,N))))
      AKVVK9(K,N)=AKVVK9(K,N)+WDT*(-4*AMIT*((H(K)*HZW*PJ(1,N))+
     &(H(K)*HZW*PJ(2,N))+(H(K)*HZW*PJ(3,N))))
C     4.191
      AKOMV(K,N)=AKOMV(K,N)+WDT*(H(K)*ZVXOM*H(N)+H(K)*ZVYOM*H(N)+
     &H(K)*ZVZOM*H(N))*GUSM
C     4.192
      AKVVOM1(K,N)=AKVVOM1(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV1*PJ(1,N))+(H(K)*ZVXV1*PJ(2,N))+(H(K)*ZVXV1*PJ(3,N))))
      AKVVOM2(K,N)=AKVVOM2(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV2*PJ(1,N))+(H(K)*ZVXV2*PJ(2,N))+(H(K)*ZVXV2*PJ(3,N))))
      AKVVOM3(K,N)=AKVVOM3(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV3*PJ(1,N))+(H(K)*ZVXV3*PJ(2,N))+(H(K)*ZVXV3*PJ(3,N))))
     
      AKVVOM4(K,N)=AKVVOM4(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV1*PJ(1,N))+(H(K)*ZVYV1*PJ(2,N))+(H(K)*ZVYV1*PJ(3,N))))
      AKVVOM5(K,N)=AKVVOM5(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV2*PJ(1,N))+(H(K)*ZVYV2*PJ(2,N))+(H(K)*ZVYV2*PJ(3,N))))
      AKVVOM6(K,N)=AKVVOM6(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV3*PJ(1,N))+(H(K)*ZVYV3*PJ(2,N))+(H(K)*ZVYV3*PJ(3,N))))
     
      AKVVOM7(K,N)=AKVVOM7(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVZV1*PJ(1,N))+(H(K)*ZVZV1*PJ(2,N))+(H(K)*ZVZV1*PJ(3,N))))
      AKVVOM8(K,N)=AKVVOM8(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVZV2*PJ(1,N))+(H(K)*ZVZV2*PJ(2,N))+(H(K)*ZVZV2*PJ(3,N))))
      AKVVOM9(K,N)=AKVVOM9(K,N)+WDT*(-4*ALFAK*AMIT*TETA21*
     &((H(K)*ZVZV3*PJ(1,N))+(H(K)*ZVZV3*PJ(2,N))+(H(K)*ZVZV3*PJ(3,N))))
C--------------------------------------------------------------------
C DESNA STRANA JEDNACINE 4.182
C------------------------------------------------------
C     4.198
      AKVV11(K,N)=AKVV11(K,N)+WDT*(2*AMIT*((H(K)*ZVXV1*PJ(1,N))+
     &(H(K)*ZVXV1*PJ(2,N))+(H(K)*ZVXV1*PJ(3,N))))
      AKVV12(K,N)=AKVV12(K,N)+WDT*(2*AMIT*((H(K)*ZVXV2*PJ(1,N))+
     &(H(K)*ZVXV2*PJ(2,N))+(H(K)*ZVXV2*PJ(3,N))))
      AKVV13(K,N)=AKVV13(K,N)+WDT*(2*AMIT*((H(K)*ZVXV3*PJ(1,N))+
     &(H(K)*ZVXV3*PJ(2,N))+(H(K)*ZVXV3*PJ(3,N))))
     
      AKVV14(K,N)=AKVV14(K,N)+WDT*(2*AMIT*((H(K)*ZVYV1*PJ(1,N))+
     &(H(K)*ZVYV1*PJ(2,N))+(H(K)*ZVYV1*PJ(3,N))))
      AKVV15(K,N)=AKVV15(K,N)+WDT*(2*AMIT*((H(K)*ZVYV2*PJ(1,N))+
     &(H(K)*ZVYV2*PJ(2,N))+(H(K)*ZVYV2*PJ(3,N))))
      AKVV16(K,N)=AKVV16(K,N)+WDT*(2*AMIT*((H(K)*ZVYV3*PJ(1,N))+
     &(H(K)*ZVYV3*PJ(2,N))+(H(K)*ZVYV3*PJ(3,N))))
     
      AKVV17(K,N)=AKVV17(K,N)+WDT*(2*AMIT*((H(K)*ZVZV1*PJ(1,N))+
     &(H(K)*ZVZV1*PJ(2,N))+(H(K)*ZVZV1*PJ(3,N))))
      AKVV18(K,N)=AKVV18(K,N)+WDT*(2*AMIT*((H(K)*ZVZV2*PJ(1,N))+
     &(H(K)*ZVZV2*PJ(2,N))+(H(K)*ZVZV2*PJ(3,N))))
      AKVV19(K,N)=AKVV19(K,N)+WDT*(2*AMIT*((H(K)*ZVZV3*PJ(1,N))+
     &(H(K)*ZVZV3*PJ(2,N))+(H(K)*ZVZV3*PJ(3,N))))
C     4.199
      AKVV21(K,N)=AKVV21(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV1*PJ(1,N))+(H(K)*ZVXV1*PJ(2,N))+
     &(H(K)*ZVXV1*PJ(3,N))))
      AKVV22(K,N)=AKVV22(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV2*PJ(1,N))+(H(K)*ZVXV2*PJ(2,N))+
     &(H(K)*ZVXV2*PJ(3,N))))
      AKVV23(K,N)=AKVV23(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVXV3*PJ(1,N))+(H(K)*ZVXV3*PJ(2,N))+
     &(H(K)*ZVXV3*PJ(3,N))))
     
      AKVV24(K,N)=AKVV24(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV1*PJ(1,N))+(H(K)*ZVYV1*PJ(2,N))+
     &(H(K)*ZVYV1*PJ(3,N))))
      AKVV25(K,N)=AKVV25(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV2*PJ(1,N))+(H(K)*ZVYV2*PJ(2,N))+
     &(H(K)*ZVYV2*PJ(3,N))))
      AKVV26(K,N)=AKVV26(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVYV3*PJ(1,N))+(H(K)*ZVYV3*PJ(2,N))+
     &(H(K)*ZVYV3*PJ(3,N))))
     
      AKVV27(K,N)=AKVV27(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVZV1*PJ(1,N))+(H(K)*ZVZV1*PJ(2,N))+
     &(H(K)*ZVZV1*PJ(3,N))))
      AKVV28(K,N)=AKVV28(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVZV2*PJ(1,N))+(H(K)*ZVZV2*PJ(2,N))+
     &(H(K)*ZVZV2*PJ(3,N))))
      AKVV29(K,N)=AKVV29(K,N)+WDT*(2*ALFAK*AMIT*TETA21*
     &((H(K)*ZVZV3*PJ(1,N))+(H(K)*ZVZV3*PJ(2,N))+
     &(H(K)*ZVZV3*PJ(3,N))))
C--------------------------------------------------------------------
      
      
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
 
       CALL JACTT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     &,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,
     &FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVZK,ZVXOM,ZVYOM,ZVZOM,
     &ZVXV1,ZVXV2,ZVXV3,ZVYV1,ZVYV2,ZVYV3,ZVZV1,ZVZV2,ZVZV3,CORD,ID,
     &PRITISAK,PENALT,PRIT,IDPRIT,GUSM)
 
        WDT=WT*DET1
  
        DO  K=1,NDIM
        DO  N=1,NDIM
          PENXX(K,N)=PENXX(K,N)+WDT*PJ(1,K)*PJ(1,N)*PENALT
          PENXY(K,N)=PENXY(K,N)+WDT*PJ(1,K)*PJ(2,N)*PENALT
          PENXZ(K,N)=PENXZ(K,N)+WDT*PJ(1,K)*PJ(3,N)*PENALT
          PENYX(K,N)=PENYX(K,N)+WDT*PJ(2,K)*PJ(1,N)*PENALT
          PENYY(K,N)=PENYY(K,N)+WDT*PJ(2,K)*PJ(2,N)*PENALT
          PENYZ(K,N)=PENYZ(K,N)+WDT*PJ(2,K)*PJ(3,N)*PENALT
          PENZX(K,N)=PENZX(K,N)+WDT*PJ(3,K)*PJ(1,N)*PENALT
          PENZY(K,N)=PENZY(K,N)+WDT*PJ(3,K)*PJ(2,N)*PENALT
          PENZZ(K,N)=PENZZ(K,N)+WDT*PJ(3,K)*PJ(3,N)*PENALT
        ENDDO
        ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDIF

C======================================================================= 
C POVRSINSKE SILE
      

      DO 250 JBRPS=1,MAXSIL
      IF (NGPSIL(1,JBRPS).EQ.NBREL) THEN 
      CALL STRANA(NEL,NDIM,NGPSIL(1,JBRPS),NPOV)
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

       CALL JACTT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     &,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,
     &FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVZK,ZVXOM,ZVYOM,ZVZOM,
     &ZVXV1,ZVXV2,ZVXV3,ZVYV1,ZVYV2,ZVYV3,ZVZV1,ZVZV2,ZVZV3,CORD,ID,
     &PRITISAK,PENALT,PRIT,IDPRIT,GUSM)
C
      WDT=WT*DET
      POVRS=POVRS+WDT

      DO K=1,NDIM
       RS1(K)=RS1(K)+(H(K)*SF1)*WDT*NGPSIL(6,JBRPS)
       RS2(K)=RS2(K)+(H(K)*SF2)*WDT*NGPSIL(7,JBRPS)
       RS3(K)=RS3(K)+(H(K)*SF3)*WDT*NGPSIL(8,JBRPS)
C------------------------------------------------------
C DESNA STRANA JEDNACINE 4.182
C------------------------------------------------------
C     4.196
      FSK(K+5*NDIM)=FSK(K+5*NDIM)+WDT*(AMI+SIGMAZV*AMIT)*H(K)*
     &FSK1*NGPSIL(4,JBRPS)
C     4.197
      FSOMEGA(K+6*NDIM)=FSOMEGA(K+6*NDIM)+WDT*(AMI+SIGMAT*AMIT)*H(K)*
     &FSOMEGA1*NGPSIL(4,JBRPS)
C------------------------------------------------------
      ENDDO
  210 CONTINUE    
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
C PAKOVANJE MATRICE C,Kk,Jk U MATRICU SKEF ZA TEMPERATURE
      DO 263 I=1,NDIM
      DO 262 J=1,NDIM
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
C=========================================================================
C K-OMEGA TURBULENTNI MODEL
C------------------------------------------------------
C------------------------------------------------------
C LEVA STRANA K JEDNACINE
C------------------------------------------------------
C
      SKEF(I+5*NDIM,J+5*NDIM)=AKVK(I,J)+AKMK(I,J)+AKBETAK(I,J) 
      SKEF(I+5*NDIM,J)=SKEF(I+5*NDIM,J)+AKKV(I,J)+AKVVK1(I,J)+
     &AKVVK2(I,J)+AKVVK3(I,J)
      SKEF(I+5*NDIM,J+NDIM)=SKEF(I+5*NDIM,J+NDIM)+AKKV(I,J)+AKVVK4(I,J)+
     &AKVVK5(I,J)+AKVVK6(I,J)
      SKEF(I+5*NDIM,J+2*NDIM)=SKEF(I+5*NDIM,J+2*NDIM)+AKKV(I,J)+
     &AKVVK7(I,J)+AKVVK8(I,J)+AKVVK9(I,J)
     
      IF (NSTAC.EQ.0) THEN
      SKEF(I+5*NDIM,J+5*NDIM)=SKEF(I+5*NDIM,J+5*NDIM)+AMK(I,J)/TIME
      ENDIF
C------------------------------------------------------
C LEVA STRANA OMEGA JEDNACINE
C------------------------------------------------------
      SKEF(I+6*NDIM,J+6*NDIM)=AKVK(I,J)+AKMOMEGA(I,J)+AKBETAOM(I,J)
      SKEF(I+6*NDIM,J)=SKEF(I+6*NDIM,J)+AKOMV(I,J)+AKVVOM1(I,J)+
     &AKVVOM2(I,J)+AKVVOM3(I,J)
      SKEF(I+6*NDIM,J+NDIM)=SKEF(I+6*NDIM,J+NDIM)+AKOMV(I,J)+
     &AKVVOM4(I,J)+AKVVOM5(I,J)+AKVVOM6(I,J)
      SKEF(I+6*NDIM,J+2*NDIM)=SKEF(I+6*NDIM,J+2*NDIM)+AKOMV(I,J)+
     &AKVVOM7(I,J)+AKVVOM8(I,J)+AKVVOM9(I,J)
     
      IF (NSTAC.EQ.0) THEN
      SKEF(I+6*NDIM,J+6*NDIM)=SKEF(I+6*NDIM,J+6*NDIM)+AMK(I,J)/TIME
      ENDIF
C========================================================================= 

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

C DESNA STRANA:
      DO K=0,2
      DO 290 I=1,NDIM
        II=I+K*NDIM
      DO 285 J=1,NDIM
        JJ=J+K*NDIM
        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
      IF (NSTAC.EQ.0) THEN
        F92(II)=F92(II)-AMV2(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
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
C DESNA STRANA USLED PRITISKA
      DO 310 I=1,NDIM
      DO 300 J=1,NDIM
      IF (J.LE.8) THEN
       F92(I)=F92(I)-AKV1P(I,J)*TT21(J+3*NDIM)
       F92(I+NDIM)=F92(I+NDIM)-AKV2P(I,J)*TT21(J+3*NDIM)
       F92(I+2*NDIM)=F92(I+2*NDIM)-AKV3P(I,J)*TT21(J+3*NDIM)
      ENDIF
      IF (I.LE.8) THEN
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV1P(J,I)*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV2P(J,I)*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV3P(J,I)*TT21(J+2*NDIM)
      ENDIF
C-------------------------------------------------------------------- 
C K-OMEGA TURBULENTNI MODEL
C------------------------------------------------------
C DESNA STRANA K JEDNACINE
C------------------------------------------------------
C     4.193
      F92(I+5*NDIM)=F92(I+5*NDIM)-(AKVK(I,J)+AKMK(I,J)+
     &0.5*AKBETAK(I,J))*TT21(J+5*NDIM)
C------------------------------------------------------
C     JEDNACINA  4.198
C------------------------------------------------------
      F92(I+5*NDIM)=F92(I+5*NDIM)+((AKVV11(I,J)+AKVV12(I,J)+AKVV13(I,J))
     &*TT21(J)+(AKVV14(I,J)+AKVV15(I,J)+AKVV16(I,J))*TT21(J+NDIM)+
     &(AKVV17(I,J)+AKVV18(I,J)+AKVV19(I,J))*TT21(J+2*NDIM))
     
      IF (NSTAC.EQ.0) THEN
      II5=I+5*NDIM
      JJ5=J+5*NDIM
      F92(II5)=F92(II5)-AMK(I,J)*(TT21(JJ5)-TT210(JJ5))/TIME
      ENDIF
      
      F92(I+5*NDIM)=F92(I+5*NDIM)+FSK(I)
C------------------------------------------------------
C DESNA STRANA OMEGA JEDNACINE
C------------------------------------------------------
C     4.194
      F92(I+6*NDIM)=F92(I+6*NDIM)-(AKVK(I,J)+AKMOMEGA(I,J)+
     &0.5*AKBETAOM(I,J))*TT21(J+6*NDIM)
C------------------------------------------------------
C     JEDNACINA 4.199
C------------------------------------------------------
      F92(I+6*NDIM)=F92(I+6*NDIM)+(AKVV21(I,J)+AKVV22(I,J)+AKVV23(I,J))*
     &TT21(J)+(AKVV24(I,J)+AKVV25(I,J)+AKVV26(I,J))*TT21(J+NDIM)+
     &(AKVV27(I,J)+AKVV28(I,J)+AKVV29(I,J))*TT21(J+2*NDIM)
     
      IF (NSTAC.EQ.0) THEN
      II6=I+6*NDIM
      JJ6=J+6*NDIM
      F92(II6)=F92(II6)-AMK(I,J)*(TT21(JJ6)-TT210(JJ6))/TIME
      ENDIF

      F92(I+6*NDIM)=F92(I+6*NDIM)+FSOMEGA(I)
C------------------------------------------------------

 300  CONTINUE
 310  CONTINUE

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
      
      IF (IMUMPS.EQ.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN 
	IF (IMUMPS.eq.1) CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,7)
#if(MUMPS_CLUSTER)
	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,7)
#else
       CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
#endif
       CALL SPAKDE (SILE,F92,LM2,NDES)
	 GOTO 400
      ENDIF



	IF(NJUTRA.EQ.1.AND.ITER.GT.0) THEN
       CALL SPAKDE (SILE,F92,LM2,NDES)
	ELSE
       IF (ISYMMS.EQ.1) THEN
	    CALL MATSTE (ALEVO,MAXA,SILE,SKE,F92,LM2,NDES,1)
	 ELSE
          CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F92,MAXA,LM2,NDES,1)
	 ENDIF
	ENDIF


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
C=====================================================================

