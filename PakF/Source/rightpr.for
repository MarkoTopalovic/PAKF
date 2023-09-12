C=========================================================================
C=========================================================================
      SUBROUTINE RIGHTPR(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,AF,INDEL,IPRESS,IPASS,PRITIS,AMASA,CPRESS,ID1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine RACU2D is used for 2D analysis
CE It is used global loop per elements
C

C SAMO PRIVREMENO UBACEN COMMON ZA AKIRIN PRIMER 
C STRUJANJA FLUIDA KROZ ELASTICNU CEV
C      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,DPRES
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)
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
      DIMENSION SPSIL(2,*),PRES(3,*),VMESH(3,*)
	DIMENSION IPRESS(*)
	DIMENSION PRITIS(2,*)
	DIMENSION AMASA(2,*)
	DIMENSION CPRESS(4,2,*)


      DIMENSION X(9),Y(9),LM2(4*9),TT21(4*9),TT210(4*9),TT21A(2*9)
      DIMENSION R(3,3),S(3,3),W(3,3)
C      DIMENSION R1(2),S1(2),W1(2)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION AKVV2(9,9)
      DIMENSION AKMIV2(9,9)
	DIMENSION AKU1(9,9),AKU_1(4,4)
      DIMENSION AMV2(9,9)
	DIMENSION F31(1),LM22(4),LMX(4)

      CALL CLEAR(DESNO,NET)


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
 100  CONTINUE

C	CALL MAXATF1(MAXA,MHT,ID,NEL,NET,NDIM,JEDN1,NWK1,5,NDIM,ID1)


      DO 400 NBREL=1,NET

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
      LM2(KLM)=ID(ID1,NEL(KLM,NBREL))
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
       AKVV2(K,N)=0.D0
       AKMIV2(K,N)=0.D0
       AMV2(K,N)=0.D0
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
C SURFACE FORCES AND BODY FORCES
  163 CONTINUE

      F31(1)=0.D0
      DO I=1,NDIM
	 DO J=1,NDIM
	  AKU_1(K,N)=0.D0
	 ENDDO
	ENDDO

C===========================================================================

C INTEGRACIJA U GAUSOVIM TACKAMA
C integration per gauss points

      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT

C subroutine FNTERP return to us interpolation functions...
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM

      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ)
       IF (IALE.EQ.1) CALL ALEHV(TT21,TT21A,HV2,HV3,NDIM,H)
C       IF (IALE.EQ.1) CALL ALEHV4(TT21,TT21A,HV2,HV3,NDIM,H,ZVXV2,
C     &ZVYV3,ZVYV2,ZVXV3,ZVHX,ZVHY)
       CALL AXISYF(INDAX,DEBLJ,X,H,NDIM)

      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

      AKVV2(K,N)=AKVV2(K,N)+WDT*((H(K)*HV2*ZVHX(N)+
     1H(K)*HV3*ZVHY(N))*GUSM)
      AKMIV2(K,N)=AKMIV2(K,N)+WDT*((ZVHX(K)*ZVHX(N)+
     1ZVHY(K)*ZVHY(N))*AMI)
	AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM
      

  164 CONTINUE
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE

C PAKOVANJE MATRICE C,Kk,Jk U MATRICU SKEF
C SKEF is local stiffness matrix

      DO K=1,NDIM
	DO N=1,NDIM
	  AKU1(K,N)=AMV2(K,N)/TIME+AKVV2(K,N)+AKMIV2(K,N)
	  AKU_1(K,K)=AKU_1(K,K)+DABS(AKU1(K,N))
	ENDDO
	ENDDO

C       CALL WRRF(AKU_1,16,'AKU_1=',IIZLAZ)
C=========================================================================
	IF (IPASS.EQ.0) THEN
	S=0.D0
	F=0.D0
	FX=0.D0
C LEFT AND RIGHT HAND SIDE
      DO 263 I=1,1
      DO 262 J=1,NDIM
      	F31(I)=F31(I)-AKU1(J,J)*TT21(J+NDIM*(ID1-1))
  262 CONTINUE
  263 CONTINUE

	   LM22(1)=IPRESS(NBREL)
	   NDES1=1
	   DESNO(NBREL)=0.D0
         CALL SPAKDE (DESNO,F31,LM22,NDES1)
       ENDIF
C=========================================================================

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C end of loop per elements
C=======================================================================
 400  CONTINUE


	DO NE1=1,NET
	  DO I=1,NDIM
	   LMX(I)=ID(ID1,NEL(I,NE1)) 
         IF (LMX(I).NE.0) THEN
	     SILE(NE1)=SILE(NE1)-CPRESS(I,ID1,NE1)*
     &(1.D0/AMASA(ID1,LMX(I)))*DESNO(LMX(I))
	   ENDIF
	  ENDDO
	ENDDO

C End of subroutine
      END
C=======================================================================
