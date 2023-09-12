C=========================================================================
C=========================================================================
      SUBROUTINE FORMPR(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,AF,INDEL,IPRESS,IPASS,PRITIS,AMASA,CPRESS,NWK1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine RACU2D is used for 2D analysis
CE It is used global loop per elements
C


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
	DIMENSION LMX(4),LMXT(4)
      
      DIMENSION AKV2P(9,4)
      DIMENSION AKV3P(9,4)



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
      ENDIF
  162 CONTINUE
C POVRSINSKE SILE I ZAPREMINSKE SILE:
C SURFACE FORCES AND BODY FORCES
  163 CONTINUE



C      IF (PENALT.GT.0.D0) THEN
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
      IF (N.LE.1) THEN
C SIGN IS CHANGED BECAUSE fidap 
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-ZVHX(K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-ZVHY(K)*HP(N))
       AKV2P(K,N)=AKV2P(K,N)+WDT*(ZVHX(K)*HP(N))
       AKV3P(K,N)=AKV3P(K,N)+WDT*(ZVHY(K)*HP(N))
      ENDIF
  185 CONTINUE
  190 CONTINUE
  195 CONTINUE
  200 CONTINUE
C      ENDIF
C===========================================================================
      DO I=1,1
	DO J=1,NDIM
	  CPRESS(J,1,NBREL)=AKV2P(J,I)
	  CPRESS(J,2,NBREL)=AKV3P(J,I)
	ENDDO
	ENDDO

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C end of loop per elements
C=======================================================================
 400  CONTINUE


C  	K=0
C      DO 10 NE1=1,NET
C       DO 20 NE2=NE1,1,-1
C	 C=0.D0
C	  DO 30 I=1,NDIM
C	  DO J=1,NDIM
C	   LMX(I)=ID(1,NEL(I,NE1)) 
C	   LMXT(J)=ID(1,NEL(J,NE2))
C	   
C         IF (LMX(I).EQ.LMXT(J).AND.LMX(I).NE.0) THEN
C	     C=C+CPRESS(I,1,NE1)*(1.D0/AMASA(1,LMX(I)))*CPRESS(J,1,NE2)	    
C	   ENDIF
C	  ENDDO
C
C	  DO J=1,NDIM
C	   LMX(I)=ID(2,NEL(I,NE1)) 
C	   LMXT(J)=ID(2,NEL(J,NE2))
C         
C	   IF (LMX(I).EQ.LMXT(J).AND.LMX(I).NE.0) THEN
C	     C=C+CPRESS(I,2,NE1)*(1.D0/AMASA(2,LMX(I)))*CPRESS(J,2,NE2)	    
C	   ENDIF
C	  ENDDO
C
C30       CONTINUE
C	   K=K+1
C	   IF(NE1.EQ.NE2) MAXA(NE1)=K
C	   ALEVO(K)=C
C 20   CONTINUE 
CC REDUCTION OF MAXA AND ALEVO
C        IK=K
C	IF (NE1.EQ.1) GOTO 10
C	 WRITE(IIZLAZ,*) 'MAXA(NE1) =',NE1,MAXA(NE1) 
C       DO I=IK,IK-MAXA(NE1),-1
C	   IF (ALEVO(I).EQ.0.D0) THEN 
C	     K=K-1
C	   ELSE
C	     GOTO 10
C	   ENDIF	  
C	 ENDDO
C
C 10   CONTINUE

C	MAXA(NET+1)=K+1
C	NWK1=K


C       RETURN




	K=0
      DO 10 NE1=1,NET
       DO 20 NE2=NE1,1,-1
	 C=0.D0
	  DO 30 I=1,NDIM
	  DO J=1,NDIM
	   LMX(I)=ID(1,NEL(I,NE1)) 
	   LMXT(J)=ID(1,NEL(J,NE2))
	   NODE=NEL(I,NE1)
         IF (LMX(I).EQ.LMXT(J).AND.LMX(I).NE.0) THEN
	     C=C+CPRESS(I,1,NE1)*(1.D0/AMASA(1,NODE))*CPRESS(J,1,NE2)	    
	   ENDIF
	  ENDDO

	  DO J=1,NDIM
	   LMX(I)=ID(2,NEL(I,NE1)) 
	   LMXT(J)=ID(2,NEL(J,NE2))
         
	   NODE=NEL(I,NE1)
	   IF (LMX(I).EQ.LMXT(J).AND.LMX(I).NE.0) THEN
	     C=C+CPRESS(I,2,NE1)*(1.D0/AMASA(2,NODE))*CPRESS(J,2,NE2)	    
	   ENDIF
	  ENDDO

30       CONTINUE
	   K=K+1
	   IF(NE1.EQ.NE2) MAXA(NE1)=K
	   ALEVO(K)=C
 20   CONTINUE 
C REDUCTION OF MAXA AND ALEVO
        IK=K
	IF (NE1.EQ.1) GOTO 10
	 WRITE(IIZLAZ,*) 'MAXA(NE1) =',NE1,MAXA(NE1) 
       DO I=IK,IK-MAXA(NE1),-1
	   IF (ALEVO(I).EQ.0.D0) THEN 
	     K=K-1
	   ELSE
	     GOTO 10
	   ENDIF	  
	 ENDDO

 10   CONTINUE

	MAXA(NET+1)=K+1
	NWK1=K

      
C End of subroutine
      END
C=======================================================================
