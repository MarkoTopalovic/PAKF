C==========================================================================
C==========================================================================
C    SUBROUTINE  MOVMSH
C                MOVE2D
C                SOL2D
C                INTERS
C                MATSTE
C                RESENF
C                ULAZS1
C                MAXATE
C                PSKEFN
C                PRNTSS
C                ZADNOS
C                BOUND3
C                BOUND4
C                SHELLS
C==========================================================================
C==========================================================================
      SUBROUTINE MOVMSH(A,CCORD,NZAD,ZADVRE,NEL,IDSS,IDENT,IDS,IDF,
     &BRZ,TT1S,CORDF,NDIM,NET,LMAX,NP,NTOT,NPTF,NETIP,NPTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IDENT(2,*),IDS(NP,*),IDF(4,*),NEL(NDIM+1,*)
      DIMENSION BRZ(*),TT1S(*)
      DIMENSION CCORD(3,*),CORDF(3,*)
      DIMENSION IDSS(3,*)
	DIMENSION V(3),D(3),NQS(3),NQF(3),DD(3)
      DIMENSION NZAD(3,*),ZADVRE(*)
 
      DIMENSION A(*)    
      REAL A

      
      K=0
      DO 10 NODEI=1,NPTI 
	    NODES=IDENT(1,NODEI)
        NODEF=IDENT(2,NODEI)

	  DO I=1,NETIP
C WE PUT NODEI INSTEAD NODES BECAUSE FREE NUMERATION NODES
C        NQS(I)=IDS(NODES,I)
        NQS(I)=IDS(NODEI+50,I)

        NQF(I)=IDF(I,NODEF)
        V(I)=0.D0
        D(I)=0.D0
        
        
        IF(NQS(I).NE.0) THEN
C          IF (NSTAC.EQ.0) V(I)=BRZ(NQS(I))
          D(I)=TT1S(NQS(I))
        ENDIF
          DD(I)=D(I)-(CCORD(I,NODEF)-CORDF(I,NODEF))
          K=K+1
          ZADVRE(K)=DD(I)
          NZAD(1,K)=NODES
          NZAD(2,K)=I
          NZAD(3,K)=1
        ENDDO
 10    CONTINUE
   

      NUMZAD=K
      EEE=2.0D6
      ANI=0.3D0

      CALL MOVE2D(A,CCORD,NZAD,ZADVRE,NEL,IDSS,NDIM,NUMZAD,NET,LMAX,
     &NPTF,EEE,ANI,NTOT)
      
      END
C==========================================================================
C=======================================================================
      SUBROUTINE MOVE2D(A,CCORD,NZAD,ZADVRE,NEL,ID,NDIM,NUMZAD,NET,LMAX,
     &NPT,EEE,ANI,NTOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      PARAMETER (NTOT = 500000)
      COMMON /BROJFK/ INDFOR,NULAZ
      COMMON /CDEBUG/ IDEBUG

      DIMENSION NZAD(3,*),NEL(NDIM+1,*),ID(3,*)
      DIMENSION CCORD(3,*),ZADVRE(*)
      DIMENSION A(*)    
      REAL A
      
      ISOLID=31
      IIZLAZ=32
      IUNV=33
      MAXVEC=NTOT
      NETIP=2
      ISRPS=1
      INDFOR=2
      IDEBUG=0
C      EEE=2.0D6
C      ANI=0.3D0

C==========================================================================
C Opening the files
      OPEN(ISOLID,FILE='SOLID.DAT')
      OPEN(IIZLAZ,FILE='SOLID.OUT')
C      OPEN(IUNV,FILE='SOLID.UNV')
C==========================================================================
C Reading numberation of NODES,ID-matrix, COORDinates of nodes
      REWIND (ISOLID)
      READ(ISOLID,1000) NPT
      KK=0
      DO I=1,NPT
       READ(ISOLID,1001)N,(ID(J,N),J=1,3)
      DO JJ=1,3
       IF (ID(JJ,N).EQ.0) THEN
        KK=KK+1
        ID(JJ,N)=KK
       ELSE
        ID(JJ,N)=0
       ENDIF
      ENDDO
      ENDDO
      NEQ=KK
C============================
1001   FORMAT(I5,1X,3I2)
C      CALL MEMORY(LID,LMAX,3*NPT,1,MAXVEC,IIZLAZ)
C      CALL MEMORY(LCORD,LMAX,3*NPT,2,MAXVEC,IIZLAZ)
C      CALL ULAZS1(A(LID),A(LCORD),NPT,ISOLID,IIZLAZ,3,NEQ)
C==========================================================================
C Reading finite elements
C      READ(ISOLID,1000) NET,NDIM
C      CALL MEMORY(LNEL,LMAX,NDIM*NET,1,MAXVEC,IIZLAZ)
C      CALL ULAZF2(A(LNEL),NET,NETIP,NDIM,ISOLID,IIZLAZ,ISRPS)
C==========================================================================
C Reading prescribed values
C      READ(ISOLID,1000) NUMZAD
C      CALL MEMORY(LNZAD,LMAX,3*NUMZAD,1,MAXVEC,IIZLAZ)
C      CALL MEMORY(LZADVR,LMAX,NUMZAD,2,MAXVEC,IIZLAZ)
C      CALL ULAZF3(A(LNZAD),A(LZADVR),NUMZAD,ISOLID,IIZLAZ)
C==========================================================================
      CALL MEMORY(LMHT,LMAX,NEQ+1,1,MAXVEC,IIZLAZ)
      CALL MEMORY(LMAXA,LMAX,NEQ+1,1,MAXVEC,IIZLAZ)

      CALL MAXATE(A(LMAXA),A(LMHT),ID,NEL,NET,NDIM,NEQ,NWK,3)
C      CALL MAXATF(A(LMAXA),A(LMHT),A(LID),A(LNEL),NET,NDIM,NEQ,NWK,3,
C     &NDIM)

      LMAX0=LMAX
      NDES=2*NDIM
      MEM=0.5*NDES*(NDES+1)
      CALL MEMORY(LSKE,LMAX,MEM,2,NTOT,IIZLAZ)
      CALL MEMORY(LSKEF,LMAX,NDES*NDES,2,NTOT,IIZLAZ)
      CALL MEMORY(LALEVO,LMAX,NWK,2,NTOT,IIZLAZ)
      LDESNO=LMAX
C      CALL MEMORY(LDESNO,LMAX,NWK,2,NTOT,IIZLAZ)
C      CALL MEMORY(LGNODE,LMAX,NPT*3,2,NTOT,IIZLAZ)
      CALL MEMORY(LSILE,LMAX,NEQ,2,NTOT,IIZLAZ)
      CALL MEMORY(LTT1,LMAX,NEQ,2,NTOT,IIZLAZ)

      CALL CLEARR(A(LMAX0),LMAX-LMAX0+1)
      
      CALL SOL2D(CCORD,ID,CCORD,NEL,A(LMAXA),A(LTT1),
     &A(LALEVO),A(LSILE),NZAD,ZADVRE,A(LSKEF),A(LSKE),NPT,NET,
     &NDIM,NDES,NEQ,NUMZAD,EEE,ANI,A(LDESNO),IIZLAZ,NWK)


C       CALL RESEN(A(LALEVO),A(LSILE),A(LMAXA),NEQ,NWK,1)
C       CALL RESEN(A(LALEVO),A(LSILE),A(LMAXA),NEQ,NWK,2)
C  
C       CALL UACTCF(A(LALEVO),A(LDESNO),A(LSILE),A(LMAXA),NEQ,1)
C       CALL UACTCF(A(LALEVO),A(LDESNO),A(LSILE),A(LMAXA),NEQ,2)

       CALL PRNTSS(CCORD,A(LTT1),ID,NPT,IIZLAZ)


1000  FORMAT (2I5)

C End of program SOLID
      END
C=======================================================================
      SUBROUTINE SOL2D(GNODE,ID,CORD,NEL,MAXA,TT1,ALEVO,SILE,NZAD,
     &ZADVRE,SKEF,SKE,NPT,NET,NDIM,NDES,NEQ,NUMZAD,EEE,ANI,DESNO,
     &IIZLAZ,NWK)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION TT1(*),ALEVO(*),SILE(*),CORD(3,*),SKEF(NDES,*),SKE(*)
      DIMENSION GNODE(3,*),ZADVRE(*),DESNO(*)
      DIMENSION NEL(NDIM+1,*),ID(3,*),MAXA(*),NZAD(3,*)

      DIMENSION X(9),Y(9),TT21(18),B(4,18),BT(18,4)
      DIMENSION F36(18)
      DIMENSION R(3,3),S(3,3),W(3,3),DT(4,4)
      DIMENSION ZVHX(9),ZVHY(9),H(9)
      DIMENSION LM2(18)
      



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

C      CALL ZADNOS(GNODE,ZADVRE,NZAD,NUMZAD)



C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIM.EQ.4) IBRGT=2


CE MAIN LOOP OVER ELEMENTS
      DO 400 NBREL=1,NET

	DO I=1,2*NDIM
         LM2(I)=0
        ENDDO

       IAXIS=4

      EMNOZI=EEE*(1.D0-ANI)/((1.D0+ANI)*(1.D0-2.D0*ANI))
      DT(1,1)=1.D0*EMNOZI
      DT(1,2)=ANI*EMNOZI/(1.D0-ANI)
      DT(1,3)=0.D0
      DT(2,1)=DT(1,2)
      DT(2,2)=1.D0*EMNOZI
      DT(2,3)=0.D0
      DT(3,1)=DT(1,3)
      DT(3,2)=DT(2,3)
      DT(3,3)=EMNOZI*(1.D0-2.D0*ANI)/(2.D0*(1.D0-ANI))
      DT(1,4)=DT(1,2)
      DT(2,4)=DT(1,2)
      DT(3,4)=0.D0
      DT(4,4)=EMNOZI*1.D0
      DT(4,1)=DT(1,4)
      DT(4,2)=DT(2,4)
      DT(4,3)=DT(3,4)

      DO 125 I=1,NDES
C      TT210(I)=0.D0
 125  TT21(I)=0.D0
C=========================================================================
      JJ=-1
       DO 130 KLM=1,NDIM
       NODE=NEL(KLM,NBREL)
       JJ=JJ+2
       DO NR=1,2
       IF (ID(NR,NEL(KLM,NBREL)) .NE. 0) THEN
C        TT21(JJ+NR-1)=GNODE(NR,NODE)
        TT21(JJ+NR-1)=TT1(ID(NR,NEL(KLM,NBREL)))
       ENDIF
      ENDDO
       X(KLM)=CORD(1,NODE)
       Y(KLM)=CORD(2,NODE)
       JX=ID(1,NODE)
       JY=ID(2,NODE)
       LM2(JJ)=JX
       LM2(JJ+1)=JY
C       TT21(KLM)=GNODE(1,NODE)        
C       TT21(KLM+NDIM)=GNODE(2,NODE)        
C       IF (JX.NE.0) TT21(KLM)=TT1(JX)
C       IF (JY.NE.0) TT21(KLM+NDIM)=TT1(JY)
 130  CONTINUE

C INCIJALIZACIJA MATRICE SKEF I F36
      DO 260 I=1,NDES
      DO 258 J=1,NDES
       SKEF(I,J)=0.D0
 258  CONTINUE
       F36(I)=0.D0
 260  CONTINUE

       DO I=1,NDES*(NDES+1)*0.5
        SKE(I)=0.D0
       ENDDO


C INTEGRACIJA U GAUSOVIM TACKAMA
      NGAUS=0
      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT
      NGAUS=NGAUS+1
      CALL INTERS(R(IBRGT,I),S(IBRGT,J),TT21,DETJ,X,Y,B,BT,NDIM,NBREL,
     &ZVHX,ZVHY,H,BT)
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

C        CALL CLEAR(STRAIN,4) 
C        CALL CLEAR(TAU,4) 

    
C       DO II=1,IAXIS
C        DO NN=1,NDIM*2
C         STRAIN(II)=STRAIN(II)+B(II,NN)*TT21(NN)
C        ENDDO
C       ENDDO
      

      

      DO 165 K=1,NDIM*2
      DO 164 N=1,NDIM*2

C Calculation of stifness matrix Kuu
      DO II=1,IAXIS
      DO JJ=1,IAXIS
       SKEF(K,N)=SKEF(K,N)+BT(K,II)*DT(II,JJ)*B(JJ,N)*WDT
      ENDDO
      ENDDO
  164 CONTINUE
  165 CONTINUE
  170 CONTINUE
  180 CONTINUE








C=========================================================================

C       DO I=1,NDES
C        DO J=1,NDES
C         F36(I)=F36(I)-SKEF(I,J)*TT21(J)
C        ENDDO
C       ENDDO       


       CALL PSKEFN(SKEF,SKE,NDES)	  
       CALL MATSTE (ALEVO,MAXA,SILE,SKE,F36,LM2,NDES,1)

C      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F36,MAXA,LM2,NDES,1)

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE


      DO 405 I=1,NUMZAD
       JJ=ID(NZAD(2,I),NZAD(1,I))       
       IF (JJ.NE.0) THEN
        SILE(JJ)=ZADVRE(I)*1.0D35
        ALEVO(MAXA(JJ))=1.0D35
       ENDIF
 405   CONTINUE

C       DO 2423 I=1,NWK
C 2423    WRITE(IIZLAZ,*)I,'ALEVO',ALEVO(I)

C      DO 424 I=1,NEQ
C 424    WRITE(IIZLAZ,*)I,'SILE',SILE(I)

      CALL RESENF(ALEVO,SILE,MAXA,NEQ,NWK,1)
      CALL RESENF(ALEVO,SILE,MAXA,NEQ,NWK,2)

C        CALL UACTCF(ALEVO,DESNO,SILE,MAXA,NEQ,1)
C        CALL UACTCF(ALEVO,DESNO,SILE,MAXA,NEQ,2)

       DO I=1,NEQ
        TT1(I)=TT1(I)+SILE(I)
       ENDDO 
      

      END

C=======================================================================
C======================================================================
      SUBROUTINE INTERS(R,S,TT21,DETJ,X,Y,B,BT,NDIM,NBREL,ZVHX,ZVHY,H,
     &BS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DIMENSION X(*),Y(*),TT21(*),B(3,*),BT(3*9,*),BS(2,*)
      DIMENSION H(*),ZVHX(*),ZVHY(*),ZVHR(9),ZVHS(9),AJ(2,2)


      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIM.GT.4) THEN
      IF (NDIM.EQ.8) THEN
        H(9)=0.D0
      ELSE
        H(9)=RR*SS
      ENDIF
      H(8)=0.5*RP*SS-0.5*H(9)
      H(7)=0.5*RR*SM-0.5*H(9)
      H(6)=0.5*RM*SS-0.5*H(9)
      H(5)=0.5*RR*SP-0.5*H(9)
      H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
      H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
      H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
      H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSE
      H(4)=0.25*RP*SM
      H(3)=0.25*RM*SM
      H(2)=0.25*RM*SP
      H(1)=0.25*RP*SP

      ENDIF



      IF (NDIM.GT.4) THEN
      IF (NDIM.EQ.8) THEN
        ZVHR(9)=0.D0
        ZVHS(9)=0.D0
      ELSE
        ZVHR(9)=-2.D0*R*SS
        ZVHS(9)=-2.D0*S*RR
      ENDIF
      ZVHR(8)=0.5*SS-0.5*ZVHR(9)
      ZVHS(8)=-RP*S-0.5*ZVHS(9)
      ZVHR(7)=-R*SM-0.5*ZVHR(9)
      ZVHS(7)=-0.5*RR-0.5*ZVHS(9)
      ZVHR(6)=-0.5*SS-0.5*ZVHR(9)
      ZVHS(6)=-RM*S-0.5*ZVHS(9)
      ZVHR(5)=-R*SP-0.5*ZVHR(9)
      ZVHS(5)=0.5*RR-0.5*ZVHS(9)
      ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))-0.25*ZVHS(9)
      ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))-0.25*ZVHR(9)
      ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))-0.25*ZVHS(9)
      ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))-0.25*ZVHR(9)
      ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))-0.25*ZVHS(9)
      ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))-0.25*ZVHS(9)

      ELSE

      ZVHR(4)=0.25*SM
      ZVHS(4)=-0.25*RP
      ZVHR(3)=-0.25*SM
      ZVHS(3)=-0.25*RM
      ZVHR(2)=-0.25*SP
      ZVHS(2)=0.25*RM
      ZVHR(1)=0.25*SP
      ZVHS(1)=0.25*RP

      ENDIF



      AJ(1,1)=DOT(ZVHR,X,NDIM)
      AJ(1,2)=DOT(ZVHR,Y,NDIM)
      AJ(2,1)=DOT(ZVHS,X,NDIM)
      AJ(2,2)=DOT(ZVHS,Y,NDIM)


      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
       WRITE(*,*)'DETERMINANT LESS THEN ZERO!!!'
       WRITE(*,*)'FOR ELEMENT NUMBER ',NBREL
      STOP
      ENDIF

      DO 20 I=1,NDIM
       ZVHX(I)=(ZVHR(I)*AJ(2,2)-ZVHS(I)*AJ(1,2))/DETJ
       ZVHY(I)=(ZVHS(I)*AJ(1,1)-ZVHR(I)*AJ(2,1))/DETJ
  20  CONTINUE  

      

      DO J=1,3*NDIM
        BS(1,J)=0.D0
        BS(2,J)=0.D0
       DO I=1,3
        B(I,J)=0.D0
       ENDDO
      ENDDO


C BATHE:
C      JG=0
C      DO I=1,NDIM
C       IG=JG+1
C       BS(1,IG)=ZVHY(I)
C       BS(2,IG)=ZVHX(I)
C       IG=IG+1
C       JG=IG+1
C       BS(2,IG)=H(I)
C       BS(1,JG)=-H(I)
C        B(1,IG)=ZVHX(I)
C        B(3,IG)=ZVHY(I)
C        B(2,JG)=-ZVHY(I)
C        B(3,JG)=-ZVHX(I)
C      ENDDO


C HINTON:
      JG=0
      DO I=1,NDIM
       IG=JG+1
       BS(1,IG)=ZVHX(I)
       BS(2,IG)=ZVHY(I)
       IG=IG+1
       JG=IG+1
       BS(1,IG)=-H(I)
       BS(2,JG)=-H(I)
        B(1,IG)=-ZVHX(I)
        B(3,IG)=-ZVHY(I)
        B(2,JG)=-ZVHY(I)
        B(3,JG)=-ZVHX(I)
      ENDDO



      DO I=1,3
       DO J=1,3*NDIM
        BT(J,I)=B(I,J)
       ENDDO
      ENDDO

C End of subroutine INTERS
      END
C======================================================================


C=======================================================================
      SUBROUTINE MATSTE(SK,MAXA,F,SKE,FE,LM,NCV,INDSK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS  PAKOVANJE MATRICA I VEKTORA ELEMENATA U MATRICE I VEKTORE SISTEMA
CE  INSERTING ELEMENT MATRIXES AND VECTORS
CE  INTO SYSTEM MATRIXES AND VECTORS
      DIMENSION SK(*),MAXA(*),F(*),SKE(*),FE(*),LM(*)
C
      K=0
      DO 200 I=1,NCV
      IVR=LM(I)
      IF(INDSK.EQ.0) GO TO 110
      DO 100 J=I,NCV
      K=K+1
      KOL=LM(J)
C     IF(IVR.EQ.0.OR.KOL.EQ.0) GO TO 100
      IF(IVR.LE.0.OR.KOL.LE.0) GO TO 100
      IF(IVR-KOL) 10,10,20
   10 KS=MAXA(KOL) + KOL - IVR
      GO TO 50
   20 KS=MAXA(IVR)+IVR-KOL
   50 SK(KS)=SK(KS)+SKE(K)
C
  100 CONTINUE
C
C 110 IF(IVR.EQ.0) GO TO 200
  110 IF(IVR.LE.0) GO TO 200
      F(IVR)=F(IVR)+FE(I)
C
  200 CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE RESENF(A,V,MAXA,NN,NWK,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
C        WRITE(3,*) 'NN,NWK',NN,NWK
C        CALL IWRR(MAXA,NN,'JDIA')
C        CALL WRR(V,NN,'V   ')
C        CALL WRR(A,NWK,'AA  ')
C
CS     L*D*L(T) FAKTORIZACIJA
CE     L*D*L(T) FACTORIZATION
C
      IF(KKK-2)40,150,150
   40 DO 140 N=1,NN
      KN=MAXA(N)
      KL=KN+1
      KU=MAXA(N+1)-1
      KH=KU-KL
      IF(KH)110,90,50
   50 K=N-KH
      IC=0
      KLT=KU
      DO 80 J=1,KH
      IC=IC+1
      KLT=KLT-1
      KI=MAXA(K)
      ND=MAXA(K+1)-KI-1
      IF(ND)80,80,60
   60 KK=MIN0(IC,ND)
      C=0.D0
      DO 70 L=1,KK
   70 C=C+A(KI+L)*A(KLT+L)
      A(KLT)=A(KLT)-C
   80 K=K+1
   90 K=N
      B=0.
      DO 100 KK=KL,KU
      K=K-1
      KI=MAXA(K)
C OVDE UBACENO:
      IF (DABS(A(KI)).LE.1.D-15) A(KI)=1.D-15
C      IF (A(KI).EQ.0.) A(KI)=1.D-15
C
      C=A(KK)/A(KI)
      B=B+C*A(KK)
  100 A(KK)=C
      A(KN)=A(KN)-B
  110 IF(A(KN).LE. 0.) THEN
        IF(ISRPS.EQ.0)
     *  WRITE(*,2000)N,A(KN)
        IF(ISRPS.EQ.1)
     *  WRITE(*,6000)N,A(KN)
C        STOP
      ENDIF
  140 CONTINUE
      RETURN
C
CS     REDUKOVANJE SLOBODNOG VEKTORA
CE     FORWARD REDUCTION
C
  150 DO 180 N=1,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF(KU-KL)180,160,160
  160 K=N
      C=0.
      DO 170 KK=KL,KU
      K=K-1
  170 C=C+A(KK)*V(K)
      V(N)=V(N)-C
  180 CONTINUE
C
CS     ZAMENA UNAZAD
CE     BACK SUBSTITUTION
C
      DO 200 N=1,NN
      K=MAXA(N)
  200 V(N)=V(N)/A(K)
      IF(NN.EQ.1)RETURN
      N=NN
      DO 230 L=2,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF (KU-KL)230,210,210
  210 K=N
      DO 220 KK=KL,KU
      K=K-1
  220 V(K)=V(K)-A(KK)*V(N)
  230 N=N-1
      RETURN
C
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'
     1//' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2D20.12)
 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITION'
     1//' ','PIVOT IS NOT POSITIVE FOR EQUATION NUM.',I4,//' ','PIVOT=',
     2D20.12)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZS1(ID,CORD,NPT,IULAZ,IIZLAZ,IDIM,NEQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CHARACTER*250ACOZ
      DIMENSION ID(IDIM,*),CORD(3,*)

      KK=0
      DO 10 I=1,NPT
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1006) N,(ID(II,N),II=1,IDIM),(CORD(J,N),J=1,3)
      DO JJ=1,IDIM
       IF (ID(JJ,N).EQ.0) THEN
        KK=KK+1
        ID(JJ,N)=KK
       ELSE
        ID(JJ,N)=0
       ENDIF
      ENDDO

   10 CONTINUE
      NEQ=KK

C 1005 FORMAT(I5,4(3X,I2),3F10.6,2I2)
 1006 FORMAT(I5,1X,3(I2),2X,3F10.6,I5)
      WRITE(IIZLAZ,2000)
 2000 FORMAT(//
     *11X,'REDNI BROJEVI,OGRANICENJA I KOORDINATE CVOROVA'/)
      DO 12 I=1,NPT
      WRITE(IIZLAZ,1006) I,(ID(II,I),II=1,IDIM),(CORD(J,I),J=1,3)
   12 CONTINUE
      END
C==========================================================================
C==========================================================================
      SUBROUTINE MAXATE(MAXA,MHT,ID,NEL,NE,NTE,JEDN,NWK,IDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    PODPROGRAM ZA FORMIRANJE VEKTORA VISINA STUBOVA I MAXA
CS    KONACNO SE SMESTAJU U ISTI PROSTOR
CE    PROGRAM TO DETERMINE COLUMN HEIGHTS VECTOR AND MAXA
C
C
      DIMENSION MAXA(*),MHT(*),NEL(NTE+1,*),LM(63),ID(IDIM,*)
C
CS    PETLJA PO ELEMENTIMA
CE    ELEMENT LOOP
C

      DO I=1,63
       LM(I)=0
      ENDDO

      DO I=1,JEDN+1
       MHT(I)=0
       MAXA(I)=0
      ENDDO

      DO 100 NLM=1,NE
         KK=0
         DO 2 I=1,NTE
            IF(NEL(I,NLM).EQ.0) GO TO 2
            N=NEL(I,NLM)
               DO 1 J=1,IDIM
                  IF(ID(J,N).LE.0) GO TO 1
                  KK=KK+1
                  LM(KK)=ID(J,N)
    1          CONTINUE
C            ENDIF
    2    CONTINUE
C
         LS=JEDN+1
         DO 10 I=1,KK
            IF (LM(I).LT.LS) LS=LM(I)
   10    CONTINUE

C
         DO 20 I=1,KK
            II=LM(I)
            ME=II-LS
            IF(ME.GT.MHT(II)) MHT(II)=ME
   20    CONTINUE
C
  100 CONTINUE
C
CS    VEKTOR MAXA
CE    VECTOR MAXA
C
      MAXA(1)=1
      MAXA(2)=2
      DO 200 I=2,JEDN
       MAXA(I+1)=MAXA(I)+MHT(I)+1
 200  CONTINUE
      NWK=MAXA(JEDN+1)-1
      LS = JEDN+1
      DO 210 I=1,LS
 210  MHT(I)=MAXA(I)

      RETURN
      END
C==========================================================================
C=========================================================================
      SUBROUTINE PSKEFN(SKEF,SKEFN,NDES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION SKEF(NDES,*),SKEFN(*)

	 K=0
       DO I=1,NDES
        DO J=I,NDES
          K=K+1
         	SKEFN(K)=SKEF(I,J)
        ENDDO
       ENDDO

      END
C==========================================================================
C=========================================================================
      SUBROUTINE PRNTSS(GNODE,TT1,ID,NPT,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(3,*),TT1(*),ID(3,*)


       WRITE(IIZLAZ,*)'RESULTS '
       DO NODE=1,NPT
         JX=ID(1,NODE)
         JY=ID(2,NODE)
         JZ=ID(3,NODE)
        IF (JX.NE.0) GNODE(1,NODE)=GNODE(1,NODE)+TT1(JX)
        IF (JY.NE.0) GNODE(2,NODE)=GNODE(2,NODE)+TT1(JY)
        IF (JZ.NE.0) GNODE(3,NODE)=GNODE(3,NODE)+TT1(JZ)
        WRITE(IIZLAZ,1000)NODE,(GNODE(J,NODE),J=1,3)
       ENDDO

1000   FORMAT(I5,3(D13.5))

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ZADNOS(GNODE,ZADVRE,NZAD,NUMZAD,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(3,*),ZADVRE(*)
       DIMENSION NZAD(3,*)

C
CE Subroutine ZADNOS is used for inclusion prescribed values
C

      
      DO 425 I=1,NUMZAD
        GNODE(NZAD(2,I),NZAD(1,I))=ZADVRE(I)
  425  CONTINUE 
      
      END
C==========================================================================
C==========================================================================
      SUBROUTINE BOUND3(CORD,NEL,NDIM,NPT,NET,NODE,NZAD,ZADVRE,NUMZAD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*),NZAD(3,*)
	DIMENSION ZADVRE(*)

C
CE Subroutine BOUND3 is used for definition of boundary constraints
CE for 3D walls
C
      
      DO I=1,NPT
        NODE(I)=0
      ENDDO

      DO 20 I=1,NET
       DO 10 J=1,NDIM
        NODE(NEL(J,I))=NODE(NEL(J,I))+1
   10  CONTINUE 
   20  CONTINUE 


      OPEN(100,FILE='OGRANIC.DAT')

      zmax=-1.d10
      do i=1,npt
	  if (cord(3,i).gt.zmax) zmax=cord(3,i)
	enddo


      DO 30 I=1,NPT
	   IX=1
	   IY=1
	   IZ=1

        IF (NODE(I).eq.8) THEN
	     ix=0
	     iy=0
	     iz=0
	  endif

c         if (cord(3,i).eq.zmax) then
c	     IX=0
c	     IY=0
c	     IZ=0
c	   endif

c	  IF (CORD(2,I).LT.1.D-8) then 
c	  ix=0
c	  IY=1
c	  endif
c	  IF (dabs(CORD(1,I)-0.1d0).LT.1.D-8) then 
c	  ix=0
c	  IY=0
c	  endif
      DO J=1,NUMZAD
	  IF (NZAD(1,J).EQ.I) THEN 
	   IF (ZADVRE(J).GT.1.D0) THEN
	     IX=0
	     IY=0
	     IZ=0
	   ENDIF
	  ENDIF
	ENDDO
C          WRITE(100,100)I,ix,iy,iz,1,1,(CORD(J,I)*0.1d0,J=1,3),NODE(I)
          WRITE(100,100)I,ix,iy,iz,1,1,(CORD(J,I),J=1,3),NODE(I)

 30   CONTINUE
 
 
      
  100  FORMAT(6I5,3F10.7,I5)

       CLOSE(100)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE BOUND3DFACE(CORD,NEL,NDIM,NPT,NET,NODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*)
	dimension ITR(6,4),ITR1(6,4),N(8),N1(8)

C
CE Subroutine BOUND3 is used for definition of boundary constraints
CE for 3D walls
C
      tol=1.e-6
      
      DO I=1,NPT
        NODE(I)=0
      ENDDO

 

      OPEN(100,FILE='OGRANIC.DAT')
      
      DO nbrel=1,(net-1)
       
	DO jj=1,NDIM
	  N(jj)=NEL(jj,NBREL)
      ENDDO

      ITR(1,1)=N(8)
      ITR(1,2)=N(4)
      ITR(1,3)=N(5)
      ITR(1,4)=N(1)

      ITR(2,1)=N(7)
      ITR(2,2)=N(3)
      ITR(2,3)=N(6)
      ITR(2,4)=N(2)

      ITR(3,1)=N(6)
      ITR(3,2)=N(2)
      ITR(3,3)=N(5)
      ITR(3,4)=N(1)

      ITR(4,1)=N(7)
      ITR(4,2)=N(3)
      ITR(4,3)=N(8)
      ITR(4,4)=N(4)

      ITR(5,1)=N(3)
      ITR(5,2)=N(2)
      ITR(5,3)=N(4)
      ITR(5,4)=N(1)

      ITR(6,1)=N(7)
      ITR(6,2)=N(6)
      ITR(6,3)=N(8)
      ITR(6,4)=N(5)


	 do 20 i=1,6
	 ind=0 

       tx=0.d0
	 ty=0.d0
	 tz=0.d0
	  do j=1,4
         tx=tx+cord(1,itr(i,j))
         ty=ty+cord(2,itr(i,j))
         tz=tz+cord(3,itr(i,j))
	  enddo


   
       DO 10 nbrel1=1,net

	if (nbrel1.eq.nbrel) goto 10

	DO jj=1,NDIM
	  N1(jj)=NEL(jj,NBREL1)
      ENDDO




      ITR1(1,1)=N1(8)
      ITR1(1,2)=N1(4)
      ITR1(1,3)=N1(5)
      ITR1(1,4)=N1(1)

      ITR1(2,1)=N1(7)
      ITR1(2,2)=N1(3)
      ITR1(2,3)=N1(6)
      ITR1(2,4)=N1(2)

      ITR1(3,1)=N1(6)
      ITR1(3,2)=N1(2)
      ITR1(3,3)=N1(5)
      ITR1(3,4)=N1(1)

      ITR1(4,1)=N1(7)
      ITR1(4,2)=N1(3)
      ITR1(4,3)=N1(8)
      ITR1(4,4)=N1(4)

      ITR1(5,1)=N1(3)
      ITR1(5,2)=N1(2)
      ITR1(5,3)=N1(4)
      ITR1(5,4)=N1(1)

      ITR1(6,1)=N1(7)
      ITR1(6,2)=N1(6)
      ITR1(6,3)=N1(8)
      ITR1(6,4)=N1(5)

       do ii=1,6 
	 
	      
	  tx1=0.d0
	  ty1=0.d0
	  tz1=0.d0
	  do jj=1,4
         tx1=tx1+cord(1,itr1(ii,jj))
         ty1=ty1+cord(2,itr1(ii,jj))
         tz1=tz1+cord(3,itr1(ii,jj))
	  enddo

	if(dabs(tx-tx1).lt.tol.and.dabs(ty-ty1).lt.tol.and.
     &dabs(tz-tz1).lt.tol) then 
         ind=1
c	   do j=1,4
c	     node(itr(i,j))=0
c	   enddo
c	write(100,*)'jeste zajednicka strana',i,nbrel
c	write(100,*) (node(itr(i,j)),j=1,4)
	   goto 20
	  endif
       enddo
10    continue
      
c	if (ind.eq.1) then

c	endif
c	if (ind.eq.0) write(100,*)'nije zajednicka strana',i,nbrel

	if (ind.eq.0) then
	   do j=1,4
	     node(itr(i,j))=1
	   enddo
	endif
      
20    continue
      enddo

      
c      do i=1,npt
c	  write(100,'(i10)') i
c	  write(100,'(2(e13.5))') dexp(-50*cord(1,i)**2),0.0
c	enddo
c     stop


      DO 30 I=1,NPT
	   IX=0
	   IY=0
	   IZ=0

        IF (NODE(I).eq.1) THEN
	     ix=1
	     iy=1
	     iz=1
	  endif

       WRITE(100,100)I,ix,iy,iz,1,1,(CORD(J,I),J=1,3),NODE(I)

 30   CONTINUE
      
  100  FORMAT(6I5,3F10.7,I5)

       CLOSE(100)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE FLUX3DFACE(CORD,NEL,NDIM,NPT,NET,ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*)
	dimension ITR(6,4),N(8)

C
CE Subroutine FLUX3DFACE is used for definition of boundary constraints
CE for 3D walls
C

      OPEN(100,FILE='FLUX.DAT')

	DO I=1,NPT
c	 WRITE(100,1000) I,ID(1,I),ID(2,I),ID(3,I)
       x=cord(1,i)*300.d0
       y=cord(2,i)*300.d0
       z=cord(3,i)*600.d0
	 WRITE(100,2000) x,y,z
	ENDDO
	STOP
      
      DO nbrel=1,net
       
	DO jj=1,NDIM
	  N(jj)=NEL(jj,NBREL)
      ENDDO

      ITR(1,1)=N(8)
      ITR(1,2)=N(4)
      ITR(1,3)=N(5)
      ITR(1,4)=N(1)

      ITR(2,1)=N(7)
      ITR(2,2)=N(3)
      ITR(2,3)=N(6)
      ITR(2,4)=N(2)

      ITR(3,1)=N(6)
      ITR(3,2)=N(2)
      ITR(3,3)=N(5)
      ITR(3,4)=N(1)

      ITR(4,1)=N(7)
      ITR(4,2)=N(3)
      ITR(4,3)=N(8)
      ITR(4,4)=N(4)

      ITR(5,1)=N(3)
      ITR(5,2)=N(2)
      ITR(5,3)=N(4)
      ITR(5,4)=N(1)

      ITR(6,1)=N(7)
      ITR(6,2)=N(6)
      ITR(6,3)=N(8)
      ITR(6,4)=N(5)


	DO 20 I=1,6
	 
	  N1=ITR(I,1) 
	  N2=ITR(I,2) 
	  N3=ITR(I,3) 
	  N4=ITR(I,4) 
	  IF (ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.AND.
     &ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.AND.
     &ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.AND.
     &ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0) THEN
	 Z=CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4)
	IF (Z.LT.1.D-5) GOTO 20
	    WRITE(100,1000) NBREL,N1,N2,N3,N4,2
	  ENDIF
   20 CONTINUE
	ENDDO
	  

1000  FORMAT(6I10)
2000  FORMAT(3(e13.5,1x))
   

       CLOSE(100)
	STOP

      END
C==========================================================================
C==========================================================================
      SUBROUTINE MERGEN(CORD,NEL,NDIM,NPT,NET,NODE,NIZ,NIZNEW,IIZLAZ,
     &NZAD,ZADVRE,NUMZAD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*)

      DIMENSION NIZ(*),NIZNEW(*)
      DIMENSION NZAD(3,*),ZADVRE(*)

C
CE Subroutine MERGE is used for merge of nodes
C

      TOL=1.D-3

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CORD(1,I)-CORD(1,J)).LT.TOL .AND.
     & DABS(CORD(2,I)-CORD(2,J)).LT.TOL.AND.
     & DABS(CORD(3,I)-CORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      

c      K=0
c	DO I=1,NPT
c	 NIZNEW(I)=0
c	IF (NIZ(I).GT.0) THEN
c	  K=K+1
c	DO J=1,3
c	  CORD(J,K)=CORD(J,I)
c	ENDDO
c	NIZNEW(I)=K
c	ENDIF
c	ENDDO



c      DO I=1,NPT
c	 IF (NIZNEW(I).EQ.0) NIZNEW(I)=NIZNEW(-NIZ(I))
c	write(iz,3000) i,niznew(i)
c	ENDDO




      DO I=1,NET
	 DO J=1,NDIM
	   NEL(J,I)=ABS(NIZNEW(NEL(J,I)))
	 ENDDO
	ENDDO
      

	WRITE(IZ,*)'C NODES' 
	DO I=1,NPT
	 IF (NIZNEW(I).GT.0) THEN
	  WRITE(IZ,1000) ABS(NIZNEW(I)),1,1,1,1,1,(CORD(J,I),J=1,3)
	 ENDIF
	ENDDO

      NPT=K
      
	WRITE(IZ,*)'C ELEMENTS' 
	DO I=1,NET
	 WRITE(IZ,2000) I,(NEL(J,I),J=1,NDIM),1
	ENDDO

	WRITE(IZ,*)'C VELOCITIES' 
	DO I=1,NUMZAD
	 WRITE(IZ,2100) ABS(NIZNEW(NZAD(1,I))),
     & NZAD(2,I),NZAD(3,I),ZADVRE(I)
	ENDDO
 
      goto 500

      WRITE(IZ,*)'9425 = ',ABS(NIZNEW(9425))
      WRITE(IZ,*)'10721 = ',ABS(NIZNEW(10721))
      WRITE(IZ,*)'19825 = ',ABS(NIZNEW(19825))
      WRITE(IZ,*)'25346 = ',ABS(NIZNEW(25346))
      WRITE(IZ,*)'30854 = ',ABS(NIZNEW(30854))
      WRITE(IZ,*)'29571 = ',ABS(NIZNEW(29571))
      WRITE(IZ,*)'36379 = ',ABS(NIZNEW(36379))
      WRITE(IZ,*)'37099 = ',ABS(NIZNEW(37099))
      WRITE(IZ,*)'41254 = ',ABS(NIZNEW(41254))
      WRITE(IZ,*)'43924 = ',ABS(NIZNEW(43924))
      WRITE(IZ,*)'44644 = ',ABS(NIZNEW(44644))
      WRITE(IZ,*)'48799 = ',ABS(NIZNEW(48799))
      WRITE(IZ,*)'59269 = ',ABS(NIZNEW(59269))
      WRITE(IZ,*)'59989 = ',ABS(NIZNEW(59989))
      WRITE(IZ,*)'73314 = ',ABS(NIZNEW(73314))

c	GOTO 500

	WRITE(IZ,*)'C PRESCRIBED VELOCITY' 
	NSLOJ=29*25*13+15*81
	WRITE(IZ,2000) 81,9,9
	DO I=NSLOJ+1,NSLOJ+81
	 WRITE(IZ,2000) ABS(NIZNEW(I))
	ENDDO

	NSLOJ=NSLOJ+81+45*25*13+35*17*9
	WRITE(IZ,2000) 153,17,9
	DO I=NSLOJ+1,NSLOJ+153
	 WRITE(IZ,2000) ABS(NIZNEW(I))
	ENDDO

	NSLOJ=NSLOJ+153+17*25*13+15*9*5
	WRITE(IZ,2000) 45,9,5
	DO I=NSLOJ+1,NSLOJ+45
	 WRITE(IZ,2000) ABS(NIZNEW(I))
	ENDDO

	NSLOJ=NSLOJ+45+21*25*13+15*45
	WRITE(IZ,2000) 45,9,5
	DO I=NSLOJ+1,NSLOJ+45
	 WRITE(IZ,2000) ABS(NIZNEW(I))
	ENDDO


	NSLOJ=NSLOJ+45+45*25*13+15*45
	WRITE(IZ,2000) 45,9,5
	DO I=NSLOJ+1,NSLOJ+45
	 WRITE(IZ,2000) ABS(NIZNEW(I))
	ENDDO

	NSLOJ=NSLOJ+45+40*25*13
	WRITE(IZ,2000) 25*13,25,13
	DO I=NSLOJ+1,NSLOJ+25*13
	 WRITE(IZ,2000) ABS(NIZNEW(I))
	ENDDO

 500     CLOSE (IZ)

c      STOP


c1000  FORMAT(I7,5I5,3F10.7)
c2000  FORMAT(10I7)
1000  FORMAT(I5,5I5,3F10.7)
2000  FORMAT(10I5)
3000  FORMAT(5I10)
2100  FORMAT(3I5,F10.7)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE BOUN33(CORD,NEL,NDIM,NPT,NET,NODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*)

C
CE Subroutine BOUND3 is used for definition of boundary constraints
CE for 3D walls
C
      
      DO I=1,NPT
        NODE(I)=0
      ENDDO

      DO 20 I=1,NET
       DO 10 J=1,NDIM
        NODE(NEL(J,I))=NODE(NEL(J,I))+1
   10  CONTINUE 
   20  CONTINUE 


      OPEN(100,FILE='OGRANIC.DAT')

      DO 30 I=1,NPT
	   IT=0
        if (i.GT.npt-1594) IT=1	  
	  IX=0
	  IY=0
	  IZ=0
	  IF (DABS(CORD(3,I)).LT.-1.D5) THEN
c	  IF (DABS(CORD(3,I)).LT.1.D-5) THEN
	     IF (NODE(I).EQ.4) THEN
              IX=0
	        IY=0
	        IZ=1
	     ELSE
              IX=1
	        IY=1
              IZ=1
    		 ENDIF
	   
	   ELSE
	     IF (NODE(I).LT.8) THEN
	       IX=1
	       IY=1
	       IZ=1
	     ENDIF
	   ENDIF

C	  X=CORD(1,I)*0.425/100.D0
C	  Y=CORD(2,I)*0.425/100.D0
C	  Z=CORD(3,I)*0.425/100.D0

	  X=CORD(1,I)
	  Y=CORD(2,I)
	  Z=CORD(3,I)

C         y0=-0.4d0

C       if (i.le.66*8) then
C	  y1=(cord(2,(mod((i-1),66)+1+66*8))-y0)/8.d0
C	  Y=y0+((i-1)/66)*y1
C	 endif

	  dod=1.D0
C	  isloj=190
C	  if (i.gt.17*91.and.i.le.25*91) dod=0.965D0
C	  if (i.gt.25*91.and.i.le.(25+8)*91) dod=0.918D0
C	  if (i.gt.(25+8)*91.and.i.le.(25+24)*91) dod=0.965D0
C	  if (i.gt.(25+24)*91.and.i.le.(74)*91) dod=0.918D0

C	  if (i.gt.21*isloj.and.i.le.27*isloj) dod=0.965D0
C	  if (i.gt.27*isloj.and.i.le.(27+10)*isloj) dod=0.918D0
C	  if (i.gt.(27+10)*isloj.and.i.le.(27+26)*isloj) dod=0.965D0
C	  if (i.gt.(27+26)*isloj.and.i.le.(80)*isloj) dod=0.918D0

C	  if (i.gt.(80)*isloj.and.i.le.(83)*isloj) dod=0.918D0

C	  X=X*dod
C	  Z=Z*dod
       
c	 WRITE(100,100)I,IX,IY,IZ,1,1,X*1.d-4,Y*1.d-4,Z*1.d-4,NODE(I)
	 WRITE(100,100)I,IX,IY,IZ,1,IT,X,Y*10.d0,Z,NODE(I)
	  
 30   CONTINUE





C 100  FORMAT(6I5,3F10.7,I5)
 100  FORMAT(6I5,3F10.5,I5)
c  100  FORMAT(i7,5I5,3(f10.7),I5)

       CLOSE(100)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE BOUNDT3(CORD,NEL,NDIM,NPT,NET,NODE,NZAD,ZADVRE,NUMZAD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*),NZAD(3,*)
	DIMENSION ZADVRE(*)

C
CE Subroutine BOUNDT3 
C
      OPEN(100,FILE='OGRANIC.DAT')
      
      DO I=1,NPT
        NODE(I)=0
      ENDDO

      DO J=1,NUMZAD
	  IF (NZAD(2,J).EQ.5) THEN 
	   N=NZAD(1,J)
	   NODE(N)=1
	   ENDIF
	ENDDO
C          WRITE(100,100)I,ix,iy,iz,1,1,(CORD(J,I)*0.1d0,J=1,3),NODE(I)
      DO I=1,NPT
       WRITE(100,100)I,0,0,0,1,NODE(I),(CORD(J,I),J=1,3)
	ENDDO

 
 
      
  100  FORMAT(6I5,3F10.7,I5)

       CLOSE(100)

	STOP

      END
C==========================================================================
C==========================================================================
      SUBROUTINE nikola2(CORD,NEL,NDIM,NPT,NET,INDEL,CCORD,KKORAK,GNODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*),CCORD(3,*)
      DIMENSION NEL(NDIM+1,*),INDEL(*)
      DIMENSION GNODE(2,6,*)

      do i=1,npt
 	 if (ccord(3,i).lt.0.d0) ccord(3,i)=cord(3,i)/2.0d0


c 	 ccord(1,i)=cord(1,i)*1.5d0
c	 ccord(2,i)=cord(2,i)*1.5d0
c 	 ccord(1,i)=cord(1,i)*1.0d0

c	 ccord(2,i)=cord(2,i)*0.2250d0
	 
	 
      if (kkorak.le.20) then	 
	 ccord(2,i)=cord(2,i)*(1.50d0-1.275d0*kkorak/20.d0)
        if (ccord(3,i).gt.0.d0) then
	   v=(-100.d0*(15.d0-cord(3,i))/15.d0)*(kkorak/20.d0)
 	   GNODE(1,3,i)=v
 	   GNODE(2,3,i)=v
	  endif
	else 
 	   GNODE(1,3,i)=0.d0
 	   GNODE(2,3,i)=0.d0
	 ccord(2,i)=cord(2,i)*(1.50d0-1.275d0)
	endif


      enddo


c      if (kkorak.gt.20) then
c	
c	return
c	endif

      return



      N1=637
C	N2=1594
	N2=1577
	NE1=476




	step=20.d0



c      DO I=1,N1
c        INDEL(I)=0
c      ENDDO

      

c      DO I=1,NE1
c       DO J=1,4
c        INDEL(NEL(J,I))=INDEL(NEL(J,I))+1
c	 ENDDO
c      ENDDO

C
      OPEN(100,FILE='DISPL.DAT')
      DO I=1,N1
      READ(100,*) N,X,Y,iv
	x=x/3.d0
	y=y/3.d0


      do k=0,5
      J2=I+N1*k
	CCORD(1,J2)=CORD(1,J2)+KKORAK*X/step
	CCORD(2,J2)=CORD(2,J2)+KKORAK*Y/step
      IF (iv.eq.1) THEN
	  GNODE(2,1,J2)=X/step
	  GNODE(2,2,J2)=Y/step
	ENDIF
	enddo

      do k=1,5
      J3=I+N1*5+N2*k
	CCORD(1,J3)=CORD(1,J3)+KKORAK*X/step
	CCORD(2,J3)=CORD(2,J3)+KKORAK*Y/step
c      IF (INDEL(I).NE.4) THEN
c	  GNODE(2,1,J3)=X/step
c	  GNODE(2,2,J3)=Y/step
c	ENDIF
	enddo


	ENDDO
      

C
      CLOSE (100)



      END


C==============================================================
C==============================================================
C==============================================================
C==============================================================
      SUBROUTINE nikola(CORD,NEL,NDIM,NPT,NET,NODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*)

C=============================================================================================================
C=============================================================================================================
C      OPEN(100,FILE='obod.txt')
C      OPEN(101,FILE='bound.txt')
C      NN1=637
c	NE1=476 
c      do i=1,306
c	 read(100,'(3i5)') iel,n1,n2
c	do j=1,5
c       write(101,'(6i5)') iel+(j-1)*NE1,
c     &n1+(j-1)*NN1,n2+(j-1)*NN1,n2+j*NN1,n1+j*NN1,1
c	enddo
c
c	enddo
c
c     stop
C=============================================================================================================
C=============================================================================================================



      OPEN(100,FILE='3D.DAT')
c      GOTO 1100

      DO I=1,NPT
        NODE(I)=0
      ENDDO

      

      DO I=1,NET
       DO J=1,NDIM
        NODE(NEL(J,I))=NODE(NEL(J,I))+1
	 ENDDO
      ENDDO


      DO I=1,5*476


	N1=NEL(1,I)
	N2=NEL(5,I)
	N3=NEL(6,I)
	N4=NEL(2,I)
       IF (NODE(N1).LT.8.AND.NODE(N2).LT.8.AND.
     &NODE(N3).LT.8.AND.NODE(N4).LT.8) THEN
	  WRITE(100,'(6I5)') I,N1,N2,N3,N4
      ENDIF	 


	N1=NEL(1,I)
	N2=NEL(4,I)
	N3=NEL(8,I)
	N4=NEL(5,I)
       IF (NODE(N1).LT.8.AND.NODE(N2).LT.8.AND.
     &NODE(N3).LT.8.AND.NODE(N4).LT.8) THEN
	  WRITE(100,'(6I5)') I,N1,N2,N3,N4
      ENDIF	 

	N1=NEL(4,I)
	N2=NEL(8,I)
	N3=NEL(7,I)
	N4=NEL(3,I)
       IF (NODE(N1).LT.8.AND.NODE(N2).LT.8.AND.
     &NODE(N3).LT.8.AND.NODE(N4).LT.8) THEN
	  WRITE(100,'(6I5)') I,N1,N2,N3,N4
      ENDIF	 

	N1=NEL(3,I)
	N2=NEL(7,I)
	N3=NEL(6,I)
	N4=NEL(2,I)
       IF (NODE(N1).LT.8.AND.NODE(N2).LT.8.AND.
     &NODE(N3).LT.8.AND.NODE(N4).LT.8) THEN
	  WRITE(100,'(6I5)') I,N1,N2,N3,N4
      ENDIF	 

	ENDDO



       CLOSE(100)

       RETURN 

1100   CONTINUE

	N1=637
	NE1=476
      

      DO I=1,NPT
	  cord(1,i)=cord(1,i)/3.d0
	  cord(2,i)=cord(2,i)/3.d0
	  cord(3,i)=cord(3,i)/3.d0
	ENDDO

      DO I=1,N1
	 WRITE (100,'(6I5,3F10.4)') I,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)+10.d0
	ENDDO

      DO I=1,N1
	 WRITE (100,'(6I5,3F10.4)') I+n1,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)+8.d0
	ENDDO


      DO I=1,N1
	 WRITE (100,'(6I5,3F10.4)') I+2*n1,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)+6.d0
	ENDDO

      DO I=1,N1
	 WRITE (100,'(6I5,3F10.4)') I+3*n1,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)+4.d0
	ENDDO

      DO I=1,N1
	 WRITE (100,'(6I5,3F10.4)') I+4*n1,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)+2.d0
	ENDDO

      DO I=1,NPT
	 WRITE (100,'(6I5,3F10.4)') I+5*n1,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)+0.d0
	ENDDO

      DO I=1,NPT
	WRITE (100,'(6I5,3F10.4)') I+5*n1+NPT,0,0,0,0,0,(CORD(J,I),J=1,2),
     &cord(3,i)-2.d0
	ENDDO

      DO I=1,NPT
	 WRITE (100,'(6I5,3F10.4)') I+5*n1+2*NPT,0,0,0,0,0,
     &(CORD(J,I),J=1,2),cord(3,i)-4.d0
	ENDDO

      DO I=1,NPT
	 WRITE (100,'(6I5,3F10.4)') I+5*n1+3*NPT,0,0,0,0,0,
     &(CORD(J,I),J=1,2),cord(3,i)-6.d0
	ENDDO

      DO I=1,NPT
	 WRITE (100,'(6I5,3F10.4)') I+5*n1+4*NPT,0,0,0,0,0,
     &(CORD(J,I),J=1,2),cord(3,i)-8.d0
	ENDDO


      DO I=1,NPT
	 WRITE (100,'(6I5,3F10.4)') I+5*n1+5*NPT,0,0,0,0,0,
     &(CORD(J,I),J=1,2),cord(3,i)-10.d0
	ENDDO
   
      WRITE(100,*)'ELEMENTI'

      DO I=1,NE1
	 WRITE (100,'(10I5)') I,(NEL(J,I),J=1,4),(NEL(J,I)+N1,J=1,4),1
	ENDDO

      DO I=1,NE1
	 WRITE (100,'(10I5)') I+NE1,(NEL(J,I)+N1,J=1,4),
     &(NEL(J,I)+2*N1,J=1,4),1
	ENDDO


      DO I=1,NE1
	 WRITE (100,'(10I5)') I+2*NE1,(NEL(J,I)+2*N1,J=1,4),
     &(NEL(J,I)+3*N1,J=1,4),1
	ENDDO

      DO I=1,NE1
	 WRITE (100,'(10I5)') I+3*NE1,(NEL(J,I)+3*N1,J=1,4),
     &(NEL(J,I)+4*N1,J=1,4),1
	ENDDO

      DO I=1,NE1
	 WRITE (100,'(10I5)') I+4*NE1,(NEL(J,I)+4*N1,J=1,4),
     &(NEL(J,I)+5*N1,J=1,4),1
	ENDDO


      DO I=1,NET
	 WRITE (100,'(10I5)') I+5*NE1,(NEL(J,I)+5*N1,J=1,4),
     &(NEL(J,I)+5*N1+NPT,J=1,4),1
	ENDDO

      DO I=1,NET
	 WRITE (100,'(10I5)') I+5*NE1+NET,(NEL(J,I)+5*N1+NPT,J=1,4),
     &(NEL(J,I)+5*N1+2*NPT,J=1,4),1
	ENDDO

      DO I=1,NET
	 WRITE (100,'(10I5)') I+5*NE1+2*NET,(NEL(J,I)+5*N1+2*NPT,J=1,4),
     &(NEL(J,I)+5*N1+3*NPT,J=1,4),1
	ENDDO

      DO I=1,NET
	 WRITE (100,'(10I5)') I+5*NE1+3*NET,(NEL(J,I)+5*N1+3*NPT,J=1,4),
     &(NEL(J,I)+5*N1+4*NPT,J=1,4),1
	ENDDO

      DO I=1,NET
	 WRITE (100,'(10I5)') I+5*NE1+4*NET,(NEL(J,I)+5*N1+4*NPT,J=1,4),
     &(NEL(J,I)+5*N1+5*NPT,J=1,4),1
	ENDDO


       CLOSE(100)

       RETURN 

      DO I=1,NPT
        NODE(I)=0
      ENDDO

      

      DO 20 I=1,NET
       DO 10 J=1,NDIM
        NODE(NEL(J,I))=NODE(NEL(J,I))+1
   10  CONTINUE 
   20  CONTINUE 



      DO 30 I=1,NPT
	  IX=0
	  IY=0
	  IZ=0
	  IF (DABS(CORD(3,I)).LT.1.D-5) THEN
	     IF (NODE(I).EQ.4) THEN
              IX=0
	        IY=0
	        IZ=1
	     ELSE
              IX=1
	        IY=1
              IZ=1
    		 ENDIF
	   
	   ELSE
	     IF (NODE(I).LT.8) THEN
	       IX=1
	       IY=1
	       IZ=1
	     ENDIF
	   ENDIF

C	  X=CORD(1,I)*0.425/100.D0
C	  Y=CORD(2,I)*0.425/100.D0
C	  Z=CORD(3,I)*0.425/100.D0

	  X=CORD(1,I)
	  Y=CORD(2,I)
	  Z=CORD(3,I)

C         y0=-0.4d0

C       if (i.le.66*8) then
C	  y1=(cord(2,(mod((i-1),66)+1+66*8))-y0)/8.d0
C	  Y=y0+((i-1)/66)*y1
C	 endif

	  dod=1.D0
C	  isloj=190
C	  if (i.gt.17*91.and.i.le.25*91) dod=0.965D0
C	  if (i.gt.25*91.and.i.le.(25+8)*91) dod=0.918D0
C	  if (i.gt.(25+8)*91.and.i.le.(25+24)*91) dod=0.965D0
C	  if (i.gt.(25+24)*91.and.i.le.(74)*91) dod=0.918D0

C	  if (i.gt.21*isloj.and.i.le.27*isloj) dod=0.965D0
C	  if (i.gt.27*isloj.and.i.le.(27+10)*isloj) dod=0.918D0
C	  if (i.gt.(27+10)*isloj.and.i.le.(27+26)*isloj) dod=0.965D0
C	  if (i.gt.(27+26)*isloj.and.i.le.(80)*isloj) dod=0.918D0

C	  if (i.gt.(80)*isloj.and.i.le.(83)*isloj) dod=0.918D0

C	  X=X*dod
C	  Z=Z*dod
       
c	 WRITE(100,100)I,IX,IY,IZ,1,1,X*1.d-4,Y*1.d-4,Z*1.d-4,NODE(I)
	 WRITE(100,100)I,IX,IY,IZ,1,1,X,Y,Z,NODE(I)
	  
 30   CONTINUE





C 100  FORMAT(6I5,3F10.7,I5)
 100  FORMAT(6I5,3F10.7,I5)
c  100  FORMAT(i7,5I5,3(f10.7),I5)


      END
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE BOUND4(CORD,NEL,NDIM,NPT,NET,ID,NODE,NZAD,ZADVRE,
     &NUMZAD)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(2,*),ID(6,*)
      DIMENSION NZAD(3,*),ZADVRE(*)
	  DIMENSION IID(5)

C
CE Subroutine BOUND4 is used for definition of symmetry plane
C
      
      DO I=1,NPT
        NODE(1,I)=0
        NODE(2,I)=0
      ENDDO

      OPEN(100,FILE='OGRANIC.DAT')

      NNODE=0
      NELEM=0


      DO 10 I=1,NPT
      IF (DABS(CORD(3,I)).LT.1.D-5.OR.CORD(3,I).GT.0.D0) THEN
       NODE(1,I)=1
	ENDIF
 10   CONTINUE

      DO 20 I=1,NET
       NN=0
       DO 30 J=1,NDIM
       IF (NODE(1,NEL(J,I)).EQ.1) NN=NN+1
30     CONTINUE 
        IF (NN.LT.8.AND.NN.GT.4) THEN
         DO J=1,NDIM
          IF (NODE(1,NEL(J,I)).EQ.0) THEN
            NODE(1,NEL(J,I))=2
            CORD(3,NEL(J,I))=0.D0
	   CORD(1,NEL(J,I))=0.5D0*(CORD(1,NEL(J-1,I))+CORD(1,NEL(J+1,I)))
	   CORD(2,NEL(J,I))=0.5D0*(CORD(2,NEL(J-1,I))+CORD(2,NEL(J+1,I)))
	    ID(1,NEL(J,I))=1
	    ID(2,NEL(J,I))=1
          ENDIF
         ENDDO         
C	      CORD(3,NEL(2,I))=0.5D0*(CORD(3,NEL(1,I))+CORD(3,NEL(3,I)))
C	      CORD(3,NEL(6,I))=0.5D0*(CORD(3,NEL(5,I))+CORD(3,NEL(7,I)))
       ENDIF
   20  CONTINUE 


      DO 35 I=1,NPT
        IF (NODE(1,I).NE.0) THEN
         NNODE=NNODE+1
		 NODE(2,I)=NNODE
		 DO K=1,5
		  IF (ID(K,I).GT.0) IID(K)=0
		  IF (ID(K,I).EQ.0) IID(K)=1
		 ENDDO
	IF (DABS(CORD(3,I)).LT.1.D-4) IID(3)=1

          WRITE(100,100)NODE(2,I),(IID(J),J=1,5),(CORD(J,I),J=1,3)
C     &	,NODE(1,I)
        ENDIF

 35   CONTINUE
      
      WRITE(100,*)'C NODES',NNODE

      DO 40 I=1,NET
       NN=0
       DO J=1,NDIM
        IF (NODE(1,NEL(J,I)).NE.0) THEN
		 NN=NN+1
		ENDIF
        ENDDO
       IF (NN.EQ.8) THEN
	       NELEM=NELEM+1
          WRITE(100,200)NELEM,(NODE(2,NEL(J,I)),J=1,8),1
	   ENDIF
 40   CONTINUE

      WRITE(100,*)'C ELEMENTS',NELEM
      NVELOC=0

      DO 50 I=1,NUMZAD
       IF (CORD(3,NZAD(1,I)).LT.0.D0) GOTO 50
       WRITE(100,300) NODE(2,NZAD(1,I)),(NZAD(J,I),J=2,3),ZADVRE(I)
       NVELOC=NVELOC+1
 50   CONTINUE

      WRITE(100,*)'C PRESCRIBED VELOCITY ',NVELOC


  100  FORMAT(6I5,3F10.7,I5)
  200  FORMAT(10I5)
  300  FORMAT(3I5,F10.7)

       CLOSE(100)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE BOUN3D(CORD,NEL,NDIM,NPT,NET,NODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*)

C
CE Subroutine BOUN3D is used for definition of boundary constraints
CE for 3D walls
C
      
      
      

      DO I=1,NPT
        NODE(I)=0
      ENDDO

      DO 20 I=1,NET
       DO 10 J=1,NDIM
        NODE(NEL(J,I))=NODE(NEL(J,I))+1
   10  CONTINUE 
   20  CONTINUE 


      OPEN(100,FILE='OGRANIC.DAT')

      DO 30 I=1,NPT
	  IX=0
	  IY=0
	  IZ=0
	     IF (NODE(I).LT.8) THEN
	       IX=1
	       IY=1
	       IZ=1
	     ENDIF
	  X=CORD(1,I)
	  Y=CORD(2,I)
	  Z=CORD(3,I)
       

C	if (x.le.-0.08d0) x=-0.8d0
	  
c	 WRITE(100,100)I,IX,IY,IZ,1,0,X,Y,Z,NODE(I)
	 WRITE(100,100)I,IX,IY,IZ,1,0,3.d0*X,3.d0*Y,3.d0*Z,NODE(I)
	  
 30   CONTINUE





C 100  FORMAT(6I5,3F10.7,I5)
 100  FORMAT(6I5,3F10.7,I5)
c  100  FORMAT(i7,5I5,3(f10.7),I5)

       CLOSE(100)

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE BOUN2D(CORD,NEL,NDIM,NPT,NET,NODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),NODE(*)

C
CE Subroutine BOUN3D is used for definition of boundary constraints
CE for 3D walls
C
      
      DO I=1,NPT
        NODE(I)=0
      ENDDO

      DO 20 I=1,NET
       DO 10 J=1,NDIM
        NODE(NEL(J,I))=NODE(NEL(J,I))+1
   10  CONTINUE 
   20  CONTINUE 


      OPEN(100,FILE='OGRANIC.DAT')

      DO 30 I=1,NPT
	  IX=0
	  IY=0
	  IZ=0
	     IF (NODE(I).LT.4) THEN
	       IX=1
	       IY=1
	       IZ=1
	     ENDIF
	  X=CORD(1,I)
	  Y=CORD(2,I)
	  Z=CORD(3,I)

	  
	 WRITE(100,100)I,IX,IY,1,1,1,X,Y,Z,NODE(I)
	  
 30   CONTINUE





C 100  FORMAT(6I5,3F10.7,I5)
 100  FORMAT(6I5,3F10.5,I5)
c  100  FORMAT(i7,5I5,3(f10.7),I5)

       CLOSE(100)

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE CEVI(CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)

      FACT0=0.1D0
      B=1.40D0
      B1=1.25D0


      FACT=FACT0

       DO I=13323,13697
         CORD(1,I)=CORD(1,I)-0.1951095D0*FACT
         CORD(2,I)=CORD(2,I)+0.4794649D0*FACT
         IF (MOD((I-13323+1),25).EQ.0) FACT=FACT*B
       ENDDO

      FACT=FACT0
         
       DO I=13698,14432
         CORD(1,I)=CORD(1,I)-0.0273843D0*FACT*1.4D0
         CORD(2,I)=CORD(2,I)+0.2642122D0*FACT*1.4D0
         IF (MOD((I-13698+1),49).EQ.0) FACT=FACT*B
       ENDDO

      FACT=FACT0

       DO I=14433,14807
         CORD(1,I)=CORD(1,I)+0.1840271D0*FACT
         CORD(2,I)=CORD(2,I)+0.4838757D0*FACT
         IF (MOD((I-14433+1),25).EQ.0) FACT=FACT*B
       ENDDO
       
	  FACT=FACT0


       DO I=14808,15307
         CORD(1,I)=CORD(1,I)+0.2518290D0*FACT
         CORD(2,I)=CORD(2,I)-0.1194137D0*FACT
         CORD(3,I)=CORD(3,I)-0.4359220D0*FACT
         IF (MOD((I-14808+1),25).EQ.0) FACT=FACT*B1
       ENDDO
      	FACT=FACT0
	
      DO I=15308,15807
         CORD(1,I)=CORD(1,I)+0.2510185D0*FACT
         CORD(2,I)=CORD(2,I)-0.1206178D0*FACT
         CORD(3,I)=CORD(3,I)+0.4364946D0*FACT
         IF (MOD((I-15308+1),25).EQ.0) FACT=FACT*B1
       ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE SHELLS(CORD,NEL,NDIM,NPT,NET,ID,NZAD,ZADVRE,NUMZAD)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NZAD(3,*)
      DIMENSION ZADVRE(*)

      DIMENSION NID(8),N(8),NWALL(6),ITR(6,4)


      IZ=49
      OPEN (IZ,FILE='SHELL.DAT')

      NNODES=0
      WRITE(IZ,*)'C NODES'

      DO 10 I=1,NPT
       IF (ID(1,I).NE.0.OR.ID(2,I).NE.0.OR.ID(3,I).NE.0) GOTO 10
       DO K=1,NUMZAD
        IF (NZAD(1,K).EQ.I) GOTO 10
       ENDDO

	  IIX=0
	  IIY=0
	  IIZ=0
	  MX=0
	  MY=0
	  MZ=0
	  IDOD=0
	IF (DABS(CORD(3,I)).LT.1.E-4) THEN
	  IIZ=1
	  MX=1
	  MY=1
	  IDOD=1
	ENDIF
       WRITE(IZ,100)I,IIX,IIY,IIZ,MX,MY,MZ,(CORD(J,I),J=1,3),IDOD
       NNODES=NNODES+1
       
       
 10   CONTINUE


 100  FORMAT(I5,1X,6I2,2X,3F10.7,I10)

      WRITE(IZ,*)'C TOTAL NODES ',NNODES



      WRITE(IZ,*)'C ELEMENTS'

      IEL=0
      THICK=5.D-2

      DO 20 NBREL=1,NET      

C	DO I=1,6
C	 NWALL(I)=0
C	ENDDO

	DO 15 I=1,NDIM
	  N(I)=NEL(I,NBREL)
        NID(I)=0
      IF(ID(1,N(I)).EQ.0.AND.ID(2,N(I)).EQ.0.AND.ID(3,N(I)).EQ.0) THEN
	 DO K=1,NUMZAD
        IF (NZAD(1,K).EQ.N(I)) GOTO 15
       ENDDO


       NID(I)=1
      ENDIF
  15  CONTINUE     




      ITR(1,1)=N(8)
      ITR(1,2)=N(4)
      ITR(1,3)=N(5)
      ITR(1,4)=N(1)

      ITR(2,1)=N(7)
      ITR(2,2)=N(3)
      ITR(2,3)=N(6)
      ITR(2,4)=N(2)

      ITR(3,1)=N(6)
      ITR(3,2)=N(2)
      ITR(3,3)=N(5)
      ITR(3,4)=N(1)

      ITR(4,1)=N(7)
      ITR(4,2)=N(3)
      ITR(4,3)=N(8)
      ITR(4,4)=N(4)

      ITR(5,1)=N(3)
      ITR(5,2)=N(2)
      ITR(5,3)=N(4)
      ITR(5,4)=N(1)

      ITR(6,1)=N(7)
      ITR(6,2)=N(6)
      ITR(6,3)=N(8)
      ITR(6,4)=N(5)


      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 30 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 30

       IEL=IEL+1
        WRITE(IZ,200) IEL,1,0,0,0,THICK,0,0.D0,0.D0
        WRITE(IZ,300) ITR(NPOV,1),ITR(NPOV,2),ITR(NPOV,4),ITR(NPOV,3)

 30   CONTINUE

 20   CONTINUE

 200  FORMAT(5I5,F10.4,I5,2F10.4)
 300  FORMAT(4I5)

      WRITE(IZ,*)'C TOTAL ELEMENTS ',IEL

      CLOSE(IZ)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE MAE(CORD,NPT,NSTEP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)



	do ii=0,5
      DO 10 I=1+nstep*ii,NPT,NSTEP*nstep
C	if (mod (i,(nstep*nstep-(nstep-1))).eq.0) goto 10

	X0=CORD(1,I)
	Y0=CORD(2,I)
	X1=CORD(1,I+NSTEP-1)
	Y1=CORD(2,I+NSTEP-1)
	  DP=1.D0/(3.D0*(NSTEP-1))
	  P=DP
        DO J=I+1,I+3
	    X=X0*(1-P)+X1*P
	    Y=Y0*(1-P)+Y1*P
	    CORD(1,J)=X
	    CORD(2,J)=Y
	    
	    X=X1*(1-P)+X0*P
	    Y=Y1*(1-P)+Y0*P
	    CORD(1,I+NSTEP-1-(J-I))=X
	    CORD(2,I+NSTEP-1-(J-I))=Y
	    P=P+DP
	  ENDDO	 
10    continue
       enddo

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE MAE1(CORD,NPT,NSTEP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)



	
      DO 10 I=1,NPT,NSTEP*nstep
       do ii=0,nstep*4-1
	   cord(3,i+ii)=cord(3,i+ii)/3.d0
       enddo 
	
10    continue
     

      END
C==========================================================================
C==========================================================================
      SUBROUTINE MAE2(CORD,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)

C      DO I=469,858
C       cord(1,i)=cord(1,i)*3.5d0
C       cord(2,i)=cord(2,i)*3.5d0
C      enddo

C      DO I=1,464

C      DO I=1,468
C       cord(1,i)=cord(1,i)-(10.D-4)*((100.D-4)-(cord(2,i)))/(100.D-4)
C      enddo

C      DO I=1239,NPT
	
C      DO I=859,NPT
C       cord(1,i)=cord(1,i)+(20.D-4)*((100.D-4)-(cord(2,i)))/(100.D-4)
C      enddo


      DO I=1,NPT
       cord(1,i)=cord(1,i)*1.d-4
       cord(2,i)=cord(2,i)*1.d-4
       cord(3,i)=cord(3,i)*1.d-4
      enddo
     

      END
C==========================================================================
C==========================================================================
      SUBROUTINE BOUND5(CORD,NPT,ID,NEL,NDIM,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION ID(6,*)
      DIMENSION NEL(NDIM+1,*)

C
CE Subroutine BOUND5 is used for definition of boundary constraints
C=========================================================================
C BOUNDARY CONDITIONS FOR IEL PORES IN MEDIAL MASS TRANSPORT EXAMPLE
C=========================================================================
      OPEN(45,FILE='INPUT.DAT')
      TOL=1.D-4


      DO II=1,NPT
C	  ID(5,II)=0 
C        IF (DABS(CORD(2,II)).LT.TOL) ID(2,II)=1
C        IF (DABS(CORD(2,II)-10.0D0).LT.TOL) ID(2,II)=1
c	  WRITE(45,1005) II,(ID(J,II),J=1,5),(CORD(J,II),J=1,3)
          ID1=0
	    ID2=0
	    IDC=0
         IF (ID(1,II).EQ.0) ID1=1
         IF (ID(2,II).EQ.0) ID2=1
         IF (ID(5,II).EQ.0) IDC=1
C	   IF (DABS(CORD(1,II)-(-2.5D-4)).LT.TOL) THEN
	   IF (DABS(CORD(1,II)-(5.062D-03)).LT.TOL) THEN
	    ID1=1
	    ID2=1
	    IDC=1
	   ENDIF
c	  WRITE(45,1015) II,(ID(J,II),J=1,5),(CORD(J,II)*1.D-4,J=1,3)
	  WRITE(45,1015) II,ID1,ID2,1,1,IDC,(CORD(J,II),J=1,3)
      ENDDO


      WRITE(45,*)'C NNODE, INDPR, ITIMF, VALUE'
      DO II=1,NPT
	  IF (DABS(CORD(1,II)-(5.062D-03)).LT.TOL) THEN
C	   WRITE(45,2015) II,1,1,1.0
	   WRITE(45,2015) II,5,2,1.0
	  ENDIF
      ENDDO

      
	goto 100


	DO I=1,1200
	  WRITE(45,1000)I,(NEL(J,I),J=1,4),1
	ENDDO
        
      DO K=1,18
	 DO I=1201+K*1200,2400+K*1200
	  DO J=1,4
         NEL(J,I)=NEL(J,I-1200)+1279         	  
	  ENDDO
	  WRITE(45,1000)I,(NEL(J,I),J=1,4),1
 	 ENDDO
	ENDDO


   
100   CLOSE (45)
	STOP
C=========================================================================

CE for 3D walls
C

 1005 FORMAT(I5,5(3X,I2),3F10.6,2I2)
 1015 FORMAT(I5,5(3X,I2),3(1PE10.3),2I2)
 2015 FORMAT(3I5,1PE10.3)
 1000 FORMAT(6(I5))

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE BOUND55(CORD,NPT,ID,NEL,NDIM,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION ID(6,*)
      DIMENSION NEL(NDIM+1,*)

	DIMENSION NODE(4)

C
CE Subroutine BOUND5 is used for definition of boundary constraints
C=========================================================================
C BOUNDARY CONDITIONS FOR IEL PORES IN MEDIAL MASS TRANSPORT EXAMPLE
C=========================================================================
      OPEN(45,FILE='INPUT.DAT')
      TOL=1.D-8

	DO 300 M=1,3
      DO 200 K=1,11
      DO 100 NE=1+(K-1)*1200,1200+(K-1)*1200
	 IBR=0
	 DO J=1,4
	  N=NEL(J,NE)
        X=CORD(1,N)	  
        Y=CORD(2,N)
	  P=0.D0+(K-1)*5.D-4
	  Q=5.D-4*M
	  R=2.D-4
	  IF(DABS((X-P)**2+(Y-Q)**2-R**2).LT.TOL) THEN
	    IBR=IBR+1
	    NODE(IBR)=N
	  ENDIF
	 ENDDO
	  IF (IBR.EQ.2) THEN
	    WRITE(45,1000) NE,NODE(1),NODE(2),3	    
	  ENDIF
100   CONTINUE
200   CONTINUE
300   CONTINUE

   
      CLOSE (45)
	STOP
C=========================================================================

 1000 FORMAT(6(I5))

      END
C==========================================================================
C==========================================================================
