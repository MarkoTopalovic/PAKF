C=========================================================================
C=========================================================================
C      SUBROUTINE EXIM3D
C=======================================================================
C=========================================================================
      SUBROUTINE EXIM3D(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &AMASA,VMESH,CCORD,VELOC,IDALE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine EXIM3D is used for 3D analysis
CE It is used global loop per elements
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*)
      DIMENSION AMASA(*),VMESH(3,*),CCORD(3,*),VELOC(3,*)
      DIMENSION IDALE(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21)
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
      DIMENSION C(21,21)
      DIMENSION RBP(8),POISS(8,8),TT21A(2*21)

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





      ISTEP=2

 100  CALL CLEAR(ALEVO,NWK)
      CALL CLEAR(DESNO,NWK)
      CALL CLEAR(SILE,JEDN)
      CALL CLEAR(VELOC,3*NPT)


C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET

      VOLUM=0.D0

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


      DO I=1,NDIM
       NODE=NEL(I,NBREL)
       TT21A(I)=VMESH(1,NODE)
       TT21A(I+NDIM)=VMESH(2,NODE)
       TT21A(I+2*NDIM)=VMESH(3,NODE)
      ENDDO


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
       POISS(K,N)=0.D0
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
      RBP(I)=0.D0
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

      CALL ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
 
      WDT=WT*DET1
      VOLUM=VOLUM+WDT
  
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
       AKV1P(K,N)=AKV1P(K,N)+WDT*(PJ(1,K)*HP(N))
       AKV2P(K,N)=AKV2P(K,N)+WDT*(PJ(2,K)*HP(N))
       AKV3P(K,N)=AKV3P(K,N)+WDT*(PJ(3,K)*HP(N))
      ENDIF
      IF (N.LE.8.AND.K.LE.8) THEN
        POISS(K,N)=POISS(K,N)+
     &(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*WDT
      ENDIF
 164  CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
       RBP(K)=RBP(K)+(PJ(1,K)*FB2+PJ(2,K)*FB3+PJ(3,K)*0.D0)*GUSM*WDT
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
C PRVI I CETVRTI KORAK:
C First and Fourth Step

	IF (ISTEP.EQ.1) THEN


       DO I=1,NDES
        LM2(I)=0
       ENDDO

       
      DO I=1,8
      LM2(I+3*NDIM)=ID(4,NEL(I,NBREL))
         K=I+3*NDIM
        F92(K)=F92(K)+RS2(I)
       DO J=1,NDIM
        F92(K)=F92(K)-(GUSM/TIME)*AKV1P(J,I)*TT210(J)
        F92(K)=F92(K)-(GUSM/TIME)*AKV2P(J,I)*TT210(J+NDIM)
        F92(K)=F92(K)-(GUSM/TIME)*AKV3P(J,I)*TT210(J+2*NDIM)
       ENDDO        
        F92(K)=F92(K)+RBP(I)
      ENDDO        

      DO I=1,8
       DO J=1,8
        SKEF(I+3*NDIM,J+3*NDIM)=POISS(I,J)
       ENDDO        
      ENDDO     

C      CALL WRRF(SKEF,NDES*NDES,'SKEF1',IIZLAZ)      
C      CALL WRRF(F92,NDES,'F92  ',IIZLAZ)      
        
      ELSE IF (ISTEP.EQ.4) THEN

       DO I=1,NDES
        LM2(I)=0
       ENDDO

       
      DO I=1,8
      LM2(I+3*NDIM)=ID(4,NEL(I,NBREL))
         K=I+3*NDIM
        F92(K)=F92(K)+RS2(I)
        DO J=1,NDIM
         F92(K)=F92(K)-(GUSM/TIME)*AKV1P(J,I)*TT210(J)
         F92(K)=F92(K)-(GUSM/TIME)*AKV2P(J,I)*TT210(J+NDIM)
         F92(K)=F92(K)-(GUSM/TIME)*AKV3P(J,I)*TT210(J+2*NDIM)
        ENDDO        
        F92(K)=F92(K)+RBP(I)
      ENDDO        

      DO I=1,8
       DO J=1,8
        SKEF(I+3*NDIM,J+3*NDIM)=POISS(I,J)
       ENDDO        
      ENDDO     

      ELSE IF (ISTEP.EQ.2) THEN
C=========================================================================
C DRUGI  KORAK
C The Second Step 
      
C       DO I=1,NDES
C        LM2(I)=0
C       ENDDO


      DO 270 I=1,NDIM
C      LM2(I)=ID(1,NEL(I,NBREL))
C      LM2(I+NDIM)=ID(2,NEL(I,NBREL))
      DO 265 J=1,NDIM
      SKEF(I,J)=-AKMIV1(I,J)
      SKEF(I+NDIM,J+NDIM)=-AKMIV1(I,J)
      SKEF(I+2*NDIM,J+2*NDIM)=-AKMIV1(I,J)
      IF (J.LE.8) THEN
       SKEF(I,J+3*NDIM)=AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AKV3P(I,J)
      ENDIF
 265  CONTINUE
 270  CONTINUE


        DELI=1.D0

      DO  I=1,NDIM
        F92(I)=F92(I)+RB2(I)*TIME/GUSM
        F92(I+NDIM)=F92(I+NDIM)+RB3(I)*TIME/GUSM
      ENDDO

      DO  I=1,3*NDIM
       DO  J=1,3*NDIM+8
        F92(I)=F92(I)+SKEF(I,J)*TT21(J)*(TIME/GUSM)
       ENDDO
      ENDDO




      ELSE IF (ISTEP.EQ.5) THEN
C=========================================================================
C   PETI KORAK
C The Fifth Step
      

      DO  I=1,NDIM
      DO  J=1,NDIM
       SKEF(I,J)=-AKMIV1(I,J)
       SKEF(I+NDIM,J+NDIM)=-AKMIV1(I,J)
       SKEF(I+2*NDIM,J+2*NDIM)=-AKMIV1(I,J)
      IF (J.LE.8) THEN
       SKEF(I,J+3*NDIM)=AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AKV3P(I,J)
      ENDIF
      ENDDO
      ENDDO

        DELI=1.D0

      DO  I=1,NDIM
        F92(I)=F92(I)+RB2(I)*TIME/GUSM
        F92(I+NDIM)=F92(I+NDIM)+RB3(I)*TIME/GUSM
      ENDDO

      DO  I=1,3*NDIM
       DO  J=1,3*NDIM+8
        F92(I)=F92(I)+SKEF(I,J)*TT21(J)*(TIME/GUSM)
       ENDDO
      ENDDO

C=========================================================================
C FAZA 2 KOD ALE FORMULACIJE :

C TRECI KORAK:
C The third step

	 ELSE IF(ISTEP.EQ.3) THEN

      DO  I=1,NDIM
       DO  J=1,NDIM
      F92(I)=F92(I)-(AKVV1(I,J)*TIME/GUSM)*TT210(J)
      F92(I+NDIM)=F92(I+NDIM)-(AKVV1(I,J)*TIME/GUSM)*TT210(J+NDIM)
      F92(I+2*NDIM)=F92(I+2*NDIM)-(AKVV1(I,J)*TIME/GUSM)*TT210(J+2*NDIM)
       ENDDO
      ENDDO


      ENDIF

C==========================================================================
C==========================================================================

C      CALL SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
C     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)

C       IF (ISTEP.EQ.1) CALL MASAM(AMASA,VOLUM,NEL,NDIM,NBREL,3)
       IF (ISTEP.EQ.2) CALL MASAM(AMASA,VOLUM,NEL,NDIM,NBREL,3)
      
      INDSK=0
      IF(ISTEP.EQ.1.OR.ISTEP.EQ.4) THEN
       INDSK=1
       CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F92,MAXA,LM2,NDES,INDSK)
      ELSE
       CALL ADDST1(VELOC,F92,NEL,NBREL,NDIM,IDALE)
      ENDIF



C RACUNANJE SILA INTERAKCIJE KOJIMA FLUID DELUJE NA ZIDOVE
c calculation of interaction forces between fluid and solid
C	 IF (NSTAC.EQ.0) THEN
C       DO K=0,2
C        DO I=1,NDIM
C          II=I+K*NDIM
C         DO J=1,NDIM
C          JJ=J+K*NDIM
C          SKEF(II,JJ)=SKEF(II,JJ)-AMV2(I,J)*TT210(JJ)/TIME
C         ENDDO
C        ENDDO
C       ENDDO
C       ENDIF


C       DO I=1,NDIM
C       N=NEL(I,NBREL)
C        DO J=1,NDES
C         P1=SKEF(I,J)*TT21(J)
C         P2=SKEF(I+NDIM,J)*TT21(J)
C         P3=SKEF(I+2*NDIM,J)*TT21(J)
C         SPSIL(1,N)=SPSIL(1,N)-P1
C         SPSIL(2,N)=SPSIL(2,N)-P2
C         SPSIL(3,N)=SPSIL(3,N)-P3
C       ENDDO
C      ENDDO



C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE

C End of subroutine

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
