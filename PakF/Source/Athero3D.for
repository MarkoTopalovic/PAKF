C=======================================================================
C=========================================================================
      SUBROUTINE ATHERO3D(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,VELOC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine RACU3D is used for 3D analysis
CE It is used global loop per elements
	COMMON /MUMPS/ IMUMPS, MUFILE,MUFILE2
      COMMON /ATHERO/ HM,ZVXM,ZVYM,HOX
      COMMON /ATHERO1/ IATHERO,AK1,GAMA,O_THRES,ALAMBDA,DD,AKleg,Rw,ALp,
     &Dpw,Pp,SigmaF,PlaqMove
C

      DIMENSION GNODE(2,11,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(11,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)
      DIMENSION VELOC(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21),TT21A(92)
      DIMENSION H(21),HP(8)
      DIMENSION AKVV1(21,21),AKMIV1(21,21),AMV2(21,21),AJV1V1(21,21)
      DIMENSION AJV1V2(21,21),AJV1V3(21,21),AJV2V1(21,21),AJV2V2(21,21)
      DIMENSION AJV2V3(21,21),AJV3V1(21,21),AJV3V2(21,21),AJV3V3(21,21)
      DIMENSION AKTV1(21,21),AKTV2(21,21),AKTV3(21,21),AKK(8,8)
      DIMENSION AKV1P(21,8),AKV2P(21,8),AKV3P(21,8),RS1(21),RS2(21)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(8,8),PENXY(8,8),PENXZ(8,8)
      DIMENSION PENYX(8,8),PENYY(8,8),PENYZ(8,8)
      DIMENSION PENZX(8,8),PENZY(8,8),PENZZ(8,8)
      DIMENSION C(21,21)
	DIMENSION SKE(46*93)
	DIMENSION ATGS1(8,8),ATGS2(8,8),ATGS3(8,8)
C	DIMENSION SKE(12*25)

      DIMENSION AKDIVV(8,8)
      DIMENSION AKMV2(8,8)
      DIMENSION AKMV3(8,8)
      DIMENSION AKTMV2(8,8)
      DIMENSION AKTMV3(8,8)
      DIMENSION AFIFI(8,8),AJK(8,8)

      DIMENSION AK1OX(8,8)
      DIMENSION AK1M(8,8)
      DIMENSION FTSS(8),FFIFI(8)
      DIMENSION FK1OXM(8)
      DIMENSION FGAMA(8)
      DIMENSION IGAUS(8)
      DIMENSION FSS(8),DD(5)

	DATA IGAUS /1,5,4,8,2,6,3,7/

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


      akapa=1.d-18 
      akapa=1.d0 

      call INDELSSTRES2(NEL,INDEL,NET,NPT,NDIM,ID)


      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
      if (kkorak.eq.0.and. iter.eq.0) 
     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)
	REWIND(MUFILE2)
	ENDIF

      POVRS=0.D0
	povrsila=0.d0
      ZAPRE=0.D0


C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET

      DO 125 I=1,92
      TT210(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0

      tez=0.d0
C=========================================================================
      DO 130 KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
	tez=tez+0.125d0*CORD(3,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      LM2(KLM+2*NDIM)=ID(3,NEL(KLM,NBREL))
      IF (KLM.LE.8) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
      LM2(KLM+3*NDIM+8)=ID(5,NEL(KLM,NBREL))
	DO J=4,8
        LM2(KLM+J*NDIM+8)=ID(J+2,NEL(KLM,NBREL))
	ENDDO
 130  CONTINUE
C====================
      mat=nel(ndim+1,nbrel)
c      mat=2
C      ALAMBDA=0.D0
c      dd(1)=1.d0
c
c	if (mat.eq.1) akt=0.075d0
c	if (mat.ge.2) akt=0.0025d0
c      akt=75.d0
c	if (tez.gt.0.d0) akt=1.8d0

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
	DO J=4,8
       TT21(KLM+J*NDIM+8)=GNODE(2,J+2,NEL(KLM,NBREL))
       TT210(KLM+J*NDIM+8)=GNODE(1,J+2,NEL(KLM,NBREL))
	ENDDO
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
c	ATGS1(K,N)=0.D0
c	ATGS2(K,N)=0.D0
c	ATGS3(K,N)=0.D0
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
      AKDIVV(K,N)=0.D0
      AKMV2(K,N)=0.D0
      AKMV3(K,N)=0.D0
      AKTMV2(K,N)=0.D0
      AKTMV3(K,N)=0.D0

      AK1OX(K,N)=0.D0
      AK1M(K,N)=0.D0
	AFIFI(K,N)=0.D0
      AJK(K,N)=0.D0
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
       FSS(I)=0.D0
       FTSS(I)=0.D0
	 FFIFI(I)=0.D0
       FK1OXM(I)=0.D0
       FGAMA(I)=0.D0
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


      IG=0
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
      IG=IG+1
 


       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)

      IF (IALE.EQ.1) CALL ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
      IF (INDAMI.EQ.1) CALL NENJ3D(PJ,TT21,AMI,NDIM,IIZLAZ)
 
      WDT=WT*DET1
      ZAPRE=ZAPRE+WDT
  
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM


C===========================================================================
C Plaque growing
C===========================================================================
	IF (MAT.EQ.2) THEN
c	ZVXW1=DOT(PJ(1,1),TT21(1+6*NDIM+8),NDIM)
c      ZVYW2=DOT(PJ(2,1),TT21(1+7*NDIM+8),NDIM)


c     DIVV=ZVXW2+ZVYW3
c      AKDIVV(K,N)=AKDIVV(K,N)+WDT*(H(K)*DIVV*H(N))
C      AKTMV2(K,N)=AKTMV2(K,N)+WDT*(H(K)*HM*ZVHX(N))
C      AKTMV3(K,N)=AKTMV3(K,N)+WDT*(H(K)*HM*ZVHY(N))
C      AKMV2(K,N)=AKMV2(K,N)+WDT*(H(K)*ZVXM*H(N))
C      AKMV3(K,N)=AKMV3(K,N)+WDT*(H(K)*ZVYM*H(N))
     
C      AK1OX(K,N)=AK1OX(K,N)+AK1*(H(K)*HOX*H(N))*WDT
C      AK1M(K,N)=AK1M(K,N)+AK1*(H(K)*HM*H(N))*WDT

      AFIFI(K,N)=AFIFI(K,N)+(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)
     &+PJ(3,K)*PJ(3,N))*WDT

	 OOX=TT21(K+4*NDIM+8)
	 AMM=TT21(K+5*NDIM+8)
	 SS=TT21(K+6*NDIM+8)
	 NN=NEL(K,NBREL)
	 TSS=DSQRT(PRES(1,NN)**2+PRES(2,NN)**2+PRES(3,NN)**2)
c	 vvv=DSQRT(gnode(2,1,NN)**2+gnode(2,2,NN)**2+gnode(2,3,NN)**2)
C       FSS(K)=FSS(K)+H(K)*(SS/(1.D0+SS))*WDT
Cc       FTSS(K)=FTSS(K)+H(K)*(TSS*OOX)*WDT
Cc       OD=1.39*TSS**(-0.118D0)
	IF (DABS(TSS).LT.1.D-10) tss=1.d0
c	IF (DABS(vvv).LT.1.D-10) vvv=1.d0
       tauOx=1.d-4
Cc	 if (kkorak.gt.1) tauOx=OOX
C       FTSS(K)=FTSS(K)+H(K)*(TSS*tauOx)*WDT
Cc       FTSS(K)=FTSS(K)+H(K)*(OD)*WDT
       
       FK1OXM(K)=FK1OXM(K)+H(K)*(AK1*OOX*AMM)*WDT
       FGAMA(K)=FGAMA(K)+H(K)*(GAMA*(OOX-O_THRES))*WDT
      ENDIF
C===========================================================================
C===========================================================================



      AKVV1(K,N)=AKVV1(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM)

c      AKVV1(K,N)=AKVV1(K,N)+WDT*((Hp(K)*HV1*PJ(1,N)+
c     1Hp(K)*HV2*PJ(2,N)+Hp(K)*HV3*PJ(3,N))*GUSM)

c      AKTV1(K,N)=AKTV1(K,N)+WDT*(H(K)*H(N)*ZVXT*GUSM*CC)
c      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVYT*GUSM*CC)
c      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVZT*GUSM*CC)
      AKTV1(K,N)=AKTV1(K,N)+WDT*(H(K)*H(N)*ZVXT*CC)
      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVYT*CC)
      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVZT*CC)
C================================================================
      AKMIV1(K,N)=AKMIV1(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
C Balansing tensor diffusivity for unsteady case
c      ATGS1(K,N)=ATGS1(K,N)+WDT*HV1*(HV1*PJ(1,K)*PJ(1,N)+
c     1PJ(1,K)*PJ(2,N)*HV2+PJ(1,K)*PJ(3,N)*HV3)*TIME/2.D0
c      ATGS2(K,N)=ATGS2(K,N)+WDT*HV2*(HV1*PJ(2,K)*PJ(1,N)+
c     1PJ(2,K)*PJ(2,N)*HV2+PJ(2,K)*PJ(3,N)*HV3)*TIME/2.D0
c      ATGS3(K,N)=ATGS3(K,N)+WDT*HV3*(HV1*PJ(3,K)*PJ(1,N)+
c     1PJ(3,K)*PJ(2,N)*HV2+PJ(3,K)*PJ(3,N)*HV3)*TIME/2.D0
CC      AKMIV1(K,N)=AKMIV1(K,N)+WDT*(HV1*H(K)*PJ(1,N)+
cC     1H(K)*PJ(2,N)*HV2+H(K)*PJ(3,N)*HV3)*TIME/2.D0
C================================================================
      AKK(K,N)=AKK(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
c     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AKT)
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N)))
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
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT

      NODE=NEL(K,NBREL) 
      NN=NEL(IGAUS(IG),NBREL)
      VELOC(1,NN)=VELOC(1,NN)+PJ(1,K)*TT21(K+7*NDIM+4)/INDEL(NODE)
      VELOC(2,NN)=VELOC(2,NN)+PJ(2,K)*TT21(K+7*NDIM+4)/INDEL(NODE)
      VELOC(3,NN)=VELOC(3,NN)+PJ(3,K)*TT21(K+7*NDIM+4)/INDEL(NODE)
     
c      VELOC(1,NN)=VELOC(1,NN)+PJ(1,K)*TT21(K+7*NDIM+4)
c      VELOC(2,NN)=VELOC(2,NN)+PJ(2,K)*TT21(K+7*NDIM+4)
c      VELOC(3,NN)=VELOC(3,NN)+PJ(3,K)*TT21(K+7*NDIM+4)


  165 CONTINUE

  170 CONTINUE
  175 CONTINUE
  180 CONTINUE


C      DO  I=1,NPT
C       WRITE(IIZLAZ,'(3E13.5)') VELOC(1,I),VELOC(2,I),VELOC(3,I)
C      ENDDO
C      STOP


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
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)

 
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
      
        INDX=0

      DO 250 JBRPS=1,MAXSIL
      IF (NGPSIL(1,JBRPS).EQ.NBREL) THEN 
      CALL STRANA(NEL,NDIM,NGPSIL(1,JBRPS),NPOV)
c      goto 250
C      IF(JBRPS.LE.4) NPOV=1
C      IF(JBRPS.GT.4) NPOV=2
C       WRITE(IIZLAZ,*)'STRANA=',NPOV
        INDX=1
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
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)
C
      WDT=WT*DET
      POVRS=POVRS+WDT


c      CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NGPSIL(6,JBRPS)),
c     &NGPSIL(6,JBRPS),NTABFT,IIZLAZ)

      
      DO 206 K=1,NDIM 
 	ivr=NGPSIL(6,JBRPS)
	if (ivr.eq.4) then
       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(ivr),ivr,
     &NTABFT,IIZLAZ)
	 RS1(K)=RS1(K)+(H(K))*WDT*SF1*FK1
	 RS2(K)=RS2(K)+(H(K))*WDT*SF2*FK1
	 RS3(K)=RS3(K)+(H(K))*WDT*SF3*FK1
c	 povrsila=povrsila+(H(K))*WDT*SF1
c	 povrsila=povrsila+(H(K))*WDT*SF2
c	 povrsila=povrsila+(H(K))*WDT*SF3
c      ENDDO
C      WRITE(IIZLAZ,*)'SF1=',SF1
C       WRITE(IIZLAZ,*)'NPOV',NPOV,'NBREL=',NBREL
C       CALL IWRRF(NGPSIL(1,JBRPS),5,'NGPSI',IIZLAZ)
C       CALL WRRF(RS3,21,'RS3=',IIZLAZ)
       else
C===========================================================================
C Plaque growing
C===========================================================================
       
       
	 ALDL=TT21(K+3*NDIM+8)
	 OOX=TT21(K+4*NDIM+8)
	 AMM=TT21(K+5*NDIM+8)
	 SS=TT21(K+6*NDIM+8)
	 NN=NEL(K,NBREL)
	 TSS=DSQRT(PRES(1,NN)**2+PRES(2,NN)**2+PRES(3,NN)**2)
	 vvv=DSQRT(gnode(2,1,NN)**2+gnode(2,2,NN)**2+gnode(2,3,NN)**2)
       FSS(K)=FSS(K)+H(K)*(SS/(1.D0+SS))*WDT
c       FTSS(K)=FTSS(K)+H(K)*(TSS*OOX)*WDT
       OD=1.39*TSS**(-0.118D0)
	IF (DABS(TSS).LT.1.D-10) TSS=1.D0
	IF (DABS(vvv).LT.1.D-10) vvv=1.D0
	IF (DABS(TSS).LT.1.D-10) OD=1.D0
       tauOx=1.d-6
	AJs=tauOx*1.d-6*TSS
c	PP=2.0d-10
c      SIGMAF=0.997D0
c	ALP=3.D-14
    
      PRESSURE=GNODE(2,4,NEL(K,NBREL))
      Pw=GNODE(2,9,NEL(K,NBREL))

c	PRESSURE=9531.D0
C      DPW=PRESSURE-ALAMBDA
C      DPW=-DPW
      DPW=PRESSURE-Pw
C	DPW=9531.D0
      AJv=ALP*DPW
c      AJv=DPW/tss/1.d9
c      write(iizlaz,*)'AJv= ',AJv
c      AJv=0.D0


	AJs=PP*(ALDL-OOX)+(1.D0-SIGMAF)*AJv*(ALDL+OOX)*0.50D0
      AJs=AJs*1.d5


c       sign=0.d0
c       if (cord(2,NGPSIL(2,JBRPS)).gt.2.4e-2.and.
c     &cord(2,NGPSIL(2,JBRPS)).lt.3.32e-2) then
       sign=1.d0
c       endif
       
       FTSS(K)=FTSS(K)+H(K)*(sign*WDT/tss/1.d19)
c       FTSS(K)=FTSS(K)+H(K)*(AJs*WDT)
       FFIFI(K)=FFIFI(K)+H(K)*(AJv)*WDT
c===============================================================
      endif
 206  continue

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
C=========================================================================      
C Plaque growing   
C=========================================================================   
	IF (MAT.EQ.1) THEN
	 SKEF(I+3*NDIM+8,J+3*NDIM+8)=AF*(AKK(I,J)*DD(1)+AJK(I,J))
	ELSEIF (MAT.EQ.2) THEN
	 SKEF(I+4*NDIM+8,J+4*NDIM+8)=AF*(AKK(I,J)*DD(2))
	 SKEF(I+5*NDIM+8,J+5*NDIM+8)=AF*(AKK(I,J)*DD(3))
	 SKEF(I+6*NDIM+8,J+6*NDIM+8)=AF*(AKK(I,J)*DD(4))
       SKEF(I+7*NDIM+8,J+7*NDIM+8)=AF*(AFIFI(I,J))*akapa
       SKEF(I+8*NDIM+8,J+8*NDIM+8)=AF*(AFIFI(I,J))
C=========================================================================      
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
	SKEF(I+6*NDIM+8,J+6*NDIM+8)=SKEF(I+6*NDIM+8,J+6*NDIM+8)+
     &ALAMBDA*AMV2(I,J)/gusm
	ENDIF
C=========================================================================      



c      AJK(I,J)=AKVV1(I,J)*CC

      C(I,J)=AMV2(I,J)*CC              

	IF(MAT.EQ.1) THEN
c      SKEF(I+3*NDIM+8,J+3*NDIM+8)=AF*(AKK(I,J)+AKVV1(I,J)*CC)
      SKEF(I+3*NDIM+8,J+3*NDIM+8)=SKEF(I+3*NDIM+8,J+3*NDIM+8)+
     &AF*(AKVV1(I,J)*CC/gusm)
      SKEF(I+3*NDIM+8,J)=AF*AKTV1(I,J)
      SKEF(I+3*NDIM+8,J+NDIM)=AF*AKTV2(I,J)
      SKEF(I+3*NDIM+8,J+2*NDIM)=AF*AKTV3(I,J)
      endif
      
      IF (NSTAC.EQ.0) THEN
C=========================================================================   
C Plaque growing   
C=========================================================================   
      SKEF(I+2*NDIM+8,J+2*NDIM+8)=SKEF(I+2*NDIM+8,J+2*NDIM+8)+
     &C(I,J)/TIME/GUSM
c	 IF(MAT.EQ.2) THEN
      SKEF(I+3*NDIM+8,J+3*NDIM+8)=SKEF(I+3*NDIM+8,J+3*NDIM+8)+
     &C(I,J)/TIME/GUSM
      SKEF(I+4*NDIM+8,J+4*NDIM+8)=SKEF(I+4*NDIM+8,J+4*NDIM+8)+
     &C(I,J)/TIME/GUSM
      SKEF(I+5*NDIM+8,J+5*NDIM+8)=SKEF(I+5*NDIM+8,J+5*NDIM+8)+
     &C(I,J)/TIME/GUSM
      SKEF(I+6*NDIM+8,J+6*NDIM+8)=SKEF(I+6*NDIM+8,J+6*NDIM+8)+
     &C(I,J)/TIME/GUSM
      SKEF(I+7*NDIM+8,J+7*NDIM+8)=SKEF(I+7*NDIM+8,J+7*NDIM+8)+
     &C(I,J)/TIME/GUSM
c      ENDIF
C=========================================================================      

c      SKEF(I+3*NDIM+8,J+3*NDIM+8)=SKEF(I+3*NDIM+8,J+3*NDIM+8)+
c     &AMV2(I,J)*CC/TIME
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
      if (MAT.EQ.1) THEN
      DO 270 I=1,NDIM
      DO 265 J=1,NDIM
C       SKEF(I,J)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J)+ATGS1(I,J))
       SKEF(I,J)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J))
       SKEF(I,J+NDIM)=AF*AJV1V2(I,J)
       SKEF(I,J+2*NDIM)=AF*AJV1V3(I,J)
      IF (J.LE.8) THEN
       SKEF(I,J+3*NDIM)=AF*AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AF*AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AF*AKV3P(I,J)
      ENDIF
      IF (I.LE.8) THEN
       SKEF(I+3*NDIM,J)=AF*AKV1P(J,I)
       SKEF(I+3*NDIM,J+NDIM)=AF*AKV2P(J,I)
       SKEF(I+3*NDIM,J+2*NDIM)=AF*AKV3P(J,I)
      ENDIF

       SKEF(I+NDIM,J)=AF*AJV2V1(I,J)
       SKEF(I+NDIM,J+NDIM)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV2V2(I,J))
c     &+ATGS2(I,J))
       SKEF(I+NDIM,J+2*NDIM)=AF*AJV2V3(I,J)

      SKEF(I+2*NDIM,J)=AF*AJV3V1(I,J)
      SKEF(I+2*NDIM,J+NDIM)=AF*AJV3V2(I,J)
      SKEF(I+2*NDIM,J+2*NDIM)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV3V3(I,J))
c    &+ATGS3(I,J))
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
        SKEF(I1,J1)=SKEF(I1,J1)+AF*PENXX(I,J)
        SKEF(I1,J2)=SKEF(I1,J2)+AF*PENXY(I,J)
        SKEF(I1,J3)=SKEF(I1,J3)+AF*PENXZ(I,J)
        SKEF(I2,J1)=SKEF(I2,J1)+AF*PENYX(I,J)
        SKEF(I2,J2)=SKEF(I2,J2)+AF*PENYY(I,J)
        SKEF(I2,J3)=SKEF(I2,J3)+AF*PENYZ(I,J)
        SKEF(I3,J1)=SKEF(I3,J1)+AF*PENZX(I,J)
        SKEF(I3,J2)=SKEF(I3,J2)+AF*PENZY(I,J)
        SKEF(I3,J3)=SKEF(I3,J3)+AF*PENZZ(I,J)
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
C DESNA STRANA:
      if (mat.eq.1) then 
      DO K=0,2
        II=I+K*NDIM
      DO 285 J=1,NDIM
        JJ=J+K*NDIM
        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
	IF (K.EQ.0) F92(II)=F92(II)!-ATGS1(I,J)*TT21(JJ)
	IF (K.EQ.1) F92(II)=F92(II)!-ATGS2(I,J)*TT21(JJ)
	IF (K.EQ.2) F92(II)=F92(II)!-ATGS3(I,J)*TT21(JJ)
      IF (NSTAC.EQ.0) THEN
        F92(II)=F92(II)-AMV2(I,J)*(TT21(JJ)-TT210(JJ))/TIME
C       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
 285  CONTINUE							      
 290  CONTINUE
      ENDDO
      endif

C=========================================================================
      DO 291 I=1,NDIM
      if (indx.eq.1) then
c	 write(iizlaz,*) ftss(i),FSS(I),FGAMA(I)
c       if (mod(kkorak,3).eq.1) then
        if (mat.eq.1) then
        F92(I+3*NDIM+8)=F92(I+3*NDIM+8)+FTSS(I)
        elseif (mat.eq.2) then
        F92(I+4*NDIM+8)=F92(I+4*NDIM+8)+FTSS(I)
c	 endif
c       if (mod(kkorak,3).eq.0) then
        F92(I+5*NDIM+8)=F92(I+5*NDIM+8)+FSS(I)
c	 endif
c       if (mod(kkorak,3).eq.2) then
        F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FGAMA(I)
        endif
c	 endif
	endif

       IF (MAT.EQ.2) THEN
       F92(I+4*NDIM+8)=F92(I+4*NDIM+8)-FK1OXM(I)
       F92(I+5*NDIM+8)=F92(I+5*NDIM+8)-FK1OXM(I)
       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FK1OXM(I)
       
c       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)-ALAMBDA*TT21(I+6*NDIM+8)
     
       F92(I+7*NDIM+8)=F92(I+7*NDIM+8)+FK1OXM(I)
	 F92(I+7*NDIM+8)=F92(I+7*NDIM+8)+FFIFI(I)
       F92(I+8*NDIM+8)=F92(I+8*NDIM+8)+FK1OXM(I)*dd(5)
c       F92(I+8*NDIM+8)=F92(I+8*NDIM+8)+FFIFI(I)
c======================================================
C       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FGAMA(I)
c======================================================

	ENDIF


      DO 286 J=1,NDIM
C=========================================================================   
C Plaque growing   
C=========================================================================   
C=========================================================================   
      IF (MAT.EQ.2) THEN 
	 DO KK=4,6
        II=I+KK*NDIM+8
        JJ=J+KK*NDIM+8
        F92(II)=F92(II)-AKK(I,J)*DD(KK-2)*TT21(JJ)
	 ENDDO
       F92(I+7*NDIM+8)=F92(I+7*NDIM+8)-akapa*AFIFI(I,J)*TT21(J+7*NDIM+8)
       F92(I+8*NDIM+8)=F92(I+8*NDIM+8)-AFIFI(I,J)*TT21(J+8*NDIM+8)
       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)-
     &(ALAMBDA*AMV2(I,J)/gusm)*TT21(J+6*NDIM+8)
      ENDIF

C=========================================================================   
C Plaque growing   
C========================================================================= 
C	II=I+6*NDIM+4
C	JJ=J+6*NDIM+4

C      F92(II)=F92(II)-(AKMIV2(I,J))*TT21(JJ)
C      F92(II+NDIM)=F92(II+NDIM)-(AKMIV3(I,J))*TT21(JJ+NDIM)
C      F92(II)=F92(II)-A12(I,J)*TT21(JJ+NDIM)
C      F92(II+NDIM)=F92(II+NDIM)-A21(I,J)*TT21(JJ)


C      F92(I+4*NDIM+4)=F92(I+4*NDIM+4)-AKMV2(I,J)*TT21(J)
C     1-AKMV3(I,J)*TT21(J+NDIM)
C      F92(I+4*NDIM+4)=F92(I+4*NDIM+4)-AKTMV2(I,J)*TT21(J)
C     1-AKTMV3(I,J)*TT21(J+NDIM)
C      F92(I+4*NDIM+4)=F92(I+4*NDIM+4)-AKDIVV(I,J)*TT21(J+4*NDIM+4)
     
C========================================================================= 
C========================================================================= 


 286  CONTINUE							      
 291  CONTINUE
    

C=========================================================================
C ovde privremeno za Nikolin primer, 27 March, 2006
      IF (MAT.EQ.1) THEN 
      DO 294 I=1,NDIM
      II=I+3*NDIM+8
C	 F92(II)=F92(II)+RS1(I) 
      DO 292 J=1,NDIM
      JJ=J+3*NDIM+8
       F92(II)=F92(II)-AKTV1(I,J)*TT21(J)-AKTV2(I,J)*TT21(J+NDIM)
     &-AKTV3(I,J)*TT21(J+2*NDIM)-AKK(I,J)*DD(1)*TT21(JJ)
c     &-AKTV3(I,J)*TT21(J+2*NDIM)-Afifi(I,J)*DD(1)*TT21(JJ)

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
      DO 315 I=1,NDIM
       F92(I)=F92(I)+RB1(I)
       F92(I+NDIM)=F92(I+NDIM)+RB2(I)
       F92(I+2*NDIM)=F92(I+2*NDIM)+RB3(I)
 315  CONTINUE
C==========================================================================
C       CALL WRR(RS1,21,'RS1= ')
C       CALL WRR(FPOM,21,'FPOM=')
C==========================================================================
C POVRSINSKE SILE SA DESNE STRANE
      DO 320 I=1,NDIM
C ZADAVANJE PRITISAKA NA DESNOJ STRANI
       F92(I)=F92(I)+RS1(I)
       F92(I+NDIM)=F92(I+NDIM)+RS2(I)
       F92(I+2*NDIM)=F92(I+2*NDIM)+RS3(I)
320   CONTINUE
C==========================================================================
      endif



      if (mat.eq.1) then 
      CALL SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)
      endif

c      write(iizlaz,*)'nbrel=',nbrel
c      CALL WRRF(F92,NDES,'F92= ',IIZLAZ)      
c      CALL WRRF(afifi,8*8,'afifi=',IIZLAZ)      
c      CALL WRRF(tt21,NDES,'tt21=',IIZLAZ)      
c      CALL WRRF(SKEF,NDES*NDES,'SKEF=',IIZLAZ)      
c      stop
      
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

      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,NETIP)
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
C
C
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
C
C

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE
C End of subroutine
      CALL MUMPSRIGHT(SILE,JEDN)
c      CALL WRRF(SILE,JEDN,'desno=',IIZLAZ)
c      stop
c      write(iizlaz,*)'povrs= ',povrs
c      write(iizlaz,*)'povrsila= ',povrsila
c      write(iizlaz,*)'zapre= ',zapre
c	stop
c      DO  I=1,NPT
c       WRITE(IIZLAZ,'(i10,i10)') i,indel(i)
c      ENDDO
c      STOP

      END
C=======================================================================
C=======================================================================
C=========================================================================
      SUBROUTINE PLAQFORM3D(CCORD,NPT,TIME,VELOC,NEL,NDIM,iizlaz) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION CCORD(3,*),VELOC(3,*)
      DIMENSION NEL(NDIM+1,*)	

c      DO  I=1,NPT
c       WRITE(IIZLAZ,'(i10,3E13.5)') i,VELOC(1,I),VELOC(2,I),VELOC(3,I)
c      ENDDO
c      STOP

      
      DO 100 I=1,NPT
	 CCORD(1,I)=CCORD(1,I)-TIME*VELOC(1,I)*1.D-9
	 CCORD(3,I)=CCORD(3,I)-TIME*VELOC(3,I)*1.D-9
100   CONTINUE	

C	DO I=1,NPT
C	 IF (NEL(NDIM+1,I).EQ.2) THEN
C	   NEND=NEL(1,I) 
C	   GOTO 200
C	  ENDIF
C	ENDDO
C	
C 200	NENDD=NEL(4,I)-NEL(1,I)



      DO K=0,9
	NSTART=1+K*58
	NEND=4+K*58
	NSTEP=4

      dp=1.d0/(NSTEP-1)
      do j=0,NSTEP*2-2
      n=NSTART+j*NSTEP
	m=NEND+j*NSTEP
        p=0.d0       
	do i=n+1,m-1 
	  p=p+dp
	  ccord(1,i)=(1.d0-p)*ccord(1,n)+p*ccord(1,m)
	  ccord(3,i)=(1.d0-p)*ccord(3,n)+p*ccord(3,m)
	 enddo 
	enddo
	ENDDO


C End of subroutine
      END
C=======================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE INDELSSTRES2(NEL,INDEL,NET,NPT,NDIM,ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
     	
      DIMENSION NEL(NDIM+1,*),INDEL(*),ID(11,*)
	DIMENSION NID(21),N(21),NWALL(6),ITR(6,4)

     
C
CE Subroutine INITIA is used for initialisation all global variable
C


C INITIALISATION:
       DO I=1,NPT
        INDEL(I)=0
       ENDDO


c======================================================
c======================================================
c       DO NBREL=1,NET
c        DO I=1,NDIM
c         NODE=NEL(I,NBREL)
c           INDEL(NODE)=INDEL(NODE)+1
c        ENDDO        
c       ENDDO

c       return
c======================================================
c======================================================

      do 100 nbrel=1,net
      mat=nel(ndim+1,nbrel)
      if (mat.ne.1) goto 100

	DO I=1,NDIM
	  N(I)=NEL(I,NBREL)
        NID(I)=0
      IF(ID(1,N(I)).EQ.0.AND.ID(2,N(I)).EQ.0.AND.ID(3,N(I)).EQ.0)
     &NID(I)=1
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


      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300
      do j=1,4
	 node=itr(npov,j)
	 indel(node)=indel(node)+1
      enddo
300   continue
100   continue



      END
C==========================================================================
C==========================================================================
C=======================================================================
      SUBROUTINE PLAQPEN3(GNODE,NEL,CCORD,ID,NET,NDIM,IIZLAZ,AMI,ISRPS,
     &INDEL,TIME,NPT,CORD,IBRGT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION GNODE(2,11,*),NEL(NDIM+1,*),CCORD(3,*),ID(11,*),INDEL(*)
      DIMENSION CORD(3,*)
      DIMENSION TT21(92),PJ(3,21),H(21),HP(10),SHEAR(3)
	dimension ck(21,3)

      DIMENSION XG(15),WGT(15),NREF(6)



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

C
CE Subroutine PLAQPEN3 is used for calculation of moving wall mesh due to plaque development
C


      DO I=1,NPT
       INDEL(I)=0
      ENDDO

       DO NBREL=1,NET
        DO I=1,NDIM
         NODE=NEL(I,NBREL)
         INDEL(NODE)=INDEL(NODE)+1
        ENDDO        
       ENDDO


   
      DO 505 NBREL=1,NET

      MAT=NEL(NDIM+1,NBREL)
      IF (MAT.NE.2) GOTO 505


      DO I=1,92
        TT21(I)=0.D0
      ENDDO


C=========================================================================
      DO  KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))

      TT21(KLM)=GNODE(2,10,NEL(KLM,NBREL))
C      LM2(KLM)=ID(1,NEL(KLM,NBREL))
C      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
C      LM2(KLM+2*NDIM)=ID(3,NEL(KLM,NBREL))
C      IF (KLM.LE.8) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
      ENDDO
C=======================================================================
C=======================================================================

      NGAUSX=IBRGT
      NGAUSY=IBRGT
      NGAUSZ=IBRGT
C INTEGRACIJA U GAUSOVIM TACKAMA
C      DO 180 I=1,IBRGT
C      DO 170 J=1,IBRGT
C      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0)
C      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ


      IG=0
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
      IG=IG+1
 


       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)


       DO 150 I=1,NDIM
        NODE=NEL(I,NBREL)
c        write(iizlaz,*)'node',node
        if (cord(1,node).ne.0.d0) then
c         CCORD(1,NODE)=CCORD(1,NODE)+PJ(1,I)*TIME*TT21(I)/INDEL(NODE)
        endif
        CCORD(2,NODE)=CCORD(2,NODE)+PJ(2,I)*TIME*TT21(I)/INDEL(NODE)
        if (cord(3,node).ne.0.d0) then
        CCORD(3,NODE)=CCORD(3,NODE)+PJ(3,I)*TIME*TT21(I)/INDEL(NODE)
        endif
150    continue

170    CONTINUE
175    CONTINUE
180    CONTINUE

505    CONTINUE

c       do i=1,npt
c        write(iizlaz,*)i,ccord(1,i)-cord(1,i)
c        write(iizlaz,*)i,ccord(2,i)-cord(2,i)
c        write(iizlaz,*)i,ccord(3,i)-cord(3,i)
c       enddo 
c       stop
      

C END OF LOOP PER ELEMENT    
      END
C=========================================================================
C=======================================================================
      SUBROUTINE MESHSOL3D(GNODE,ID,IDS,CORD,CCORD,NEL,MAXAS,TT1,ALEVO,
     &SILE,NPT,NET,NDIM,DESNO,MHTS,IIZLAZ)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION TT1(*),ALEVO(*),SILE(*),CORD(3,*),CCORD(3,*)
      DIMENSION GNODE(2,11,*),DESNO(*)
      DIMENSION NEL(NDIM+1,*),MAXAS(*),MHTS(*)
      DIMENSION IDS(4,*),ID(11,*)

      DIMENSION SKEF(24,24)
      DIMENSION SKE(12*25)
      DIMENSION CK(21,3),TT21(142),TT210(142),TT21W(63),PJ(3,21)
      DIMENSION H(21),HP(8),PJP(3,8),F36(142)
      DIMENSION B(6,63)
      DIMENSION LM2(142)
      DIMENSION DT(6,6),DTM(6)
      DIMENSION L(24),M(24),BB(24)
      DIMENSION STRAIN(6),HE(21,4),XJ(3,3),TAU(6)
      DIMENSION XG(15),WGT(15),NREF(6)
      

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



   
       EEE=1.D14
       ANI=0.3d0

C ODREDJIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIM.EQ.8) IBRGT=2

      NGAUSX=IBRGT
      NGAUSY=IBRGT
      NGAUSZ=IBRGT

      NDES=NDIM*3

      Ymax=-1.d10
      DO I=1,NPT
       DO J=1,4
        IDS(J,I)=0
       ENDDO
       if (cord(2,i).gt.Ymax) Ymax=cord(2,i)
      ENDDO


      do nbrel=1,net
       MAT=NEL(NDIM+1,NBREL)
       IF (MAT.EQ.1) THEN
        DO I=1,NDIM
         NODE=NEL(I,NBREL)
         IDS(4,NODE)=1
         IDS(1,NODE)=IDS(1,NODE)+1
        ENDDO
       endif
      enddo


       
      do nbrel=1,net
       MAT=NEL(NDIM+1,NBREL)
       IF (MAT.EQ.2) THEN
        DO I=1,NDIM
         NODE=NEL(I,NBREL)
         if (ABS(IDS(4,NODE)).eq.1) then
           ids(4,node)=-1
         else 
           ids(4,node)=2
         endif
        ENDDO 
       endif
      enddo

  
      do i=1,npt
       if (ids(4,i).eq.1.and.ids(1,i).le.4.and.
     &(cord(2,i).lt.1.e-7).or.dabs(cord(2,i)-Ymax)
     &.lt.1.e-7) ids(4,i)=3
      enddo

c      do i=1,npt
c       write(iizlaz,*) i,ids(4,i)
c      enddo
c      stop


     
      DO I=1,NPT
       IF (IDS(4,I).LE.1) THEN
         IDS(1,I)=0
         IDS(2,I)=0
         IDS(3,I)=0
       ELSE
         IDS(1,I)=1
         IDS(2,I)=1
         IDS(3,I)=1
       ENDIF
       if (dabs(cord(1,i)).lt.1.e-7) then
         ids(1,i)=1
       endif
       if (dabs(cord(3,i)).lt.1.e-7) then
         ids(3,i)=1
       endif
      ENDDO
      
      KK=0
      DO N=1,NPT
      DO JJ=1,3
       IF (IDS(JJ,N).EQ.0) THEN
        KK=KK+1
        IDS(JJ,N)=KK
       ELSE
        IDS(JJ,N)=0
       ENDIF
      ENDDO
      ENDDO
      
      NEQ=KK
      
      CALL MAXATE2(MAXAS,MHTS,IDS,NEL,NET,NDIM,NEQ,NWK,3)      
      CALL CLEARR(ALEVO,NWK)
      CALL CLEARR(SILE,NEQ)
      CALL DT6(DT,EEE,ANI)

CE MAIN LOOP OVER ELEMENTS
      DO 400 NBREL=1,NET
      MAT=NEL(NDIM+1,NBREL)
      IF (MAT.NE.1) GOTO 400

	DO I=1,3*NDIM
         LM2(I)=0
      ENDDO

       IAXIS=4


      DO 125 I=1,NDES
C      TT210(I)=0.D0
 125  TT21(I)=0.D0
C=========================================================================
      JJ=-2
       DO 130 KLM=1,NDIM
       NODE=NEL(KLM,NBREL)
       JJ=JJ+3
C       DO NR=1,3
C       IF (IDS(NR,NEL(KLM,NBREL)) .NE. 0) THEN
C        TT21(JJ+NR-1)=TT1(IDS(NR,NEL(KLM,NBREL)))
C       ENDIF
C      ENDDO
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
       JX=IDS(1,NODE)
       JY=IDS(2,NODE)
       JZ=IDS(3,NODE)
       LM2(JJ)=JX
       LM2(JJ+1)=JY
       LM2(JJ+2)=JZ
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

      IG=0
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
      IG=IG+1
 


       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)


      WDT=WT*DET1

      DO I=1,6
        DO J=1,63
         B(I,J)=0.D0
        ENDDO
      ENDDO

      JJ=-2
      DO I=1,NDIM
       JJ=JJ+3
        B(1,JJ)=PJ(1,I)
        B(2,JJ+1)=PJ(2,I)
        B(3,JJ+2)=PJ(3,I)
        B(4,JJ)=PJ(2,I)
        B(4,JJ+1)=PJ(1,I)
        B(5,JJ+1)=PJ(3,I)
        B(5,JJ+2)=PJ(2,I)
        B(6,JJ)=PJ(3,I)
        B(6,JJ+2)=PJ(1,I)
      ENDDO

C        CALL CLEAR(STRAIN,6) 
C        CALL CLEAR(TAU,6) 
C      IF (IATYP.EQ.1) THEN
C       DO II=1,6
C        DO NN=1,NDIM*3
C         STRAIN(II)=STRAIN(II)+B(II,NN)*TT21(NN)
C        ENDDO
C       ENDDO
C      ENDIF

      
      CALL INTEGK(SKE,B,DT,LM2,WDT,NDES,6,6)


c       k=0
c       do i=1,ndes
c        do j=1,ndes
c        k=k+1
c        if (i.eq.j) WRITE(IIZLAZ,*) 'SKE',i,ig,nbrel,SKE(k)
c       enddo
c       enddo
      

c       do ii=1,6
c        write(iizlaz,'(i5,x,6e13.5)')nbrel,(b(ii,ii))
c       enddo

c      IF (NMODM.EQ.0) THEN
c      DO II=1,6
c       TAU(II)=0.D0
c        DO JJ=1,6
c         TAU(II)=TAU(II)+DT(II,JJ)*STRAIN(JJ)
c        ENDDO
c      ENDDO
c      
c      ENDIF

c      IF (NDIN.EQ.0) CALL INTEGF(SILE,B,TAU,LM2,-WDT,NDESS,6)
     


  170 CONTINUE
  175 CONTINUE
  180 CONTINUE









C=========================================================================

C       DO I=1,NDES
C        DO J=1,NDES
C         F36(I)=F36(I)-SKEF(I,J)*TT21(J)
C        ENDDO
C       ENDDO       


c       CALL PSKEFN(SKEF,SKE,NDES)	  
C       DO I=1,NDES*(NDES+1)*0.5
C        WRITE(IIZLAZ,*) 'SKE',NBREL,SKE(I)
C       ENDDO
       CALL MATSTE (ALEVO,MAXAS,SILE,SKE,F36,LM2,NDES,1)

C      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F36,MAXA,LM2,NDES,1)

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE
c      STOP

      DO 405 I=1,NPT
      if (ids(4,i).ne.-1) goto 405
       DO J=1,3
       JJ=IDS(J,I)       
       IF (JJ.ne.0) THEN
        SILE(JJ)=(CCORD(J,I)-CORD(J,I))*1.0D35
        ALEVO(MAXAS(JJ))=1.0D35
c        WRITE(IIZLAZ,*)'SILE',I,JJ,(CCORD(J,I)-CORD(J,I))
c        WRITE(IIZLAZ,*)JJ,'ALEVO',ALEVO(MAXAS(JJ))
       ENDIF
       ENDDO
 405   CONTINUE
c       STOP
C       DO 2423 I=1,NWK
C 2423    WRITE(IIZLAZ,*)I,'ALEVO',ALEVO(I)

c      DO 424 I=1,NEQ
c        WRITE(IIZLAZ,*)I,'ALEVO',ALEVO(MAXAS(I))
c       WRITE(IIZLAZ,*)I,'SILE',SILE(I)
c 424   continue       

c     STOP
      CALL RESENF(ALEVO,SILE,MAXAS,NEQ,NWK,1)
      CALL RESENF(ALEVO,SILE,MAXAS,NEQ,NWK,2)


c===================================================================
c      DO I=1,NPT
c       DO J=1,3
c        JJ=IDS(J,I)       
c        IF (JJ.ne.0) THEN
c         WRITE(IIZLAZ,*)I,JJ,'SILE',SILE(jj)
c        ENDIF
c       ENDDO
c       ENDDO

C      DO I=1,NEQ
C        WRITE(IIZLAZ,*)'MAXA',MAXAS(I)
C      ENDDO

c      STOP
c===================================================================

      DO I=1,NPT
       DO J=1,3
        JJ=IDS(J,I)
        IF (JJ.NE.0) THEN
         CCORD(J,I)=CORD(J,I)+SILE(JJ)
        ENDIF
       ENDDO
      ENDDO
      

      END
C=======================================================================
C======================================================================
C=======================================================================
      SUBROUTINE DT6(DT,EEE,ANI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DT(6,*)

      DO I=1,6
       DO J=1,6
        DT(I,J)=0.D0
       ENDDO
      ENDDO

      EMNOZI=EEE*(1.D0-ANI)/((1.D0+ANI)*(1.D0-2.D0*ANI))
      DT(1,1)=1.D0*EMNOZI
      DT(1,2)=ANI*EMNOZI/(1.D0-ANI)
      DT(1,3)=DT(1,2)
      DT(2,1)=DT(1,2)
      DT(2,2)=DT(1,1)
      DT(2,3)=DT(2,1)
      DT(3,1)=DT(1,3)
      DT(3,2)=DT(2,3)
      DT(3,3)=DT(1,1)
      DT(4,4)=EMNOZI*(1.D0-2.D0*ANI)/(2.D0*(1.D0-ANI))
      DT(5,5)=DT(4,4)
      DT(6,6)=DT(4,4)

      END
C=======================================================================
C=======================================================================
      SUBROUTINE INTEGK(A,B,C,LM,W,II,JJ,NBDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE MATRIX - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU MATRICE - A
C .
CE.       A(NWE)   - MATRIX
CE.       B(JJ,II) - MATRIX
CE.       C(JJ,JJ) - MATRIX
CE.       LM(II)   - EQUATION NUMBER
CE.       W        - WEIGHTS KOEFICIENTS
CS.       A(NWE)   - MATRICA
CS.       B(JJ,II) - MATRICA
CS.       C(JJ,JJ) - MATRICA
CS.       LM(II)   - BROJEVI JEDNACINA
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (BT * C * B) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(NBDIM,*),C(NBDIM,*),P(63),LM(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEGK'
      IJ=0
      DO 10 I=1,II
C         IF(LM(I).GT.0) THEN
         IF(LM(I).GT.0.OR.LM(I).LT.0) THEN
            DO 20 J=1,JJ
               P(J)=0.D0
            DO 20 K=1,JJ
   20       P(J)=P(J)+B(K,I)*C(K,J)
         ENDIF
      DO 10 J=I,II
         IJ=IJ+1
         IF(LM(I).EQ.0.OR.LM(J).EQ.0) GO TO 10
         X=0.D0
         DO 30 K=1,JJ
   30    X=X+P(K)*B(K,J)
         A(IJ)=A(IJ)+X*W
   10 CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE SKEKUU(SKE,AKUU,M,N,NKUU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
     
      DIMENSION SKE(*),AKUU(NKUU,*)

      NN=0
      DO I=1,M
       DO J=I,N
        NN=NN+1
        AKUU(I,J)=SKE(NN)
        AKUU(J,I)=SKE(NN)
       ENDDO
      ENDDO      

      END
C=======================================================================
C=======================================================================
C==========================================================================
      SUBROUTINE MAXATE2(MAXA,MHT,ID,NEL,NE,NTE,JEDN,NWK,IDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    PODPROGRAM ZA FORMIRANJE VEKTORA VISINA STUBOVA I MAXA
CS    KONACNO SE SMESTAJU U ISTI PROSTOR
CE    PROGRAM TO DETERMINE COLUMN HEIGHTS VECTOR AND MAXA
C
C
      DIMENSION MAXA(*),MHT(*),NEL(NTE+1,*),LM(63),ID(IDIM+1,*)
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
c      if (nel(nte+1,nlm).ne.1) goto 100
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
