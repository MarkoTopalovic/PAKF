C=======================================================================
C=========================================================================
      SUBROUTINE ATHERO3DFitting(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,VELOC,WW,YYY)
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
      DIMENSION WW(6,*),YYY(*)


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
      DIMENSION FK1OXM(8),FK1OXM1(8)
      DIMENSION FGAMA(8),FGAMA1(8)
      DIMENSION IGAUS(8)
      DIMENSION FSS(8),DD(5)
      
      DIMENSION WWW(5,5),YYYY(6)
      DIMENSION LE(6),ME(6)

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


      NDFIT=5

      call INDELSSTRES2(NEL,INDEL,NET,NPT,NDIM,ID)


      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
      if (kkorak.eq.0.and. iter.eq.0) 
     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)
	REWIND(MUFILE2)
	ENDIF

      POVRS=0.D0
	povrsila=0.d0
      ZAPRE=0.D0


      DO I=1,JEDN
        YYY(I)=0.D0
       DO J=1,NDFIT
        WW(J,I)=0.D0
       ENDDO
      ENDDO

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
c      ALAMBDA=1.D-5
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
       FK1OXM1(I)=0.D0
       FGAMA(I)=0.D0
       FGAMA1(I)=0.D0
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

      AFIFI(K,N)=AFIFI(K,N)+DD(1)*(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)
     &+PJ(3,K)*PJ(3,N))*WDT

	 OOX=TT21(K+4*NDIM+8)
	 AMM=TT21(K+5*NDIM+8)
	 SS=TT21(K+6*NDIM+8)
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
       
       FK1OXM(K)=FK1OXM(K)+H(K)*(AK1*OOX*AMM)*WDT
       FK1OXM1(K)=FK1OXM1(K)+H(K)*(OOX*AMM)*WDT
       FGAMA(K)=FGAMA(K)+H(K)*(GAMA*(OOX-O_THRES))*WDT
       FGAMA1(K)=FGAMA1(K)+H(K)*((OOX-O_THRES))*WDT
      ENDIF
C===========================================================================
C===========================================================================



      AKVV1(K,N)=AKVV1(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM)

c      AKVV1(K,N)=AKVV1(K,N)+WDT*((Hp(K)*HV1*PJ(1,N)+
c     1Hp(K)*HV2*PJ(2,N)+Hp(K)*HV3*PJ(3,N))*GUSM)

      AKTV1(K,N)=AKTV1(K,N)+WDT*(H(K)*H(N)*ZVXT*GUSM*CC)
      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVYT*GUSM*CC)
      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVZT*GUSM*CC)
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
       FSS(K)=FSS(K)+H(K)*(SS/(1.D0+SS))*WDT
c       FTSS(K)=FTSS(K)+H(K)*(TSS*OOX)*WDT
       OD=1.39*TSS**(-0.118D0)
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
c      AJv=DPW
c      write(iizlaz,*)'AJv= ',AJv
c      AJv=0.D0
	AJs=PP*(ALDL-OOX)+(1.D0-SIGMAF)*AJv*(ALDL+OOX)*0.50D0
c      AJs=0.d0




       FTSS(K)=FTSS(K)+H(K)*(AJs)*WDT
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
       SKEF(I+7*NDIM+8,J+7*NDIM+8)=AF*(AFIFI(I,J))
       SKEF(I+8*NDIM+8,J+8*NDIM+8)=AF*(AFIFI(I,J))
C=========================================================================      
C FITTING PROCEDURE
C=========================================================================      
       J1=LM2(I+4*NDIM+8)
       IF (J1.NE.0) WW(1,J1)=WW(1,J1)+AF*(AKK(I,J))*TT21(J+4*NDIM+8)
       
       J1=LM2(I+5*NDIM+8)
       IF (J1.NE.0) WW(2,J1)=WW(2,J1)+AF*(AKK(I,J))*TT21(J+5*NDIM+8)
       
       J1=LM2(I+6*NDIM+8)
       IF (J1.NE.0) WW(3,J1)=WW(3,J1)+AF*(AKK(I,J))*TT21(J+6*NDIM+8)
       
       J1=LM2(I+7*NDIM+8)
       IF (J1.NE.0) YYY(J1)=YYY(J1)-AF*(AFIFI(I,J))*TT21(J+7*NDIM+8)
       
       J1=LM2(I+8*NDIM+8)
      IF (J1.NE.0) YYY(J1)=YYY(J1)-AF*(AFIFI(I,J))*TT21(J+8*NDIM+8)
       
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
     &ALAMBDA*AMV2(I,J)
C=========================================================================      
C FITTING PROCEDURE
C=========================================================================      
       J1=LM2(I+6*NDIM+8)
C       IF (J1.NE.0) YYY(J1)=YYY(J1)-ALAMBDA*(AMV2(I,J))*TT21(J+6*NDIM+8)
c       IF (J1.NE.0)WW(5,J1)=WW(5,J1)+(AMV2(I,J))*TT21(J+6*NDIM+8)
C=========================================================================      
       
	ENDIF
C=========================================================================      



C       AJK(I,J)=AKVV1(I,J)*CC

      C(I,J)=AMV2(I,J)*CC              

	IF(MAT.EQ.1) THEN
c      SKEF(I+3*NDIM+8,J+3*NDIM+8)=AF*(AKK(I,J)+AKVV1(I,J)*CC)
      SKEF(I+3*NDIM+8,J+3*NDIM+8)=SKEF(I+3*NDIM+8,J+3*NDIM+8)+
     &AF*(AKVV1(I,J)*CC)
      SKEF(I+3*NDIM+8,J)=AF*AKTV1(I,J)
      SKEF(I+3*NDIM+8,J+NDIM)=AF*AKTV2(I,J)
      SKEF(I+3*NDIM+8,J+2*NDIM)=AF*AKTV3(I,J)
      endif
      
      IF (NSTAC.EQ.0) THEN
C=========================================================================   
C Plaque growing   
C=========================================================================   
      SKEF(I+2*NDIM+4,J+2*NDIM+4)=SKEF(I+2*NDIM+4,J+2*NDIM+4)+
     &C(I,J)/TIME/GUSM
c	 IF(MAT.EQ.2) THEN
      SKEF(I+3*NDIM+4,J+3*NDIM+4)=SKEF(I+3*NDIM+4,J+3*NDIM+4)+
     &C(I,J)/TIME/GUSM
      SKEF(I+4*NDIM+4,J+4*NDIM+4)=SKEF(I+4*NDIM+4,J+4*NDIM+4)+
     &C(I,J)/TIME/GUSM
      SKEF(I+5*NDIM+4,J+5*NDIM+4)=SKEF(I+5*NDIM+4,J+5*NDIM+4)+
     &C(I,J)/TIME/GUSM
      SKEF(I+6*NDIM+4,J+6*NDIM+4)=SKEF(I+6*NDIM+4,J+6*NDIM+4)+
     &C(I,J)/TIME/GUSM
      SKEF(I+7*NDIM+4,J+7*NDIM+4)=SKEF(I+7*NDIM+4,J+7*NDIM+4)+
     &C(I,J)/TIME/GUSM
c      ENDIF
C=========================================================================      

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
C       if (mod(kkorak,3).eq.0) then
        F92(I+5*NDIM+8)=F92(I+5*NDIM+8)+FSS(I)
c	 endif
c       if (mod(kkorak,3).eq.2) then
C        F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FGAMA(I)
c        F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FGAMA1(I)
C=========================================================================      
C FITTING PROCEDURE
C=========================================================================      
       J1=LM2(I+4*NDIM+8)
       IF (J1.NE.0) YYY(J1)=YYY(J1)+FTSS(I)

       J1=LM2(I+5*NDIM+8)
       IF (J1.NE.0) YYY(J1)=YYY(J1)+FSS(I)

       J1=LM2(I+6*NDIM+8)
       IF (J1.NE.0) WW(5,J1)=WW(5,J1)-FGAMA1(I)
       IF (J1.NE.0) YYY(J1)=YYY(J1)-ALAMBDA*TT21(I+6*NDIM+8)
C       IF (J1.NE.0) YYY(J1)=YYY(J1)+FGAMA(I)
      
C=========================================================================      
        
C        endif
	 endif
	endif

       IF (MAT.EQ.2) THEN
       F92(I+4*NDIM+8)=F92(I+4*NDIM+8)-FK1OXM(I)
       F92(I+5*NDIM+8)=F92(I+5*NDIM+8)-FK1OXM(I)
       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FK1OXM(I)
       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)-ALAMBDA*TT21(I+6*NDIM+8)
       F92(I+7*NDIM+8)=F92(I+7*NDIM+8)+FK1OXM(I)
	 F92(I+7*NDIM+8)=F92(I+7*NDIM+8)+FFIFI(I)
       F92(I+8*NDIM+8)=F92(I+8*NDIM+8)+FK1OXM(I)
c       F92(I+8*NDIM+8)=F92(I+8*NDIM+8)+FFIFI(I)
c======================================================
C       F92(I+6*NDIM+8)=F92(I+6*NDIM+8)+FGAMA(I)
c======================================================
C=========================================================================      
C FITTING PROCEDURE
C=========================================================================      
       J1=LM2(I+4*NDIM+8)
       IF (J1.NE.0) WW(4,J1)=WW(4,J1)+FK1OXM1(I)
       J1=LM2(I+5*NDIM+8)
       IF (J1.NE.0) WW(4,J1)=WW(4,J1)+FK1OXM1(I)
       J1=LM2(I+6*NDIM+8)
       IF (J1.NE.0) WW(4,J1)=WW(4,J1)-FK1OXM1(I)
       J1=LM2(I+7*NDIM+8)
       IF (J1.NE.0) WW(4,J1)=WW(4,J1)-FK1OXM1(I)
       J1=LM2(I+7*NDIM+8)
       IF (J1.NE.0) YYY(J1)=YYY(J1)+FFIFI(I)
       J1=LM2(I+8*NDIM+8)
       IF (J1.NE.0) WW(4,J1)=WW(4,J1)-FK1OXM1(I)
C=========================================================================      

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
       F92(I+7*NDIM+8)=F92(I+7*NDIM+8)-AFIFI(I,J)*TT21(J+7*NDIM+8)
       F92(I+8*NDIM+8)=F92(I+8*NDIM+8)-AFIFI(I,J)*TT21(J+8*NDIM+8)
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

c       write(iizlaz,*)'left side'
c       do i=1,jedn
c         write(iizlaz,'(5e13.5)')(ww(j,i),j=1,ndfit)
c       enddo
c       write(iizlaz,*)'right side'
c       do i=1,jedn
c         write(iizlaz,'(e13.5)') yyy(i)
c       enddo
c       stop

       DO I=1,NDFIT
       DO J=1,NDFIT
        WWW(I,J)=0.D0
        DO K=1,JEDN
         WWW(I,J)=WWW(I,J)+WW(I,K)*WW(J,K)
        ENDDO
       ENDDO
       ENDDO

       DO I=1,NDFIT
        YYYY(I)=0.D0
        DO K=1,JEDN
         YYYY(I)=YYYY(I)+WW(I,K)*YYY(K)
        ENDDO
       ENDDO
       

      call minv(WWW,NDFIT,det1,le,me);

       DO I=1,NDFIT
        YYY(I)=0.D0
        DO J=1,NDFIT
         YYY(I)=YYY(I)+WWW(I,J)*YYYY(J)
        ENDDO
       ENDDO


       WRITE(IIZLAZ,*)'RESENJA SU: '
       DO I=1,NDFIT
        write(iizlaz,*),YYY(I)
       enddo

       write(iizlaz,*)DD(2),DD(3),DD(4),ALAMBDA,GAMA,AK1
C       stop

C End of subroutine
C      CALL MUMPSRIGHT(SILE,JEDN)
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
