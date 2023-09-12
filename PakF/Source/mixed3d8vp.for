#define MUMPS_CLUSTER .FALSE.
C=========================================================================
      SUBROUTINE Mixed3D8VP(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,KKORAK,NUMPASS,METOD,id1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine Mixed3D8VP is used for 3D analysis in mixed formulation
CE It is used global loop per elements
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
	COMMON /TRANSP/ INDFL,LID1
	COMMON /ALPHA_SEG/ ALPHAU,ALPHAV,ALPHAW,ALPHAP
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*),id1(6,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21),TT21A(3*21)
      DIMENSION H(21),HP(8)
      DIMENSION AKVV1(21,21),AKMIV1(21,21),AMV2(21,21),AJV1V1(21,21)
      DIMENSION AJV1V2(21,21),AJV1V3(21,21),AJV2V1(21,21),AJV2V2(21,21)
      DIMENSION AJV2V3(21,21),AJV3V1(21,21),AJV3V2(21,21),AJV3V3(21,21)
      DIMENSION AKTV1(21,21),AKTV2(21,21),AKTV3(21,21),AKK(21,21)
	DIMENSION AKVV1P(21,21),AKVV2P(21,21),AKVV3P(21,21)
	DIMENSION AKMIV1P(21,21),AKMIV2P(21,21),AKMIV3P(21,21)
C	DIMENSION AJK(21,21)
      DIMENSION AKV1P(21,8),AKV2P(21,8),AKV3P(21,8),RS1(21),RS2(21)
	DIMENSION RSVP(21)
	DIMENSION AKPP(8,8)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(21,21),PENXY(21,21),PENXZ(21,21)
      DIMENSION PENYX(21,21),PENYY(21,21),PENYZ(21,21)
      DIMENSION PENZX(21,21),PENZY(21,21),PENZZ(21,21)
      DIMENSION C(21,21)
	DIMENSION SKE(46*93),TT21P(92)
	DIMENSION ATGS1(8,8),ATGS2(8,8),ATGS3(8,8)
	DIMENSION AKIIX(8,8),AKIIY(8,8),AKIIZ(8,8)
	DIMENSION RP1(8),RP2(8),RP3(8)
	DIMENSION NELL(5)
C	DIMENSION SKE(12*25)



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
      NDIMP=8
      


C=========================================================================
c      if (ipass.gt.0) then
      DO NODE=1,NPT
       DO NZDT=1,NUMZAD
        IF(NODE.EQ.NZAD(1,NZDT)) THEN
         CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
     &NTABFT,IIZLAZ)
         MESTO=NZAD(2,NZDT)
c	if (mesto.eq.3) fk1=0.d0
	   GNODE(2,MESTO,NODE)=ZADVRE(NZDT)*FK1
c	   GNODE(1,MESTO,NODE)=ZADVRE(NZDT)*FK1
        ENDIF  
	 ENDDO
	ENDDO
c=======================================================================
c      endif
  



      call INDELSSTRES(NEL,INDEL,NET,NPT,NDIM,ID)


      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
#if(MUMPS_CLUSTER)
      if ((kkorak.eq.1.and.iter.eq.0) .or. METOD.eq.4)
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)
     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,5)
	!REWIND(MUFILE2) bogdan
#else
        CALL sparseassembler_init(0)
#endif
	ENDIF

      POVRS=0.D0
	povrsila=0.d0
      ZAPRE=0.D0


C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET


      if (nel(5,nbrel).eq.0) THEN
      NELL(1)=NEL(1,NBREL)
      NELL(2)=NEL(2,NBREL)
      NELL(3)=NEL(3,NBREL)
      NELL(4)=NEL(4,NBREL)
      call tetra_seg_elem(GNODE,ALEVO,DESNO,SILE,
     &NEL,ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,METOD,NUMPASS,id1,nbrel,4,
     &NELL)
      GOTO 400
      ENDIF
      
c==========================================================================      
      if (nel(5,nbrel).NE.0) THEN
      NELL(1)=NEL(1,NBREL)
      NELL(2)=NEL(2,NBREL)
      NELL(3)=NEL(3,NBREL)
      NELL(4)=NEL(7,NBREL)
      call tetra_seg_elem(GNODE,ALEVO,DESNO,SILE,
     &NEL,ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,METOD,NUMPASS,id1,nbrel,4,
     &NELL)
      NELL(1)=NEL(1,NBREL)
      NELL(2)=NEL(5,NBREL)
      NELL(3)=NEL(6,NBREL)
      NELL(4)=NEL(7,NBREL)
      call tetra_seg_elem(GNODE,ALEVO,DESNO,SILE,
     &NEL,ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,METOD,NUMPASS,id1,nbrel,4,
     &NELL)
      NELL(1)=NEL(2,NBREL)
      NELL(2)=NEL(5,NBREL)
      NELL(3)=NEL(6,NBREL)
      NELL(4)=NEL(7,NBREL)
      call tetra_seg_elem(GNODE,ALEVO,DESNO,SILE,
     &NEL,ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,METOD,NUMPASS,id1,nbrel,4,
     &NELL)
      GOTO 399
      ENDIF

c==========================================================================
      DO 125 I=1,92
      TT210(I)=0.D0
      TT21P(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0

C      tez=0.d0
C=========================================================================
      DO 130 KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
C	tez=tez+0.125d0*CORD(3,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      LM2(KLM+2*NDIM)=ID(3,NEL(KLM,NBREL))
      IF (KLM.LE.NDIMP) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
C      IF (KLM.LE.NDIMP) LM2(KLM+3*NDIM)=ID(4,NEL(ndim+1,NBREL))
      LM2(KLM+3*NDIM+NDIMP)=ID(5,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
c      mat=nel(ndim+1,nbrel)
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
        TT21P(3*NDIM+NR)=ZADVRE(NZDT)*FK1
        ELSE
        TT21P((MESTO-1)*NDIM+NR)=ZADVRE(NZDT)*FK1
        ENDIF
      ENDIF  
 155  CONTINUE
 160  CONTINUE
C=======================================================================
C=========================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
C       if (ipass.eq.0.and.kkorak.eq.1) goto 135
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE

C      IF (KLM.LE.8.AND.ID(4,NEL(KLM,NBREL)).NE.0) THEN
      IF (KLM.LE.NDIMP) THEN
C        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(ndim+1,NBREL))
C        TT210(KLM+3*NDIM)=GNODE(1,4,NEL(ndim+1,NBREL))
        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
        TT21(KLM+3*NDIM+NDIMP)=GNODE(2,5,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM+NDIMP)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE
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
  
       PENXX(K,N)=0.D0
       PENXY(K,N)=0.D0
       PENXZ(K,N)=0.D0
       PENYX(K,N)=0.D0
       PENYY(K,N)=0.D0
       PENYZ(K,N)=0.D0
       PENZX(K,N)=0.D0
       PENZY(K,N)=0.D0
       PENZZ(K,N)=0.D0
      IF ((K.LE.NDIMP).AND.(N.LE.NDIMP)) THEN
       AKPP(K,N)=0.D0
      ENDIF
      IF (N.LE.NDIMP) THEN
       AKV1P(K,N)=0.D0
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
       AKVV1P(K,N)=0.D0
       AKVV2P(K,N)=0.D0
       AKVV3P(K,N)=0.D0
       AKMIV1P(K,N)=0.D0
       AKMIV2P(K,N)=0.D0
       AKMIV3P(K,N)=0.D0
      ENDIF
      AKVV1(K,N)=0.D0
      AKTV1(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKMIV1(K,N)=0.D0
	AKIIX(K,N)=0.D0
	AKIIY(K,N)=0.D0
	AKIIZ(K,N)=0.D0
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
C      AJK(K,N)=0.D0
  162 CONTINUE
  163 CONTINUE


C POVRSINSKE SILE I ZAPREMINSKE SILE:
      DO 199 I=1,NDIM
      RS1(I)=0.D0
      RS2(I)=0.D0
      RS3(I)=0.D0
      RSVP(I)=0.D0
      RB1(I)=0.D0
      RB2(I)=0.D0
      RB3(I)=0.D0
      RP1(I)=0.D0
      RP2(I)=0.D0
      RP3(I)=0.D0
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
 


       CALL JACTNP(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)

      IF (IALE.EQ.1) CALL ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
c      IF (INDAMI.EQ.1) CALL NENJ3D(PJ,TT21,AMI,NDIM,IIZLAZ)
 
      WDT=WT*DET1
      ZAPRE=ZAPRE+WDT
  
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM
      AKVV1(K,N)=AKVV1(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM)
      
c=================================================================
c=================================================================
	AKVV1P(K,N)=AKVV1P(K,N)+WDT*PJ(1,K)*
     1(HV1*PJ(1,N)+HV2*PJ(2,N)+HV3*PJ(3,N))*GUSM
	AKVV2P(K,N)=AKVV2P(K,N)+WDT*PJ(2,K)*
     1(HV1*PJ(1,N)+HV2*PJ(2,N)+HV3*PJ(3,N))*GUSM
	AKVV3P(K,N)=AKVV3P(K,N)+WDT*PJ(3,K)*
     1(HV1*PJ(1,N)+HV2*PJ(2,N)+HV3*PJ(3,N))*GUSM
     
c=================================================================

	AKMIV1P(K,N)=AKMIV1P(K,N)+WDT*PJ(1,K)*(
     1(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
      AKMIV2P(K,N)=AKMIV2P(K,N)+WDT*PJ(2,K)*((PJ(1,K)*PJ(1,N)+
     1+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
      AKMIV3P(K,N)=AKMIV3P(K,N)+WDT*PJ(3,K)*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
c=================================================================
c=================================================================








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
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AKT)
C      AKMIV3(K,N)=AKMIV1(K,N)
      AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM

C      AJV1V1(K,N)=AJV1V1(K,N)+H(K)*HXU*H(N)*GUSM*WDT
C      AJV1V2(K,N)=AJV1V2(K,N)+H(K)*HYU*H(N)*GUSM*WDT
C      AJV1V3(K,N)=AJV1V3(K,N)+H(K)*HZU*H(N)*GUSM*WDT
C      AJV2V1(K,N)=AJV2V1(K,N)+H(K)*HXV*H(N)*GUSM*WDT
C      AJV2V2(K,N)=AJV2V2(K,N)+H(K)*HYV*H(N)*GUSM*WDT
C      AJV2V3(K,N)=AJV2V3(K,N)+H(K)*HZV*H(N)*GUSM*WDT
C      AJV3V1(K,N)=AJV3V1(K,N)+H(K)*HXW*H(N)*GUSM*WDT
C      AJV3V2(K,N)=AJV3V2(K,N)+H(K)*HYW*H(N)*GUSM*WDT
C      AJV3V3(K,N)=AJV3V3(K,N)+H(K)*HZW*H(N)*GUSM*WDT

      if (K.LE.NDIMP.AND.N.LE.NDIMP) THEN
C	AKPP(K,N)=AKPP(K,N)+
C     &WDT*(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))
      ENDIF
       IF (N.LE.NDIMP.AND.PENALT.LT.1.D0) THEN
C	 HP(N)=1.d0
       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
      ENDIF
 164  CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  175 CONTINUE
  180 CONTINUE





         DO I=1,NDIM
	   DO J=1,NDIM
	    AKIIX(I,I)=AKIIX(I,I)+
     &DABS(AKMIV1(I,J))+DABS(AKVV1(I,J))+DABS(AJV1V1(I,J))
	    AKIIY(I,I)=AKIIY(I,I)+
     &DABS(AKMIV1(I,J))+DABS(AKVV1(I,J))+DABS(AJV2V2(I,J))
	    AKIIZ(I,I)=AKIIZ(I,I)+
     &DABS(AKMIV1(I,J))+DABS(AKVV1(I,J))+DABS(AJV3V3(I,J))
	   ENDDO
C	    AKIIX(I,I)=1.D0/AKIIX(I,I)
C	    AKIIY(I,I)=1.D0/AKIIY(I,I)
C	    AKIIZ(I,I)=1.D0/AKIIZ(I,I)
	  ENDDO

      IF (PENALT.LT.1.D0) THEN

       NGAUSX=IBRGT !-1
       NGAUSY=IBRGT !-1
       NGAUSZ=IBRGT !-1
       FPP=0.D0

       DO  NGX=1,NGAUSX
       JR=NREF(NGAUSX) + NGX
       R = XG(JR)
 
       DO  NGY=1,NGAUSY
       JS=NREF(NGAUSY) + NGY
       S = XG(JS)
c 
       DO  NGZ=1,NGAUSZ
       JT=NREF(NGAUSZ) + NGZ
       T = XG(JT)
 
       WT=WGT(JR)*WGT(JS)*WGT(JT)
 
       CALL JACTNP(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,ndimp,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)
 
        WDT=WT*DET1
  

        
        DO  N=1,NDIMP
        DO  K=1,NDIMP
C        AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
C        AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
C        AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))

	  DO  K1=1,NDIM
C
	AKPP(N,K)=AKPP(N,K)+
     &AKV1P(N,K1)*(1.D0/AKIIX(K1,K1))*AKV1P(K,K1)+
     &AKV2P(N,K1)*(1.D0/AKIIY(K1,K1))*AKV2P(K,K1)+
     &AKV3P(N,K1)*(1.D0/AKIIZ(K1,K1))*AKV3P(K,K1)

C	AKPP(N,K)=AKPP(N,K)+WDT*(PJ(1,N)*(1.D0/AKIIX(K,K))*PJ(1,K)
C     &+PJ(2,N)*(1.D0/AKIIY(K,K))*PJ(2,K)
C     &+PJ(3,N)*(1.D0/AKIIZ(K,K))*PJ(3,K))


       
c=================================================================
c	AKVV1P(K,N)=AKVV1P(K,N)+WDT*
c     1(HV1*PJ(1,K)+HV2*PJ(2,K)+HV3*PJ(3,K))*GUSM
c	AKVV2P(K,N)=AKVV2P(K,N)+WDT*
c     1(HV1*PJ(1,K)+HV2*PJ(2,K)+HV3*PJ(3,K))*GUSM
c	AKVV3P(K,N)=AKVV3P(K,N)+WDT*
c     1(HV1*PJ(1,K)+HV2*PJ(2,K)+HV3*PJ(3,K))*GUSM
     
c=================================================================

c	AKMIV1P(K,N)=AKMIV1P(K,N)+WDT*(
c     1(PJ(1,K)*PJ(1,K)+PJ(2,K)*PJ(2,K)+PJ(3,K)*PJ(3,K))*AMI)
c      AKMIV2P(K,N)=AKMIV2P(K,N)+WDT*((PJ(1,K)*PJ(1,K)+
c     1+PJ(2,K)*PJ(2,K)+PJ(3,K)*PJ(3,K))*AMI)
c      AKMIV3P(K,N)=AKMIV3P(K,N)+WDT*((PJ(1,K)*PJ(1,K)+
c     1PJ(2,K)*PJ(2,K)+PJ(3,K)*PJ(3,K))*AMI)
c=================================================================
c=================================================================
       ENDDO
        ENDDO
        ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDIF






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
 
       CALL JACTNP(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)

 
        WDT=WT*DET1
  
        DO  K=1,NDIM
        DO  N=1,NDIM
C      IF (N.LE.NDIMP.AND.PENALT.LT.1.D0) THEN
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
       CALL JACTNP(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)
C
      WDT=WT*DET
      POVRS=POVRS+WDT


C      CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NGPSIL(6,JBRPS)),
C     &NGPSIL(6,JBRPS),NTABFT,IIZLAZ)

C      flux=0.d0
C	xx=0.d0
C       TEMP=DOT(H,TT21(3*NDIM+NDIMP+1),NDIM)
      DO K=1,NDIM
      
	if (indfl.eq.1) then
C	 RS1(K)=RS1(K)+(H(K))*WDT*SF1*FK1
C	 RS2(K)=RS2(K)+(H(K))*WDT*SF2*FK1
C	 RS3(K)=RS3(K)+(H(K))*WDT*SF3*FK1
       FK1=1.D0
	 RS1(K)=RS1(K)+(H(K))*WDT*SF1*FK1
	 RS2(K)=RS2(K)+(H(K))*WDT*SF2*FK1
	 RS3(K)=RS3(K)+(H(K))*WDT*SF3*FK1

	 RSVP(K)=RSVP(K)+(PJ(1,K))*WDT*SF1*FK1
	 RSVP(K)=RSVP(K)+(PJ(2,K))*WDT*SF2*FK1
	 RSVP(K)=RSVP(K)+(PJ(3,K))*WDT*SF3*FK1
C	 RSVP(K)=RSVP(K)+(H(K))*WDT*SF1*FK1
C	 RSVP(K)=RSVP(K)+(H(K))*WDT*SF2*FK1
C	 RSVP(K)=RSVP(K)+(H(K))*WDT*SF3*FK1
c	 povrsila=povrsila+(H(K))*WDT*SF1
c	 povrsila=povrsila+(H(K))*WDT*SF2
c	 povrsila=povrsila+(H(K))*WDT*SF3
	 

      elseif (indfl.eq.0) then
C	AK1=2.D-8
C	VW=4.D-6

C	AK1=2.D-9
C	VW=4.D-7


      	
	VW=FK1
	AK1=VW/1.D2/2.D0

c      flux=flux+DTEMPDN*AKT*WDT
c	xx=xx+cord(3,nel(k,nbrel))
C	 RS1(K)=RS1(K)+H(K)*(1.d-6)*wdt
c	RS1(K)=RS1(K)-1.d-6
c	RS1(K)=RS1(K)-(DTEMPDN*AKT-TEMP*(VW-AK1))*WDT
	RS1(K)=RS1(K)-(-H(K)*TEMP*(VW-AK1))*WDT
	AKK(K,K)=AKK(K,K)-(VW-AK1)*H(K)*WDT
	endif



      ENDDO
C      WRITE(IIZLAZ,*)'SF1=',SF1
C       WRITE(IIZLAZ,*)'NPOV',NPOV,'NBREL=',NBREL
C       CALL IWRRF(NGPSIL(1,JBRPS),5,'NGPSI',IIZLAZ)
C       CALL WRRF(RS3,21,'RS3=',IIZLAZ)


  210 CONTINUE    
C       WRITE(IIZLAZ,*)'POVRSINA',JBRPS,'=',POVRS
c      WRITE(IIZLAZ,*)'flux ',xx,flux
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
      C(I,J)=AMV2(I,J)*CC              

      SKEF(I+3*NDIM+NDIMP,J+3*NDIM+NDIMP)=AF*(AKK(I,J)+AKVV1(I,J)*CC)
      SKEF(I+3*NDIM+NDIMP,J)=AF*AKTV1(I,J)
      SKEF(I+3*NDIM+NDIMP,J+NDIM)=AF*AKTV2(I,J)
      SKEF(I+3*NDIM+NDIMP,J+2*NDIM)=AF*AKTV3(I,J)
      IF (NSTAC.EQ.0) THEN
      SKEF(I+3*NDIM+NDIMP,J+3*NDIM+NDIMP)=
     &SKEF(I+3*NDIM+NDIMP,J+3*NDIM+NDIMP)+AMV2(I,J)*CC/TIME
      ENDIF
 262  CONTINUE
 263  CONTINUE
C=========================================================================
C     CALL NUL(FPOM,21)
C
C     DO I=1,NDIM
C     DO J=1,NDIM
C       FPOM(I)=FPOM(I)+AKMIV1(I,J)*TT21(J)
C      IF (J.LE.NDIMP) FPOM(I)=FPOM(I)+AKV1P(I,J)*TT21(J+3*NDIM)
C     ENDDO
C       FPOM(I)=FPOM(I)-RS1(I)
C     ENDDO
C
C     CALL WRR(FPOM,21,'FPOM=')

C LEVA STRANA:
      DO 270 I=1,NDIM
      DO 265 J=1,NDIM
C       SKEF(I,J)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J)+ATGS1(I,J))
       SKEF(I,J)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J))
       SKEF(I,J+NDIM)=AF*AJV1V2(I,J)
       SKEF(I,J+2*NDIM)=AF*AJV1V3(I,J)
      IF (J.LE.NDIMP) THEN
       SKEF(I,J+3*NDIM)=AF*AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AF*AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AF*AKV3P(I,J)
      ENDIF
      IF (I.LE.NDIMP) THEN
       SKEF(I+3*NDIM,J)=AF*AKV1P(J,I)
       SKEF(I+3*NDIM,J+NDIM)=AF*AKV2P(J,I)
       SKEF(I+3*NDIM,J+2*NDIM)=AF*AKV3P(J,I)
C==========================================================================
C Poisson equation
C==========================================================================
      IF (I.LE.NDIMP.AND.I.LE.NDIMP) THEN
	 IF (IPASS.EQ.0.AND.METOD.EQ.4) THEN
        SKEF(I+3*NDIM,J+3*NDIM)=SKEF(I+3*NDIM,J+3*NDIM)+AKPP(I,J)
	 ENDIF
	 IF (IPASS.EQ.4.AND.METOD.EQ.4) THEN
        SKEF(I+3*NDIM,J+3*NDIM)=SKEF(I+3*NDIM,J+3*NDIM)+AKPP(I,J)
	 ENDIF
	ENDIF
C==========================================================================
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
C        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
C        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21P(JJ)
        F92(II)=F92(II)-SKEF(II,JJ)*TT21P(JJ)
      IF (NSTAC.EQ.0) THEN
        F92(II)=F92(II)-AMV2(I,J)*(TT21(JJ)-TT210(JJ))/TIME
C       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
C     F92(I)=F92(I)-A12(I,J)*TT21(J+NDIM)
 285  CONTINUE							      
 290  CONTINUE
      ENDDO

C=========================================================================
C ovde privremeno za Nikolin primer, 27 March, 2006
      DO 294 I=1,NDIM
      II=I+3*NDIM+NDIMP
c	 F92(II)=F92(II)+RS1(I) 
      DO 292 J=1,NDIM
      JJ=J+3*NDIM+NDIMP
       F92(II)=F92(II)-AKTV1(I,J)*TT21(J)-AKTV2(I,J)*TT21(J+NDIM)
     &-AKTV3(I,J)*TT21(J+2*NDIM)-AKK(I,J)*TT21(JJ)

      IF (NSTAC.EQ.0) THEN
       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
 292  CONTINUE
 294  CONTINUE


C==========================================================================
C==========================================================================
      IF (IPASS.GE.1.AND.IPASS.LE.3.AND.METOD.EQ.4) THEN
      DO I=1,NDIM
	  II=I+NDIM
	  III=I+2*NDIM
       DO J=1,NDIM
	   JJ=J+NDIM
	   JJJ=J+2*NDIM

        SKEF(I,J)=SKEF(I,J)+AKIIX(I,J)*ALPHAU/(1.D0-ALPHAU)
        SKEF(II,JJ)=SKEF(II,JJ)+AKIIY(I,J)*ALPHAV/(1.D0-ALPHAV)
        SKEF(III,JJJ)=SKEF(III,JJJ)+AKIIZ(I,J)*ALPHAW/(1.D0-ALPHAW)
	  F92(I)=F92(I)+AKIIX(I,J)*TT21(J)*ALPHAU/(1.D0-ALPHAU)
	  F92(II)=F92(II)+AKIIY(I,J)*TT21(JJ)*ALPHAV/(1.D0-ALPHAV)
	  F92(III)=F92(III)+AKIIZ(I,J)*TT21(JJJ)*ALPHAW/(1.D0-ALPHAW)

	 ENDDO
	ENDDO
	ENDIF
C==========================================================================
C==========================================================================



	IF (IPASS.NE.0.AND.IPASS.NE.4) THEN
     
      DO 310 I=1,NDIM
      DO 300 J=1,NDIM
      IF (J.LE.NDIMP) THEN
C       F92(I)=F92(I)+AKV1P(J,I)*TT21(J+3*NDIM)
C       F92(I+NDIM)=F92(I+NDIM)+AKV2P(J,I)*TT21(J+3*NDIM)
C       F92(I+2*NDIM)=F92(I+2*NDIM)+AKV3P(J,I)*TT21(J+3*NDIM)
       F92(I)=F92(I)-AKV1P(I,J)*TT21(J+3*NDIM)
       F92(I+NDIM)=F92(I+NDIM)-AKV2P(I,J)*TT21(J+3*NDIM)
       F92(I+2*NDIM)=F92(I+2*NDIM)-AKV3P(I,J)*TT21(J+3*NDIM)
      ENDIF
      IF (I.LE.NDIMP) THEN
	IF (METOD.EQ.1) THEN
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV1P(J,I)*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV2P(J,I)*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV3P(J,I)*TT21(J+2*NDIM)
	ENDIF
      ENDIF
 300  CONTINUE
 310  CONTINUE
      ENDIF
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
       if (indfl.eq.1) then
        F92(I)=F92(I)+RS1(I)
	 elseif (indfl.eq.0) then
        II=I+3*NDIM+NDIMP
        F92(II)=F92(II)+RS1(I)
	 endif
       F92(I+NDIM)=F92(I+NDIM)+RS2(I)
       F92(I+2*NDIM)=F92(I+2*NDIM)+RS3(I)
320   CONTINUE
C==========================================================================


C==========================================================================
C Poisson equation
C==========================================================================
	IF (IPASS.EQ.0.AND.METOD.EQ.4) THEN
  


	DO J=1,NDIM
	 DO K=1,NDIM
	  RP1(J)=RP1(J)-(AKVV1(J,K)+AKMIV1(J,K))*TT21(K)
 	  RP2(J)=RP2(J)-(AKVV1(J,K)+AKMIV1(J,K))*TT21(K+NDIM)
 	  RP3(J)=RP3(J)-(AKVV1(J,K)+AKMIV1(J,K))*TT21(K+2*NDIM)
	 ENDDO
	ENDDO

      DO I=1,NDIMP
	II=I+3*NDIM
	DO J=1,NDIM
        F92(II)=F92(II)-AKV1P(J,I)*(1.D0/AKIIX(J,J))*RP1(J)
        F92(II)=F92(II)-AKV2P(J,I)*(1.D0/AKIIY(J,J))*RP2(J)
        F92(II)=F92(II)-AKV3P(J,I)*(1.D0/AKIIZ(J,J))*RP3(J)
	ENDDO
	ENDDO

      

      DO 311 I=1,NDIMP
      DO 301 J=1,NDIM
C       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV1P(J,I)*F92(J)/WDT
C       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV2P(J,I)*F92(J+NDIM)/WDT
C       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV3P(J,I)*F92(J+2*NDIM)/WDT

C       F92(I+3*NDIM)=F92(I+3*NDIM)-(AKVV1P(I,J)+AKMIV1P(I,J))*TT21(J)
C       F92(I+3*NDIM)=F92(I+3*NDIM)-(AKVV2P(I,J)+AKMIV2P(I,J))
C     &*TT21(J+NDIM)
C       F92(I+3*NDIM)=F92(I+3*NDIM)-(AKVV3P(I,J)+AKMIV3P(I,J))
C     &*TT21(J+2*NDIM)

C       ALAMBDA=0.D0
C	 IF (NUMPASS.EQ.4) ALAMBDA=1.D0
C       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV1P(J,I)*TT21(J)*ALAMBDA
C       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV2P(J,I)*TT21(J+NDIM)*ALAMBDA
C       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV3P(J,I)*TT21(J+2*NDIM)*ALAMBDA

       IF (I.LE.NDIMP.AND.J.LE.NDIMP) THEN
        F92(I+3*NDIM)=F92(I+3*NDIM)-AKPP(I,J)*TT21P(J+3*NDIM)
       ENDIF
C==========================================================================
 301  CONTINUE
c       F92(I+3*NDIM)=F92(I+3*NDIM)+RSVP(I)
 311  CONTINUE
      ENDIF
C==========================================================================

	IF (IPASS.EQ.4.AND.METOD.EQ.4) THEN
  
      DO I=1,NDIMP
      DO J=1,NDIM
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV1P(J,I)*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV2P(J,I)*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV3P(J,I)*TT21(J+2*NDIM)
c      F92(I+3*NDIM)=F92(I+3*NDIM)-
c     &AKPP(I,J)*(TT21(J+3*NDIM)-TT210(J+3*NDIM))

C==========================================================================
      ENDDO
      ENDDO
      ENDIF

	IF (IPASS.EQ.5.AND.METOD.EQ.4) THEN
  
      DO I=1,NDIM
	NODE=NEL(I,NBREL)
      DO J=1,NDIMP
C	DO K=1,NDIM
	IF (LM2(I).NE.0) THEN
       GNODE(2,1,NODE)=GNODE(2,1,NODE)+
     &(AKV1P(J,I))*(1.D0/AKIIX(I,I))*TT21(J+4*NDIM)
	ENDIF
	IF (LM2(I+NDIM).NE.0) THEN
       GNODE(2,2,NODE)=GNODE(2,2,NODE)+
     &(AKV2P(J,I))*(1.D0/AKIIY(I,I))*TT21(J+4*NDIM)
	ENDIF
	IF (LM2(I+2*NDIM).NE.0) THEN
       GNODE(2,3,NODE)=GNODE(2,3,NODE)+
     &(AKV3P(J,I))*(1.D0/AKIIZ(I,I))*TT21(J+4*NDIM)
	ENDIF
     
C==========================================================================
C      ENDDO
      ENDDO
      ENDDO
      ENDIF
C==========================================================================



C      CALL WRRF(F92,NDES,'F92= ',IIZLAZ)      
C      CALL WRRF(SKEF,NDES*NDES,'SKEF=',IIZLAZ)      
      
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
c	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,NETIP)
#if(MUMPS_CLUSTER)
	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,5)
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
c      if (iter.gt.0) then
c      DO I=1,NDIM
c        I1=I
c        I2=I+NDIM
c        I3=I+2*NDIM
c        PPP=PPP+
c     &(PJ(1,I)*TT21(I1)+PJ(2,I)*TT21(I2)+PJ(3,I)*TT21(I3))
c       ENDDO
c      WRITE(IIZLAZ,*)'PPP= ',PPP,nbrel
c	endif
399   CALL SSTRES8(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)
 400  CONTINUE

C End of subroutine
cc      CALL MUMPSRIGHT(SILE,JEDN)
C      CALL WRRF(SILE,JEDN,'desno=',IIZLAZ)
c      if (iter.gt.0) stop
c      write(iizlaz,*)'povrs= ',povrs
c      write(iizlaz,*)'povrsila= ',povrsila
c      write(iizlaz,*)'zapre= ',zapre
c	stop
      END
C=======================================================================
C=======================================================================
C==========================================================================
      SUBROUTINE NODEADD(CORD,NPT,ID,NEL,NDIM,NET,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      DIMENSION ID(6,*)
      DIMENSION NEL(NDIM+1,*)

C
C=========================================================================
      OPEN(45,FILE='INPUT.DAT')
      
	WRITE(IIZLAZ,*)'CVOROVI'
      NNEL=NPT
      DO 100 I=1,101
      
	DO J=(I-1)*116+1,(I-1)*116+116
	 IDX=1
	 IF (ID(1,J).NE.0) IDX=0
	 IDY=1
	 IF (ID(2,J).NE.0) IDY=0
	 IDZ=1
	 IF (ID(3,J).NE.0) IDZ=0
       WRITE(IIZLAZ,'(I10,5I5,3E13.5)')J,IDX,IDY,IDZ,1,1,
     &(CORD(II,J),II=1,3)
	ENDDO
	
	IF (I.EQ.101) GOTO 100
	  DO K=1,100
         WRITE(IIZLAZ,'(I10,5I5,3E13.5)')NNEL+K,1,1,1,0,1,0.D0,0.D0,0.D0
	  ENDDO
	  NNEL=NNEL+100
100   CONTINUE      
    	
	WRITE(IIZLAZ,*)'ELEMENTS'
	NNEL=NPT
	DO II=1,NET
	  NNEL=NNEL+1
	  WRITE(IIZLAZ,'(10I10)')II,(NEL(J,II),J=1,8),NNEL
	ENDDO
	
	

	STOP
C=========================================================================

      END
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
       SUBROUTINE JACTNP(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION H(21),P(3,21),NEL(NDIM+1,*)
       DIMENSION XJ(3,3),LE(3),ME(3),CK(21,3),PJ(3,21),HP(8)
       DIMENSION TT21(*),V1(21),V2(21),V3(21),PP(8),XJJ(3,3),TEMP(21)
       DIMENSION VTANG(3,21),AN(3),VSHX(3),VSHY(3),VSHZ(3),SHEAR(*)
       DIMENSION IPERM(8),NOD9(13)
C
CE  Subroutine JACTNP is used for integration 3D finite elements, 8 node - 1 pressure element
C
C       COMMON /VISKOZ/ AMI,INDAMI
C       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C       COMMON /SRPSKI/ ISRPS
C       COMMON /ULAZNI/ IULAZ,IIZLAZ
C       COMMON /AJVVV/ HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW
C       COMMON /UPWIND/ IUPWIN
	COMMON /TRANSP/ INDFL,LID1

       DATA IPERM/2,3,4,1,6,7,8,5/
       DATA NOD9/9,10,11,12,13,14,15,16,17,18,19,20,21/


       IELX=NDIM
       INDUP=0
       NND9=13
       CALL CLEAR(SHEAR,3)

      RP=1.0+R
      SP=1.0+S
      TP=1.0+T
      RM=1.0-R
      TM=1.0-T
      SM=1.0-S
      RR=1.0-R*R
      SS=1.0-S*S
      TT=1.0-T*T
      DO 82 I=1,21
      H(I)=0.D0
      DO 82 J=1,3
   82 P(J,I)=0.D0
C
C     INTERPOLACIJSKE FINKCIJE I NJIHOVI IZVODI
CE    INTERPOLATION FUNCTIONS AND DERIVATIVES
C
C
CS    PRVIH 8 CVOROVA
CE    FIRST 4 NODES
C
      H(1)=0.125*RP*SP*TP
      H(2)=0.125*RM*SP*TP
      H(3)=0.125*RM*SM*TP
      H(4)=0.125*RP*SM*TP
      H(5)=0.125*RP*SP*TM
      H(6)=0.125*RM*SP*TM
      H(7)=0.125*RM*SM*TM
      H(8)=0.125*RP*SM*TM
      
100   DO JJ=1,8
       HP(JJ)=H(JJ)
      ENDDO

	if (ndimp.eq.1) hp(1)=1.d0
C
      P(1,1)=0.125*SP*TP
      P(1,2)=-P(1,1)
      P(1,3)=-0.125*SM*TP
      P(1,4)=-P(1,3)
      P(1,5)=0.125*SP*TM
      P(1,6)=-P(1,5)
      P(1,7)=-0.125*SM*TM
      P(1,8)=-P(1,7)
C
      P(2,1)=0.125*RP*TP
      P(2,2)=0.125*RM*TP
      P(2,3)=-P(2,2)
      P(2,4)=-P(2,1)
      P(2,5)=0.125*RP*TM
      P(2,6)=0.125*RM*TM
      P(2,7)=-P(2,6)
      P(2,8)=-P(2,5)
C
      P(3,1)=0.125*RP*SP
      P(3,2)=0.125*RM*SP
      P(3,3)=0.125*RM*SM
      P(3,4)=0.125*RP*SM
      P(3,5)=-P(3,1)
      P(3,6)=-P(3,2)
      P(3,7)=-P(3,3)
      P(3,8)=-P(3,4)
C
      IF (IELX.EQ.8) GO TO 50
C
CS    STEPENE SLOBODE ZA CVOROVE PREKO 8
CE     DEGREES OF FREADOM FOR NODES OVER 8
C
      I=0
    2 I=I+1
      IF (I.GT.NND9) GO TO 30
      NN=NOD9(I)-8
      GO TO (9,10,11,12,13,14,15,16,17,18,19,20,21),NN
C
    9 H(9)=0.25*RR*SP*TP
      P(1,9)=-0.50*R*SP*TP
      P(2,9)=0.25*RR*TP
      P(3,9)=0.25*RR*SP
      GO TO 2
   10 H(10)=0.25*RM*SS*TP
      P(1,10)=-0.25*SS*TP
      P(2,10)=-0.50*RM*S*TP
      P(3,10)=0.25*RM*SS
      GO TO 2
   11 H(11)=0.25*RR*SM*TP
      P(1,11)=-0.50*R*SM*TP
      P(2,11)=-0.25*RR*TP
      P(3,11)=0.25*RR*SM
      GO TO 2
   12 H(12)=0.25*RP*SS*TP
      P(1,12)=0.25*SS*TP
      P(2,12)=-0.50*RP*S*TP
      P(3,12)=0.25*RP*SS
      GO TO 2
   13 H(13)=0.25*RR*SP*TM
      P(1,13)=-0.50*R*SP*TM
      P(2,13)=0.25*RR*TM
      P(3,13)=-0.25*RR*SP
      GO TO 2
   14 H(14)=0.25*RM*SS*TM
      P(1,14)=-0.25*SS*TM
      P(2,14)=-0.50*RM*S*TM
      P(3,14)=-0.25*RM*SS
      GO TO 2
   15 H(15)=0.25*RR*SM*TM
      P(1,15)=-0.50*R*SM*TM
      P(2,15)=-0.25*RR*TM
      P(3,15)=-0.25*RR*SM
      GO TO 2
   16 H(16)=0.25*RP*SS*TM
      P(1,16)=0.25*SS*TM
      P(2,16)=-0.50*RP*S*TM
      P(3,16)=-0.25*RP*SS
      GO TO 2
   17 H(17)=0.25*RP*SP*TT
      P(1,17)=0.25*SP*TT
      P(2,17)=0.25*RP*TT
      P(3,17)=-0.50*RP*SP*T
      GO TO 2
   18 H(18)=0.25*RM*SP*TT
      P(1,18)=-0.25*SP*TT
      P(2,18)=0.25*RM*TT
      P(3,18)=-0.50*RM*SP*T
      GO TO 2
   19 H(19)=0.25*RM*SM*TT
      P(1,19)=-0.25*SM*TT
      P(2,19)=-0.25*RM*TT
      P(3,19)=-0.50*RM*SM*T
      GO TO 2
   20 H(20)=0.25*RP*SM*TT
      P(1,20)=0.25*SM*TT
      P(2,20)=-0.25*RP*TT
      P(3,20)=-0.50*RP*SM*T
      GO TO 2
   21 H(21)=RR*SS*TT
      P(1,21)=-2.0*R*SS*TT
      P(2,21)=-2.0*S*RR*TT
      P(3,21)=-2.0*T*RR*SS
      GO TO 2
C
CS    KOREKCIJE PRVIH 20 FINKCIJA AKO JE UPOTREBLJEN CVOR 21
CE    CORECTION OF FIRST 20 FUNCTIONS IF NODE 21 EXISTS
C
   30 IN=NOD9(NND9)
      IF(IN.NE.21) GO TO 40
      DO 36 I=1,8
      H(I)=H(I)-0.125*H(21)
      DO 36 J=1,3
   36 P(J,I)=P(J,I)-0.125*P(J,21)
      IF(NND9.EQ.1) GO TO 51
      DO 37 I=1,NND9-1
      IN=NOD9(I)
      H(IN)=H(IN)-0.25*H(21)
      DO 37 J=1,3
   37 P(J,IN)=P(J,IN)-0.25*P(J,21)
C
CS    KOREKCIJE PRVIH 8 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 8
CE    CORECTION OF FIRST 8 FUNCTIONS IF NODES OVER 8 EXISTS
C
   40 IH=0
   41 IH=IH+1
      IF(IH.GT.NND9) GO TO 50
      IN=NOD9(IH)
      IF(IN.GT.16) GO TO 46
      I1=IN-8
      I2=IPERM(I1)
C
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 45 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   45 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   46 IF(IN.EQ.21) GO TO 51
      I1=IN-16
      I2=I1+4
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 47 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   47 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   51 H(NND9+8)=H(21)
      DO 39 J=1,3
   39 P(J,NND9+8)=P(J,21)

      
   50 HH=0.D0
      DO I=1,NDIM
       HH=HH+H(I)
      ENDDO
C      WRITE(IIZLAZ,*)'ZBIR H=',HH

      DO I=1,NDIM
      V1(I)=TT21(I)
      V2(I)=TT21(I+NDIM)
      V3(I)=TT21(I+2*NDIM)
      TEMP(I)=TT21(I+3*NDIM+8)
      IF (I.LE.NDIMP) THEN
        PP(I)=TT21(I+3*NDIM)
      ENDIF
      ENDDO
      
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C

      DO I=1,3
      DO J=1,3
       XJ(I,J)=0.D0
       XJJ(I,J)=0.D0
        DO  KK=1,NDIM
         XJ(I,J)=XJ(I,J)+P(I,KK)*CK(KK,J)
         XJJ(I,J)=XJJ(I,J)+P(I,KK)*CK(KK,J)
        ENDDO    
C      WRITE(IIZLAZ,*)'XJ=',XJ(I,J)
      ENDDO    
      ENDDO    

      HV1=DOT(H,V1,NDIM)
      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)
  
      





      CALL MINV(XJJ,3,DET1,LE,ME)

c      if (dabs(r).lt.1.d-12.and.dabs(s).lt.1.d-12.and.
c     &dabs(t).lt.1.d-12) then
c      if (det1.gt.1.d-15) then
c	  write(iizlaz,'(10i5)') nbrel,(nel(i,nbrel),i=1,ndim),1
c      else
c	  write(iizlaz,'(10i5)') nbrel,
c     &nel(5,nbrel),nel(6,nbrel),nel(7,nbrel),nel(8,nbrel),
c     &nel(1,nbrel),nel(2,nbrel),nel(3,nbrel),nel(4,nbrel),1
c 	endif
c	endif



      IF (DET1.LT.1.D-15) THEN
c       WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
c       WRITE(IIZLAZ,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
c       WRITE(IIZLAZ,*)'DETERMINANTE= ',DET1
c       WRITE(IIZLAZ,*)'NODES  COORDINATES'
c       DO I=1,NDIM
c        WRITE(IIZLAZ,1000) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
c        WRITE(IIZLAZ,1001) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
c	  ENDDO
c       return 
      ENDIF
 1000 FORMAT(I5,3(D13.5))
 1001 FORMAT(I5,3(f10.6))

      DO 85 I=1,3
      DO 85 JJ=1,NDIM
      PJ(I,JJ)=0.D0
      DO 85 K=1,3
      PJ(I,JJ)=PJ(I,JJ) + XJJ(I,K)*P(K,JJ)
   85 CONTINUE


      IF (IUPWIN.EQ.1.AND.INDUP.EQ.0) THEN
c       CALL INTER1(CK,V1,V2,V3,H,PJ,NDIM,AMI)
c      	AK=2.D0/AMI 
      	AK=2.D0/2.d0
       VV=sqrt(HV1**2+HV2**2+HV3**2)
	if (dabs(vv).gt.1.d-15) then
      DO I=1,NDIM
c        PP=AKK*(HV1*PJ(1,I)+HV2*PJ(2,I)+HV3*PJ(3,I))/VV
        P1=AK*(HV1*PJ(1,I)+HV2*PJ(2,I)+HV3*PJ(3,I))/VV
        Hp(I)=Hp(I)+P1
      ENDDO
	endif

       INDUP=1
c      GOTO 100
      ENDIF
      

      HXU=0.D0
      HYU=0.D0
      HZU=0.D0
      HXV=0.D0
      HYV=0.D0
      HZV=0.D0
      HXW=0.D0
      HYW=0.D0
      HZW=0.D0
      ZVXT=0.D0
      ZVYT=0.D0
      ZVZT=0.D0
      DO L =1,NDIM
        HXU=HXU+PJ(1,L)*V1(L)
        HYU=HYU+PJ(2,L)*V1(L)
        HZU=HZU+PJ(3,L)*V1(L)
        HXV=HXV+PJ(1,L)*V2(L)
        HYV=HYV+PJ(2,L)*V2(L)
        HZV=HZV+PJ(3,L)*V2(L)
        HXW=HXW+PJ(1,L)*V3(L)
        HYW=HYW+PJ(2,L)*V3(L)
        HZW=HZW+PJ(3,L)*V3(L)
        ZVXT=ZVXT+PJ(1,L)*TEMP(L)
        ZVYT=ZVYT+PJ(2,L)*TEMP(L)
        ZVZT=ZVZT+PJ(3,L)*TEMP(L)
      ENDDO




       

      IF(KFIX.GT.0) GO TO 70
      RETURN



C
CS     DETERMINATA POVRSINSKOG JAKOBIJANA
CE     SURFACE JACOBIAN DETERMINANT
C
  70   GO TO (71,72,73),KFIX
CS     KONSTANTNO KSI
CE     CONSTANT KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      ANX=R*(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))
      ANY=R*(XJ(3,1)*XJ(2,3)-XJ(3,3)*XJ(2,1))
      ANZ=R*(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ETA
CE     CONSTANT ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      ANX=-S*(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))
      ANY=S*(XJ(1,1)*XJ(3,3)-XJ(1,3)*XJ(3,1))
      ANZ=-S*(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ZETA
CE     CONSTANT ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
      ANX=T*(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))
      ANY=-T*(XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1))
      ANZ=T*(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))
   74 DET=DSQRT(DET)
      ANX=ANX/DET
      ANY=ANY/DET
      ANZ=ANZ/DET
      AN(1)=ANX
      AN(2)=ANY
      AN(3)=ANZ
C      WRITE(IIZLAZ,*)'ANX=',ANX
C      WRITE(IIZLAZ,*)'ANY=',ANY
C      WRITE(IIZLAZ,*)'ANZ=',ANZ
      CALL SHEARS(AN,TT21,VTANG,NDIM)
      
      V1X=0.D0
      V1Y=0.D0
      V1Z=0.D0
      V2X=0.D0
      V2Y=0.D0
      V2Z=0.D0
      V3X=0.D0
      V3Y=0.D0
      V3Z=0.D0
	DTDX=0.D0
	DTDY=0.D0
	DTDZ=0.D0
      VSHX(1)=0.D0
      VSHX(2)=0.D0
      VSHX(3)=0.D0
      VSHY(1)=0.D0
      VSHY(2)=0.D0
      VSHY(3)=0.D0
      VSHZ(1)=0.D0
      VSHZ(2)=0.D0
      VSHZ(3)=0.D0
      DO I=1,NDIM
       V1X=V1X+PJ(1,I)*V1(I)
       V1Y=V1Y+PJ(2,I)*V1(I)
       V1Z=V1Z+PJ(3,I)*V1(I)
       V2X=V2X+PJ(1,I)*V2(I)
       V2Y=V2Y+PJ(2,I)*V2(I)
       V2Z=V2Z+PJ(3,I)*V2(I)
       V3X=V3X+PJ(1,I)*V3(I)
       V3Y=V3Y+PJ(2,I)*V3(I)
       V3Z=V3Z+PJ(3,I)*V3(I)
       VSHX(1)=VSHX(1)+PJ(1,I)*VTANG(1,I)
       VSHX(2)=VSHX(2)+PJ(2,I)*VTANG(1,I)
       VSHX(3)=VSHX(3)+PJ(3,I)*VTANG(1,I)
       VSHY(1)=VSHY(1)+PJ(1,I)*VTANG(2,I)
       VSHY(2)=VSHY(2)+PJ(2,I)*VTANG(2,I)
       VSHY(3)=VSHY(3)+PJ(3,I)*VTANG(2,I)
       VSHZ(1)=VSHZ(1)+PJ(1,I)*VTANG(3,I)
       VSHZ(2)=VSHZ(2)+PJ(2,I)*VTANG(3,I)
       VSHZ(3)=VSHZ(3)+PJ(3,I)*VTANG(3,I)
	 DTDX=DTDX+PJ(1,I)*TEMP(I)
	 DTDY=DTDY+PJ(2,I)*TEMP(I)
	 DTDZ=DTDZ+PJ(3,I)*TEMP(I)
      ENDDO
      SHEAR(1)=DOT(AN,VSHX,3)
      SHEAR(2)=DOT(AN,VSHY,3)
      SHEAR(3)=DOT(AN,VSHZ,3)
	

      PRIT=DOT(HP,PP,NDIMP)
      
	if (indfl.eq.1) then
C       SF1=ANX
C       SF2=ANY
C       SF3=ANZ
C	else if (indfl.eq.0) then
       SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
       SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
       SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
	endif
	
	
	DTEMPDN=DTDX*ANX+DTDY*ANY+DTDZ*ANZ



      IF ( DET.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000) NBREL,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     *WRITE(3,6000) NBREL,KFIX,R,S,T,DET
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
C

      END
C=======================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE SSTRES8(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NEL(NDIM+1,*),ID(6,*),INDEL(*),NZAD(3,*)
      DIMENSION CK(21,*),PJ(3,*),H(*),HP(*),TT21(*),PRES(3,*),ZADVRE(*)
	DIMENSION NID(21),N(21),NWALL(6),ITR(6,4)
      DIMENSION XG(*),WGT(*),NREF(*),SHEAR(3)

C
CE Subroutine SSTRES is used for calculation shear stresses
C

      NGAUSX=IBRGT
      NGAUSY=IBRGT
      NGAUSZ=IBRGT
      ndimp=8


C      WRITE(IIZLAZ,*)'ELEMENT= ',NBREL
      
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
      if (n(3).eq.n(4).and.n(7).eq.n(8).and.npov.eq.4) goto 300
       IF (NWALL(NPOV).EQ.0) GOTO 300

C
CS  PETLJA PO GAUSOVIM TACKAMA
CE  GAUSS POINTS LOOP
C
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
   60 KK=0
      DO 210 NGX=1,NGXP
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

      KK=KK+1
      NODE=ITR(NPOV,KK)

C      do jjj=1,4
C       NODE=ITR(NPOV,jjj)
C      DO JJ=1,NUMZAD
C       IF(NODE.EQ.NZAD(1,JJ)) GOTO 300
C      ENDDO
C      enddo
       
C	   IF(NODE.le.36) GOTO 300

C      DO JJ=1,NUMZAD
CC       IF(NODE.EQ.NZAD(1,JJ).AND.DABS(ZADVRE(JJ)).LT.1.D-12) GOTO 150
C       IF(NODE.EQ.NZAD(1,JJ)) GOTO 300
C      ENDDO
C      GOTO 300

 150   CALL JACTNP(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &ndimp)
c      WRITE(IIZLAZ,*)'SHEAR STRESS= ',SHEAR(1),SHEAR(2),SHEAR(3)
c	write(iizlaz,*)'indel(',node,')= ',indel(node)
c	write(iizlaz,*)'npov= ',npov
c	write(iizlaz,*)'nbrel= ',nbrel

C      write(iizlaz,*)'Pres(2,node)= ',Pres(2,node)
C      if (mod(node,121).eq.111) indel(node)=2
C      if (mod(node,441).eq.421) indel(node)=2

	PRES(1,NODE)=PRES(1,NODE)-AMI*SHEAR(1)/INDEL(NODE)
      PRES(2,NODE)=PRES(2,NODE)-AMI*SHEAR(2)/INDEL(NODE)
      PRES(3,NODE)=PRES(3,NODE)-AMI*SHEAR(3)/INDEL(NODE)

C      write(iizlaz,*)'node, shear (1)',node, shear(1)

C	PRES(1,NODE)=1.d0*INDEL(NODE)
C	PRES(2,NODE)=1.d0*INDEL(NODE)
C	PRES(3,NODE)=1.d0*INDEL(NODE)
  
C      INDEL(NODE)=INDEL(NODE)+1
      

  210 CONTINUE    

  300 CONTINUE    

      END
C==========================================================================
C==========================================================================
      SUBROUTINE FILLN2P(GNODE,SILE,ID,NPT,IPASS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(2,6,*),SILE(*),ID(6,*)

	COMMON /ALPHA_SEG/ ALPHAU,ALPHAV,ALPHAW,ALPHAP


c      ALPHA=2.D0/3.D0
C
      ALPHA=0.0D0
c      ALPHAPP=0.999D0
      
c      ALPHAPP=0.0D0

C      ALPHAP=0.999D0
      ALPHAPP=ALPHAP

c      ALPHAP=0.0D0
      
c      ALPHA=0.999D0
C      ALPHAP=ALPHA

      ALPHAU=ALPHA
      ALPHAV=ALPHA
      ALPHAW=ALPHA


CE Subroutine FILLN is used for re-writing values at the end of time step
C
       DO NODE=1,NPT 
        DO I=1,5
C        DO I=1,6
c This should be change for acoustics problem
         JJ=ID(I,NODE)

         
	   IF (JJ.NE.0) THEN
          IF (IPASS.NE.4) THEN
	      IF (IPASS.EQ.0) THEN
        GNODE(2,I,NODE)=ALPHAP*GNODE(2,I,NODE)+(1.D0-ALPHAP)*SILE(JJ)
c		   GNODE(2,I,NODE)=GNODE(2,I,NODE)+SILE(JJ)
c		   GNODE(2,I,NODE)=SILE(JJ)
	      ELSE
C		   GNODE(2,I,NODE)=GNODE(2,I,NODE)+SILE(JJ)
		   GNODE(2,I,NODE)=SILE(JJ)
	      ENDIF
	    ELSE IF (IPASS.EQ.4) THEN
C		   GNODE(2,5,NODE)=SILE(JJ)
		   GNODE(2,5,NODE)=(1.D0-ALPHAPP)*SILE(JJ)
		   GNODE(2,4,NODE)=GNODE(2,4,NODE)+(1.D0-ALPHAPP)*SILE(JJ)
C		   GNODE(2,4,NODE)=GNODE(2,4,NODE)+SILE(JJ)
C           GNODE(2,I,NODE)=ALPHAP*GNODE(2,I,NODE)+(1.D0-ALPHAP)*SILE(JJ)
C		   GNODE(2,4,NODE)=SILE(JJ)
   
C		   GNODE(2,I,NODE)=GNODE(2,I,NODE)+SILE(JJ)
	    ENDIF
	   ENDIF


        ENDDO
       ENDDO
     
 
            
      END
C==========================================================================
