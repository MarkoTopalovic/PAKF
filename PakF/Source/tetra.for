C=========================================================================
      SUBROUTINE RACU3Dtethra(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,METOD,NUMPASS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine RACU3D is used for 3D analysis
CE It is used global loop per elements
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21),TT21A(3*21)
      DIMENSION H(21),HP(15),hpv(15),PJP(3,21)
      DIMENSION AKVV1(21,21),AKMIV1(21,21),AMV2(21,21),AJV1V1(21,21)
      DIMENSION AJV1V2(21,21),AJV1V3(21,21),AJV2V1(21,21),AJV2V2(21,21)
      DIMENSION AJV2V3(21,21),AJV3V1(21,21),AJV3V2(21,21),AJV3V3(21,21)
      DIMENSION AKTV1(21,21),AKTV2(21,21),AKTV3(21,21),AKK(21,21)
      DIMENSION AKV1P(21,15),AKV2P(21,15),AKV3P(21,15),RS1(21),RS2(21)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(15,15),PENXY(15,15),PENXZ(15,15)
      DIMENSION PENYX(15,15),PENYY(15,15),PENYZ(15,15)
      DIMENSION PENZX(15,15),PENZY(15,15),PENZZ(15,15)
	DIMENSION AKPP(15,15),AKVV1P(15,15),AKVV2P(15,15),AKVV3P(15,15)
	DIMENSION AKMIV1P(15,15),AKMIV2P(15,15),AKMIV3P(15,15)
	DIMENSION AKIIX(15,15),AKIIY(15,15),AKIIZ(15,15)
      DIMENSION C(21,21)
	DIMENSION SKE(46*93)
c	DIMENSION ATGS1(8,8),ATGS2(8,8),ATGS3(8,8)
C	DIMENSION SKE(12*25)
      dimension nizel(4)



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



c       if (ipass.gt.0) penalt=1.d9
c       if (ndim.eq.8) ipass=-1

C=========================================================================
      DO NODE=1,NPT
       DO NZDT=1,NUMZAD
        IF(NODE.EQ.NZAD(1,NZDT)) THEN
         CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
     &NTABFT,IIZLAZ)
         MESTO=NZAD(2,NZDT)
	   GNODE(2,MESTO,NODE)=ZADVRE(NZDT)*FK1
c	   GNODE(1,MESTO,NODE)=ZADVRE(NZDT)*FK1
        ENDIF  
	 ENDDO
	ENDDO
c=======================================================================

c      R0=1.d-3
c      do i=1,npt
c	rr=dsqrt(cord(1,i)**2+cord(2,i)**2)
c	v=1.d0*(1.d0-(rr/R0)**2)
c	if(id(1,i).ne.0.and.id(3,i).ne.0.and.
c     &dabs(cord(3,i)).lt.1.d-7) then
c	  write(iizlaz,'(3i5,f10.5)')i,3,1,v
c	endif
c	enddo
c     stop


      if (ndim.eq.10) ndimp=4
      if (ndim.eq.11) ndimp=11
c      if (ndim.eq.4) ndimp=1
      if (ndim.eq.4) ndimp=4
      if (ndim.eq.5) ndimp=4
      if (ndim.eq.8) ndimp=4

      PEN1=1.D0
C      if (ndim.eq.4) PEN1=1.D9


c      call INDELSSTRES(NEL,INDEL,NET,NPT,NDIM,ID)


      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
c      if (kkorak.eq.0.and. iter.eq.0) 
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)
c	REWIND(MUFILE2)
c        CALL sparseassembler_init(0)
	ENDIF

      povrs=0.0d0
      povrsila=0.0d0
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
      if (KLM.LE.ndimp) then
	if (nel(ndim+1,nbrel).ne.0) then
	   LM2(klm+3*NDIM)=ID(4,NEL(ndim+1,NBREL))
      else
	   LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
	endif
	endif
c      IF (KLM.LE.1) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
      LM2(KLM+3*NDIM+ndimp)=ID(5,NEL(KLM,NBREL))
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
C=========================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
C      IF (KLM.LE.8.AND.ID(4,NEL(KLM,NBREL)).NE.0) THEN
      IF (KLM.LE.ndimp) THEN
	if (nel(ndim+1,nbrel).ne.0) then
        TT21(klm+3*NDIM)=GNODE(2,4,NEL(ndim+1,NBREL))
        TT210(klm+3*NDIM)=GNODE(1,4,NEL(ndim+1,NBREL))
	else
        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
	endif
      ENDIF
        TT21(KLM+3*NDIM+ndimp)=GNODE(2,5,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM+ndimp)=GNODE(1,5,NEL(KLM,NBREL))
        
 140  CONTINUE
      
c        tt21(3*ndim+1)=2.d2-tez
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
  
c      IF ((K.LE.8).AND.(N.LE.8)) THEN
       PENXX(K,N)=0.D0
       PENXY(K,N)=0.D0
       PENXZ(K,N)=0.D0
       PENYX(K,N)=0.D0
       PENYY(K,N)=0.D0
       PENYZ(K,N)=0.D0
       PENZX(K,N)=0.D0
       PENZY(K,N)=0.D0
       PENZZ(K,N)=0.D0
c      ENDIF
      IF (N.LE.ndimp) THEN
       AKV1P(K,N)=0.D0
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
      ENDIF
       AKPP(K,N)=0.D0
       AKVV1P(K,N)=0.D0
       AKVV2P(K,N)=0.D0
       AKVV3P(K,N)=0.D0
       AKMIV1P(K,N)=0.D0
       AKMIV2P(K,N)=0.D0
       AKMIV3P(K,N)=0.D0
      AKIIX(K,N)=0.D0
      AKIIY(K,N)=0.D0
      AKIIZ(K,N)=0.D0


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


c      alfa=0.5854102966249685
c	beta=0.138196601125015
c	wt=0.25



      nkk=4
      if (ndim.eq.10)nkk=4
      if (ndim.eq.11)nkk=4
      if (ndim.eq.5) nkk=5
      if (ndim.eq.8) nkk=4
      if (ndim.eq.4) nkk=4


	do 180 kk=1,nkk

      if (nkk.eq.4) then
        call Integ4points(r,s,t,wt,kk)      
	elseif (nkk.eq.5) then
        call Integ5points(r,s,t,wt,kk)      
	elseif (nkk.eq.11) then
        call Integ10points(r,s,t,wt,kk)      
	elseif (nkk.eq.1) then
 	  wt=1.0d0
	  r=0.25d0
	  s=0.25d0
	  t=0.25d0
	endif
      


c      DO 180 NGX=1,NGAUSX
c      JR=NREF(NGAUSX) + NGX
c      R = XG(JR)
 
c      DO 175 NGY=1,NGAUSY
c      JS=NREF(NGAUSY) + NGY
c      S = XG(JS)
 
c      DO 170 NGZ=1,NGAUSZ
c      JT=NREF(NGAUSZ) + NGZ
c      T = XG(JT)
 
c      WT=WGT(JR)*WGT(JS)*WGT(JT)
     

       
	if (ndim.eq.4) then
       call jact4node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,ndimp,hpv)
c	 if (penalt.lt.1.d0) hp(1)=1.d0
	elseif(ndim.eq.10) then
       call jact10node1 (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL)
	elseif(ndim.eq.11) then
       call jact11node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL)
	elseif(ndim.eq.8) then
c       call jact4node (r,s,t,ck,h,tt21,nbrel,4,det1,pj,hp,hv1,
c     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,4,hpv)
       call jact8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,iizlaz,ndimp)
	elseif(ndim.eq.5) then
       call jact5node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL)
	endif

c       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
c     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
c     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)

c      IF (IALE.EQ.1) CALL ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
c      IF (INDAMI.EQ.1) CALL NENJ3D(PJ,TT21,AMI,NDIM,IIZLAZ)
 
      WDT=WT*DET1
      ZAPRE=ZAPRE+WDT
  
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM

C      AKVV1(K,N)=AKVV1(K,N)+WDT*((Hpv(K)*HV1*PJ(1,N)+
C     1Hpv(K)*HV2*PJ(2,N)+Hpv(K)*HV3*PJ(3,N))*GUSM)

      AKVV1(K,N)=AKVV1(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM)

c      AKTV1(K,N)=AKTV1(K,N)+WDT*(H(K)*H(N)*ZVXT*GUSM*CC)
c      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVYT*GUSM*CC)
c      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVZT*GUSM*CC)
C================================================================
      AKMIV1(K,N)=AKMIV1(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
C================================================================
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


c       IF (N.LE.ndimp.AND.PENALT.LT.1.D0) THEN
c       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
c       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
c       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
c      ENDIF
 164  CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
  165 CONTINUE

c  170 CONTINUE
c  175 CONTINUE
  180 CONTINUE


       NGAUSX=IBRGT-1
       NGAUSY=IBRGT-1
       NGAUSZ=IBRGT-1


        DO I=1,NDIM
	   DO J=1,NDIM
	    AKIIX(I,I)=AKIIX(I,I)+
     &DABS(AKMIV1(I,J))+DABS(AKVV1(I,J))+DABS(AJV1V1(I,J))
	    AKIIY(I,I)=AKIIY(I,I)+
     &DABS(AKMIV1(I,J))+DABS(AKVV1(I,J))+DABS(AJV2V2(I,J))
	    AKIIZ(I,I)=AKIIZ(I,I)+
     &DABS(AKMIV1(I,J))+DABS(AKVV1(I,J))+DABS(AJV3V3(I,J))
	   ENDDO
	    AKIIX(I,I)=1.D0/AKIIX(I,I)
	    AKIIY(I,I)=1.D0/AKIIY(I,I)
	    AKIIZ(I,I)=1.D0/AKIIZ(I,I)
	  ENDDO

      

c       DO  NGX=1,NGAUSX
c       JR=NREF(NGAUSX) + NGX
c       R = XG(JR)
 
c       DO  NGY=1,NGAUSY
c       JS=NREF(NGAUSY) + NGY
c       S = XG(JS)
 
c       DO  NGZ=1,NGAUSZ
c       JT=NREF(NGAUSZ) + NGZ
c       T = XG(JT)
c       WT=WGT(JR)*WGT(JS)*WGT(JT)


c      nkk=1
      nkk=ndimp
	if (ndim.eq.11) nkk=5
	if (ndim.eq.10) nkk=4
	if (ndim.eq.8) nkk=1
	if (ndim.eq.4) nkk=4


	do kk=1,nkk

      if (nkk.eq.4) then
        call Integ4points(r,s,t,wt,kk)      
	elseif (nkk.eq.5) then
        call Integ5points(r,s,t,wt,kk)      
	elseif (nkk.eq.11) then
        call Integ10points(r,s,t,wt,kk)      
	elseif (nkk.eq.1) then
 	  wt=1.0d0
	  r=0.25d0
	  s=0.25d0
	  t=0.25d0
	endif
      
	if (ndim.eq.4) then
       call jact4node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,ndimp,hpv)
	elseif(ndim.eq.8) then
       call jact8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,iizlaz,ndimp)
	elseif(ndim.eq.10) then
       call jact4node1 (r,s,t,ck,h,tt21,nbrel,4,det1,pjp,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,4,hpv)
       call jact10node1 (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL)
	elseif(ndim.eq.11) then
       call jact4node (r,s,t,ck,h,tt21,nbrel,4,det1,pjp,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,4,hpv)
       call jact11node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL)
	endif
c       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
c     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
c     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)

 
        WDT=WT*DET1
  
        DO  K=1,NDIM
        DO  N=1,NDIM
C      IF (N.LE.8.AND.PENALT.LT.1.D0) THEN
C       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
C      ENDIF
         IF (PENALT.GT.0.D0) THEN
c          PEN=1.d-2
          PEN=penalt
          PENXX(K,N)=PENXX(K,N)+WDT*PJ(1,K)*PJ(1,N)*(PEN)
          PENXY(K,N)=PENXY(K,N)+WDT*PJ(1,K)*PJ(2,N)*(PEN)
          PENXZ(K,N)=PENXZ(K,N)+WDT*PJ(1,K)*PJ(3,N)*(PEN)
          PENYX(K,N)=PENYX(K,N)+WDT*PJ(2,K)*PJ(1,N)*(PEN)
          PENYY(K,N)=PENYY(K,N)+WDT*PJ(2,K)*PJ(2,N)*(PEN)
          PENYZ(K,N)=PENYZ(K,N)+WDT*PJ(2,K)*PJ(3,N)*(PEN)
          PENZX(K,N)=PENZX(K,N)+WDT*PJ(3,K)*PJ(1,N)*(PEN)
          PENZY(K,N)=PENZY(K,N)+WDT*PJ(3,K)*PJ(2,N)*(PEN)
          PENZZ(K,N)=PENZZ(K,N)+WDT*PJ(3,K)*PJ(3,N)*(PEN)



	   else
       IF (N.LE.ndimp.AND.PENALT.LT.1.D0) THEN


C	AKPP(N,N)=AKPP(N,N)+WDT*
C     &(PJ(1,K)*AKIIX(K,K)*PJ(1,K)+PJ(2,K)*AKIIY(K,K)*PJ(2,K)+
C     &PJ(3,K)*AKIIZ(K,K)*PJ(3,K))

      if (k.le.ndimp) then
c	AKPP(K,N)=AKPP(K,N)+WDT*
c     &(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))
	AKPP(K,N)=AKPP(K,N)+WDT*
     &(PJp(1,K)*PJp(1,N)+PJp(2,K)*PJp(2,N)+PJp(3,K)*PJp(3,N))
      endif


C	FPP=FPP+
C     &(PJ(1,K)*AKIIX(K,K)*F92X(K)+PJ(2,K)*AKIIY(K,K)*F92Y(K)+
C     &PJ(3,K)*AKIIZ(K,K)*F92Z(K))


        AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
        AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
        AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
       
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



C       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
      ENDIF

         ENDIF
        ENDDO
        ENDDO
        ENDDO
c       ENDDO
c       ENDDO



C======================================================================= 
C POVRSINSKE SILE
      

      DO 250 JBRPS=1,MAXSIL
      IF (NGPSIL(1,JBRPS).EQ.NBREL) THEN 
      
      nn=0

      do i=1,4
	 do j=2,4
	   if (nel(i,nbrel).eq.ngpsil(j,jbrps)) then
	     nn=nn+1
	     nizel(nn)=nel(i,nbrel)
	   endif
	 enddo
	enddo


      if (nn.eq.3) then
c	write(*,*) nbrel,nizel(1),nizel(2),nizel(3)

	do 210 kk=1,3

      call strana_tetra(nizel,nel,nbrel,r,s,t,wt,kk,ndim)
	call Integ3points(s,t,wt,kk)
c	call Integ3points(r,s,wt,kk)
c      call surfjact10 (r,s,nizel,h,gnode,sf1,sf2,sf3)
c      call Integ4points(r,s,t,wt,kk)      
c	t=1.d0-r-s
	r=0.d0
      if (ndim.eq.10) then
	call jactSurf10node (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami)
      elseif (ndim.eq.8) then
	call jactSurf8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami,ndimp)
	elseif(ndim.eq.5) then
	call jactSurf5node (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami)
	elseif(ndim.eq.4) then
	call jactSurf4node (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami,ndimp)
	endif

      det=SURF1(nizel,CORD)
C
      WDT=WT*DET

      povrs=povrs+wdt  
c      do 209 ii=1,3
      DO 209 K=1,NDIM
c	 if (nel(k,nbrel).ne.nizel(ii)) goto 209
	 RS1(K)=RS1(K)+(H(K))*WDT*SF1
	 RS2(K)=RS2(K)+(H(K))*WDT*SF2
	 RS3(K)=RS3(K)+(H(K))*WDT*SF3
	povrsila=povrsila+(H(K))*WDT*SF1
	povrsila=povrsila+(H(K))*WDT*SF2
	povrsila=povrsila+(H(K))*WDT*SF3
  209 continue

  210 CONTINUE    
      endif 
c       WRITE(IIZLAZ,*)'POVRSINA',JBRPS,'=',POVRS
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

      SKEF(I+3*NDIM+ndimp,J+3*NDIM+ndimp)=AF*(AKK(I,J)+AKVV1(I,J)*CC)
      SKEF(I+3*NDIM+ndimp,J)=AF*AKTV1(I,J)
      SKEF(I+3*NDIM+ndimp,J+NDIM)=AF*AKTV2(I,J)
      SKEF(I+3*NDIM+ndimp,J+2*NDIM)=AF*AKTV3(I,J)
      IF (NSTAC.EQ.0) THEN
      SKEF(I+3*NDIM+ndimp,J+3*NDIM+ndimp)=
     &SKEF(I+3*NDIM+ndimp,J+3*NDIM+ndimp)+AMV2(I,J)*CC/TIME
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
      DO 270 I=1,NDIM
      DO 265 J=1,NDIM
       SKEF(I,J)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J))!+ATGS1(I,J))
C       SKEF(I,J)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J))
       SKEF(I,J+NDIM)=AF*AJV1V2(I,J)
       SKEF(I,J+2*NDIM)=AF*AJV1V3(I,J)
      IF (J.LE.ndimp) THEN
       SKEF(I,J+3*NDIM)=AF*AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AF*AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AF*AKV3P(I,J)
      ENDIF
      IF (I.LE.ndimp) THEN
       SKEF(I+3*NDIM,J)=AF*AKV1P(J,I)
       SKEF(I+3*NDIM,J+NDIM)=AF*AKV2P(J,I)
       SKEF(I+3*NDIM,J+2*NDIM)=AF*AKV3P(J,I)
      ENDIF

	 IF (IPASS.EQ.0.AND.METOD.EQ.4.AND.I.LE.NDIMP.AND.J.LE.NDIMP) THEN
         SKEF(I+3*NDIM,J+3*NDIM)=SKEF(I+3*NDIM,J+3*NDIM)+AKPP(I,J)
	 ENDIF

       SKEF(I+NDIM,J)=AF*AJV2V1(I,J)
       SKEF(I+NDIM,J+NDIM)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV2V2(I,J))
c     &+ATGS2(I,J))
       SKEF(I+NDIM,J+2*NDIM)=AF*AJV2V3(I,J)

      SKEF(I+2*NDIM,J)=AF*AJV3V1(I,J)
      SKEF(I+2*NDIM,J+NDIM)=AF*AJV3V2(I,J)
      SKEF(I+2*NDIM,J+2*NDIM)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV3V3(I,J))
c     &+ATGS3(I,J))
      IF (NSTAC.EQ.0) THEN
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       SKEF(I+2*NDIM,J+2*NDIM)=SKEF(I+2*NDIM,J+2*NDIM)+AMV2(I,J)/TIME
      ENDIF
 265  CONTINUE
 270  CONTINUE
       
      IF (PENALT.GT.0.D0) THEN
c      IF (PEN.GT.0.D0) THEN

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
        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
c	IF (K.EQ.0) F92(II)=F92(II)-ATGS1(I,J)*TT21(JJ)
c	IF (K.EQ.1) F92(II)=F92(II)-ATGS2(I,J)*TT21(JJ)
c	IF (K.EQ.2) F92(II)=F92(II)-ATGS3(I,J)*TT21(JJ)
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
      II=I+3*NDIM+ndimp
C	 F92(II)=F92(II)+RS1(I) 
      DO 292 J=1,NDIM
      JJ=J+3*NDIM+ndimp
       F92(II)=F92(II)-AKTV1(I,J)*TT21(J)-AKTV2(I,J)*TT21(J+NDIM)
     &-AKTV3(I,J)*TT21(J+2*NDIM)-AKK(I,J)*TT21(JJ)

      IF (NSTAC.EQ.0) THEN
       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
 292  CONTINUE
 294  CONTINUE


C==========================================================================
C       CALL NUL(FPOM,21)
C DESNA STRANA USLED PRITISKA
C	IF (IPASS.NE.0.AND.IPASS.NE.4) THEN

      DO 310 I=1,NDIM
      DO 300 J=1,NDIM
      IF (J.LE.ndimp) THEN
       F92(I)=F92(I)-AKV1P(I,J)*TT21(J+3*NDIM)
       F92(I+NDIM)=F92(I+NDIM)-AKV2P(I,J)*TT21(J+3*NDIM)
       F92(I+2*NDIM)=F92(I+2*NDIM)-AKV3P(I,J)*TT21(J+3*NDIM)

      ENDIF
      IF (I.LE.ndimp) THEN
C      IF (I.LE.ndimp) THEN
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV1P(J,I)*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV2P(J,I)*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)-AKV3P(J,I)*TT21(J+2*NDIM)
      ENDIF
 300  CONTINUE
 310  CONTINUE
C      endif

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



C==========================================================================
C Poisson equation
C==========================================================================
	IF (IPASS.EQ.0.AND.METOD.EQ.4) THEN
C  
      DO 311 I=1,NDIMP
      DO 301 J=1,NDIM

       F92(I+3*NDIM)=F92(I+3*NDIM)-(AKVV1P(J,I)+AKMIV1P(J,I))*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)-(AKVV2P(J,I)+AKMIV2P(J,I))
     &*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)-(AKVV3P(J,I)+AKMIV3P(J,I))
     &*TT21(J+2*NDIM)

       ALAMBDA=0.D0
	 IF (NUMPASS.EQ.4) ALAMBDA=1.D0
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV1P(J,I)*TT21(J)*ALAMBDA
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV2P(J,I)*TT21(J+NDIM)*ALAMBDA
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV3P(J,I)*TT21(J+2*NDIM)*ALAMBDA

C       F92(I+3*NDIM)=F92(I+3*NDIM)+FPP

       F92(I+3*NDIM)=F92(I+3*NDIM)-AKPP(I,J)*TT21(J+3*NDIM)
C==========================================================================
 301  CONTINUE
C       F92(I+3*NDIM)=F92(I+3*NDIM)+RSVP(I)
 311  CONTINUE
      ENDIF
C==========================================================================
C
	IF (IPASS.EQ.4) THEN
  
      DO I=1,NDIMP
      DO J=1,NDIM
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV1P(J,I)*TT21(J)
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV2P(J,I)*TT21(J+NDIM)
       F92(I+3*NDIM)=F92(I+3*NDIM)+AKV3P(J,I)*TT21(J+2*NDIM)
C
C==========================================================================
      ENDDO
      ENDDO
      ENDIF

	IF (IPASS.EQ.5) THEN
  
      DO J=1,NDIMP
      DO I=1,NDIM
	F92(I)=0.D0
	NODE=NEL(I,NBREL)
       GNODE(2,1,NODE)=GNODE(2,1,NODE)+
     &(AKV1P(I,J))*AKIIX(I,I)*TT21(J+3*NDIM+NDIMP)
       GNODE(2,2,NODE)=GNODE(2,2,NODE)+
     &(AKV2P(I,J))*AKIIY(I,I)*TT21(J+3*NDIM+NDIMP)
       GNODE(2,3,NODE)=GNODE(2,3,NODE)+
     &(AKV3P(I,J))*AKIIZ(I,I)*TT21(J+3*NDIM+NDIMP)
     
C==========================================================================
      ENDDO
      ENDDO
      ENDIF
C==========================================================================


c      CALL SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
c     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)

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
       CALL SPAKDE (SILE,F92,LM2,NDES)
c       CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
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
c        PPP=0.D0

c      if (iter.gt.0) then      
c	write(iizlaz,*)'nbrel= ',nbrel
c	pp=0.d0
c      DO I=1,NDIM
c	 do j=1,ndimp
c      pp=pp-AKV3P(I,j)*TT21(j+3*NDIM)
c	enddo
c	write(iizlaz,*)'pp= ',nbrel,pp
c      ENDdo

c       ppp=0.d0
c       DO I=1,NDIM
c        I1=I
c        I2=I+NDIM
c        I3=I+2*NDIM
c        PPP=PPP+
c     &(PJp(1,I)*TT21(I1)+PJp(2,I)*TT21(I2)+PJp(3,I)*TT21(I3))
c       ENDDO
c      WRITE(IIZLAZ,*)'PPP= ',PPP,nbrel
c	endif
 400  CONTINUE

c      CALL WRRF(SILE,JEDN,'desno=',IIZLAZ)
c      if (iter.eq.1) stop
     
      write(iizlaz,*)'povrs= ',povrs
      write(iizlaz,*)'povrsila= ',povrsila
      write(iizlaz,*)'zapre= ',zapre
c 	stop
       
C End of subroutine
      if (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) then
c       CALL MUMPSRIGHT(SILE,JEDN)
      endif
C      CALL WRRF(SILE,JEDN,'desno=',IIZLAZ)
c      penalt=0.d0
      END
C=======================================================================
C=======================================================================
