C=========================================================================
      SUBROUTINE Tetra3D4_10nodes(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASS,KKORAK,METOD,NUMPASS,ID1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine RACU3D is used for 3D analysis
CE It is used global loop per elements
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
	COMMON /TRANSP/ INDFL,LID1
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*),ID1(6,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)

      

      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92)
      DIMENSION F92(92)
      DIMENSION ME(46),LE(46)
	DIMENSION IDL(45)


      PARAMETER (ND=16)
      PARAMETER (NK=24)
C      PARAMETER (ND=12)
C      PARAMETER (NK=18)
      
      DIMENSION AKIII(NK,NK)
      DIMENSION AKII(NK,NK),AKBI(ND,NK),AKIB(NK,ND),AKBIII(ND,NK)
	DIMENSION SKEF(ND+NK,*)


      NDIMP=4

	NDES=ND+NK
	NV=4
C	NV=3
	NDIM10=10


      


C=========================================================================
c      if (ipass.gt.0) then
      DO NODE=1,NPT
       DO NZDT=1,NUMZAD
        IF(NODE.EQ.NZAD(1,NZDT)) THEN
         CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
     &NTABFT,IIZLAZ)
         MESTO=NZAD(2,NZDT)
C	if (mesto.eq.3) GNODE(2,MESTO,NODE)=5.D0
	   GNODE(2,MESTO,NODE)=ZADVRE(NZDT)*FK1
c	   GNODE(1,MESTO,NODE)=ZADVRE(NZDT)*FK1
        ENDIF  
	 ENDDO
	ENDDO
c=======================================================================
c      endif
  
C      DO I=1,116
C	 GNODE(2,3,I)=5.D0
C	 GNODE(1,3,I)=5.D0
C	ENDDO 


C      call INDELSSTRES(NEL,INDEL,NET,NPT,NDIM,ID)


      IF (IMUMPS.eq.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) THEN
c      if (kkorak.eq.0.and. iter.eq.0) 
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)
c     &CALL MUMPSINIT(NEL,ID,JEDN,NET,NDIM,5)
c	REWIND(MUFILE2)
c        CALL sparseassembler_init(0)
	ENDIF


C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET

      DO 125 I=1,92
      TT210(I)=0.D0
      LM2(I)=0
	F92(I)=0.D0
 125  TT21(I)=0.D0

      DO I=1,NDES
	 IDL(I)=0
	DO J=1,NDES
	 SKEF(I,J)=0.D0
	ENDDO
	ENDDO 

C      tez=0.d0
C=========================================================================
      DO 130 KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
C	tez=tez+0.125d0*CORD(3,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
      DO KLM=1,NDIM
	 DO J=1,NV
        LM2(J+(KLM-1)*NV)=ID(J,NEL(KLM,NBREL))
	  IDL(J+(KLM-1)*NV)=ID1(J,NEL(KLM,NBREL))
        TT21(J+(KLM-1)*NV)=GNODE(2,J,NEL(KLM,NBREL))
        TT210(J+(KLM-1)*NV)=GNODE(1,J,NEL(KLM,NBREL))
	 ENDDO
	ENDDO

c fill of additional 10 nodes of tetrahedral element
	  do j=1,NV
	  tt21((5-1)*NV+j)=0.5d0*(tt21((1-1)*NV+j)+tt21((2-1)*NV+j))
	  tt21((6-1)*NV+j)=0.5d0*(tt21((2-1)*NV+j)+tt21((3-1)*NV+j))
	  tt21((7-1)*NV+j)=0.5d0*(tt21((1-1)*NV+j)+tt21((3-1)*NV+j))
	  tt21((8-1)*NV+j)=0.5d0*(tt21((1-1)*NV+j)+tt21((4-1)*NV+j))
	  tt21((9-1)*NV+j)=0.5d0*(tt21((2-1)*NV+j)+tt21((4-1)*NV+j))
	  tt21((10-1)*NV+j)=0.5d0*(tt21((3-1)*NV+j)+tt21((4-1)*NV+j))
	  enddo

c fill of additional 10 nodes of tetrahedral element
	  do j=1,NV
	  IDL((5-1)*NV+j)=IDL((1-1)*NV+j).AND.IDL((2-1)*NV+j)
	  IDL((6-1)*NV+j)=IDL((2-1)*NV+j).AND.IDL((3-1)*NV+j)
	  IDL((7-1)*NV+j)=IDL((1-1)*NV+j).AND.IDL((3-1)*NV+j)
	  IDL((8-1)*NV+j)=IDL((1-1)*NV+j).AND.IDL((4-1)*NV+j)
	  IDL((9-1)*NV+j)=IDL((2-1)*NV+j).AND.IDL((4-1)*NV+j)
	  IDL((10-1)*NV+j)=IDL((3-1)*NV+j).AND.IDL((4-1)*NV+j)
	  enddo



C      WRITE(IIZLAZ,'(I6,45I5)')NBREL,(IDL(I),I=1,45)

      

      DO J=1,3
       CK(5,J)=0.5D0*(CK(1,J)+CK(2,J))
       CK(6,J)=0.5D0*(CK(2,J)+CK(3,J))
       CK(7,J)=0.5D0*(CK(1,J)+CK(3,J))
       CK(8,J)=0.5D0*(CK(1,J)+CK(4,J))
       CK(9,J)=0.5D0*(CK(2,J)+CK(4,J))
       CK(10,J)=0.5D0*(CK(3,J)+CK(4,J))
	ENDDO
	
      CALL Tetra4_10(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM10,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,TT21,TT210,CK,LM2,
     &NBREL,F92,NV,IDL)

C==========================================================================
C==========================================================================
C STATIC CONDENSATION ON THE ELEMENT LEVEL
C      ND=18
C	DO I=13,ndes
C	  IF (IDL(I).EQ.1) ND=ND-1
C	ENDDO

      DO I=1,NK
	DO J=1,NK
	 AKII(I,J)=0.D0
	 AKIII(I,J)=0.D0
	ENDDO
	ENDDO

      M=0
	DO 200 I=ND+1,NDES
C	IF (MOD(I-16,4).EQ.0) GOTO 200
	M=M+1
	N=0
	DO 100 J=ND+1,NDES
C	IF (MOD(J-16,4).EQ.0) GOTO 100
	N=N+1
	 AKII(M,N)=SKEF(I,J)
100   CONTINUE	
200   CONTINUE	

      CALL MINV(AKII,NK,DET1,LE,ME)

      M=0
	N=0
	DO I=1,ND
	DO J=1,NK
	 AKBI(I,J)=0.D0
	 AKIB(J,I)=0.D0
	ENDDO
	ENDDO

      DO 201 I=1,ND
	N=0
	DO 101 J=ND+1,NDES
C	IF (MOD(J-16,4).EQ.0) GOTO 101
	N=N+1
	 AKBI(I,N)=SKEF(I,J)
	 AKIB(N,I)=SKEF(J,I)
101   CONTINUE	
201   CONTINUE	

      DO I=1,ND
	 DO J=1,NK
	  AKBIII(I,J)=0.D0
	  DO K=1,NK
	   AKBIII(I,J)=AKBIII(I,J)+AKBI(I,K)*AKII(K,J)
	  ENDDO
	 ENDDO
	ENDDO
       

            

C      WRITE(IIZLAZ,*)'NBREL'
c      goto 301

      DO I=1,ND
	M=0
	 DO 293 J=1,NK
C	  IF (MOD(J,4).EQ.0) GOTO 293
	    M=M+1
          F92(I)=F92(I)-AKBIII(I,M)*F92(J+ND)
293	 CONTINUE
	ENDDO


      DO 300 I=1,ND
	 DO 299 J=1,ND
	  DO 298 K=1,NK
	    SKEF(I,J)=SKEF(I,J)-AKBIII(I,K)*AKIB(K,J) 
  298 CONTINUE
  299 CONTINUE
  300 CONTINUE

301    CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F92,MAXA,LM2,NDES,1)


C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE


C      STOP
C End of subroutine
cc      CALL MUMPSRIGHT(SILE,JEDN)
C      CALL WRRF(SILE,JEDN,'desno=',IIZLAZ)
C      stop
c      write(iizlaz,*)'povrs= ',povrs
c      write(iizlaz,*)'povrsila= ',povrsila
c      write(iizlaz,*)'zapre= ',zapre
c	stop
      END
C=======================================================================
C=======================================================================
C=========================================================================
      SUBROUTINE Tetra4_10(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,TT21,TT210,CK,LM2,NBREL,
     &F92,NV,IDL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine RACU3D is used for 3D analysis
CE It is used global loop per elements
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
	COMMON /TRANSP/ INDFL,LID1
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(4+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)


      DIMENSION CK(21,*),LM2(*),TT21(*),TT210(*),PJ(3,21),TT21A(3*21)
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
	DIMENSION AKPP(10,10)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(*),SHEAR(3)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(21,21),PENXY(21,21),PENXZ(21,21)
      DIMENSION PENYX(21,21),PENYY(21,21),PENYZ(21,21)
      DIMENSION PENZX(21,21),PENZY(21,21),PENZZ(21,21)
      DIMENSION C(21,21)
	DIMENSION SKE(46*93)
C	DIMENSION ATGS1(8,8),ATGS2(8,8),ATGS3(8,8)
      dimension nizel(4)
	DIMENSION IDL(*)


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


      NDIMP=4
	IBRGT=2
	

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
       call jact10node10 (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,NV)
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
       call jact8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,pjp,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,iizlaz,ndimp)
	elseif(ndim.eq.10) then
c       call jact4node1 (r,s,t,ck,h,tt21,nbrel,4,det1,pjp,hp,hv1,
c     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,4,hpv)
       call jact10node10 (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,NV)
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

C      if (k.le.ndimp) then
C	AKPP(K,N)=AKPP(K,N)+WDT*
C     &(PJP(1,K)*PJP(1,N)+PJP(2,K)*PJP(2,N)+PJP(3,K)*PJP(3,N))
C      endif


C	FPP=FPP+
C     &(PJ(1,K)*AKIIX(K,K)*F92X(K)+PJ(2,K)*AKIIY(K,K)*F92Y(K)+
C     &PJ(3,K)*AKIIZ(K,K)*F92Z(K))


        AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
        AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
        AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
       



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
C      DO 263 I=1,NDIM
C      DO 262 J=1,NDIM
C       AJK(I,J)=AKVV1(I,J)*CC
C      C(I,J)=AMV2(I,J)*CC              

C      SKEF(I+3*NDIM+ndimp,J+3*NDIM+ndimp)=AF*(AKK(I,J)+AKVV1(I,J)*CC)
C      SKEF(I+3*NDIM+ndimp,J)=AF*AKTV1(I,J)
C      SKEF(I+3*NDIM+ndimp,J+NDIM)=AF*AKTV2(I,J)
C      SKEF(I+3*NDIM+ndimp,J+2*NDIM)=AF*AKTV3(I,J)
C      IF (NSTAC.EQ.0) THEN
C      SKEF(I+3*NDIM+ndimp,J+3*NDIM+ndimp)=
C     &SKEF(I+3*NDIM+ndimp,J+3*NDIM+ndimp)+AMV2(I,J)*CC/TIME
C      ENDIF
C 262  CONTINUE
C 263  CONTINUE
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

      NI=NV
	NJ=NV 

      DO 270 I=1,NDIM
      DO 265 J=1,NDIM
       SKEF(1+(I-1)*NI,1+(J-1)*NJ)=
     &AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J))
       SKEF(1+(I-1)*NI,2+(J-1)*NJ)=AF*AJV1V2(I,J)
       SKEF(1+(I-1)*NI,3+(J-1)*NJ)=AF*AJV1V3(I,J)
      SKEF(2+(I-1)*NI,1+(J-1)*NJ)=AF*AJV2V1(I,J)
      SKEF(2+(I-1)*NI,2+(J-1)*NJ)=
     &AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV2V2(I,J))
      SKEF(2+(I-1)*NI,3+(J-1)*NJ)=AF*AJV2V3(I,J)

      SKEF(3+(I-1)*NI,1+(J-1)*NJ)=AF*AJV3V1(I,J)
      SKEF(3+(I-1)*NI,2+(J-1)*NJ)=AF*AJV3V2(I,J)
      SKEF(3+(I-1)*NI,3+(J-1)*NJ)=
     &AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV3V3(I,J))

      IF (NSTAC.EQ.0) THEN
      SKEF(1+(I-1)*NI,1+(J-1)*NJ)=
     &SKEF(1+(I-1)*NI,1+(J-1)*NJ)+AMV2(I,J)/TIME
      SKEF(2+(I-1)*NI,2+(J-1)*NJ)=
     &SKEF(2+(I-1)*NI,2+(J-1)*NJ)+AMV2(I,J)/TIME
      SKEF(3+(I-1)*NI,3+(J-1)*NJ)=
     &SKEF(3+(I-1)*NI,3+(J-1)*NJ)+AMV2(I,J)/TIME
      ENDIF
 265  CONTINUE
 270  CONTINUE
       
      IF (PENALT.GT.0.D0) THEN
      DO I=1,NDIM
        I1=1+(I-1)*NI
        I2=2+(I-1)*NI
        I3=3+(I-1)*NI
       DO J=1,NDIM
        J1=1+(J-1)*NJ
        J2=2+(J-1)*NJ
        J3=3+(J-1)*NJ
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


C DESNA STRANA:
      DO K=1,3
      DO 290 I=1,NDIM
        II=K+(I-1)*NI
      DO 285 J=1,NDIM
        JJ=K+(J-1)*NJ
        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
      IF (NSTAC.EQ.0) THEN
        F92(II)=F92(II)-AMV2(I,J)*(TT21(JJ)-TT210(JJ))/TIME
      ENDIF
 285  CONTINUE							      
 290  CONTINUE
      ENDDO

C ZAPREMINSKE SILE SA DESNE STRANE
      DO 315 I=1,NDIM
       F92(1+(I-1)*NI)=F92(1+(I-1)*NI)+RB1(I)
       F92(2+(I-1)*NI)=F92(2+(I-1)*NI)+RB2(I)
       F92(3+(I-1)*NI)=F92(3+(I-1)*NI)+RB3(I)
 315  CONTINUE
C==========================================================================
C       CALL WRR(RS1,21,'RS1= ')
C       CALL WRR(FPOM,21,'FPOM=')
C==========================================================================
C POVRSINSKE SILE SA DESNE STRANE
      DO 320 I=1,NDIM
C ZADAVANJE PRITISAKA NA DESNOJ STRANI
       if (indfl.eq.1) then
        F92(1+(I-1)*NI)=F92(1+(I-1)*NI)+RS1(I)
	 elseif (indfl.eq.0) then
        II=5+(I-1)*NI
        F92(II)=F92(II)+RS1(I)
	 endif
       F92(2+(I-1)*NI)=F92(2+(I-1)*NI)+RS2(I)
       F92(3+(I-1)*NI)=F92(3+(I-1)*NI)+RS3(I)
320   CONTINUE
C==========================================================================
C      WRITE(IIZLAZ,'(3E13.5,3I6)') (SKEF(I,I),I=1,3),(LM2(J),J=1,3)
C      WRITE(IIZLAZ,'(3E13.5)') (F92(I),I=1,3)

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

C      CALL SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
C     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)

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
c	 CALL MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,5)
c       CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
c       CALL SPAKDE (SILE,F92,LM2,NDES)
	 RETURN
      ENDIF



	IF(NJUTRA.EQ.1.AND.ITER.GT.0) THEN
C       CALL SPAKDE (SILE,F92,LM2,NDES)
	ELSE
       IF (ISYMMS.EQ.1) THEN
C	    CALL MATSTE (ALEVO,MAXA,SILE,SKE,F92,LM2,NDES,1)
	 ELSE
C          CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F92,MAXA,LM2,NDES,1)
	 ENDIF
	ENDIF

      DO I=1,NDES
	  IF (IDL(I).EQ.1) THEN
	   DO J=1,NDES
	    SKEF(I,J)=0.D0
	    SKEF(J,I)=0.D0
	   ENDDO
	  ENDIF
	ENDDO
	

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
C 400  CONTINUE

      END
C=======================================================================
c==============================================================================
      subroutine jact10node10 (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,
     &hp,hv1,hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,nel,NV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(3,10),pj(3,*),xjj(3,3),xj(3,3),ck(21,*)
      dimension me(16),le(16)
      dimension xjja(20),Xjjj(17)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension V1(10),V2(10),V3(10),pp(10),temp(10),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(3),vtang(3,10)
	dimension xxjjj(4,4)

                 

       do i=1,ndim
		 h(i)=0.d0 
         do j=1,3
		   p(j,i)=0.d0 
	     enddo
	   enddo

c       r=0.0d0
c       s=0.5d0
c       t=0.5d0
  
       e1=r
       e2=s
       e3=t
       e4=1.0d0-r-s-t

       hp(1)=e1
       hp(2)=e2
       hp(3)=e3
       hp(4)=e4


      h(5)=4.0*e1*e2
      h(6)=4.0*e2*e3
      h(7)=4.0*e1*e3
      h(8)=4.0*e1*e4
      h(9)=4.0*e2*e4
      h(10)=4.0*e3*e4

      h(1)=e1*(2.0*e1-1)
      h(2)=e2*(2.0*e2-1)
      h(3)=e3*(2.0*e3-1)
      h(4)=e4*(2.0*e4-1)
	
c      open (13,file='res.txt')
c      do iii=1,10
c       write(13,*)'h( ',iii,')=',h(iii)   
c      enddo
c      stop

      p(1,1)=4.0*e1-1.0
      p(2,2)=4.0*e2-1.0
      p(3,3)=4.0*e3-1.0

      p(1,4)=-(4.0*e4-1.0)
      p(2,4)=-(4.0*e4-1.0)
      p(3,4)=-(4.0*e4-1.0)

      p(1,5)=4.0*e2
      p(2,5)=4.0*e1

      p(2,6)=4.0*e3
      p(3,6)=4.0*e2

      p(1,7)=4.0*e3
      p(3,7)=4.0*e1


      p(1,8)=4.0*(e4-e1)
      p(2,8)=-4.0*e1
      p(3,8)=-4.0*e1

      p(1,9)=-4.0*e2
      p(2,9)=4.0*(e4-e2)
      p(3,9)=-4.0*e2


      p(1,10)=-4.0*e3
      p(2,10)=-4.0*e3
      p(3,10)=4.0*(e4-e3)


c      k=1
      do i=1,3
      do j=1,3
       xj(i,j)=0.d0
       xjj(i,j)=0.d0
	   do kk=1,ndim
           xj(i,j)=xj(i,j)+p(i,kk)*ck(kk,j)
           xjj(i,j)=xjj(i,j)+p(i,kk)*ck(kk,j)
	   enddo
c		   xjja(k)=xjj(i,j)
c	    k=k+1
	enddo
	enddo
		   
     
c========================================================== 
c      call minv(xjja,3,det1,le,me);
      call minv(xjj,3,det1,le,me);
c========================================================== 
      if (dabs(det1)<1.d-15) then
       write (*,*) 'determinante less then zero for element',nbrel
       return
      endif
	  
  

      aJ_1=(1.0d0/6.0d0)*det1

      det1=aJ_1



c      if (det1<0.d0) stop

      do i=1,3
       do ijj=1,ndim
         pj(i,ijj)=0.d0
	    do k=1,3
            pj(i,ijj)=pj(i,ijj) + xjj(i,k)*p(k,ijj)
	    enddo
	 enddo
	enddo






      DO I=1,NDIM
      V1(I)=TT21(1+(I-1)*NV)
      V2(I)=TT21(2+(I-1)*NV)
      V3(I)=TT21(3+(I-1)*NV)
C      TEMP(I)=TT21(I+3*NDIM+8)
      IF (I.LE.NDIMP) THEN
        PP(I)=TT21(4+(I-1)*NV)
      ENDIF
      ENDDO
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C

      HV1=DOT(H,V1,NDIM)
      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)
  
  


c      IF (DET1.LT.1.D-15) THEN
c       WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
c       WRITE(IIZLAZ,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
c       WRITE(IIZLAZ,*)'DETERMINANTE= ',DET1
c       WRITE(IIZLAZ,*)'NODES  COORDINATES'
c       DO I=1,NDIM
c        WRITE(IIZLAZ,1000) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
c        WRITE(IIZLAZ,1001) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
c	  ENDDO
c       return 
c      ENDIF
 1000 FORMAT(I5,3(D13.5))
 1001 FORMAT(I5,3(f10.6))


      iupwin=0
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




       
      
c      IF(KFIX.GT.0) GO TO 70
      RETURN

       end
c==============================================================================
