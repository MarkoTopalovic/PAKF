C=========================================================================
      SUBROUTINE Tetra3D4_15nodes(GNODE,ALEVO,DESNO,SILE,NEL,
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

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(45,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92)
      DIMENSION CKP(21,3),LM2P(92),TT21P(92),TT210P(92),F92(92)
      DIMENSION IBRICK1(8),IBRICK2(8),IBRICK3(8),IBRICK4(8),ITETRA(4)
	DIMENSION IBRICK(4,8),IBRICK0(8),AKIII(33,33)
      DIMENSION ME(46),LE(46)
      DIMENSION AKII(33,33),AKBI(12,33),AKIB(33,12),AKBIII(12,33)
	DIMENSION IDL(45)

C      DATA IBRICK1/1,5,12,8,7,11,15,13/
C      DATA IBRICK2/2,9,12,5,6,14,15,11/
C      DATA IBRICK3/3,6,11,7,10,14,15,13/
C      DATA IBRICK4/4,8,12,9,10,13,15,14/

	DATA IBRICK1/1,8,12,5,7,13,15,11/
      DATA IBRICK2/2,5,12,9,6,11,15,14/
      DATA IBRICK3/3,7,11,6,10,13,15,14/
      DATA IBRICK4/4,9,12,8,10,14,15,13/
      DATA ITETRA/1,2,3,4/





      NDIMP=8
	NDES=45

      


      DO I=1,8
	 IBRICK(1,I)=IBRICK1(I)
	 IBRICK(2,I)=IBRICK2(I)
	 IBRICK(3,I)=IBRICK3(I)
	 IBRICK(4,I)=IBRICK4(I)
	ENDDO

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
	 DO J=1,3
        LM2(J+(KLM-1)*3)=ID(J,NEL(KLM,NBREL))
	  IDL(J+(KLM-1)*3)=ID1(J,NEL(KLM,NBREL))
        TT21(J+(KLM-1)*3)=GNODE(2,J,NEL(KLM,NBREL))
        TT210(J+(KLM-1)*3)=GNODE(1,J,NEL(KLM,NBREL))
	 ENDDO
	ENDDO

c fill of additional 15 nodes of tetrahedral element
	  do j=1,3
	  tt21((5-1)*3+j)=0.5d0*(tt21((1-1)*3+j)+tt21((2-1)*3+j))
	  tt21((6-1)*3+j)=0.5d0*(tt21((2-1)*3+j)+tt21((3-1)*3+j))
	  tt21((7-1)*3+j)=0.5d0*(tt21((1-1)*3+j)+tt21((3-1)*3+j))
	  tt21((8-1)*3+j)=0.5d0*(tt21((1-1)*3+j)+tt21((4-1)*3+j))
	  tt21((9-1)*3+j)=0.5d0*(tt21((2-1)*3+j)+tt21((4-1)*3+j))
	 tt21((10-1)*3+j)=0.5d0*(tt21((3-1)*3+j)+tt21((4-1)*3+j))
	 tt21((11-1)*3+j)=(1.d0/3.d0)*
     &(tt21((1-1)*3+j)+tt21((2-1)*3+j)+tt21((3-1)*3+j))
	   tt21((12-1)*3+j)=(1.d0/3.d0)*
     &(tt21((1-1)*3+j)+tt21((2-1)*3+j)+tt21((4-1)*3+j))
	   tt21((13-1)*3+j)=(1.d0/3.d0)*
     &(tt21((1-1)*3+j)+tt21((3-1)*3+j)+tt21((4-1)*3+j))
	   tt21((14-1)*3+j)=(1.d0/3.d0)*
     &(tt21((2-1)*3+j)+tt21((3-1)*3+j)+tt21((4-1)*3+j))
	   tt21((15-1)*3+j)=(1.d0/4.d0)*(tt21((1-1)*3+j)+
     &tt21((2-1)*3+j)+tt21((3-1)*3+j)+tt21((4-1)*3+j))
	  enddo

c fill of additional 15 nodes of tetrahedral element
	  do j=1,3
	  IDL((5-1)*3+j)=IDL((1-1)*3+j).AND.IDL((2-1)*3+j)
	  IDL((6-1)*3+j)=IDL((2-1)*3+j).AND.IDL((3-1)*3+j)
	  IDL((7-1)*3+j)=IDL((1-1)*3+j).AND.IDL((3-1)*3+j)
	  IDL((8-1)*3+j)=IDL((1-1)*3+j).AND.IDL((4-1)*3+j)
	  IDL((9-1)*3+j)=IDL((2-1)*3+j).AND.IDL((4-1)*3+j)
	  IDL((10-1)*3+j)=IDL((3-1)*3+j).AND.IDL((4-1)*3+j)
	  IDL((11-1)*3+j)=IDL((1-1)*3+j).AND.IDL((2-1)*3+j).
     &AND.IDL((3-1)*3+j)
	  IDL((12-1)*3+j)=IDL((1-1)*3+j).AND.IDL((2-1)*3+j).
     &AND.IDL((4-1)*3+j)
	  IDL((13-1)*3+j)=IDL((1-1)*3+j).AND.IDL((3-1)*3+j).
     &AND.IDL((4-1)*3+j)
	  IDL((14-1)*3+j)=IDL((2-1)*3+j).AND.IDL((3-1)*3+j).
     &AND.IDL((4-1)*3+j)
	  IDL((15-1)*3+j)=IDL((1-1)*3+j).AND.IDL((2-1)*3+j).
     &AND.IDL((3-1)*3+j).AND.IDL((4-1)*3+j)
	  enddo

C      WRITE(IIZLAZ,'(I6,45I5)')NBREL,(IDL(I),I=1,45)

      

      DO J=1,3
       CK(5,J)=0.5D0*(CK(1,J)+CK(2,J))
       CK(6,J)=0.5D0*(CK(2,J)+CK(3,J))
       CK(7,J)=0.5D0*(CK(1,J)+CK(3,J))
       CK(8,J)=0.5D0*(CK(1,J)+CK(4,J))
       CK(9,J)=0.5D0*(CK(2,J)+CK(4,J))
       CK(10,J)=0.5D0*(CK(3,J)+CK(4,J))
       CK(11,J)=(1.D0/3.D0)*(CK(1,J)+CK(2,J)+CK(3,J))
       CK(12,J)=(1.D0/3.D0)*(CK(1,J)+CK(2,J)+CK(4,J))
       CK(13,J)=(1.D0/3.D0)*(CK(1,J)+CK(3,J)+CK(4,J))
       CK(14,J)=(1.D0/3.D0)*(CK(2,J)+CK(3,J)+CK(4,J))
       CK(15,J)=(1.D0/4.D0)*(CK(1,J)+CK(2,J)+CK(3,J)+CK(4,J))
	ENDDO
	
C MAIN LOOP OVER 4 BRICKS INSIDE TETRAHEDRAL ELEMENT
	DO K=1,4


      DO KLM=1,NDIMP
	 DO J=1,3
        LM2P(J+(KLM-1)*3)=LM2(J+(IBRICK(K,KLM)-1)*3)
        TT21P(J+(KLM-1)*3)=TT21(J+(IBRICK(K,KLM)-1)*3)
        TT210P(J+(KLM-1)*3)=TT210(J+(IBRICK(K,KLM)-1)*3)
	 ENDDO
	ENDDO
      

	DO I=1,8
	  IBRICK0(I)=IBRICK(K,I)
	 DO J=1,3
	  CKP(I,J)=CK(IBRICK(K,I),J)
	 ENDDO
	ENDDO
C=========================================================================
C      WRITE(IIZLAZ,'(I6,8I5)')NBREL,(LM2P(I),I=1,8)
C      WRITE(IIZLAZ,'(I6,8E13.5)')NBREL,(CKP(I,3),I=1,8)
C      WRITE(IIZLAZ,'(I6,24E13.5)')NBREL,(TT21P(I),I=1,24)

      CALL Tetra4Bricks(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,8,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,TT21P,TT210P,CKP,LM2P,
     &NBREL,IBRICK0,F92)
	ENDDO

      CALL Tetra4Bricks(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,4,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,TT21P,TT210P,CK,LM2P,
     &NBREL,ITETRA,F92)

C==========================================================================
C==========================================================================
C STATIC CONDENSATION ON THE ELEMENT LEVEL
      ND=33
	DO I=13,45
	  IF (IDL(I).EQ.1) ND=ND-1
	ENDDO

      DO I=1,33
	DO J=1,33
	 AKII(I,J)=0.D0
	 AKIII(I,J)=0.D0
	ENDDO
	ENDDO

	DO 200 I=13,45
C	IF (IDL(I).EQ.1) GOTO 200
	DO 100 J=13,45
C	IF (IDL(J).EQ.1) GOTO 100
	 AKII(I-12,J-12)=SKEF(I,J)
C	IF (DABS(SKEF(I,J)).GT.(0.01*PENALT)) 
C     &AKII(I-12,J-12)=SKEF(I,J)/PENALT
100   CONTINUE	
200   CONTINUE	

      CALL MINV(AKII,33,DET1,LE,ME)


C	DO I=13,45
C	DO J=13,45
C	IF (IDL(I).EQ.1.AND.IDL(J).EQ.1) AKII(I,J)=0.D0
C	IF (DABS(SKEF(I,J)).GT.(0.01*PENALT)) 
C    &AKII(I-12,J-12)=AKII(I-12,J-12)*PENALT
C      ENDDO
C      ENDDO

C      M=0
	
C      DO 202 I=13,45
C	IF (IDL(I).EQ.1) GOTO 200
C	M=M+1
C	N=0
C	DO 102 J=13,45
C	IF (IDL(J).EQ.1) GOTO 100
C	N=N+1
C	 AKII(I,J)=AKIII(M,N)
C102   CONTINUE	
C202   CONTINUE	

 
      M=0
	N=0
	DO I=1,12
	DO J=1,33
	 AKBI(I,J)=0.D0
	 AKIB(J,I)=0.D0
	ENDDO
	ENDDO

      DO 201 I=1,12
C	IF (IDL(I).EQ.1) GOTO 201
C	M=M+1
	DO 101 J=13,45
C	IF (IDL(J).EQ.1) GOTO 101
C	N=N+1
C	 AKBI(M,N)=SKEF(I,J)
C	 AKIB(N,M)=SKEF(J,I)
	 AKBI(I,J-12)=SKEF(I,J)
	 AKIB(J-12,I)=SKEF(J,I)
101   CONTINUE	
201   CONTINUE	

      DO I=1,12
	 DO J=1,33
	  AKBIII(I,J)=0.D0
	  DO K=1,33
	   AKBIII(I,J)=AKBIII(I,J)+AKBI(I,K)*AKII(K,J)
	  ENDDO
	 ENDDO
	ENDDO
       

            

C      WRITE(IIZLAZ,*)'NBREL'

      DO I=1,12
	 DO J=1,33
C	   IF (IDL(J+12).EQ.0) THEN
          F92(I)=F92(I)-AKBIII(I,J)*F92(J+12)
C	   ENDIF
	 ENDDO
	ENDDO


      DO 300 I=1,12
C	   IF (IDL(I).EQ.1) GOTO 300
	 DO 299 J=1,12
C	   IF (IDL(J).EQ.1) GOTO 299
	  DO 298 K=1,33
C	   IF (IDL(K+12).EQ.1) GOTO 298
	    SKEF(I,J)=SKEF(I,J)-AKBIII(I,K)*AKIB(K,J) 
  298 CONTINUE
  299 CONTINUE
C	 WRITE(IIZLAZ,'(12E13.5)')(SKEF(I,J),J=1,12)
  300 CONTINUE

      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F92,MAXA,LM2,45,1)


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
      SUBROUTINE Tetra4Bricks(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF1,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,TT21,TT210,CK,LM2,NBREL,
     &IBRICK,F921)
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

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(45,45),SKEF1(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)
	DIMENSION IBRICK(*)


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
	DIMENSION AKPP(8,8)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3),F921(*)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(21,21),PENXY(21,21),PENXZ(21,21)
      DIMENSION PENYX(21,21),PENYY(21,21),PENYZ(21,21)
      DIMENSION PENZX(21,21),PENZY(21,21),PENZZ(21,21)
      DIMENSION C(21,21)
	DIMENSION SKE(46*93)
	DIMENSION ATGS1(8,8),ATGS2(8,8),ATGS3(8,8)

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


      NDIMP=8
	IBRGT=2
	

      
C GLAVNA PETLJA PO ELEMENTIMA
C      DO 400 NBREL=1,NET

      
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
 199  CONTINUE

C===========================================================================
      IF (NDIM.EQ.8) THEN
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
 


       CALL JACTNPP(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
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

      AJV1V1(K,N)=AJV1V1(K,N)+H(K)*HXU*H(N)*GUSM*WDT
      AJV1V2(K,N)=AJV1V2(K,N)+H(K)*HYU*H(N)*GUSM*WDT
      AJV1V3(K,N)=AJV1V3(K,N)+H(K)*HZU*H(N)*GUSM*WDT
      AJV2V1(K,N)=AJV2V1(K,N)+H(K)*HXV*H(N)*GUSM*WDT
      AJV2V2(K,N)=AJV2V2(K,N)+H(K)*HYV*H(N)*GUSM*WDT
      AJV2V3(K,N)=AJV2V3(K,N)+H(K)*HZV*H(N)*GUSM*WDT
      AJV3V1(K,N)=AJV3V1(K,N)+H(K)*HXW*H(N)*GUSM*WDT
      AJV3V2(K,N)=AJV3V2(K,N)+H(K)*HYW*H(N)*GUSM*WDT
      AJV3V3(K,N)=AJV3V3(K,N)+H(K)*HZW*H(N)*GUSM*WDT

      if (K.LE.NDIMP.AND.N.LE.NDIMP) THEN
c	AKPP(K,N)=AKPP(K,N)+
c     &WDT*(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))
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
      ENDIF
C Za Borisa trebalo nesto:

C      DO  K=1,NDIM
C      DO  n=1,NDIM
C     	   write(iizlaz,8787) k,n,amv2(k,n)
C      enddo
C	enddo
C8787  format ('Mass(',i2,',',i2,')= ',(1pd13.5))	
      IF (PENALT.LT.1.D0.AND.NDIM.EQ.8) THEN

       NGAUSX=IBRGT-1
       NGAUSY=IBRGT-1
       NGAUSZ=IBRGT-1
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
 
       CALL JACTNPP(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,ndimp,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)
 
        WDT=WT*DET1
  

        
        DO  N=1,NDIMP
        DO  K=1,NDIMP

	AKPP(K,N)=AKPP(K,N)+
     &WDT*(PJ(1,K)*PJ(1,N)+PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))

c	AKPP(N,N)=AKPP(N,N)+WDT*
c     &(PJ(1,K)*AKIIX(K,K)*PJ(1,K)+PJ(2,K)*AKIIY(K,K)*PJ(2,K)+
c     &PJ(3,K)*AKIIZ(K,K)*PJ(3,K))


c	FPP=FPP+
c     &(PJ(1,K)*AKIIX(K,K)*F92X(K)+PJ(2,K)*AKIIY(K,K)*F92Y(K)+
c     &PJ(3,K)*AKIIZ(K,K)*F92Z(K))


c        AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
c        AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
c        AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
       
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
      ENDIF

     

      IF (PENALT.GT.0.D0.AND.NDIM.EQ.4) THEN


      NKK=1
	do kk=1,nkk

      if (nkk.eq.4) then
        call Integ4points(r,s,t,wt,kk)      
	elseif (nkk.eq.5) then
        call Integ5points(r,s,t,wt,kk)      
	elseif (nkk.eq.11) then
        call Integ10points(r,s,t,wt,kk)      
	elseif (nkk.eq.1) then
 	  wt=-0.8d0
	  r=0.25d0
	  s=0.25d0
	  t=0.25d0
	endif
      

       call jact4node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,4,hpv)

 
        WDT=WT*DET1
  
        DO  K=1,NDIM
        DO  N=1,NDIM
C      IF (N.LE.8.AND.PENALT.LT.1.D0) THEN
C       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
C      ENDIF
         
          PENXX(K,N)=PENXX(K,N)+WDT*PJ(1,K)*PJ(1,N)*(PENALT)*1.D5
          PENXY(K,N)=PENXY(K,N)+WDT*PJ(1,K)*PJ(2,N)*(PENALT)*1.D5
          PENXZ(K,N)=PENXZ(K,N)+WDT*PJ(1,K)*PJ(3,N)*(PENALT)*1.D5
          PENYX(K,N)=PENYX(K,N)+WDT*PJ(2,K)*PJ(1,N)*(PENALT)*1.D5
          PENYY(K,N)=PENYY(K,N)+WDT*PJ(2,K)*PJ(2,N)*(PENALT)*1.D5
          PENYZ(K,N)=PENYZ(K,N)+WDT*PJ(2,K)*PJ(3,N)*(PENALT)*1.D5
          PENZX(K,N)=PENZX(K,N)+WDT*PJ(3,K)*PJ(1,N)*(PENALT)*1.D5
          PENZY(K,N)=PENZY(K,N)+WDT*PJ(3,K)*PJ(2,N)*(PENALT)*1.D5
          PENZZ(K,N)=PENZZ(K,N)+WDT*PJ(3,K)*PJ(3,N)*(PENALT)*1.D5
        ENDDO
	  ENDDO
	  ENDDO

        ENDIF


      IF (PENALT.GT.0.D0.AND.NDIM.EQ.8) THEN
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
 
       CALL JACTNPP(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
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
          
		PEN=1.D0/PENALT
c	    IF (K.EQ.1.or.N.EQ.1) PEN=penalt
          PENXX(K,N)=PENXX(K,N)+WDT*PJ(1,K)*PJ(1,N)*PEN
          PENXY(K,N)=PENXY(K,N)+WDT*PJ(1,K)*PJ(2,N)*PEN
          PENXZ(K,N)=PENXZ(K,N)+WDT*PJ(1,K)*PJ(3,N)*PEN
          PENYX(K,N)=PENYX(K,N)+WDT*PJ(2,K)*PJ(1,N)*PEN
          PENYY(K,N)=PENYY(K,N)+WDT*PJ(2,K)*PJ(2,N)*PEN
          PENYZ(K,N)=PENYZ(K,N)+WDT*PJ(2,K)*PJ(3,N)*PEN
          PENZX(K,N)=PENZX(K,N)+WDT*PJ(3,K)*PJ(1,N)*PEN
          PENZY(K,N)=PENZY(K,N)+WDT*PJ(3,K)*PJ(2,N)*PEN
          PENZZ(K,N)=PENZZ(K,N)+WDT*PJ(3,K)*PJ(3,N)*PEN
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
       CALL JACTNPP(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
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
C      DO 263 I=1,NDIM
C      DO 262 J=1,NDIM
CC       AJK(I,J)=AKVV1(I,J)*CC
C      C(I,J)=AMV2(I,J)*CC              
C   
C
C      SKEF(5+(I-1)*3,5+(J-1)*3)=AF*(AKK(I,J)+AKVV1(I,J)*CC)
C      SKEF(5+(I-1)*3,1+(J-1)*3)=AF*AKTV1(I,J)
C      SKEF(5+(I-1)*3,2+(J-1)*3)=AF*AKTV2(I,J)
C      SKEF(5+(I-1)*3,3+(J-1)*3)=AF*AKTV3(I,J)
C      IF (NSTAC.EQ.0) THEN
C      SKEF(5+(I-1)*3,5+(J-1)*3)=
C     &SKEF(5+(I-1)*3,5+(J-1)*3)+AMV2(I,J)*CC/TIME
C      ENDIF
C 262  CONTINUE
C 263  CONTINUE
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
       SKEF(1+(I-1)*3,1+(J-1)*3)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J))
       SKEF(1+(I-1)*3,2+(J-1)*3)=AF*AJV1V2(I,J)
       SKEF(1+(I-1)*3,3+(J-1)*3)=AF*AJV1V3(I,J)
C      IF (J.LE.NDIMP) THEN
C       SKEF(1+(I-1)*3,4+(J-1)*3)=AF*AKV1P(I,J)
C       SKEF(2+(I-1)*3,4+(J-1)*3)=AF*AKV2P(I,J)
C       SKEF(3+(I-1)*3,4+(J-1)*3)=AF*AKV3P(I,J)
C      ENDIF
C      IF (I.LE.NDIMP) THEN
C       SKEF(4+(I-1)*3,1+(J-1)*3)=AF*AKV1P(J,I)
C       SKEF(4+(I-1)*3,2+(J-1)*3)=AF*AKV2P(J,I)
C       SKEF(4+(I-1)*3,3+(J-1)*3)=AF*AKV3P(J,I)
C==========================================================================
C==========================================================================
C      ENDIF


       SKEF(2+(I-1)*3,1+(J-1)*3)=AF*AJV2V1(I,J)
       SKEF(2+(I-1)*3,2+(J-1)*3)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV2V2(I,J))
       SKEF(2+(I-1)*3,3+(J-1)*3)=AF*AJV2V3(I,J)

      SKEF(3+(I-1)*3,1+(J-1)*3)=AF*AJV3V1(I,J)
      SKEF(3+(I-1)*3,2+(J-1)*3)=AF*AJV3V2(I,J)
      SKEF(3+(I-1)*3,3+(J-1)*3)=AF*(AKMIV1(I,J)+AKVV1(I,J)+AJV3V3(I,J))
c    &+ATGS3(I,J))
      IF (NSTAC.EQ.0) THEN
      SKEF(1+(I-1)*3,1+(J-1)*3)=SKEF(1+(I-1)*3,1+(J-1)*3)+AMV2(I,J)/TIME
      SKEF(2+(I-1)*3,2+(J-1)*3)=SKEF(2+(I-1)*3,2+(J-1)*3)+AMV2(I,J)/TIME
      SKEF(3+(I-1)*3,3+(J-1)*3)=SKEF(3+(I-1)*3,3+(J-1)*3)+AMV2(I,J)/TIME
      ENDIF
 265  CONTINUE
 270  CONTINUE
       
      IF (PENALT.GT.0.D0) THEN
      DO I=1,NDIM
        I1=1+(I-1)*3
        I2=2+(I-1)*3
        I3=3+(I-1)*3
       DO J=1,NDIM
        J1=1+(J-1)*3
        J2=2+(J-1)*3
        J3=3+(J-1)*3
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
      DO K=1,3
      DO 290 I=1,NDIM
        II=K+(I-1)*3
      DO 285 J=1,NDIM
        JJ=K+(J-1)*3
        F92(II)=F92(II)-(AKMIV1(I,J)+AKVV1(I,J))*TT21(JJ)
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
C      DO 294 I=1,NDIM
C      II=5+(I-1)*3
Cc	 F92(II)=F92(II)+RS1(I) 
C      DO 292 J=1,NDIM
C      JJ=5+(J-1)*3
C       F92(II)=F92(II)-AKTV1(I,J)*TT21(J)-AKTV2(I,J)*TT21(J+NDIM)
C     &-AKTV3(I,J)*TT21(J+2*NDIM)-AKK(I,J)*TT21(JJ)
C
C      IF (NSTAC.EQ.0) THEN
C       F92(II)=F92(II)-C(I,J)*(TT21(JJ)-TT210(JJ))/TIME
C      ENDIF
C 292  CONTINUE
C 294  CONTINUE


C==========================================================================
C       CALL NUL(FPOM,21)
C DESNA STRANA USLED PRITISKA
C	IF (IPASS.NE.0.AND.IPASS.NE.4) THEN
     
C      DO 310 I=1,NDIM
C      DO 300 J=1,NDIM
C      IF (J.LE.NDIMP) THEN
C       F92(1+(I-1)*3)=F92(1+(I-1)*3)-AKV1P(I,J)*TT21(4+(J-1)*3)
C       F92(2+(I-1)*3)=F92(2+(I-1)*3)-AKV2P(I,J)*TT21(4+(J-1)*3)
C       F92(3+(I-1)*3)=F92(3+(I-1)*3)-AKV3P(I,J)*TT21(4+(J-1)*3)
C      ENDIF
C      IF (I.LE.NDIMP) THEN
C       F92(4+(I-1)*3)=F92(4+(I-1)*3)-AKV1P(J,I)*TT21(1+(J-1)*3)
C       F92(4+(I-1)*3)=F92(4+(I-1)*3)-AKV2P(J,I)*TT21(2+(J-1)*3)
C       F92(4+(I-1)*3)=F92(4+(I-1)*3)-AKV3P(J,I)*TT21(3+(J-1)*3)
C      ENDIF
C 300  CONTINUE
C 310  CONTINUE
C      ENDIF
C==========================================================================
C ZAPREMINSKE SILE SA DESNE STRANE
      DO 315 I=1,NDIM
       F92(1+(I-1)*3)=F92(1+(I-1)*3)+RB1(I)
       F92(2+(I-1)*3)=F92(2+(I-1)*3)+RB2(I)
       F92(3+(I-1)*3)=F92(3+(I-1)*3)+RB3(I)
 315  CONTINUE
C==========================================================================
C       CALL WRR(RS1,21,'RS1= ')
C       CALL WRR(FPOM,21,'FPOM=')
C==========================================================================
C POVRSINSKE SILE SA DESNE STRANE
      DO 320 I=1,NDIM
C ZADAVANJE PRITISAKA NA DESNOJ STRANI
       if (indfl.eq.1) then
        F92(1+(I-1)*3)=F92(1+(I-1)*3)+RS1(I)
	 elseif (indfl.eq.0) then
        II=5+(I-1)*3
        F92(II)=F92(II)+RS1(I)
	 endif
       F92(2+(I-1)*3)=F92(2+(I-1)*3)+RS2(I)
       F92(3+(I-1)*3)=F92(3+(I-1)*3)+RS3(I)
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
       CALL SPAKDE (SILE,F92,LM2,NDES)
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


      DO K=1,3
	 DO I=1,NDIM
        M=IBRICK(I) 
        DO J=1,NDIM
	   
	   N=IBRICK(J) 
	   SKEF1((M-1)*3+K,(N-1)*3+K)=SKEF1((M-1)*3+K,(N-1)*3+K)+
     &SKEF((I-1)*3+K,(J-1)*3+K)
	  ENDDO      
	   F921((M-1)*3+K)=F921((M-1)*3+K)+F92((I-1)*3+K)

	 ENDDO      
	ENDDO      


C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
C 400  CONTINUE

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
C==========================================================================
       SUBROUTINE JACTNPP(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN,
     &NDIMP)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION H(21),P(3,21),NEL(4+1,*)
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
      V1(I)=TT21(1+(I-1)*3)
      V2(I)=TT21(2+(I-1)*3)
      V3(I)=TT21(3+(I-1)*3)
C      TEMP(I)=TT21(I+3*NDIM+8)
      IF (I.LE.NDIMP) THEN
        PP(I)=TT21(4+(I-1)*3)
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



      IF (DET1.LT.1.D-17) THEN
       WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
       WRITE(IIZLAZ,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
       WRITE(IIZLAZ,*)'DETERMINANTE= ',DET1
       WRITE(IIZLAZ,*)'NODES  COORDINATES'
       DO I=1,NDIM
        WRITE(IIZLAZ,1000) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
	  ENDDO
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
