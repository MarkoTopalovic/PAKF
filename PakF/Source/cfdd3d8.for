#define MUMPS_CLUSTER .FALSE.
C=========================================================================
      SUBROUTINE CFDD3d8(A,NASLOV,VREME,NTABF,ID,CORD,NEL,NZAD,ZADVRE,
     &TABF,ITFMAX,NGPSIL,MAXA,MHT,TT1,TT10,CCORD,SPAR1,AKc,
     &IKc,ROW,INDEL,SPSIL,MAXIT,NETIP,IULAZ,IIZLAZ,INDAMI,NUMZAD,NPT,
     &NDIM,MAXSIL,JEDN,NWK,NET,NPER,NTABFT,NDES,IDPRIT,IFORM,NSTAC,INDAX
     &,KKORAK,GNODE,GUSM,CC,AKT,EPSTR,AMI,BETA,TETAO,FB2,FB3,PENALT,
     &VVREME,SKEF,PRES,PRIT,ALEVO,DESNO,SILE,NNPER,IF,IS,LMAX,LMAX2,
     &NTOTF,AMASA,VMESH,VELOC,IDALE,GNOD0,CCORD0,VMESH0,IPRESS,PRITIS,
     &CPRESS,IAGAIN,ID1,PRES1,TAU,VOSI)
     
      use servis
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#if(!MUMPS_CLUSTER)
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
#endif
C
CE Subroutine CFDD makes 2-D and 3-D analysis with fluid-structure coupling
C
      COMMON /INTERA/ IINTER,NPTI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE,POS(3)
      COMMON /AKIRA1/ VOLT1,PBALL1,VSR1,PTUBE1
      COMMON /ALE/ IALE,METOD,IALFA
      COMMON /PRIKAZ/ INDSC
      COMMON /STAUNV/ NPRINT
	COMMON /NJUTN/ NJUTRA,ISYMMS
 	COMMON /SPBCON/ IBOUND
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
	COMMON /RESTART/ IRESTA,RESTVR
	COMMON /TRANSP/ INDFL,LID1
	COMMON /ALPHA_SEG/ ALPHAU,ALPHAV,ALPHAW,ALPHAP
#if(!MUMPS_CLUSTER)
	COMMON /MUMPS/ mumps_par
      TYPE (DMUMPS_STRUC) mumps_par
#endif
C------------------------------------------
C ZA TURBULENCIJIU
C------------------------------------------
      COMMON /TURB/ ITURB
      COMMON /POCK/ POCK
      COMMON /POCO/ POCO
      COMMON /DELTAL/ DELTAL
      COMMON /IZID/ IZID
      COMMON /NASLOVF/ NASLOVF
      CHARACTER(*)NASLOV
      CHARACTER*250 NASLOVF
      CHARACTER*1 IMEF*20
      DIMENSION TT1(*),ALEVO(*),DESNO(*),SILE(*),ZADVRE(*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(7,*),NZAD(3,*),ID1(7,*)
      DIMENSION NGPSIL(8,*),MAXA(*)
      DIMENSION SKEF(NDES,*),TT10(*),MHT(*),TABF(2,NTABFT,*)
      DIMENSION ITFMAX(*),INDEL(*),NTABF(*)
      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),CCORD(3,*),PRES(3,*)
      DIMENSION AMASA(*),VMESH(3,*),VELOC(3,*),IDALE(3,*)
      DIMENSION GNOD0(3,*),CCORD0(3,*),VMESH0(3,*)
      DIMENSION AKc(*),IKc(2,*),ROW(*)
      DIMENSION LM2(92)
      DIMENSION VREME(*),GNODE(2,7,*)
      DIMENSION IF(*),IS(*)
	DIMENSION IPRESS(*)
	DIMENSION PRITIS(2,*)
	DIMENSION CPRESS(4,2,*)
C------------------------------------------
C ZA OSI
C------------------------------------------
	DIMENSION PRES1(100,3,*)
	DIMENSION TAU(100,2,*)
	DIMENSION VOSI(100,*)
C------------------------------------------
      DIMENSION A(*)
      REAL A
 
      dimension ipassn(8)
      data ipassn /4,5,1,2,3/
      numpass=5
      ITER=0
      INDSK=1
      IF (IRESTA.EQ.1) THEN
	  VVREME=VVREME+RESTVR
	ENDIF

C VVREME is global time
C KKORAK is global step

      KORAK=0

C KORAK is local step at each period NNPER

      KORAK=KORAK+1
      TIME=VREME(NNPER)
      IPROLAZ=0

111   CONTINUE

       ITER=0
c TT10-vector of unknowns values at start of every time step
      IF (INDFL.EQ.20) THEN
	IF (IPROLAZ.EQ.1.AND.NSTAC.EQ.0) NSTAC=1
	IF (IPROLAZ.EQ.0) NSTAC=0
      NEQF=0 
	DO 122 I=1,NPT
	IF (IPROLAZ.EQ.0) THEN
         ID(1,I)=ID1(1,I)
         ID(2,I)=ID1(2,I)
         ID(3,I)=ID1(3,I)
	   ID(4,I)=1
	   ID(5,I)=1
	   ID(6,I)=1
	   ID(7,I)=1
	ELSE
         ID(1,I)=1
         ID(2,I)=1
         ID(3,I)=1
	   ID(4,I)=1
	   ID(5,I)=ID1(5,I)
	   ID(6,I)=1
	   ID(7,I)=1
	ENDIF
        DO J=1,7 
         IF (ID(J,I).EQ.0) THEN
	     NEQF=NEQF+1
	     ID(J,I)=NEQF	     
	   ELSE
	     ID(J,I)=0
	   ENDIF
	  ENDDO
122	CONTINUE
       CALL MAXATF(MAXA,MHT,ID,NEL,NET,NDIM,NEQF,NWKF,7,NDIM+1,iizlaz)
c       stop
	 JEDN1=NEQF
	 NWK1=NWKF
	ELSE
       NWK1=NWK
	 JEDN1=JEDN
	ENDIF     
      DO 20 I=1,NDES
 20   LM2(I)=0
      IF (IALE.EQ.2) THEN
      DO NODE=1,NPT
       DO I=1,3
        GNOD0(I,NODE)=GNODE(1,I,NODE)
       ENDDO
      ENDDO
      ENDIF

      IF (NPTI.GT.0.AND.NSTAC.EQ.1) CALL CLEAR(GNODE,2*7*NPT)

      CALL ZADNOD(GNODE,ZADVRE,NZAD,TABF,VVREME,ITFMAX,NTABFT,
     &IIZLAZ,NUMZAD,NPT,CORD)

      DO NODE=1,NPT
       DO I=1,7
        GNODE(1,I,NODE)=GNODE(2,I,NODE)
       ENDDO
      ENDDO
      IF (IALE.EQ.2) CALL TCORDC(CCORD,CCORD0,NPT,NETIP)

      IF (NPTI.GT.0.AND.KKORAK.GT.1) THEN
        CALL WALLPS(A,IF,LMAX,GNODE,A(1),NETIP,NSTAC,TIME,CCORD,NPT,
     &NEL,NDIM,NET)
        CALL GREZON(VMESH,CCORD,CCORD0,NPT,NETIP,TIME)
      ENDIF
C==========================================================================
C INITIAL AND BOUNDARY CONDITION FOR SEGREGATED PROCEDURE
      IPASS=1
      if (kkorak.eq.1) CALL CLEAR(prit,idprit*net)
C=========================================================================	
99    CONTINUE
      IF (METOD.EQ.4) THEN
      NEQF=0 
	DO 123 I=1,NPT
	    ID(1,I)=ID1(1,I)
          ID(2,I)=ID1(2,I)
          ID(3,I)=ID1(3,I)
          ID(4,I)=ID1(4,I)
          ID(5,I)=1
          ID(6,I)=1
	   IF (IPASSN(IPASS).EQ.0.OR.IPASSN(IPASS).EQ.4) THEN
	    ID(1,I)=1
	    ID(2,I)=1
	    ID(3,I)=1
C	    ID(4,I)=1
  	   ELSEIF (IPASSN(IPASS).EQ.1) THEN
c	    ID(1,I)=1
          ID(2,I)=1
	    ID(3,I)=1
	    ID(4,I)=1
  	   ELSEIF (IPASSN(IPASS).EQ.2) THEN
	    ID(1,I)=1
C	    ID(2,I)=1
	    ID(3,I)=1
	    ID(4,I)=1
  	   ELSEIF (IPASSN(IPASS).EQ.3) THEN
	    ID(1,I)=1
	    ID(2,I)=1
C	    ID(3,I)=1
          ID(4,I)=1
  	   ELSEIF (IPASSN(IPASS).EQ.5) THEN
c	    ID(1,I)=1
c	    ID(2,I)=1
c	    ID(3,I)=1
	    ID(4,I)=1
	   ENDIF
        DO J=1,7 
         IF (ID(J,I).EQ.0) THEN
	     NEQF=NEQF+1
	     ID(J,I)=NEQF
	   ELSE
	     ID(J,I)=0
	   ENDIF
        ENDDO
123	CONTINUE
       NN=0
       IF(IMUMPS.EQ.0) THEN
       IOSA=3
	IF (IPASSN(IPASS).EQ.3) THEN
	ENDIF
	CALL MAXATF(MAXA,MHT,ID,NEL,NET,NDIM+NN,NEQF,NWKF,7,NDIM+1,iizlaz)
       ENDIF
       JEDN1=NEQF
       NWK1=NWKF
	ELSE
       NWK1=NWK
	 JEDN1=JEDN
	ENDIF
C==========================================================================
      DO I=1,NPT
	 DO J=1,7
        IF (ID(J,I).NE.0) TT1(ID(J,I))=GNODE(2,J,I)
	 ENDDO
	ENDDO  
C==========================================================================
C==========================================================================
C LABEL 100:
C==========================================================================
C MAXIT is maximal number of iteration

 100  IF (ITER.GT.MAXIT) THEN

        WRITE(IIZLAZ,3001) MAXIT
 3001 FORMAT(//
     1,'IT IS ACHIEVED MAXIMUM ',I5,' ITERATIONS WITHOUT MAKING
     1 CONVERGENCE !!'//)
       STOP
       ENDIF

C NWK is number of terms of stiffness matrix in skyline procedure
C ALEVO is upper profile of skyline for asymmetric system
C DESNO is lower profile of skyline for asymmetric system
C SILE is right-hand side of system

      CALL CLEAR(SPSIL,NETIP*NPT)
      IF (ITER.EQ.0.OR.NJUTRA.NE.1) THEN 
       CALL CLEAR(ALEVO,NWK1)
       IF (ISYMMS.EQ.0) CALL CLEAR(DESNO,NWK1)
	ENDIF

      CALL CLEAR(SILE,JEDN1)
      CALL CLEAR(pres,3*npt)

CE DEFINITION OF NUMBER OF GAUSS POINT FOR INTEGRATION
CS ODREDJIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
C gauss point integration for 9-node and 8-node element is 3x3
C gauss point integration for 4-node element is 2x2
C NDIM is number of nodes per element
C NPT is global number of nodes

	IF (NETIP.EQ.2.AND.NDIM.EQ.4) IBRGT=2
	IF (NETIP.EQ.2.AND.NDIM.EQ.9) IBRGT=3
	IF (NETIP.EQ.3.AND.NDIM.EQ.8) IBRGT=2
	IF (NETIP.EQ.3.AND.NDIM.EQ.21) IBRGT=3

C SPSIL are forces from fluid dynamic analysis to solid
C================================================================
CE ONLY TEMPORARILY
CS SAMO PRIVREMENO
      IUPWIN=0
C================================================================
C================================================================
C METOD VREMENSKE INTEGRACIJE
       AF=1.D0
       IF (IALFA.EQ.2) AF=1.D0/2.D0
       IF (IALFA.EQ.3) AF=2.D0/3.D0
C================================================================
C================================================================
C GLAVNA PETLJA PO ELEMENTIMA
C GLOBAL LOOP PER ELEMENTS
C NBREL is counter of elements
	IF (NETIP.EQ.3 ) THEN
      CALL EXPAN3D(GNODE,CORD,CCORD,ID,VVREME,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI/GUSM,KKORAK)
	endif

	IF (NETIP.EQ.2 ) THEN
      IF (ITURB.EQ.1) THEN
       IF(ITER.EQ.0) THEN
       DO I=1,NPT
        GNODE(1,6,I)=POCK
        GNODE(1,7,I)=POCO
       ENDDO
       ENDIF
      CALL RACU2DT(GNODE,NEL,ID,CORD,SKEF,TABF,ITFMAX,
     1NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     1NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,TETAO,FB2,FB3,
     1NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,INDAMI,BETA,NDES,IDPRIT,
     1IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,NPER,NTABFT,PRES,VMESH,
     1IALE,AF,NJUTRA,ISYMMS,DELTAL,PRES1,TAU,VOSI,DT,IBKOR)
      ELSEIF (ITURB.EQ.2) THEN
      CALL RACU2DL(GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,DELTAL)
      ENDIF
          
      IF (IALE.EQ.2) THEN
      CALL ALE2D1 (GNODE,NEL,ID,CORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,CCORD,CCORD0,VMESH0,GNOD0)
      ELSEIF (METOD.EQ.4) THEN
      CALL RACU2D(GNODE,NEL,ID,CCORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,AF,INDEL,NJUTRA,ISYMMS,NZAD,ZADVRE)
     
     	ENDIF
     	
      CALL RACU2D(GNODE,NEL,ID,CCORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE,AF,INDEL,NJUTRA,ISYMMS,NZAD,ZADVRE)

      ENDIF

	IF(NETIP.EQ.3) THEN
       IF (METOD.EQ.2) THEN
      CALL EXIM3D(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &AMASA,VMESH,CCORD,VELOC,IDALE)
      GOTO 480
      ELSEIF (IALE.EQ.2) THEN
      CALL ALE3DN(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE)
      ELSE

	IF (NETIP.EQ.3) THEN
	IF (PENALT.LT.1.D0.AND.NETIP.EQ.3) THEN
	
      if(NDIM.eq.10.or.ndim.eq.4.or.ndim.eq.11) then
       CALL RACU3Dtetra_SEG(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),METOD,NUMPASS,id1)

      else !if (ndim.eq.6) then

      CALL Mixed3D6VP(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),KKORAK,NUMPASS,
     *METOD,id1)
      
      endif

      ELSE

      if (ndim.eq.4) then
      CALL RACU3Dtethra(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),METOD,NUMPASS)
	ENDIF
      if (ndim.eq.10.or.ndim.eq.11) then
      CALL RACU3Dtethra(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),METOD,NUMPASS)
       elseif (ndim.eq.8) then
       if (ITURB.EQ.2) then 
      CALL RACU3DLES(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,DELTAL)
       ELSEIF (ITURB.EQ.1) THEN
       
       IF(ITER.EQ.0) THEN
       DO I=1,NPT
        GNODE(1,6,I)=POCK
        GNODE(1,7,I)=POCO
       ENDDO
       ENDIF
        
      CALL RACU3DT(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK,DELTAL,PRES1,TAU,VOSI,DT,
     &IBKOR)
     
      ELSEif (metod.eq.1) then 
      CALL RACU3D(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,KKORAK)
      elseif (metod.eq.4) then
        CALL Mixed3D8VP(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),KKORAK,NUMPASS,
     &METOD)
      elseif (ndim.eq.21) then
      CALL Mixed3D8VP(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),KKORAK,NUMPASS,
     &METOD)
      endif
      endif

      ENDIF
	else
      CALL RACU3Dtethra(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE,AF,NJUTRA,ISYMMS,ITER,IPASSN(IPASS),METOD,NUMPASS)
      endif
       ENDIF
      ENDIF

      if (jedn1.eq.0) goto 422
      IF (NJUTRA.EQ.1.AND.ITER.GT.0) GOTO 420
       IF(METOD.EQ.4.AND.IPASSN(IPASS).EQ.5) GOTO 480
      IF (IMUMPS.EQ.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) GOTO 422 
      IF (ISYMMS.EQ.1) THEN 
        CALL RESENF(ALEVO,SILE,MAXA,JEDN1,NWK1,1)
	ELSE
        CALL UACTCF(ALEVO,DESNO,SILE,MAXA,JEDN1,1)
	ENDIF
420   IF (ISYMMS.EQ.1) THEN 
        CALL RESENF(ALEVO,SILE,MAXA,JEDN1,NWK1,2)
    	ELSE 
        CALL UACTCF(ALEVO,DESNO,SILE,MAXA,JEDN1,2)
	ENDIF
 422    continue

C Poziv MUMPS paralelnog solvera sa silama i brojem jednacina
      IF (IMUMPS.eq.1) THEN
#if(MUMPS_CLUSTER)
        CALL MUMPSRIGHT(SILE,JEDN1)
        CALL SOLVER(SILE, JEDN1)
#else
	  ISPARSE_N = JEDN1
        
        CALL sparseassembler_getnonzero(ISPARSE_NZ)
        CALL sparseassembler_savesparsefile()
        CALL sparseassembler_kill()
        
        if(ISPARSE_N.gt.0.and.ISPARSE_NZ.gt.0)then
         CALL MUMPS_INIT(ISPARSE_N, ISPARSE_NZ)
         CALL MUMPS_READ_SPARSE_FILE(ISPARSE_NZ)        
         CALL MUMPS_FACTOR(mumps_par)
         CALL MUMPS_SOLUTION(SILE)
         CALL MUMPS_END
        endif
#endif
      ENDIF
      IF(IMUMPS.GE. 5 .and. IMUMPS.LE.9) THEN
        phase = 1
        CALL SEND_INT_ARRAY(1, phase)
      	CALL SEND_INT_ARRAY(1, IMUMPS)
        CALL SOLVER_REMOTE(SILE, JEDN1, 1)
      ENDIF
        DO 440 I=1,jedn1
           TT1(I)=TT1(I)+AF*SILE(I)
440     CONTINUE
C===========================================================================
CE TRANSFER DATA FROM VECTOR TT1 TO MATRIX GNODE 
CE AT THE END OF CURRENT TIME STEP
C===========================================================================
        if (metod.eq.4) then
         call FILLN2p(GNODE,SILE,ID,NPT,IPASSN(IPASS))
	  elseif (metod.eq.1) then
         call FILLN2(GNODE,SILE,ID,NPT)
        elseif (metod.eq.2) then
         call FILLN2(GNODE,SILE,ID,NPT)
	  endif
	WRITE(IIZLAZ,*)'METOD', METOD
  480  CONTINUE

      WRITE(*,*)'ITER= ',ITER
	IF (METOD.EQ.4) WRITE(*,*)'IPASS= ',IPASS
      WRITE(*,*)'PERIOD= ',NNPER
      WRITE(*,*)'STEP= ',KKORAK
      WRITE(*,*)'TIME= ',VVREME
      if (metod.eq.4.and.ipass.eq.1)
     &CALL KONVTF2(TT1,SILE,KONVP,4,ID,ITER,NPT,EPSTR,ISRPS,6,alphap)

      IF (METOD.EQ.4.AND.IPASS.LT.numpass) THEN
		IPASS=IPASS+1
	    GOTO 99
	ENDIF

       
	 if (metod.eq.1) then
	 CALL KONVTF(TT1,SILE,KONVV1,1,ID,ITER,NPT,EPSTR,ISRPS,7)
       CALL KONVTF(TT1,SILE,KONVV2,2,ID,ITER,NPT,EPSTR,ISRPS,7)
       CALL KONVTF(TT1,SILE,KONVV3,3,ID,ITER,NPT,EPSTR,ISRPS,7)
       CALL KONVTF(TT1,SILE,KONVP,4,ID,ITER,NPT,EPSTR,ISRPS,7)
       CALL KONVTF(TT1,SILE,KONVT,5,ID,ITER,NPT,EPSTR,ISRPS,7)
       CALL KONVTF(TT1,SILE,KONVK,6,ID,ITER,NPT,EPSTR,ISRPS,7)
       CALL KONVTF(TT1,SILE,KONVO,7,ID,ITER,NPT,EPSTR,ISRPS,7)
       IF (KONVV1*KONVV2*KONVV3*KONVP*KONVT*KONVK*KONVO.EQ.0) THEN
         ITER=ITER+1
         GO TO 100
       endif       
       endif
       if (metod.eq.4.and.konvp.eq.0) then
         iter=iter+1
         ipass=1
         goto 99
	 ENDIF
 
C FOR CALCULATION FORCES FROM FLUID ON THE WALLS
      IF (NETIP.EQ.2)
     &CALL FORC2D(GNODE,NEL,ID,CCORD,SKEF,
     &TABF,ITFMAX,NBREL,TIME,KKORAK,VVREME,SPSIL,ALEVO,DESNO,SILE,ITER,
     &NGPSIL,MAXA,IBRGT,NASLOV,GUSM,CC,AKT,IIZLAZ,AMI,
     &INDAMI,BETA,TETAO,FB2,FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,
     &NPER,NTABFT,NDES,IDPRIT,IFORM,PENALT,PRESS,NSTAC,INDAX,IUPWIN,
     &PRES,VMESH,IALE)
     
      IF (PENALT.GT.1.D0) THEN
      IF(NETIP.EQ.2) THEN 
       
       IF (ITURB.EQ.1) THEN
      CALL PRIT2DT(GNODE,NEL,CCORD,ID,TT10,PRIT,INDAX,NDIM,
     &IDPRIT,PENALT,NUMZAD,NPT,MAXSIL,JEDN1,NWK1,NET,AMI,INDAMI,IIZLAZ,
     &IUPWIN,GUSM,AKT)
      ELSE
      CALL PENTPR(GNODE,NEL,CCORD,ID,TT10,PRIT,INDAX,NDIM,
     &IDPRIT,PENALT,NUMZAD,NPT,MAXSIL,JEDN1,NWK1,NET,AMI,INDAMI,IIZLAZ,
     &IUPWIN,GUSM,AKT)
      ENDIF
     
      ELSE IF(NETIP.EQ.3) THEN
      if (ndim.eq.8) then
      
      IF (ITURB.EQ.1) THEN
      CALL PENTP3T(GNODE,NEL,CCORD,ID,PRIT,NET,NDIM,IDPRIT,PENALT,
     &IIZLAZ,AMI,ISRPS)
      ELSE
      CALL PENTP3(GNODE,NEL,CCORD,ID,PRIT,NET,NDIM,IDPRIT,PENALT,IIZLAZ,
     &AMI,ISRPS)
      ENDIF
      
	else
      CALL PENTP3Tetra(GNODE,NEL,CCORD,ID,PRIT,NET,NDIM,IDPRIT,PENALT,
     &IIZLAZ,AMI,ISRPS)
	endif
       ENDIF
      ENDIF  

      CALL RETPRI(GNODE,INDEL,NPT,TT1,ID,NDIM,PENALT,PRIT,
     &IDPRIT,NEL,NET,NETIP,CORD,IIZLAZ)

      IF (NETIP.EQ.3)
     &CALL FORC3D(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CCORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE)

 777  CALL IZLLSF(GNODE,IIZLAZ,IDPRIT,ISRPS,NPT,PRES,IFORM)

      IF (NPRINT.NE.0.AND.(MOD(KKORAK,NPRINT).NE.0)) RETURN 
C======================================================================
C ZA STAMPU NEUTRAL FILE, DZIGA 2018
c======================================================================
      CALL STAU09F(GNODE,NPT,69,1,KKORAK)
      CALL STAU09F(GNODE,NPT,69,11,KKORAK)
C------------------------------------------------
      CALL STAU09F(GNODE,NPT,69,21,KKORAK)
      CALL STAU09F(GNODE,NPT,69,22,KKORAK)
C------------------------------------------------
      CALL STAU09F(GNODE,NPT,69,41,KKORAK)
      CALL STAU10F(PRES,NPT,69,51,KKORAK)
C=========================================================================
C WRITING *.VTK FILE IN EVERY STEP    
      CALL VTKSTP(IMEF,PAKUSR,IDUZIF,ISRPS,NASLOV,GNODE,PRES,
     1  KKORAK,NPT,CORD,NETIP,NEL,NDIM,NET)
C=========================================================================
C=========================================================================
C WRITING OSI CALCULATION FILE
C      IF(NETIP.EQ.3.AND.KKORAK.EQ.IBKOR) THEN
      OPEN (79,FILE='OSI_INDEX.vtk',
     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
      CALL VTKHD(79,NASLOV,SKORAK)
      CALL VTKMSH(79,NPT,CORD,NETIP,NEL,NDIM,NET)
      CALL VTKPHD(79,NPT)
       WRITE(79,878)
      WRITE(79,879)
 878  FORMAT('SCALARS OSI_INDEX double')
 879  FORMAT('LOOKUP_TABLE default')
      CALL OSICALC(KKORAK,IBKOR,PRES,PRES1,VOSI,TAU,DT,VVREME,NPT,
     &NUMZAD,NZAD,ZADVRE)
C     ENDIF

C=========================================================================
      IF (INDSC.EQ.1) THEN
       CALL ASCIIF(GNODE,NPT,CCORD,CORD,KKORAK,59,NEL,NET,NDIM,VVREME,
     &PRES)
      ELSEIF (netip.eq.2) THEN
      	III=59
      CALL STAGPF(GNODE,NASLOV,VVREME,KKORAK,1,NPT,III,1,NET,NEL,PRIT,
     &CCORD,CORD,NDIM,NETIP,PENALT,IDPRIT,ISRPS,SPSIL,PRES,ID)
    
      ELSEIF(netip.eq.3.and.(ndim.eq.4.or.ndim.eq.8.or.ndim.eq.10.or.
     &ndim.eq.21.or.ndim.eq.11))then
       call printpos(gnode,npt,nel,net,ndim,cord,penalt,id,pres,spsil)
	III=59
      CALL STAGPF(GNODE,NASLOV,VVREME,KKORAK,1,NPT,III,1,NET,NEL,PRIT,
     &CCORD,CORD,NDIM,NETIP,PENALT,IDPRIT,ISRPS,SPSIL,PRES,ID)
	endif


      
      IF (PENALT.GT.1.D0) THEN
       DO I=1,NPT
        GNODE(2,4,I)=0.D0
       ENDDO
      ENDIF
C==========================================================================
C WRITING FORCES FOR PAK-S IN FILE 'ZFLUID'
      IF (NPTI.GT.0) 
     &CALL ZFLUID(A,SPSIL,NPT,NETIP,LMAX,NTOTF,IIZLAZ)
C==========================================================================
      RETURN
      END
C==========================================================================
C=========================================================================
