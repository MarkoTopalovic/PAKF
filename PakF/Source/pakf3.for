C$DEBUG      
C===========================================================================
C===========================================================================
C     SUBROUTINE  INPAKF
C                 DATTIM
C                 READFL
C                 NAMEF
C                 ULAZF1
C                 ULAZF2
C                 ULAZF4
C                 INPU06
C                 INPU07
C                 INPU08
C                 TIDPRI
C                 TNDIM
C                 INPF03
C                SPSOUT
C                IZLLSF
C                NENJUT
C                SSPSIL
C                IWRIT
C                WRIT
C===========================================================================
C===========================================================================
      SUBROUTINE INPAKF(A,NTOT,LMAX,NKDT,DTDT,NPER,IBKOR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /INTERA/ IINTER,NPTI
C
CE Subroutine INPAKF is used for prepare data for reading input 
CE data for fluid analysis
C
      PARAMETER (IIDIM = 46)
      PARAMETER (IRDIM = 11)

      CHARACTER*250NASLOV
      CHARACTER*250 NASLOVF
      COMMON /NASLOVF/ NASLOVF
      DIMENSION A(*),NKDT(*),DTDT(*)
C      EQUIVALENCE(A(1),IA(1))
      REAL A

      NPTI=0


       CALL MEMORY(LMAX0,LMAX,0,1,NTOT,51)
       LMAX0=LMAX
       LMAX=LMAX+4
       CALL MEMORY(LHEAD,LMAX,80,1,NTOT,51)

       CALL MEMORY(LIDAT,LMAX,IIDIM,1,NTOT,51)
       CALL MEMORY(LDAT,LMAX,IRDIM,2,NTOT,51)

       CALL REAINT (A(LMAX0),IIDIM)
       CALL REAINT (A(LMAX0+1),IRDIM)
       CALL REAINT (A(LMAX0+2),LIDAT)
       CALL REAINT (A(LMAX0+3),LDAT)

       CALL READFL(NASLOV,A,A(LIDAT),A(LDAT),NTOT,LMAX,IBKOR)
       CALL DATTIM(A,A(LIDAT),NKDT,DTDT,NPER)
       CALL HEADIN(A(LHEAD),NASLOV) 
              
       DO I=1,80        
       NASLOVF(I:I)=NASLOV(I:I)
       ENDDO
      
      END
C===========================================================================
C===========================================================================
       SUBROUTINE DATTIM(A,IF,NKDT,DTDT,NPERF)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION A(*),NKDT(*),IF(*),DTDT(*)
       REAL A
C
CE Subroutine DATTIM is used for definition vector of time periods
C
       NPERF=IF(36)
       LNTABF=IF(2)     
       LVREME=IF(1)     

  	   CALL JEDNA1(DTDT,A(LVREME),NPERF)

       DO NNPER=1,NPERF
  	   NKDT(NNPER)=NENPER(A(LNTABF),NNPER)
C  	   DTDT(NNPER)=A(LVREME+2*(NNPER-1))
       ENDDO
       

      END
C===========================================================================
      SUBROUTINE READFL(NASLOV,A,IDATF,DATF,NTOT,LMAX,IBKOR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine READFL is used for reading input data for fluid analysis
C

      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /BROJFK/ INDFOR,NULAZ
      COMMON /ALE/ IALE,METOD,IALFA
      COMMON /PRIKAZ/ INDSC
	COMMON /NJUTN/ NJUTRA,ISYMMS
	COMMON /SPBCON/ IBOUND
      COMMON /CHECKJACOB3D/ JACOB3D
	COMMON /TRANSP/ INDFL,LID1
C------------------------------------------
C ZA TURBULENCIJIU
C------------------------------------------
      COMMON /TURB/ ITURB
      COMMON /POCK/ POCK
      COMMON /POCO/ POCO
      COMMON /DELTAL/ DELTAL
      COMMON /IZID/ IZID
C------------------------------------------

      DIMENSION A(*),IDATF(*),DATF(*)
      REAL A
      
      CHARACTER*250NASLOV
      CHARACTER*1 IMEF*20
      CHARACTER *24 LSTFL,UNVFL,PAKNEU

C NASLOV-heading of problem

C physical properties:
C GUSM-density of fluid
C CC-specific heat coefficient when pressure is constant
C AKT-conduction coefficient

C EPSTR-relative tolerance
C MAXIT-maximal number of iterations
C NETIP-type of finite element (2-2D or 3-3D)
C IULAZ-input device number
C IIZLAZ-output device number
C AMI-dynamic viscosity
C INDAMI-inicator for nonlinear model of dynamic viscosity
C BETA-thermic expansion coefficient
C TETAO-referential temperature when thermic expansion is using
C FB2-gravity acceleration at x direction
C FB3-gravity acceleration at y direction
C NUMZAD-number of prescribed values
C NPT-total number of nodes
C NDIM-number of nodes per finite element
C MAXSIL-number of presribed surface forces
C NEQF-number of equation for fluid system
C NWKF-number of terms in system of equations for fluid
C NET-total number of element
C NPER-total number of periods
C NTABFT-total number of time functions						 
C NDES-number of unknowns per finite element
C IDPRIT-number of unknowns pressure per finite element
C IFORM-indicator for format condition
C IALFA-indicator for integration method
C IPENAL-indicator for using penalty method
C MBAND-indicator for definition of equations
C PENALT-penalty factor
C NSTAC-indicator for steady analysis 1-steady 0-unsteady
C INDAX-indicator for axisymmetric condition
C NUMZE-number prescribed values for moving fluid mesh (elasticity)
C NEQE-number of equations for moving fluid mesh (elasticity)
C NWKE-number of terms in system equations for moving fluid mesh (elasticity)





C Initialisation constants
      MAXVEC=NTOT

      IULAZ=57
      IIZLAZ=56
       INDFL=0
       IPENAL=0
       IALE=0
       MBAND=0
      ISRPS=1
      IDEBUG=0
      INDFOR=2
      INDIZL=0
      INDGRA=0
      IFORM=0
      IPRBR=0
      NGE=1

C==========================================================================
C==========================================================================
C Opening input file
      CALL OTVORF(ISRPS,IULAZ,IIZLAZ,IMEF,LSTFL,UNVFL,PAKNEU,IDUZIF)
C Opening output file
      CALL OTVIZF(ISRPS,IIZLAZ,IMEF,LSTFL,UNVFL,IDUZIF)
C==========================================================================
      CALL ZAGLAF(ISRPS,IIZLAZ)
C Reading heading of problem
      CALL INPU01(IULAZ,NASLOV,IIZLAZ,ISRPS)
C==========================================================================
C Reading indicator for format of data
      CALL INPU02(IULAZ,INDFOR,IIZLAZ,ISRPS)
C==========================================================================
      CALL OTVGRF(IMEF,LSTFL,UNVFL,IDUZIF,59)
      CALL OTVNEUF(IMEF,PAKNEU,IDUZIF,ISRPS,69)
C
C Reading basic data of problem
      CALL INPF03 (IULAZ,NPT,NSTAC,NPER,INDFL,IPENAL,IIZLAZ,ISRPS,IALE,
     &ITURB,IZID)
C==========================================================================
C Reading basic data of problem
      CALL INPU04(IULAZ,METOD,IFORM,MAXIT,EPSTA,EPSTR,NJRAP,
     &MBAND,IIZLAZ,ISRPS,IALFA)
C==========================================================================
C Reading indicator for executing of problem
      CALL INPU05 (IULAZ,IREST,IIZLAZ,ISRPS)
C==========================================================================
C Reading time steps
      CALL MEMORY(LVREME,LMAX,NPER,2,MAXVEC,IIZLAZ)
      CALL MEMORY(LNTABF,LMAX,NPER,1,MAXVEC,IIZLAZ)
      CALL PERIOF(A(LVREME),NPER,A(LNTABF),IULAZ,IIZLAZ,ISRPS,IBKOR)
C==========================================================================
C==========================================================================
C Reading numberation of NODES,ID-matrix, COORDinates of nodes
      CALL MEMORY(LID,LMAX,7*NPT,1,MAXVEC,IIZLAZ)
Cxxx
C     CALL CZINIT
C     CALL CZPROV(NPT,1,34)
Cxxx
C==========================================================================
C==========================================================================
c      CALL MEMORY(LCORD,LMAX,3*NPT,2,MAXVEC,IIZLAZ)
      CALL MEMORY(LCORD,LMAX,3*NPT,2,MAXVEC,IIZLAZ)

C privremeno izbaceno:
C      CALL MEMORY(LIDE,LMAX,2*NPT,1,MAXVEC,IIZLAZ)
      CALL MEMORY(LIDE,LMAX,0,2,MAXVEC,IIZLAZ)

      LID1=1
C       IF (INDFL.EQ.0) THEN
	   CALL MEMORY(LID1,LMAX,7*NPT,1,MAXVEC,IZLAZ)
C       ENDIF
      CALL ULAZF1(A(LID),A(LCORD),NPT,IULAZ,IIZLAZ,A(LIDE),INDFL,NEQF,
     &ISRPS,IPENAL,IALE,METOD,A(LID1),ITURB)

C	 call mae(A(LCORD),npt,16)
C	 call mae1(A(LCORD),npt,16)
C	 call mae2(A(LCORD),npt)
Cxxx
C      CALL ZASTIT
Cxxx
C==========================================================================
C Reading basic data about finite elements
      CALL INPU06(IULAZ,NETIP,NET,INDAX,NGET,NP2DMX,PENALT,IIZLAZ,ISRPS,
     +          IPENAL)
C=========================================================================
      CALL TNDIM(NDIM,NP2DMX,NDES,NETIP)
C
C==========================================================================
C Reading finite elements
      CALL MEMORY(LNEL,LMAX,(NDIM+1)*NET,1,MAXVEC,IIZLAZ)
	
C Ucitavanje elemenata i smestanje u matricu nel(9,*)
      CALL ULAZF2(A(LNEL),NET,NETIP,NDIM,IULAZ,IIZLAZ,ISRPS)
C	CALL CHECKJACOBIAN3D (A(LCORD),A(LNEL),NETIP,NET,NDIM,JACOB3D,
C	&IIZLAZ)

Cxxx
C      CALL BOUND5(A(LCORD),NPT,A(LID),A(LNEL),NDIM,NET)
c      CALL NODEADD(A(LCORD),NPT,A(LID),A(LNEL),NDIM,NET,IIZLAZ)
C      CALL BOUND55(A(LCORD),NPT,A(LID),A(LNEL),NDIM,NET)
C     CALL CZPROV(NET,1,34)
Cxxx      
C==========================================================================
C Reading material constants: 
C Density of fluid and Indicator for nonlinear dynamic viscosity
C Reading basic data about heat transfer if heat analysis is required
      CALL INPU07(IULAZ,INDFL,NANLK,NTABK,MAXTK,NANLC,NTABC,MAXTC,
     &GUSM,INDAMI,AKT,CC,FB2,FB3,AMI,BETA,TETAO,IIZLAZ,ISRPS)
C==========================================================================
C Reading number of prescribed values
      CALL INPU08(IULAZ,NUMZAD,MAXSIL,IIZLAZ,ISRPS)
C==========================================================================
      CALL MEMORY(LNZAD,LMAX,3*NUMZAD,1,MAXVEC,IIZLAZ)
      CALL MEMORY(LZADVR,LMAX,NUMZAD,2,MAXVEC,IIZLAZ)
C MANJI BAND
C      CALL NEWID1(A(LID),NEQF,NPT,PENALT,NETIP)
C VECI BAND-TACNIJI REZULTAT:
       IF (MBAND.EQ.1)
     &  CALL NEWIDF(A(LID),NEQF,NPT,PENALT,NETIP)
C Reading prescribed values
      CALL ULAZF3(A(LNZAD),A(LZADVR),NUMZAD,IULAZ,IIZLAZ,ISRPS)
Cxxx
C     CALL CZPROV(NUMZAD,1,34)
Cxxx      

C==========================================================================
C Reading initial values
      CALL INPU09(IULAZ,POCU,POCV,POCW,POCP,POCT,POCK,POCO,DELTAL,
     &IIZLAZ,ISRPS)
C========================================================================
CS Ucitavanje vremenskih funkcija
      CALL INPU10(IULAZ,NTABFT,MAXTFT,IIZLAZ,ISRPS)
CE Reading time functions
      CALL MEMORY(LTABF,LMAX,NTABFT*MAXTFT*2,2,MAXVEC,IIZLAZ)
      CALL MEMORY(LITFMA,LMAX,NTABFT,1,MAXVEC,IIZLAZ)
      CALL ULTAFF(A(LTABF),A(LITFMA),NPER,NTABFT,IULAZ,IIZLAZ,ISRPS)
C========================================================================
C Reading boundary conditions
      CALL MEMORY(LNGPSI,LMAX,MAXSIL*8,1,MAXVEC,IIZLAZ)
      CALL ULAZF4(A(LNEL),A(LNGPSI),MAXSIL,NDIM,IULAZ,IIZLAZ,ISRPS,
     &NETIP)
Cxxx
C     CALL CZPROV(MAXSIL,1,34)
Cxxx      
C==========================================================================
C Reading end of input data
      CALL INPU11(IULAZ)
C==========================================================================
C Closing output file
      CALL ZATVOF(IPRBR,INDIZL,INDGRA,IIZLAZ)
C==========================================================================
C       CALL CEVI(A(LCORD))
C==========================================================================
C ZA NEUTRAL FAJL, DZIGA
       CALL TGRMATF(1,NMATT,NETIP,69,NDIM,NASLOV,GAMA,1)
       CALL TGRAUKF(A(LCORD),A(LID),NPT,69)
C ZA ELEMENTE
       CALL TGRAF1(A(LNEL),NDIM,NET,NGE,59,NETIP,NDIM)
C     IF (NETIP.EQ.2) THEN 
C       CALL TGRAU2F(NEL,NDIM,NET,69,1)
C      ELSEIF (NETIP.EQ.3) THEN
C       CALL TGRAU3F(NEL,NDIM,NET,69,1)
C      ENDIF
C==========================================================================
       CALL TGRAFC(A(LCORD),NPT,59)
       CALL TGRBCF(A(LCORD),NPT,59,A(LID))
       CALL TGRAFF(A(LNEL),NDIM,NET,NGE,59,NETIP,NDIM)
C only for akira example, Sept 13, 2007
C       CALL SHELLWRITE(A(LNEL),NDIM,NET,NGE,59,NETIP,NDIM,A(LID),
C     &A(LCORD),NPT)
C      CALL SHELL2(A(LNEL),NET,2,NDIM,A(LID),A(LCORD),NPT,A(LCORD))
C	 STOP

C==========================================================================
      CALL RESTAF (IREST,IIZLAZ,ISRPS)

C OVDE SAMO PRIVREMENO ZA GRANICNE USLOVE OGRANICENJA ZIDOVA:
C for Schima example, merge of nodes
c       CALL genertube(IIZLAZ)
c       CALL generCURZIO(IIZLAZ)
C      CALL genercoating(IIZLAZ,A(LMAX))
c      CALL MEMORY(LNIZ,LMAX,NPT,1,MAXVEC,IIZLAZ)
c      CALL MEMORY(LNIZNEW,LMAX,NPT*50*10,1,MAXVEC,IIZLAZ)
c      CALL MEMORY(LCCCORD,LMAX,NPT*50*10*3,2,MAXVEC,IIZLAZ)
c      CALL MEMORY(LNEL1,LMAX,NET*50,1,MAXVEC,IIZLAZ)
c      CALL MEMORY(LIID,LMAX,NPT*50*10*6,1,MAXVEC,IIZLAZ)
c      CALL MEMORY(LNIZF,LMAX,NPT*50*10,1,MAXVEC,IIZLAZ)
c      CALL TETGEN41(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN4(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN4_8(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN4_6(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN41_6(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN10(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN10_6(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN11_6(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TetGen8FluxFace_6(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),
c     &A(LNIZNEW),A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN10_B5(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL TETGEN10_B55(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))
c      CALL OctGen21(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNIZNEW),
c     &A(LNEL1),A(LCCCORD),A(LIID),A(LNIZF))



c      CALL MERGEN(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),
C     *A(LNIZ),A(LNIZNEW),IIZLAZ,A(LNZAD),A(LZADVR),NUMZAD)
c       CALL BOUN33(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID))

C       CALL NIKOLA(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID))
c	stop
c       CALL BOUND3(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),
c     &A(LNZAD),A(LZADVR),NUMZAD)
c       CALL BOUNDT3(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),
c     &A(LNZAD),A(LZADVR),NUMZAD)
c       CALL BOUND3DFACE(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID))
c       CALL FLUX3DFACE(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID))
C       CALL BOUN3D(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID))
c       CALL BOUN2D(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID))
c 	stop
C       CALL BOUND4(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LMAX),
C     &A(LNZAD),A(LZADVR),NUMZAD)
C      CALL SHELLS(A(LCORD),A(LNEL),NDIM,NPT,NET,A(LID),A(LNZAD),
C     &A(LZADVR),NUMZAD)
C==========================================================================
C==========================================================================
C Writing transform ID-matrix at output LST file
C==========================================================================
C Definition of MAXA-vector
      CALL MEMORY(LMAXA,LMAX,NEQF+1,1,MAXVEC,IIZLAZ)
      CALL MEMORY(LMHT,LMAX,NEQF+1,1,MAXVEC,IIZLAZ)
      IF (IALE.EQ.1.AND.METOD.EQ.2) THEN
       CALL MEMORY(LIDALE,LMAX,NPT*3,1,MAXVEC,IZLAZ)
       CALL IDALF(A(LIDALE),A(LID),NPT,NEQF)
      ELSE
       LIDALE=LMAX
      ENDIF
	IF (ISYMMS.EQ.1) THEN
       CALL MAXATE(A(LMAXA),A(LMHT),A(LID),A(LNEL),NET,NDIM,NEQF,NWKF,7)
	ELSE
c	IF (IPENAL.EQ.0.AND.NETIP.EQ.3.AND.NDIM.EQ.8.and.metod.eq.4) THEN
c       CALL MAXATF(A(LMAXA),A(LMHT),A(LID),A(LNEL),NET,NDIM+1,NEQF,
c     &NWKF,6,NDIM+1)
c	ELSE
       CALL MAXATF(A(LMAXA),A(LMHT),A(LID),A(LNEL),NET,NDIM,NEQF,
     &NWKF,7,NDIM+1,iizlaz)
c	 ENDIF
	 ENDIF
C==========================================================================
C PRIVREMENO IZBACENO:
C      CALL ELASTO(A,LMAX,NPT,LIDE,LNZADE,NUMZE,NEQE,NWKE,LMAXAE,LTT1E,
C     &LNEL,NET,NDIM,LM2,LZADVE,LVECTJ)
C==========================================================================
      CALL MEMORY(LTT1,LMAX,NEQF,2,MAXVEC,IIZLAZ)
      CALL MEMORY(LTT10,LMAX,NEQF,2,MAXVEC,IIZLAZ)
C      CALL MEMORY(LSKEF,LMAX,NDES*NDES,2,MAXVEC,IIZLAZ)

      CALL TIDPRI(PENALT,NETIP,IDPRIT,NDIM)

      CALL MEMORY(LCCORD,LMAX,3*NPT,2,MAXVEC,IIZLAZ)

C MEMORY FOR SPARSE MATRIX:  
	 ISPARC=0
       IF (ISPARC.EQ.1) THEN
        CALL MEMORY(LSPAR1,LMAX,9*NEQF,2,MAXVEC,IIZLAZ)
        CALL SPARS1(A(LSPAR1),NEQF,A(LNEL),A(LID),NDIM,NPT,NET)
        CALL MSPARC(A(LSPAR1),NKC,NEQF,NDIM)
        CALL MEMORY(LAKc,LMAX,NKC,2,MAXVEC,IIZLAZ)
        CALL MEMORY(LIKc,LMAX,NKC,1,MAXVEC,IIZLAZ)
        CALL MEMORY(LROW,LMAX,NEQF,2,MAXVEC,IIZLAZ)
       ELSE
        LSPAR1=LMAX
        LAKc=LMAX
        LIKc=LMAX
        LROW=LMAX
       ENDIF

       CALL MEMORY(LINDEL,LMAX,NPT,1,MAXVEC,IIZLAZ)
       CALL MEMORY(LSPSIL,LMAX,NPT*NETIP,2,MAXVEC,IIZLAZ)
       CALL MEMORY(LGNODE,LMAX,NPT*7*2,2,MAXVEC,IIZLAZ)


C Initialisation of values
      CALL INITIA(A(LCORD),A(LCCORD),A(LNEL),A(LTT1),A(LINDEL),
     &A(LSPSIL),NET,NPT,NEQF,NDIM,KKORAK,VVREME,NETIP,A(LGNODE),A(LID))



C==========================================================================

      CALL NAMEF(IDATF,DATF,LMAX,NTOT,LTT1,LNEL,LID,LNZAD,LZADVR,
     &LNGPSI,LMAXA,LCORD,LTT10,LVREME,LTABF,LITFMA,LIDALE,
     &LCCORD,LAKc,LIKc,LROW,LSPAR1,LMHT,LINDEL,LNTABF,LSPSIL,LGNODE,
     &MAXIT,NETIP,IULAZ,IIZLAZ,INDAMI,NUMZAD,NPT,NDIM,MAXSIL,NEQF,
     &NWKF,NET,NPER,NTABFT,NDES,IDPRIT,IFORM,NSTAC,INDAX,KKORAK,
     &IPENAL,GUSM,CC,AKT,EPSTR,AMI,BETA,TETAO,FB2,FB3,PENALT,VVREME,
     &ITURB,IZID)

      END

C==========================================================================

      SUBROUTINE NAMEF(I,R,LMAX,NTOT,LTT1,LNEL,LID,LNZAD,LZADVR,
     &LNGPSI,LMAXA,LCORD,LTT10,LVREME,LTABF,LITFMA,LIDALE,
     &LCCORD,LAKc,LIKc,LROW,LSPAR1,LMHT,LINDEL,LNTABF,LSPSIL,LGNODE,
     &MAXIT,NETIP,IULAZ,IIZLAZ,INDAMI,NUMZAD,NPT,NDIM,MAXSIL,NEQF,
     &NWKF,NET,NPER,NTABFT,NDES,IDPRIT,IFORM,NSTAC,INDAX,KKORAK,
     &IPENAL,GUSM,CC,AKT,EPSTR,AMI,BETA,TETAO,FB2,FB3,PENALT,VVREME,
     &ITURB,IZID)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION I(*),R(*)
C
CE Subroutine NAMEF is used for memory all global variables
C

C========================================================================
C DEFINITION INTEGER DATA
      I(1)=LVREME
      I(2)=LNTABF
      I(3)=LID
      I(4)=LCORD
      I(5)=LNEL
      I(6)=LNZAD
      I(7)=LZADVR
      I(8)=LTABF
      I(9)=LITFMA
      I(10)=LNGPSI
      I(11)=LMAXA
      I(12)=LMHT
      I(13)=LTT1
      I(14)=LTT10
      I(15)=LCCORD
      I(16)=LSPAR1
      I(17)=LAKc
      I(18)=LIKc
      I(19)=LROW
      I(20)=LINDEL
      I(21)=LSPSIL
      I(22)=LMAX
      I(23)=NTOT
      I(24)=MAXIT
      I(25)=NETIP
      I(26)=IULAZ
      I(27)=IIZLAZ
      I(28)=INDAMI
      I(29)=NUMZAD
      I(30)=NPT
      I(31)=NDIM
      I(32)=MAXSIL
      I(33)=NEQF
      I(34)=NWKF
      I(35)=NET
      I(36)=NPER
      I(37)=NTABFT
      I(38)=NDES
      I(39)=IDPRIT
      I(40)=IFORM
      I(41)=NSTAC
      I(42)=INDAX
      I(43)=KKORAK
      I(44)=IPENAL
      I(45)=LIDALE
      I(46)=LGNODE
      I(47)=ITURB
      I(48)=IZID
C========================================================================
C DEFINITION REAL DATA
      R(1)=GUSM
      R(2)=CC
      R(3)=AKT
      R(4)=EPSTR
      R(5)=AMI
      R(6)=BETA
      R(7)=TETAO
      R(8)=FB2
      R(9)=FB3
      R(10)=PENALT
      R(11)=VVREME
C==========================================================================

      END
C=======================================================================
C=======================================================================
      SUBROUTINE ULAZF1(ID,CORD,NPT,IULAZ,IIZLAZ,IDE,INDFL,NEQF,ISRPS,
     &IPENAL,IALE,METOD,ID1,ITURB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C      COMMON /SAMOFL/ INDFL
C      COMMON /ULAZNI/ IULAZ,IIZLAZ

C
CE Subroutine ULAZF1 is used for reading input data for nodes,constraints 
CE and their coordinates
C
      CHARACTER*250 ACOZ
      DIMENSION ID(7,*),IDE(2,*),CORD(3,*),ID1(7,*)

       KK=0
       KJ=0

      DO 10 I=1,NPT
      CALL ISPITF(ACOZ,IULAZ)
c      READ(ACOZ,1005) N,(ID(II,N),II=1,5),(CORD(J,N),J=1,3)
C      READ(ACOZ,1015) N,(ID(II,N),II=1,6),(CORD(J,N),J=1,3)
C      READ(ACOZ,1006) N,(ID(II,N),II=1,5),(CORD(J,N),J=1,3)
C        READ(ACOZ,*) N,(ID(II,N),II=1,5),(CORD(J,N),J=1,3)
C=========================================================================
      IF (ITURB.EQ.1) THEN
      READ(ACOZ,1005) N,(ID(II,N),II=1,7),(CORD(J,N),J=1,3)
       ELSEIF (ITURB.EQ.0.OR.ITURB.EQ.2) THEN
      READ(ACOZ,1015) N,(ID(II,N),II=1,5),(CORD(J,N),J=1,3)
      ENDIF
C=========================================================================
C 1005 FORMAT(I10,7(I5),F14.6,2F13.6,2I2)
C  1015 FORMAT(I5,6(3X,I2),3F10.6,2I2)
       IF (IPENAL.EQ.1) ID(4,N)=1

      IF (INDFL.EQ.1) THEN
        ID(5,N)=1 
      ELSEIF (INDFL.EQ.2) THEN
        DO JJ=1,4  
         ID(JJ,N)=1 
        ENDDO
      ENDIF
 
       ID1(1,N)=ID(1,N)
       ID1(2,N)=ID(2,N)
       ID1(3,N)=ID(3,N)
       ID1(4,N)=ID(4,N)
       ID1(5,N)=ID(5,N)
       ID1(6,N)=ID(6,N)
       ID1(7,N)=ID(7,N)

      IF (ITURB.EQ.1) THEN
 
 
      DO JJ=1,7
       IF (ID(JJ,N).EQ.0) THEN
        KJ=KJ+1
        ID(JJ,N)=KJ
       ELSE
        ID(JJ,N)=0
       ENDIF
      ENDDO
      
      ELSE
 
      DO JJ=1,5
       IF (ID(JJ,N).EQ.0) THEN
        KK=KK+1
        ID(JJ,N)=KK
       ELSE
        ID(JJ,N)=0
       ENDIF
      ENDDO
      
      ENDIF
      
   10 CONTINUE


      IF (ITURB.EQ.1) THEN
       NEQF=KJ
      ELSE
       NEQF=KK
      ENDIF

C      IF (ITURB.EQ.1) THEN
 1005 FORMAT(I10,7(3X,I2),3F13.6,2I2)
 6005 FORMAT(I10,7I5,2F13.6,F14.6,2I2)
 1010 FORMAT(I5,7(I5),3F10.6,2I2)
 1012 FORMAT(I10,7(I7),3(x,F13.5),2I2)
C      ELSEIF (ITURB.EQ.0.OR.ITURB.EQ.2) THEN
 1015 FORMAT(I10,5(I5),3F13.6,2I2)
 6015 FORMAT(I10,5I5,2F13.6,F14.6,2I2)
 1011 FORMAT(I5,5(I5),3F10.6,2I2)
 1013 FORMAT(I10,5(I5),3(x,F13.5),2I2)
C      ENDIF
C 1015 FORMAT(I5,6(3X,I2),3F10.6,2I2)
C 1005 FORMAT(I5,5(3X,I2),3F10.6,2I2)
C 1006 FORMAT(I5,5(3X,I2),3F11.6,2I2)
C 1010 FORMAT(I5,5(I5),3F10.6,2I2)
C 1012 FORMAT(5X,I5,6(I5),3F10.6,2I2)
C 1012 FORMAT(5X,I10,6(I5),3(x,1pe10.3),2I2)
C 1013 FORMAT(I10,5(I5),3(x,1pe13.5),2I2)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,2002)
 2000 FORMAT(6X,'P O D A C I    O    C V O R O V I M A'
     1/6X,37('-')///
     16X,'REDNI BROJEVI,   BROJEVI JEDNACINA,   KOORDINATE CVOROVA'/)
 2002 FORMAT(6X,'D A T A   F O R   N O D A L   P O I N T   D A T A'
     1/6X,49('-')///
     16X,'NUM. OF NODE,   NUM. OF EQUATIONS,  COORD. OF NODES'/)

C      DO 12 I=1,NPT
C      WRITE(IIZLAZ,1013) I,(ID(II,I),II=1,5),(CORD(J,I),J=1,3)
C    &,IDE(1,I),IDE(2,I)
C   12 CONTINUE
      IF (ITURB.EQ.1) THEN
      DO 12 I=1,NPT
      WRITE(IIZLAZ,1012) I,(ID(II,I),II=1,7),(CORD(J,I),J=1,3)
   12 CONTINUE
      ELSEIF (ITURB.EQ.0.OR.ITURB.EQ.2) THEN
      DO 13 I=1,NPT
      WRITE(IIZLAZ,1013) I,(ID(II,I),II=1,5),(CORD(J,I),J=1,3)
   13 CONTINUE
   
      ENDIF
  
      END
C==========================================================================
C==========================================================================
      SUBROUTINE ULAZF2(NEL,NET,NETIP,NDIM,IULAZ,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /TIPEL/ NP2DMX
C      COMMON /TIPELM/ NETIP
C
CE Subroutine ULAZF2 is used for reading input data for finite element
C

      CHARACTER*250ACOZ
      DIMENSION NEL(NDIM+1,*)


      DO 20 J=1,NET
      IF (NETIP.EQ.2) THEN
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1007)NN,(NEL(I,NN),I=1,5)
c      READ(ACOZ,1017)NN,(NEL(I,NN),I=1,5)
C      READ(ACOZ,*) NN,(NEL(I,NN),I=1,5)
       IF (NDIM.GT.4) THEN
      CALL ISPITF(ACOZ,IULAZ)
         READ(ACOZ,1007)  (NEL(I,NN),I=5,9)
       ENDIF
      ELSEIF (NETIP.EQ.3) THEN
      CALL ISPITF(ACOZ,IULAZ)
c for tethra elements
c      if (ndim.eq.4) READ(ACOZ,1007)NN,(NEL(I,NN),I=1,5)
      if (ndim.eq.4) READ(ACOZ,1017)NN,(NEL(I,NN),I=1,5)
      if (ndim.eq.5) READ(ACOZ,1007)NN,(NEL(I,NN),I=1,5)
      if (ndim.eq.6) READ(ACOZ,1017)NN,(NEL(I,NN),I=1,6)
 
c      if (ndim.eq.8) READ(ACOZ,1007)NN,(NEL(I,NN),I=1,9)
      if (ndim.eq.8) READ(ACOZ,1017)NN,(NEL(I,NN),I=1,9)

      if (ndim.eq.10) READ(ACOZ,1017)NN,(NEL(I,NN),I=1,10)
      if (ndim.eq.11) READ(ACOZ,1017)NN,(NEL(I,NN),I=1,11)
c      if (ndim.eq.10) READ(ACOZ,1007)NN,(NEL(I,NN),I=1,10)

C      READ(ACOZ,1017)NN,(NEL(I,NN),I=1,9)
       IF (NDIM.eq.21) THEN
	   READ(ACOZ,1017)NN,(NEL(I,NN),I=1,8)
         CALL ISPITF(ACOZ,IULAZ)
         READ(ACOZ,1017)  (NEL(I,NN),I=9,21)
       ENDIF
      ENDIF

c      n1=nel(2,nn)
c	nel(2,nn)=nel(3,nn)
c	nel(3,nn)=n1

   20 CONTINUE
 1007 FORMAT(13I5)
 1017 FORMAT(13I10)


      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,2002)
 2000 FORMAT(6X,'P O D A C I    O    E L E M E N T I M A'
     1/6X,39('-')///
     16X,'REDNI BROJEVI ELEMENATA,   NUMERICIJA CVOROVA PO ELEMENTU'/)
 2002 FORMAT(6X,'D A T A   F O R   E L E M E N T S'
     1/6X,49('-')///
     16X,'ELEMENT NUMBER, NUMBER OF NODES PER ELEMENT'/)



      DO 22 J=1,NET
      IF (NETIP.EQ.2) THEN
      WRITE(IIZLAZ,1007) J,(NEL(I,J),I=1,4)
        IF (NDIM.GT.4) THEN
          WRITE(IIZLAZ,1007)  (NEL(I,J),I=5,9)
        ENDIF
      ELSEIF (NETIP.EQ.3) THEN
      WRITE(IIZLAZ,1017) J,(NEL(I,J),I=1,ndim)
c        IF (NDIM.GT.8) THEN
c          WRITE(IIZLAZ,1007)  (NEL(I,J),I=9,ndim)
c        ENDIF
      ENDIF
c      WRITE(IIZLAZ,100)
   22 CONTINUE
 100  FORMAT(80(' '))

    
c      write(iizlaz,*)'elements'
c      do i=1,net
c        if (nel(5,i).ne.0) then
c	    write(iizlaz,'(7i10)')i,(nel(j,i),j=1,3),(nel(j,i),j=5,7)
c         else
c	    write(iizlaz,'(7i10)')i,(nel(j,i),j=1,4)
c	   endif
c	enddo
c	stop

      RETURN
      END

C==========================================================================


C==========================================================================
C==========================================================================
      SUBROUTINE ULAZF4(NEL,NGPSIL,MAXSIL,NDIM,IULAZ,IIZLAZ,ISRPS,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)


      CHARACTER*250ACOZ
      DIMENSION NGPSIL(8,*),NEL(NDIM+1,*)

C
CE Subroutine ULAZF4 is used for reading input data for surface 
CE boundary tractions
C

      IF (MAXSIL.EQ.0) RETURN

      DO 100 I=1,MAXSIL
c	if (i.gt.49500) stop 'stao ovde'
      CALL ISPITF(ACOZ,IULAZ)
      IF (NETIP.EQ.2) THEN
      READ(ACOZ,1015) NGPSIL(1,I),NGPSIL(2,I),NGPSIL(3,I)
     1,NGPSIL(4,I),NGPSIL(5,I)
      ELSE IF (NETIP.EQ.3) THEN
      READ(ACOZ,1015) NGPSIL(1,I),NGPSIL(2,I),NGPSIL(3,I)
     1,NGPSIL(4,I),NGPSIL(5,I),NGPSIL(6,I),NGPSIL(7,I),NGPSIL(8,I)
      ENDIF
C      IEL=NGPSIL(1,I)
C      N1=NGPSIL(2,I)
C      N2=NGPSIL(3,I)
C     IF ((NEL(1,IEL).EQ.N1 .AND. NEL(2,IEL).EQ.N2) .OR.
C    1    (NEL(2,IEL).EQ.N1 .AND. NEL(1,IEL).EQ.N2) .OR.
C    1    (NEL(2,IEL).EQ.N1 .AND. NEL(3,IEL).EQ.N2) .OR.
C    1    (NEL(3,IEL).EQ.N1 .AND. NEL(2,IEL).EQ.N2) .OR.
C    1    (NEL(3,IEL).EQ.N1 .AND. NEL(4,IEL).EQ.N2) .OR.
C    1    (NEL(4,IEL).EQ.N1 .AND. NEL(3,IEL).EQ.N2) .OR.
C    1    (NEL(1,IEL).EQ.N1 .AND. NEL(4,IEL).EQ.N2) .OR.
C    1    (NEL(4,IEL).EQ.N1 .AND. NEL(1,IEL).EQ.N2)) THEN 
C         ELSE
C            WRITE(IIZLAZ,2000) IEL,N1,N2
C            STOP
C         ENDIF
  100 CONTINUE

 1015 FORMAT(8I5)
 1017 FORMAT(8I10)
c 1017 FORMAT(8I5)
C 2000 FORMAT (//
C     111X,'GRESKA U ULAZNIM PODACIMA KOD UCITAVANJA GRANICNIH USLOVA ZA'
C     1/11X,'2D-ELEMENTE'// 
C     111X,'IVICU ELEMENTA BROJ= ',I5,' NE CINE CVOROVI ',2I5)
C
      IF (ISRPS.EQ.0.AND.NETIP.EQ.2)
     *WRITE(IIZLAZ,1103)
      IF (ISRPS.EQ.1.AND.NETIP.EQ.2)
     *WRITE(IIZLAZ,3103)

      IF (ISRPS.EQ.0.AND.NETIP.EQ.3)
     *WRITE(IIZLAZ,1303)
      IF (ISRPS.EQ.1.AND.NETIP.EQ.3)
     *WRITE(IIZLAZ,3303)

 1103 FORMAT(6X,'P O D A C I   O   G R A N I C N I M   U S L O V I M A'/
     16X,53('-')//
     16X,'P O V R S I N S K E   S I L E'/6X,20('-')//     
     111X,'ELEMENT,CVOR1,CVOR2,IND, PRAVCA DELOVANJA POVRSINSKE SILE'//
     116X,'EQ.0; NE POSTOJI POVRSINSKA SILA U DATOM PRAVCU'/ 
     116X,'EQ.1; POSTOJI POVRSINSKA SILA U DATOM PRAVCU'///) 

 3103 FORMAT(6X,'D A T A   A B O U T   B O U N D A R Y   C O N D I T I O
     1 N S'/6X,59('-')//
     16X,'S U R F A C E   T R A C T I O N S'/6X,33('-')//
     111X,'ELEM., NODE1, NODE2,IND. FOR SURFACE TRACTION DIRECTION'/
     116X,'EQ.0; NO SURFACE TRACTION IN PRESCRIBED DIRECTION'/
     116X,'EQ.1; SURFACE TRACTION IN PRESCRIBED DIRECTION'///)

 1303 FORMAT(6X,'P O D A C I   O   G R A N I C N I M   U S L O V I M A'/
     16X,53('-')//
     16X,'P O V R S I N S K E   S I L E'/6X,20('-')//     
     111X,'ELEMENT,CV,CV2,CV3,CV4,IND. FOR SURFACE TRACTION DIRECTION'//
     116X,'EQ.0; NE POSTOJI POVRSINSKA SILA U DATOM PRAVCU'/ 
     116X,'EQ.1; POSTOJI POVRSINSKA SILA U DATOM PRAVCU'///) 

 3303 FORMAT(6X,'D A T A   A B O U T   B O U N D A R Y   C O N D I T I O
     1 N S'/6X,59('-')//
     16X,'S U R F A C E   T R A C T I O N S'/6X,33('-')//
     111X,'ELEM., N1, N2,N3,N4,IND. FOR SURFACE TRACTION DIRECTION'/
     116X,'EQ.0; NO SURFACE TRACTION IN PRESCRIBED DIRECTION'/
     116X,'EQ.1; SURFACE TRACTION IN PRESCRIBED DIRECTION'///)

      DO 72 I=1,MAXSIL
      IF (NETIP.EQ.2) THEN
      WRITE(IIZLAZ,1016) NGPSIL(1,I),NGPSIL(2,I),NGPSIL(3,I)
     1,NGPSIL(4,I),NGPSIL(5,I)
      ELSEIF (NETIP.EQ.3) THEN
      WRITE(IIZLAZ,1016) NGPSIL(1,I),NGPSIL(2,I),NGPSIL(3,I)
     1,NGPSIL(4,I),NGPSIL(5,I),NGPSIL(6,I),NGPSIL(7,I),NGPSIL(8,I)
      ENDIF
  72  CONTINUE
 1016 FORMAT(I10,7I10)
      END
C==========================================================================
C==========================================================================
C======================================================================
C======================================================================
C==========================================================================
C=======================================================================
      SUBROUTINE INPU06(IULAZ,NETIP,NET,INDAX,NGET,NP2DMX,PENALT,
     &IIZLAZ,ISRPS,IPENAL)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA
      COMMON /BIFURC/ DIAM,ALFA,BETA,VEL1,VEL2,PHASE1,PHASE2,IDIV1,
     &IDIV2,PENAL,AL,AL1,AL2
C      COMMON /STORER/ ISINCR,AITIME,AFTIME
      COMMON /CHECKJACOB3D/ JACOB3D

      CHARACTER*250ACOZ

C
CE Subroutine INPU06 is used for reading basic data for finite elements
C

      NMATT=1
      CALL ISPITF(ACOZ,IULAZ)
C      READ(ACOZ,1002) NETIP,NET,INDAX,IAKIRA,JACOB3D
      READ(ACOZ,1012) NETIP,NET,INDAX,IAKIRA,JACOB3D


      IF (IAKIRA.EQ.3) THEN
       CALL ISPITF(ACOZ,IULAZ)
        READ(ACOZ,2004) PERIO1,NUMST,VALST,TOL,PENAL
       CALL ISPITF(ACOZ,IULAZ)
        READ(ACOZ,2005) DIAM,ALFA,BETA,IDIV1,IDIV2,AL,AL1,AL2
       CALL ISPITF(ACOZ,IULAZ)
        READ(ACOZ,2003) VEL1,VEL2,PHASE1,PHASE2
      ENDIF

 2004  FORMAT(F10.3,I5,2F10.7,E10.3)
 2005  FORMAT(3F10.3,2I5,3F10.3)
 2003  FORMAT(4F10.3)
     
      IF (IAKIRA.EQ.1.OR.IAKIRA.EQ.2.OR.IAKIRA.EQ.4.OR.IAKIRA.EQ.5) THEN
       CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1005) PERIO1,SVEXC,RD,RA,ALEN,GAMA,DCEXP,MOVEW
       CALL ISPITF(ACOZ,IULAZ)
C       READ(ACOZ,1003) MESH,NUMALV,DISTAL,DIS1AL,RE,NUMST,VALST,TOL,
       READ(ACOZ,1033) MESH,NUMALV,DISTAL,DIS1AL,RE,NUMST,VALST,TOL,
     &NBSTAC,NPRVEL
      ENDIF

C      IF (IAKIRA.NE.0) THEN
C        IF (DABS(AFTIME).LE.1.D-10) AFTIME=NUMST*VALST*1.001D0
C      ENDIF

      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2011) NETIP,NET,INDAX
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6011) NETIP,NET,INDAX
 2011 FORMAT(///6X,'O S N O V N I   P O D A C I   O   E L E M E N T I M
     1 A'/6X,53('-')///
     111X,'TIP KONACNOG ELEMENTA .......................... NETIP =',I5/
     116X,'EQ.2; IZOPARAM. 2D ELEMENT U RAVNI X-Y'/
     116X,'EQ.3; IZOPARAM. 3D ELEMENT U PROSTORU'//
     211X,'BROJ KONACNIH ELEMENATA .......................... NET =',
     1I10/
     311X,'INDIKATOR ZA OSNOSIMETRICNE ELEMENTE ........... INDAX =',I5/
     116X,'EQ.0; 2D ELEMENT U RAVNI X-Y ILI 3D ELEMENT U XYZ PROSTORU'//
     116X,'EQ.1; 2D OSNOSIMETRICNI ELEMENT U XY RAVNI'///)

 6011 FORMAT(///6X,'B A S I C   D A T A   F O R    E L E M E N T S'
     1/6X,53('-')///
     111X,'ELEMENT TYPE ................................... NETIP =',I5/
     116X,'EQ.2; ISOPARAM. 2D ELEMENT IN PLANE X-Y'/
     116X,'EQ.3; ISOPARAM. 3D ELEMENT IN XYZ SPACE'//
     211X,'NUMBER OF ELEMENTS ............................... NET =',
     1I10/
     3/11X,'INDICATOR FOR AXISYMMETRIC ELEMENTS ............ INDAX =',I5
     1/16X,'EQ.0; 2D IN PLANE X-Y OR 3D IN XYZ SPACE'/
     116X,'EQ.1; AXISYMMETRIC IN PLANE XY'///)



      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1006) NMAT2D,NMAT2D,NP2DMX,PENALT
 1006 FORMAT(3I5,D10.3)
C Only for Burger example this line is with C comment
C      IF (IPENAL.EQ.1.AND.DABS(PENALT).LT.1.D-5) PENALT=1.D9

      IF(NMAT2D.GT.1) MAT2D=0
      IF(NMATT.EQ.1) MAT2D=1
      IF(NMATT.EQ.1) NMAT2D=1
      IF(NP2DMX.EQ.0) NP2DMX=4
         IF(ISRPS.EQ.0)
     *   WRITE(IIZLAZ,3002)NMAT2D,NMAT2D,NP2DMX,PENALT
         IF(ISRPS.EQ.1)
     *   WRITE(IIZLAZ,6002)NMAT2D,NMAT2D,NP2DMX,PENALT
 3002 FORMAT(6X,'U L A Z N I  P O D A C I  Z A   E L E M E N T E'
     1/6X,47('-')///
     111X,'UKUPAN BROJ RAZLICITIH MATERIJALA ............. NMAT2D =',I5/
     116X,'EQ.0; NMAT2D.EQ.1'//
     211X,'MATERIJAL BROJ ................................ NMAT2D =',I5/
     216X,'EQ.0; NMAT2D.EQ.1'//
     311X,'MAKSIMALAN BROJ CVOROVA PO ELEMENTU ........... NP2DMX =',I5/
     316X,'EQ.0; POSTAJE "4" (ZA 4 CVORA)'/
     316X,'EQ.4; CETVOROCVORNI 2D ELEMENT'/
     316X,'EQ.9; DEVETOCVORNI 2D ELEMENT'//
     311X,'PENALTI FAKTOR ................................. PENALT =',
     3  D10.3/16X,'EQ.0; POSTAJE 10^9'///)
 6002 FORMAT(6X,'I N P U T   D A T A   F O R   E L E M E N T S'
     1/6X,45('-')///
     111X,'TOTAL NUMBER OF MATERIAL ...................... NMAT2D =',I5/
     116X,'EQ.0; NMAT2D.EQ.1'//
     211X,'NUMBER OF MATERIAL ............................ NMAT2D =',I5/
     216X,'EQ.0; NMAT2D.EQ.1'//
     311X,'MAXIMUM NUMBER OF NODES PER ELEMENT ........... NP2DMX =',I5/
     316X,'EQ.0; NP2DMX.EQ.4'/
     316X,'EQ.4; FOUR NODES PER 2D ELEMENT'/
     316X,'EQ.9; NINE NODES PER 2D ELEMENT'/
     316X,'EQ.8; EIGHT NODES PER 3D ELEMENT'/
     316X,'EQ.21; TWENTY-ONE NODES PER 3D ELEMENT'//
     311X,'PENALTY FACTOR ................................ PENALT =',
     3D10.3/16X,'EQ.0; SET BY DEFAULT TO 10^9 '///)

 1002 FORMAT(8I5)
 1012 FORMAT(I5,I10,3I5)
 1005 FORMAT(F10.3,6F10.7,I5)
 1003 FORMAT(2I5,3F10.7,I5,2F10.7,2I5)
 1033 FORMAT(2I5,2F10.7,F10.3,I5,2F10.7,2I5)


      END
C=======================================================================
C=======================================================================
      SUBROUTINE INPU07(IULAZ,INDFL,NANLK,NTABK,MAXTK,NANLC,NTABC,MAXTC,
     &GUSM,INDAMI,AKT,CC,FB2,FB3,AMI,BETA,TETAO,IIZLAZ,ISRPS)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250ACOZ
      common /orbital/ ORBRAD,ORBFRE
      common /difuzija/ dd(4)
C
CE Subroutine INPU07 is used for reading material constants
C

C INITIALISATION VALUES
      NANLK=0 
      NTABK=0 
      MAXTK=0 
      NANLC=0 
      NTABC=0 
      MAXTC=0 
      INDAMI=0

      GUSM=0.D0
      AKT=0.D0
      CC=0.D0
      FB2=0.D0
      FB3=0.D0
      AMI=0.D0
      TETAO=0.D0
     
      

C     IF (INDFL.EQ.0.OR.INDFL.EQ.2) THEN
C Reading basic data about heat transfer if heat analysis is required
C     CALL ISPITF(ACOZ,IULAZ)
C     READ(ACOZ,1002) NANLK,NTABK,MAXTK,NANLC,NTABC,MAXTC
C     IF(ISRPS.EQ.0)
C    *WRITE(IIZLAZ,2005) NANLK,NTABK,MAXTK
C     IF(ISRPS.EQ.1)
C    *WRITE(IIZLAZ,6005) NANLK,NTABK,MAXTK
C2005 FORMAT(//6X,' OSNOVNI PODACI O ZAVISNOSTI MATERIJALNIH KONSTANTI
C    1 OD TEMPERATURE'/7X,58('-')///
C    111X,'ZAVISNOST KOEFICIJENATA PROVODJENJA OD TEMPERATURE:'/
C    111X,'BROJ ANALITICKI ZADATIH FUNKCIJA ............ NANLK =',I5/
C    116X,'EQ.0; SVE FUNKCIJE ZADATE TABELARNO'///
C    211X,'BROJ TABELARNO ZADATIH FUNKCIJA ............. NTABK =',I5/
C    216X,'EQ.0; SVE FUNKCIJE ZADATE ANALITICKI'///
C    311X,'MAX. BROJ TACAKA ZA DEFINISANJE KRIVE AK(T) . MAXTK =',I5/
C    316X,'EQ.0; SVE FUNKCIJE ZADATE ANALITICKI'/
C    316X,'EQ.1; KOEFICIJENTI PROVODJENJA SU LINEARNI'/
C    316X,'GT.1; KOEFICIJENTI PROVODJENJA SU NELINEARNI')
C6005 FORMAT(6X,' BASIC DATA ABOUT TEMPERATURE DEPENDENCE MATERIAL',
C    1' CONSTANTS'/7X,58('-')///
C    116X,'TEMPERATURE DEPENDENCE CONDUCTION COEFICIENT:'///
C    111X,'NUMBER OF FUNCTION GIVEN IN ANALITICAL FORM .... NANLK =',I5/
C    116X,'EQ.0; ALL FUNCTION ARE GIVEN AS TABLE'///
C    211X,'NUMBER OF FUNCTION GIVEN AS TABLE .............. NTABK =',I5/
C    216X,'EQ.0; ALL FUNCTION ARE GIVEN IN ANALITICAL FORM'///
C    311X,'MAX. NO. OF POINTS USED FOR DEFIN. FUNCTION .... MAXTK =',I5/
C    316X,'EQ.0; ALL FUNCT. ARE GIVEN IN ANALITICAL FORM'/
C    316X,'EQ.1; CONDUCTION COEFICIENTS ARE LINEAR'/
C    316X,'GT.1; CONDUCTION COEFICIENTS ARE NONLINEAR')
C     IF(ISRPS.EQ.0)
C    *WRITE(IIZLAZ,2225) NANLC,NTABC,MAXTC
C     IF(ISRPS.EQ.1)
C    *WRITE(IIZLAZ,6225) NANLC,NTABC,MAXTC
C2225 FORMAT(////
C    411X,'ZAVISNOST KOEFICIJENATA SPECIFICNE TOPLOTE OD TEMPERATURE:'/
C    411X,'BROJ ANALITICKI ZADATIH FUNKCIJA ............ NANLC =',I5/
C    416X,'EQ.0; SVE FUNKCIJE ZADATE TABELARNO'///
C    511X,'BROJ TABELARNO ZADATIH FUNKCIJA ............. NTABC =',I5/
C    516X,'EQ.0; SVE FUNKCIJE ZADATE ANALITICKI'///
C    611X,'MAX. BROJ TACAKA ZA DEFINISANJE KRIVE CP(T) . MAXTC =',I5/
C    616X,'EQ.0; SVE FUNKCIJE ZADATE ANALITICKI'/
C    616X,'EQ.1; KOEFICIJENTI SPEC. TOPLOTE SU LINEARNI'/
C    616X,'GT.1; KOEFICIJENTI SPEC. TOPLOTE SU NELINEARNI')
C6225 FORMAT(////
C    416X,'TEMPERATURE DEPENDENCE SPECIFIC HEAT COEFICIENT:'/
C    416X,'(FOR NSTAC.EQ.0 FUNCTIONS OF SPECIFIC HEAT ARE NOT GIVEN)'/
C    416X,'(FOR NSTAC.EQ.1 NUMBER OF GIVEN FUNCTIONS ARE NMATT)'///
C    411X,'NUMBER OF FUNCTION GIVEN IN ANALITICAL FORM .... NANLC =',I5/
C    416X,'EQ.0; ALL FUNCTIONS ARE GIVEN AS TABLE'///
C    511X,'NUMBER OF FUNCTION GIVEN AS TABLE .............. NTABC =',I5/
C    516X,'EQ.0; ALL FUNCTION ARE GIVEN IN ANALITICAL FORM'///
C    611X,'MAX. NO. OF POINTS USED FOR DEFIN. FUNCTION .... MAXTC =',I5/
C    616X,'EQ.0; ALL FUNCTION ARE GIVEN IN ANALITICAL FORM'/
C    616X,'EQ.1; SPECIFIC HEAT COEFICIENTS ARE LINEAR'/
C    616X,'GT.1; SPECIFIC HEAT COEFICIENTS ARE NONLINEAR')
C     ENDIF

C Reading material constants: 
C Density of fluid and Indicator for nonlinear dynamic viscosity


      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1008) GUSM,INDAMI
C=======================================================================      
C For orbital shaker example:
c      CALL ISPITF(ACOZ,IULAZ)
c      READ(ACOZ,*) ORBRAD,ORBFRE
C=======================================================================      
      
      
 1008 FORMAT(F10.2,I5)
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,3004) GUSM,INDAMI
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6004) GUSM,INDAMI

 3004 FORMAT(6X,'P O D A C I   O   M A T E R I J A L N I M   K O N S T A 
     1 N T A M A'/6X,65('-')///
     111X,'GUSTINA FLUIDA IZNOSI....................... DENS = ',F10.5//
     211X,'INDIKATOR ZA KORISCENJE NE-NJUTNOVOG FLUIDA..INDAMI = ',I5//
     216X,'EQ.0; DINAMICKA VISKOZNOST JE KONSTANTNA'/
     216X,'EQ.1; KORISTI SE GENERALISANA CASSON-OVA RELACIJA'///)

 6004 FORMAT(6X,'D A T A   F O R   M A T E R I A L   C O N S T A N T S'/
     16X,53('-')///
     111X,'DENSITY OF FLUID ........................... DENS = ',F10.5//
     211X,'INDICATOR FOR NON-NEWTONIAN FLUID ........ INDAMI = ',I5/
     216X,'EQ.0; DYNAMIC VISCOSITY IS CONSTANTS'/
     216X,'EQ.1; GENERALIZED CASSON RELATION IS USED'///)
C==========================================================================
      IF (INDFL.EQ.0.OR.INDFL.EQ.2) THEN
C Reading heat conduction coefficient if heat analysis is required
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1010)AKT
c      READ(ACOZ,1010)dd(1),dd(2),dd(3),dd(4)
 1010 FORMAT(5F10.2)

      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,3005) AKT
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6005) AKT
 3005 FORMAT(
     111X,'KOEFICIJENT KONDUKCIJE ................... CONDUC= ',D13.5//)
 6005 FORMAT(
     111X,'CONDUCTION COEFFICIENT ................... CONDUC= ',D13.5//)

      ENDIF
C==========================================================================
C==========================================================================
      IF (INDFL.EQ.0.OR.INDFL.EQ.2) THEN
C Reading specific heat (pressure-constant) if heat analysis is required
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1010) CC
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,3006) CC
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6006) CC
 3006 FORMAT(
     111X,'SPEC. TOPLOTA PRI KONST. PRITISKU ........ SPECH = ',F10.2//)
 6006 FORMAT(
     111X,'SPEC. HEAT AT CONST. PRESSURE ............ SPECH = ',F10.2//)
      ENDIF
C==========================================================================
C Reading gravity accelerations,dynamics viscosity,
C thermic expansion coefficient,referential temperature 
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1011) FB2,FB3,FB33,AMI,BETA,TETAO
 1011 FORMAT(6D10.2)
C==========================================================================
C   FOR ORBITAL SHAKER EXAMPLE
c      AMI=AMI*5.D2 
c      open (32,file='Shear stress.txt')
c      open (32,file='Percentage.txt')
C==========================================================================
 
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,1106) FB2,FB3,FB33,AMI,BETA,TETAO
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,3106) FB2,FB3,FB33,AMI,BETA,TETAO
 1106 FORMAT(
     111X,'UBRZANJE ZEMLJINE TEZE U PRAVCU OSE X ....... GX =',D10.2//
     111X,'UBRZANJE ZEMLJINE TEZE U PRAVCU OSE Y ....... GY =',D10.2//
     111X,'UBRZANJE ZEMLJINE TEZE U PRAVCU OSE Z ....... GZ =',D10.2//
     111X,'DINAMICKA VISKOZNOST ....................... AMI =',D10.2//
     111X,'KOEFICIJENT TERMICKE EKSPANZIJE ........... BETA =',D10.2//
     111X,'REFERENTNA TEMPERATURA ................... TETA0 =',D10.2///)
 3106 FORMAT(
     111X,'GRAVITY ACCELERATION IN X DIRECTION ......... GX =',D10.2//
     111X,'GRAVITY ACCELERATION IN Y DIRECTION ......... GY =',D10.2//
     111X,'GRAVITY ACCELERATION IN Z DIRECTION ......... GZ =',D10.2//
     111X,'DYNAMIC VISCOSITY .......................... AMI =',D10.2//
     111X,'COEFFICIENT OF THERMAL EXPANSION .......... BETA =',D10.2//
     111X,'REFERENCE TEMPERATURE .................... TETA0 =',D10.2///)
C==========================================================================
 1002 FORMAT(15I5)


      END
C=======================================================================
      SUBROUTINE INPU08(IULAZ,NUMZAD,MAXSIL,IIZLAZ,ISRPS)
      CHARACTER*250ACOZ
C
CE Subroutine INPU08 is used for reading basic data for finite elements
C

      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1001) NUMZAD,MAXSIL
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,3007) NUMZAD,MAXSIL
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6007) NUMZAD,MAXSIL
 3007 FORMAT(6X,'P O D A C I   O   Z A D A T I M   V R E D N O S T I M A
     1'/6X,55('-')///
     111X,'UKUPAN BROJ ZADATIH VREDNOSTI .............. NUMZAD = ',I5//
     211X,'BROJ GRANICNIH POVRSINA ZA RACUNANJE'/
     211X,'POVRSINSKIH SILA ............................ NUMST = ',I5//)
 6007 FORMAT(6X,'D A T A   F O R   P R E S C R I B E D   V A L U E S'/
     16X,51('-')///
     111X,'GLOBAL NUM. OF PRESC. VALUES AT NODES ...... NUMZAD = ',I5//
     211X,'GLOBAL NUMBER OF BOUNDARY SURFACES ON WHICH'/ 
     211X,'SURFACE TRACTION ARE CALCULATED ............. NUMST = ',I5//)
 1001 FORMAT(5I5)

      END
C=======================================================================
C==========================================================================
      SUBROUTINE TIDPRI(PENALT,NETIP,IDPRIT,NDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine TIDPRI is used for definition 
CE dimension IDPRIT for pressure matrix 
C

      IF (PENALT.GT.0.D0) THEN
       IF (NDIM.EQ.4.AND.NETIP.EQ.2) IDPRIT=1
       IF (NDIM.EQ.9.AND.NETIP.EQ.2) IDPRIT=4
       IF (NDIM.EQ.4.AND.NETIP.EQ.3) IDPRIT=1
       IF (NDIM.EQ.8.AND.NETIP.EQ.3) IDPRIT=1
       IF (NDIM.EQ.10.AND.NETIP.EQ.3) IDPRIT=1
       IF (NDIM.EQ.21.AND.NETIP.EQ.3) IDPRIT=8
      ELSE
        IDPRIT=2
      ENDIF

      END
C=======================================================================
C==========================================================================
      SUBROUTINE TNDIM(NDIM,NP2DMX,NDES,NETIP)

C
CE Subroutine TNDIM is used for definition NDES-dimension matrix per element
C


C DODELJIVANJE ZA NDIM,NDES
      NDIM=NP2DMX
      IF (NETIP.EQ.2) THEN
c This should be change for acoustics problem
c         NDES=4*NDIM+4
         NDES=6*NDIM
      ELSEIF (NETIP.EQ.3) THEN
c         NDES=4*NDIM+21
         NDES=10*NDIM
      ENDIF
 

      END
C=======================================================================
C=======================================================================
      SUBROUTINE INPF03(IULAZ,NPT,NSTAC,NPER,INDFL,IPENAL,IIZLAZ,ISRPS,
     &IALE,ITURB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /STAUNV/ NPRINT
	COMMON /SPBCON/ IBOUND
      COMMON /IZID/ IZID
C
CE Subroutine INPF03 is used for reading basic data for problem
C

C Reading basic data of problem
      CHARACTER*250ACOZ

      CALL ISPITF(ACOZ,IULAZ)
      READ (ACOZ,1022) NPT,NGET,NMATT,NSTAC,NPER,NPRINT,INDFL,IPENAL,
     &IALE,IBOUND,ITURB,IZID

 1002 FORMAT(15I5)
 1022 FORMAT (I10,15I5)
      IF(NPER.EQ.0) NPER = 1
      IF(NGET.EQ.0) NGET = 1
      IF(NMATT.EQ.0) NMATT=1
C      IF(NPRINT.EQ.0) NPRINT=1

      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2001) NPT,NGET,NMATT,NSTAC,NPER,NPRINT,INDFL,IPENAL,
     &IALE,IBOUND,ITURB,IZID
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6001) NPT,NGET,NMATT,NSTAC,NPER,NPRINT,INDFL,IPENAL,
     &IALE,IBOUND,ITURB,IZID
 2001 FORMAT(6X,'O S N O V N I    P O D A C I    O    P R O B L E M U'
     1/6X,51('-')///
     111X,'UKUPAN BROJ CVORNIH TACAKA ........................ NP =',
     1I10/
     116X,'EQ.0; PREKIDA SE IZVRSAVANJE PROGRAMA'///
     211X,'BROJ GRUPA ELEMENATA ............................ NGET =',I5/
     216X,'EQ.0; POSTAJE "1"; (MAX. 10 GRUPA)'///
     311X,'BROJ RAZLICITIH MATERIJALA ..................... NMATT =',I5/
     316X,'EQ.0; POSTAJE "1"'///
     411X,'INDIKATOR STACIONARNOSTI ....................... NSTAC =',I5/
     416X,'EQ.1; STACIONARAN PROBLEM'/
     416X,'EQ.0; NESTACIONARAN PROBLEM'///
     511X,'BROJ PERIODA SA KONSTANTNIM VREMENSKIM KORACIMA . NPER =',I5/
     516X,'EQ.0; POSTAJE "1"'///
     611X,'DEFINISANJE STAMPARSKOG KORAKA ................ NPRINT =',I5/
     616X,'EQ.0; POSTAJE "1"'///
     711X,'INDIKATOR TIPA ANALIZE  ........................ INDFL =',I5/
     716X,'EQ.2; ANALIZA PROVODJENJA TOPLOTE BEZ DINAMIKE FLUIDA'/
     716X,'EQ.1; ANALIZA DINAMIKE FLUIDA BEZ PROVODJENJA TOPLOTE'/
     716X,'EQ.0; SPREGNUTO RESAVANJE SA PROVODJENJEM TOPLOTE'///
     711X,'INDIKATOR ZA KORISCENJE PENALTY ANALIZE ....... IPENAL =',I5/
     716X,'EQ.1; KORISTI SE PENALTY METOD'/
     716X,'EQ.0; KORISTI SE MESOVITA FORMULACIJA'///
     711X,'INDIKATOR ZA KORISCENJE ALE FORMULACIJE .....,,.. IALE =',I5/
     716X,'EQ.1; KORISTI SE ALE FORMULACIJA'/
     716X,'EQ.0; NE KORISTI SE ALE FORMULACIJA'///
     711X,'INDICATOR FOR SPECIFIC BOUNDARY CONDITION.....  IBOUND =',I5/
     716X,'EQ.1; SPECIFIC BOUNDARY CONDITION ARE USED'/
     716X,'EQ.0; SPECIFIC BOUNDARY CONDITION ARE NOT USED'///
     711X,'INDIKATOR ZA TURBULENTNO STRUJANJE ............. ITURB =',I5/
     716X,'EQ.0; NE KORISTI SE TURBULENTNO STRUJANJE'/
     716X,'EQ.1; KORISTI SE K-OMEGA TURBULENTNI MODEL'/
     716X,'EQ.2; KORISTI SE LES TURBULENTNI MODEL'///
     711X,'INDIKATOR ZA ZIDNE FUNKCIJE .................... IZID =',I5/
     716X,'EQ.0; NE KORISTE SE ZIDNE FUNKCIJE'/
     716X,'EQ.1; KORISTE SE ZIDNE FUNKCIJE'///)
 6001 FORMAT(6X,'B A S I C    D A T A    F O R   T H E   P R O B L E M'
     1/6X,53('-')///
     111X,'TOTAL NUMBER OF NODAL POINTS ...................... NP =',
     1I10/
     116X,'EQ.0; PROGRAM STOP'///
     211X,'NUMBER OF ELEMENT GROUPS ........................ NGET =',I5/
     216X,'EQ.0; DEFAULT SET "1"'///
     311X,'NUMBER OF DIFERENT MATERIALS ................... NMATT =',I5/
     316X,'EQ.0; DEFAULT SET "1"'///
     411X,'STEADINES INDICATOR ............................ NSTAC =',I5/
     416X,'EQ.1; STEADY STATE'/
     416X,'EQ.0; TRANSIENT'///
     511X,'NUMBER OF CONSTANT TIME STEP PERIODS ............ NPER =',I5/
     516X,'EQ.0; DEFAULT SET "1";'///
     611X,'OUTPUT PRINTING INTERVAL ...................... NPRINT =',I5/
     616X,'EQ.0; DEFAULT SET "1"'///
     711X,'INDICATOR FOR ANALYSIS ........................  INDFL =',I5/
     716X,'EQ.0; ANALYSIS FLUID WITH COUPLED HEAT TRANSFER'/
     716X,'EQ.1; FLUID FLOW ANALYSIS ONLY'/
     716X,'EQ.2; HEAT TRANSFER ONLY'///
     711X,'INDICATOR FOR PENALTY ANALYSIS................  IPENAL =',I5/
     716X,'EQ.0; MIXED FORMULATION'/
     716X,'EQ.1; PENALTY METHOD USED'///
     711X,'INDICATOR FOR ALE FORMULATION...................  IALE =',I5/
     716X,'EQ.1; ALE FORMULATION IS USED'/
     716X,'EQ.0; ALE FORMULATION IS NOT USED'///
     711X,'INDICATOR FOR SPECIFIC BOUNDARY CONDITION.....  IBOUND =',I5/
     716X,'EQ.1; SPECIFIC BOUNDARY CONDITION ARE USED'/
     716X,'EQ.0; SPECIFIC BOUNDARY CONDITION ARE NOT USED'///
     711X,'INDICATOR FOR TURBULENT FLOW ................... ITURB =',I5/
     716X,'EQ.0; TURBULENT FLOW IS NOT USED'/
     716X,'EQ.1; K-OMEGA TURBULENT MODEL IS USED'/
     716X,'EQ.2; LES TURBULENT MODEL IS USED'///
     711X,'INDICATOR FOR WALL FUNCTIONS ................... IZID =',I5/
     716X,'EQ.0; WALL FUNCTIONS ARE NOT USED'/
     716X,'EQ.1; WALL FUNCTIONS ARE USED'///)


      END
C==========================================================================
C=========================================================================
C=========================================================================
      SUBROUTINE SPSOUT(SPSIL,NETIP,IIZLAZ,NPT,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SPSIL(NETIP,*)
C
CE Subroutine SPSOUT is used for printout forces from fluid calculations
C
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,1001)
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,1000)
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,2001)
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,2000)

      DO NODE=1,NPT
       WRITE(IIZLAZ,100) NODE,(SPSIL(J,NODE),J=1,NETIP)
C      WRITE(IIZLAZ,200) NODE,1,2,SPSIL(1,NODE)
C      WRITE(IIZLAZ,200) NODE,2,2,SPSIL(2,NODE)
C      WRITE(IIZLAZ,200) NODE,3,2,SPSIL(3,NODE)
      ENDDO           


 100  FORMAT(I5,3(D13.5))
 200  FORMAT(3I5,D13.5)
1000  FORMAT(///6X,'F O R C E S  F R O M  F L U I D  C A L C U L A T I O
     & N'/6X,53('-')//)
1001  FORMAT(///6X,'S I L E  K O J E   S E  D O B I J A J U   P R O R A 
     &C U N O M   F L U I D A'/6X,/75('-'))
2000  FORMAT('NODE    FORCES (FX,FY,FZ)'/)
2001  FORMAT('CVOR    SILE (FX,FY,FZ)'/)
      END
C=========================================================================
      SUBROUTINE IZLLSF(GNODE,IIZLAZ,IDPRIT,ISRPS,NPT,PRES,IFORM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM
C
CE  Subroutine IZLLSF is used for printing output data in "name.LST" file
C
C       COMMON /SRPSKI/ ISRPS

C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      
       DIMENSION GNODE(2,7,*),N4(4),D4(4),PRES(3,*)
C
CS  BRZINE CVOROVA U PRAVCU X1     
CE  NODE VELOCITIES IN DIRECTION X1
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)

      KG = 0
   20 DO 50 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 100
      N4(I) = KG
C      JJEDN=ID(1,KG)
C      IF(JJEDN.EQ.0) THEN
C         D4(I) = 0.
C      ELSE
C         D4(I) = TT1(JJEDN)
         D4(I) = GNODE(2,1,KG)
C      ENDIF
   50 CONTINUE
C
  100 I = I - 1
      IF(IFORM.EQ.1) GO TO 88
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 77
C
   88 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
   77 IF(KG.LT.NPT) THEN
      GO TO 20
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2005)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6005)
C
CS  BRZINE CVOROVA U PRAVCU X2     
CE  NODE VELOCITIES IN DIRECTION X2
C
      KG = 0
  120 DO 150 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1100
      N4(I) = KG
         D4(I) = GNODE(2,2,KG)
C      ENDIF
  150 CONTINUE
C
 1100 I = I - 1
      IF(IFORM.EQ.1) GO TO 188
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 177
C
  188 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
  177 IF(KG.LT.NPT) THEN
      GO TO 120
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6010)
C
CS  BRZINE CVOROVA U PRAVCU X3     
CE  NODE VELOCITIES IN DIRECTION X3
C
      KG = 0
  420 DO 450 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1400
      N4(I) = KG
         D4(I) = GNODE(2,3,KG)
  450 CONTINUE
C
 1400 I = I - 1
      IF(IFORM.EQ.1) GO TO 488
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 477
C
  488 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
  477 IF(KG.LT.NPT) THEN
      GO TO 420
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2015)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6015)
C
CS  PRITISCI U CVOROVIMA                                
CE  NODE PRESSURES 
C
      KG = 0
  220 DO 250 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1200
C      IF(INDSC.EQ.0) GO TO 256
C      IF(ISKC(KG).EQ.0) GO TO 255
      N4(I) = KG
C      JJEDN=ID(3,KG)
C      IF(JJEDN.EQ.0) THEN
C         D4(I)=PRES(KG)
C      ELSE
C         D4(I) = TT1(JJEDN)
         D4(I) = GNODE(2,4,KG)
C      ENDIF
  250 CONTINUE
C
 1200 I = I - 1
      IF(IFORM.EQ.1) GO TO 288
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 277
C
  288 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
  277 IF(KG.LT.NPT) THEN
      GO TO 220
      ENDIF
C
C
CS  TEMPERATURE U CVOROVIMA                                
CE  NODE TEMPERATURES 
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2025)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6025)



      KG = 0
  320 DO 350 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1300
      N4(I) = KG
         D4(I) = GNODE(2,5,KG)
  350 CONTINUE
C
 1300 I = I - 1
      IF(IFORM.EQ.1) GO TO 388
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 377
C
  388 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
  377 IF(KG.LT.NPT) THEN
      GO TO 320
      ENDIF





C
C
CS  ELEKTRICNI POTENCIJALI U CVOROVIMA                                
CE  NODE ELECTRICAL POTENTIAL 
C
C      IF(ISRPS.EQ.0)
C     1WRITE(IIZLAZ,2026)
C      IF(ISRPS.EQ.1)
C     1WRITE(IIZLAZ,6026)
C      KG = 0
C 3201 DO 3501 I=1,4
C      KG = KG +1
C      IF(KG.GT.NPT) GO TO 13001
C      N4(I) = KG
C         D4(I) = GNODE(2,6,KG)
C 3501 CONTINUE
C
C13001 I = I - 1
C      IF(IFORM.EQ.1) GO TO 3881
C      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
C      GO TO 3771
C
C 3881 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
C 3771 IF(KG.LT.NPT) THEN
C      GO TO 3201
C      ENDIF
C------------------------------------------------------------------
CS  TURBULENTNA KINETICKA ENERGIJA K                                
CE  TURBULENT KINETIC ENERGY K  
C------------------------------------------------------------------
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2035)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6035)
C
      KG = 0
  421 DO 451 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1401
      N4(I) = KG
         D4(I) = GNODE(2,6,KG)
  451 CONTINUE
C
 1401 I = I - 1
      IF(IFORM.EQ.1) GO TO 489
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 478
C
  489 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
  478 IF(KG.LT.NPT) THEN
      GO TO 421
      ENDIF
C------------------------------------------------------------------
C------------------------------------------------------------------
CS  SPECIFICNA DISIPACIJA KINETICKE ENERGIJE TURBULENCIJE OMEGA                                
CE  TURBULENT SPECIFIC DISSIPATION RATE OMEGA  
C------------------------------------------------------------------
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2045)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6045)
C
      KG = 0
  521 DO 551 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1501
      N4(I) = KG
         D4(I) = GNODE(2,7,KG)
  551 CONTINUE
C
 1501 I = I - 1
      IF(IFORM.EQ.1) GO TO 589
      WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      GO TO 578
C
  589 WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
  578 IF(KG.LT.NPT) THEN
      GO TO 521
      ENDIF
C------------------------------------------------------------------


      RETURN
 5002 FORMAT(4(I5,D13.5))
 5003 FORMAT(4(I5,3X,F10.3))
 2001 FORMAT(/'    C V O R N E    B R Z I N E    U    X 1    P R A V C '
     1,'U            '//
     1' CVOR   BRZINA     CVOR   BRZINA     CVOR   BRZINA     CVOR  BRZI
     1NA    '/' BROJ',3(14X,'BROJ'))
 2005 FORMAT(/'    C V O R N E    B R Z I N E    U    X 2    P R A V C '
     1,'U            '//
     1' CVOR   BRZINA     CVOR   BRZINA     CVOR   BRZINA     CVOR  BRZI
     1NA    '/' BROJ',3(14X,'BROJ'))
 2010 FORMAT(/'    C V O R N E    B R Z I N E    U    X 3    P R A V C '
     1,'U            '//
     1' CVOR   BRZINA     CVOR   BRZINA     CVOR   BRZINA     CVOR  BRZI
     1NA    '/' BROJ',3(14X,'BROJ'))
 2015 FORMAT(/'                  C V O R N I    P R I T I S C I  '//
     1' CVOR  PRITISAK    CVOR  PRITISAK    CVOR  PRITISAK    CVOR  PRIT
     1ISAK  '/' BROJ',3(14X,'BROJ'))
 2025 FORMAT(/'                  C V O R N E    T E M P E R A T U R E'//
     1' CVOR  TEMPERATURA CVOR  TEMPERATURA CVOR  TEMPERATURA CVOR TEMPE       
     1RATURA'/' BROJ',3(14X,'BROJ'))
 2035 FORMAT(/'                  SPEC.  DISIPACIJA  KIN.  EN.  OMEGA'//
     1' CVOR  OMEGA       CVOR  OMEGA       CVOR  OMEGA       CVOR OMEGA       
     1'/' BROJ',3(14X,'BROJ'))
 2045 FORMAT(/'                  TURBULENTNA KINETICKA ENERGIJA K'//
     1' CVOR  TURBULENTK  CVOR  TURBULENTK  CVOR  TURBULENTK  CVOR TURBU       
     1LENTK'/' BROJ',3(14X,'BROJ'))
 6001 FORMAT(/'    N O D A L    V E L O C I T I E S   I N   D I R E C T'
     1,' I O N   X 1 '//
     1' NODE  VELOCITY    NODE  VELOCITY    NODE  VELOCITY    NODE VELOC
     1ITY   '/'  No.',3(15X,'No.'))
 6005 FORMAT(/'    N O D A L    V E L O C I T I E S   I N   D I R E C T'
     1,' I O N   X 2 '//
     1' NODE  VELOCITY    NODE  VELOCITY    NODE  VELOCITY    NODE VELOC
     1ITY   '/'  No.',3(15X,'No.'))
 6010 FORMAT(/'    N O D A L    V E L O C I T I E S   I N   D I R E C T'
     1,' I O N   X 3 '//
     1' NODE  VELOCITY    NODE  VELOCITY    NODE  VELOCITY    NODE VELOC
     1ITY   '/'  No.',3(15X,'No.'))
 6015 FORMAT(/'                    N O D A L    P R E S S U R E S'//
     1' NODE  PRESSURE    NODE  PRESSURE    NODE  PRESSURE    NODE PRESS
     1URE   '/'  No.',3(15X,'No.'))
 6025 FORMAT(/'                  N O D A L    T E M P E R A T U R E S'//
     1' NODE  TEMPERATURE NODE  TEMPERATURE NODE  TEMPERATURE NODE TEMPE       
     1RATURE'/'  No.',3(15X,'No.'))
 6035 FORMAT(/'                  T U R B U L E N T     K'//
     1' NODE  TURBULENTK  NODE  TURBULENTK  NODE  TURBULENTK  NODE TURBU       
     1LENTK'/'  No.',3(15X,'No.'))
 6045 FORMAT(/'  SPECIFIC  DISSIPATION  OF  KINETIC  ENERGY  OMEGA'//
     1' NODE  OMEGA       NODE  OMEGA       NODE  OMEGA       NODE OMEGA       
     1'/'  No.',3(15X,'No.'))

      END
C======================================================================
C==========================================================================
      SUBROUTINE NENJUT(ZVHX,ZVHY,TT210)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ULAZNI/ IULAZ,IIZLAZ
      COMMON /VISKOZ/ AMI,INDAMI
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

      DIMENSION ZVHX(*),ZVHY(*),TT210(*),DUDX(2,2),D(2,2)
C
CE Subroutine NENJUT is used for non-Newtonian fluid calculation
C

C INITIALISATION
C==========================================================================
      DO I=1,2
       DO J=1,2
        DUDX(I,J)=0.D0
       ENDDO
      ENDDO

      DII=0.D0
      AK0=0.6125D0
      AK1=0.174D0
C==========================================================================


      DO I=1,NDIM
       DUDX(1,1)=DUDX(1,1)+ZVHX(I)*TT210(I)
       DUDX(1,2)=DUDX(1,2)+ZVHY(I)*TT210(I)
       DUDX(2,1)=DUDX(2,1)+ZVHX(I)*TT210(I+NDIM)
       DUDX(2,2)=DUDX(2,2)+ZVHY(I)*TT210(I+NDIM)
      ENDDO

    
      DO I=1,2
       DO J=1,2
        D(I,J)=0.5D0*(DUDX(I,J)+DUDX(J,I))
       ENDDO
      ENDDO

      DO I=1,2
       DO J=1,2
        DII=DII+0.5D0*(D(I,J)*D(I,J))
       ENDDO
      ENDDO

	GAMA=2.D0*DSQRT(DII)
      IF(DABS(GAMA).LT.1.D-8) RETURN
      AMI=(1.D0/(GAMA))*((AK0+AK1*DSQRT(GAMA))**2)
C     P1=1.D0/((1.D0+(0.5D0*GAMA)**1.7D0)**0.3D0)
C     AMI=0.03D0+(0.1315D0-0.03D0)*P1
      WRITE(IIZLAZ,*)'GAMA= ',GAMA,' AMI= ',AMI
      END
C==========================================================================
C==========================================================================
      SUBROUTINE NENJ3D(PJ,TT210,AMI,NDIM,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

c      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /VISKOZ/ AMI,INDAMI
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

      DIMENSION PJ(3,*),TT210(*),DUDX(3,3),D(3,3)
C
CE Subroutine NENJUT is used for non-Newtonian fluid calculation
C

C INITIALISATION
C==========================================================================
      DO I=1,3
       DO J=1,3
        DUDX(I,J)=0.D0
       ENDDO
      ENDDO

      DII=0.D0
      AK0=0.6125D0
      AK1=0.174D0
C==========================================================================


      DO I=1,NDIM
       DUDX(1,1)=DUDX(1,1)+PJ(1,I)*TT210(I)
       DUDX(1,2)=DUDX(1,2)+PJ(2,I)*TT210(I)
       DUDX(1,3)=DUDX(1,3)+PJ(3,I)*TT210(I)
       DUDX(2,1)=DUDX(2,1)+PJ(1,I)*TT210(I+NDIM)
       DUDX(2,2)=DUDX(2,2)+PJ(2,I)*TT210(I+NDIM)
       DUDX(2,3)=DUDX(2,3)+PJ(3,I)*TT210(I+NDIM)
       DUDX(3,1)=DUDX(3,1)+PJ(1,I)*TT210(I+2*NDIM)
       DUDX(3,2)=DUDX(3,2)+PJ(2,I)*TT210(I+2*NDIM)
       DUDX(3,3)=DUDX(3,3)+PJ(3,I)*TT210(I+2*NDIM)
      ENDDO

    
      DO I=1,3
       DO J=1,3
        D(I,J)=0.5D0*(DUDX(I,J)+DUDX(J,I))
       ENDDO
      ENDDO

      DO I=1,3
       DO J=1,3
        DII=DII+0.5D0*(D(I,J)*D(I,J))
       ENDDO
      ENDDO

	GAMA=2.D0*DSQRT(DII)
      IF(DABS(GAMA).LT.1.D-8) RETURN
      AMI=(1.D0/(GAMA))*((AK0+AK1*DSQRT(GAMA))**2)
C     P1=1.D0/((1.D0+(0.5D0*GAMA)**1.7D0)**0.3D0)
C     AMI=0.03D0+(0.1315D0-0.03D0)*P1
C     WRITE(IIZLAZ,*)'GAMA= ',GAMA,' AMI= ',AMI
      END
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE SSPSIL(SPSIL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

C
CE   Subroutine SPSIL is used for printing data for interactive forces 
CE   in output file *.LST
C
      DIMENSION SPSIL(2,*)

        WRITE(IIZLAZ,*)'INTERACTIVE FORCES FROM FLUID CALCULATION'
        WRITE(IIZLAZ,*)'CVOR         Fx         Fy'
 
       DO I=1,NPT
        WRITE(IIZLAZ,100) I,SPSIL(1,I),SPSIL(2,I)
       ENDDO
      
 100    FORMAT (I5,2(D13.5))

      END
C==========================================================================
        
        
C==========================================================================
      SUBROUTINE IWRIT(IA,LMAX,IZLAZ)
       DIMENSION IA(*)
C
CE  Subroutine IWRIT is used for printing integer data in output file IZLAZ
C
       DO I=1,LMAX
         WRITE (IZLAZ,100) IA(I)
       ENDDO
100   FORMAT (I10)
      END
C==========================================================================
C==========================================================================
      SUBROUTINE WRIT(A,LMAX,IZLAZ)
       DIMENSION A(*)
       REAL A
C
CE  Subroutine WRIT is used for printing real data in output file IZLAZ
C
       DO I=1,LMAX
         WRITE (IZLAZ,100) A(I)
       ENDDO
100   FORMAT (D13.5)
      END
C==========================================================================
C==========================================================================
      SUBROUTINE ROTATE(X,Y,X0,Y0,ANGLE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
	X_OLD=X
	Y_OLD=Y

      X_OLD=X_OLD-X0 
      Y_OLD=Y_OLD-Y0 
	X=X_OLD*DCOS(ANGLE)-Y_OLD*DSIN(ANGLE)
	Y=X_OLD*DSIN(ANGLE)+Y_OLD*DCOS(ANGLE)

	X=X+X0
	Y=Y+Y0

      END
C==========================================================================
