C=======================================================================
C
C  MAIN PROGRAM PAKF
C
C=======================================================================
#define MUMPS_CLUSTER .FALSE.
      PROGRAM PAKF
      use servis
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      PARAMETER (NTOT = 14000000)
C      PARAMETER (NTOT = 7000000)
C for Akira's problem 256 MB-RAM
C       PARAMETER (NTOT = 63000000)
c       PARAMETER (NTOT = 120000000)
C       PARAMETER (NTOT = 126000000)
C       PARAMETER (NTOT = 197000000)
C      PARAMETER (NTOT = 185000000)
c      PARAMETER (NTOT = 320000000)
      PARAMETER (NTOT = 250000000)
C
C ......................................................................
C .
C .                     P A K - F 
C .
CE.        FINITE ELEMENT PROGRAM FOR COMPUTATIONAL FLUID DYNAMICS
CE.        PROGRAM MAIN
C .
C ......................................................................
C
CE    NTOT - MAXIMUM TOTAL STORAGE AVAILABLE
CS    NTOT - MAKSIMALAN PROSTOR U VEKTORU A
C
      COMMON A(NTOT)
      REAL A
      DIMENSION NKDT(100),DTDT(100)
      INTEGER IBKOR
C      DIMENSION PRES1(100,3,100000)
C      DIMENSION TAU(100,2,100000)
C      DIMENSION VOSI(100,100000)
C
C     INDIKATOR ZA DVOSTRUKU PRECIZNOST
C     POCETNI REPER U VEKTORU A(LMAX)
      LMAX = 1

#define MUMPS_CLUSTER .FALSE.

C
CE    PAKF - INPUT DATA 
CS    PAKF - ULAZNI PODACI
C
C     POCETNI REPER ZA PAKF
      LPAKF=LMAX
C==========================================================================
C READING INPUT DATA FROM INPUT FILE
        IBKOR=0
        CALL INPAKF(A,NTOT,LMAX,NKDT,DTDT,NPER,IBKOR)
C==========================================================================
         CALL DELJIV(LMAX,2,INDL)
         IF(INDL.EQ.0) LMAX=LMAX+1
         LSKF=LMAX
CC
CE    TIME PERIOD LOOP
CS    OSNOVNA PETLJA PO VREMENSKIM PERIODIMA
C     -----------
      call measure_time_start()
      CALL NEUTRAF(NKDT,DTDT,NPER,1,69)     
      CALL PERIOD(A,NKDT,DTDT,NPER,LPAKF,LSKF,IBKOR)
      CALL NEUTRAF(NKDT,DTDT,NPER,2,69)
      STOP
      END
C=======================================================================
