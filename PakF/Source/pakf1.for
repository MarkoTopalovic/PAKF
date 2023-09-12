#define MUMPS_CLUSTER .FALSE.
C=======================================================================
C
C      SUBROUTINE PERIOD
C
C=======================================================================
      SUBROUTINE PERIOD(A,NKDT,DTDT,NPER,LPAKF,LSKF,IBKOR)
      use servis
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
C ......................................................................
C .
CE.    P R O G R A M
CE.        WITH LOOP OVER TIME PERIODS AND STEPS
CS.    P R O G R A M
CS.        SA PETLJOM PO VREMENSKIM PERIODIMA I KORACIMA
C .
C ......................................................................
C
      DIMENSION A(*)
      REAL A
      DIMENSION NKDT(*),DTDT(*)
      DIMENSION IS(1)

      IS(1)=0
      VREM0=0.D0
      INDT=0
      if(IMUMPS.GE. 5 .and. IMUMPS.LE.9)call servis_start()
C
CE    BASIC LOOP OVER TIME PERIODS
CS    OSNOVNA PETLJA PO VREMENSKIM PERIODIMA
C
      DO 100 IPER=1,NPER
         IINDT=NKDT(IPER)
         DT=DTDT(IPER)
         IPDT=INDT+1
         INDT=INDT+IINDT
C
CE    BASIC LOOP OVER TIME STEPS
CS    OSNOVNA PETLJA PO VREMENSKIM KORACIMA
C
      DO 500 KORBR=IPDT,INDT
         VREM0 = VREM0 + DT
            CALL PPAKF(A,LPAKF,IPER,LSKF,IS,VREM0,KORBR,IBKOR)
  500 CONTINUE 
C
  100 CONTINUE
C
	if (IMUMPS.EQ.1 .or. (IMUMPS.GE. 5 .and. IMUMPS.LE.9)) then 
        if(imumps .gt. 1 ) call servis_done()
        call measure_time_end() 	
	  close (mufile)
	  close (mufile2)
	endif
      CALL ZAGEND(56)
      
      RETURN
      END
C=======================================================================
C=======================================================================


