C=======================================================================
C=======================================================================
C   SUBROUTINE MINV
C              READDD
C              IREADD
C              CLEAR
C              ICLEAR
C              DELJIV
C              WRR
C              IWRR
C              JEDNA1
C              WRITED
C              READD
C=======================================================================
C=======================================================================
      SUBROUTINE MINV(A,N,D,L,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine MINV is used for inverting a matrix stored by columns vector
C

C
C ......................................................................
C .
CE.  P R O G R A M
CE.      TO INVERT A MATRIX STORED BY COLUMNS IN VECTOR
CS.   P R O G R A M
CS.      ZA INVERTOVANJE MATRICE A SLOZENE PO STUPCIMA U VEKTORU
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),L(*),M(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MINV'
CE    SEARCH FOR LARGEST ELEMENT
CS    NACI NAJVECI ELEMENT
      D=1.0d0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
CE    INTERCHANGE ROWS
CS    IZMENITI VRSTE
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C     INTERCHANGE COLUMNS
C     IZMENITI KOLONE
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
CE    DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
CE    CONTAINED IN BIGA)
CS    PODELITI KOLONU SA NEGATIVNIM STOZEROM-
CS    VREDNOST STOZERNOG ELEMENTA JE U BIGA.
   45 IF(BIGA) 48,46,48
   46 D=0.0D0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
CE    REDUCE MATRIX
CS    REDUKOVATI MATRICU
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
CE    DIVIDE ROW BY PIVOT
CS    PODELITI VRSTU STOZEROM
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
CE    PRODUCT OF PIVOTS
CS    PROIZVOD STOZERA
      D=D*BIGA
CE    REPLACE PIVOT BY RECIPROCAL
CS    ZAMENITI STOZER RECIPROCNOM VREDNOSCU
      A(KK)=1.0d0/BIGA
   80 CONTINUE
CE    FINAL ROW AND COLUMN INTERCHANGE
CS    POSLEDNJA IZMENA VRSTA I KOLONA
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE READDD(A,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine READDD is used for reading real vector on direct access file
C
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ REAL VECTOR ON DIRECT ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE REALNOG VEKTORA SA FILE SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READDD'
C      IF(ISKDSK.EQ.0.AND.N.EQ.NWK)RETURN
      IF(N.EQ.0) RETURN
      NK=N*8/LD
      NN=NK*LD/8
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/8+1
         NN=K*LD/8
         IF(K.EQ.NK) NN=N
         LS=LS+1
         READ(II,REC=LS) (A(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'RN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL WRR(A,N,'RN  ')
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE IREADD(IA,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine IREADD is used for reading integer vector on direct access file
C
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ INTEGER VECTOR ON DIRECT ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE CELOBROJNOG VEKTORA SA FILE SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IREADD'
      IF(N.EQ.0) RETURN
      NK=N*4/LD
      NN=NK*LD/4
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/4+1
         NN=K*LD/4
         IF(K.EQ.NK) NN=N
         LS=LS+1
         READ(II,REC=LS) (IA(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'IRN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL IWRR(IA,N,'IRN ')
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CLEAR(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine CLEAR is used for clearing a floating-point array
C
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR A FLOATING-POINT ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE REALNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR -  A)
CE.         A(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA -  A)
CS.         A(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CLEAR'
      DO 10 I=1,N
   10 A(I)=0.0D0
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ICLEAR(IA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine ICLEAR is used for clearing integer array
C
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR INTEGER ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE CELOBROJNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR - IA)
CE.        IA(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA - IA)
CS.        IA(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ICLEAR'
      DO 10 I=1,N
   10 IA(I)=0
      RETURN
      END
C=======================================================================
C
C======================================================================
      SUBROUTINE DELJIV(IBROJ,IDELI,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine DELJIV is used for testing divide of number without residue
C

C ......................................................................
C .
CE.    P R O G R A M
CE.        TO TEST DIVIDE OF NUMBER WITHOUT RESIDUE
CS.    P R O G R A M
CS.        ZA ISPITIVANJE DELJIVOSTI BROJA BEZ OSTATKA
C .
CE.            IBROJ  - BROJ KOJI DELIMO
CS.            IBROJ  - BROJ KOJI DELIMO
C .
CE.            IDELI  - BROJ KOJIM DELIMO
CS.            IDELI  - BROJ KOJIM DELIMO
C .
CE.            IND    - INDIKATOR OF DIVIDE
CE.                     =0; WITHOUT RESIDUE
CE.                     =0; WITH RESIDUE
CS.            IND    - INDIKATOR DELJIVOSTI
CS.                     =0; BEZ OSTATKA
CS.                     =0; SA OSTATKOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
C
C      IF(IDEBUG.GT.0) PRINT *, ' DELJIV'
      BROJ=IBROJ
      DELI=IDELI
      REZ=BROJ/DELI
      IREZ=REZ
      OST=REZ-IREZ
      TOL=0.0D00
      IND=0
      IF(DABS(OST).GT.TOL) IND=1
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRR(A,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine WRR is used writing real vector to output file
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRR'
      WRITE(IZLAZ,5010) CHAR
      WRITE(IZLAZ,5000) (A(I),I=1,N)
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(4(1PD18.9))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE IWRR(M,N,CHAR)
C
C
CE Subroutine IWRR is used writing integer vector to output file
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE INTEGER VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE CELOBROJNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION M(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IWRR'
      WRITE(IZLAZ,5010) CHAR
      WRITE(IZLAZ,5000) (M(I),I=1,N)
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(7I10)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE JEDNA1(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine JEDNA1 is used to equalizing two real vectors 
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO EQUALIZING 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA IZJEDNACAVANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDNA1'
      DO 10 I=1,N
   10 A(I)=B(I)
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE WRITED(A,N,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine WRITED is used 
CE for writing real vector in sequential access file
C
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR IN SEQUENTIAL ACCESS FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA U SEKVENCIJALNI FILE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRITED'
      LDUZZ=30000
      N1=1
      N2=N
      IF(N2.GT.LDUZZ) N2=LDUZZ
   10 WRITE(II) (A(I),I=N1,N2)
      N1=N1+LDUZZ
      IF(N1.GT.N) GO TO 20
      N2=N2+LDUZZ
      IF(N2.GT.N) N2=N
      GO TO 10
   20 CONTINUE
C      WRITE(3,*) 'WNS,II',N,II
C      CALL WRR(A,N,'WNS ')
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE READD(A,N,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
CE Subroutine READD is used 
CE for reading real vector from sequential access file
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ REAL VECTOR ON SEQUENTIAL ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE REALNOG VEKTORA IZ SEKVENCIJALNOG FILE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READD'
      LDUZZ=30000
      N1=1
      N2=N
      IF(N2.GT.LDUZZ) N2=LDUZZ
   10 READ(II) (A(I),I=N1,N2)
      N1=N1+LDUZZ
      IF(N1.GT.N) GO TO 20
      N2=N2+LDUZZ
      IF(N2.GT.N) N2=N
      GO TO 10
   20 CONTINUE
C      WRITE(3,*) 'RNS,II',N,II
C      CALL WRR(A,N,'RNS ')
      RETURN
      END
C=======================================================================
C
C======================================================================
      DOUBLE PRECISION FUNCTION DOT(A,B,N)
C....  SKALARNI PROIZVOD VEKTORA
      DOUBLE PRECISION A,B
      DIMENSION A(*),B(*)
      DOT=0.D0
      DO 10 I=1,N
   10 DOT=DOT+A(I)*B(I)
      RETURN
      END
C=======================================================================
C

