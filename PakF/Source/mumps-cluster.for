C=========================================================================
      SUBROUTINE MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2

      DIMENSION NEL(NDIM+1,*),ID(6,*)
	DIMENSION LM2(NETIP*NDIM)

      IF (IMUMPS.EQ.0) RETURN

       REWIND(MUFILE)

	NN=NETIP*NDIM

	WRITE(MUFILE) JEDN
	WRITE(MUFILE) NET

      n=1
	N1=0
      DO NBREL=1,NET
	N0=0
C	  WRITE(mufile,*) n
      DO II=1,NETIP*NDIM
        LM2(II)=0
      ENDDO
C=========================================================================
      DO KLM=1,NDIM
	  DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
 	  ENDDO
	ENDDO

      DO JJ=1,NN
	  IF (LM2(JJ).NE.0) THEN
	     N=N+1
	   N0=N0+1
        ENDIF
	ENDDO
	N1=N1+N0*N0  
	enddo
C	  WRITE(mufile,*) n

	WRITE(MUFILE) N-1
	WRITE(MUFILE) N1


      n=1
      DO NBREL=1,NET
	  WRITE(mufile) n
      DO II=1,NETIP*NDIM
        LM2(II)=0
      ENDDO
C=========================================================================
      DO KLM=1,NDIM
	  DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
 	  ENDDO
      ENDDO

      DO JJ=1,NN
	  IF (LM2(JJ).NE.0) N=N+1
	ENDDO
	  
	enddo
	  WRITE(mufile) n

C	WRITE(MUFILE,*) N-1

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,NETIP*NDIM
        LM2(I)=0
 125  CONTINUE
C=========================================================================
      DO 130 KLM=1,NDIM
	 DO J=1,NETIP
        LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
	 ENDDO
 130  CONTINUE

      do i=1,nn
C	      WRITE(MUFILE,*) (LM2(I),I=1,NN)
	 IF (LM2(I).NE.0) WRITE(MUFILE) LM2(I)
      enddo
      
400   CONTINUE


	END
C=======================================================================
C=========================================================================
C=========================================================================
      SUBROUTINE MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,NETIP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2
      DIMENSION NEL(NDIM+1,*),ID(6,*)
	DIMENSION LM2(NETIP*NDIM)
      DOUBLE PRECISION elembuffer(67108864)
      INTEGER BUFFERSIZE

      DIMENSION SKEF(NDES,*)



      IF (IMUMPS.EQ.0) RETURN


      DO II=1,NETIP*NDIM
        LM2(II)=0
      ENDDO
C=========================================================================
      DO KLM=1,NDIM
	  DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
 	  ENDDO
	ENDDO

C      DO I=1,NETIP*ndim
C      DO J=1,NETIP*ndim
C	 IF (LM2(J).NE.0.AND.lm2(i).ne.0) then
C	   WRITE(MUFILE2) SKEF(J,I)
C       endif
C	ENDDO
C	ENDDO

C ========================================================================
C Milos dodao baferisanje od 256MB zbog performansi, Avgust 2009.
      BUFFERSIZE = 67108864
      DO KLM=1,NDIM
          DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
          ENDDO
        ENDDO

      buffercounter = 0
      DO I=1,NETIP*ndim
      DO J=1,NETIP*ndim
         if (LM2(J).NE.0.AND.lm2(i).ne.0) then
           buffercounter = buffercounter+1
           elembuffer(buffercounter) = SKEF(J,I)
         endif
         if (buffercounter.eq.BUFFERSIZE) then
           write(mufile2) (elembuffer(k), k=1,BUFFERSIZE)
           buffercounter = 0
         endif
      ENDDO
      ENDDO

      write(mufile2) (elembuffer(k), k=1, buffercounter)


	END
C=========================================================================
C=========================================================================
C=========================================================================
C=========================================================================
      SUBROUTINE MUMPSRIGHT(A,N)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON /MUMPS_PAK/ IMUMPS, MUFILE,MUFILE2


      DIMENSION A(*)

      IF (IMUMPS.EQ.0) RETURN
	
C      do i=1,n
c	WRITE(MUFILE,*) (A(I),I=1,N)
 	WRITE(MUFILE2) (A(I), I=1,N)
C	enddo
C        CLOSE(MUFILE2)
	END
C=========================================================================
C=========================================================================
C=======================================================================
      SUBROUTINE checksolver(sile,jedn,right)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	DIMENSION SILE(*), right(*)


      iKc=0
      mufile2=107
      IIZLAZ=78
      OPEN (MUFILE2, FILE='FULLMATR.BIN',FORM='BINARY')
      OPEN (iizlaz, FILE='comparison.txt')

       
	read(mufile2) Jedn
	read(mufile2) iKc
	do i=1,jedn
	 right(i)=-right(i)
	enddo
      do k=1,iKc
       read(mufile2)i,j,a
       right(i)=right(i)+a*sile(j)
	enddo

 200  FORMAT (I5,1PD20.12)

      DO I=1,JEDN
        WRITE(IIZLAZ,200) i,right(i)
	ENDDO

      close (mufile2)
      close (iizlaz)

      END
C=======================================================================
C==========================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE FULLMA(ALEVO,DESNO,MaxKs,JEDN,SILE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION ALEVO(*),DESNO(*),sile(*)
	dimension MaxKs(*)


      iKc=0
      mufile2=107
	tol=1.0D-20
c      OPEN(iizlaz,FILE='PUNAMATR.TXT')
      OPEN (MUFILE2,FILE='/tmp/PAKF_MUMPS.OUT',FORM='BINARY')

      iKc=0

C	WRITE(mufile2) JEDN
C	WRITE(mufile2) iKc




	DO 20 I=1,JEDN
       DO 10 J=1,JEDN
        if(i.eq.1.and.j.eq.1) then
	     iKc=iKc+1
		 goto 10 
	  endif
            IF(J.GE.I) THEN
            iKs=MaxKs(j)-j+i
            if (iKs.gt.MaxKs(j-1)) then
              if (DABS(ALEVO(iKs)).GT.tol) then
                  iKc=iKc+1
C                  WRITE(mufile2,*)I,J,ALEVO(IKS)
C                  WRITE(IIZLAZ,*)I,J,ALEVO(IKS)
C                  Kc(iKc)=ALEVO(iKs)
               endif
            endif
            ELSE
             IKS=MAXKS(I)-I+J
             if (iKs.gt.MaxKs(I-1)) then
               if (DABS(DESNO(iKs)).GT.tol) then
                  iKc=iKc+1
C                  WRITE(mufile2,*)I,J,DESNO(IKS)
C                  WRITE(IIZLAZ,*)I,J,DESNO(IKS)
C                  Kc(iKc)=DESNO(iKs)
               endif
            endif
           ENDIF
C         write(iizlaz,200)i,j,vredn
  10  continue
  20  continue


C	WRITE(IIZLAZ,*) JEDN
C	WRITE(IIZLAZ,*) iKc
c      rewind (mufile2)        
	WRITE(mufile2) JEDN
	WRITE(mufile2) iKc
c	goto 50

      iKc=0

C
CE Subroutine FULLMATR IS USED FOR PRINTING FULL MATRIX
C
 
 
 
      
	DO 55 I=1,JEDN
       DO 45 J=1,JEDN
        if(i.eq.1.and.j.eq.1) then
	     iKc=iKc+1
          WRITE(mufile2)I,J,ALEVO(1)
		 goto 45 
	  endif

            IF(J.GE.I) THEN
            iKs=MaxKs(j)-j+i
            if (iKs.gt.MaxKs(j-1)) then
               if (DABS(ALEVO(iKs)).GT.tol) then
                  iKc=iKc+1
c                  WRITE(IIZLAZ,*)I,J,ALEVO(IKS)
                  WRITE(mufile2)I,J,ALEVO(IKS)
C                  Kc(iKc)=ALEVO(iKs)
               endif
            endif
            ELSE
             IKS=MAXKS(I)-I+J
             if (iKs.gt.MaxKs(I-1)) then
               if (DABS(DESNO(iKs)).GT.tol) then
                  iKc=iKc+1
c                  WRITE(IIZLAZ,*)I,J,DESNO(IKS)
                  WRITE(mufile2)I,J,DESNO(IKS)
C                  Kc(iKc)=DESNO(iKs)
               endif
            endif
           ENDIF
C         write(iizlaz,200)i,j,vredn
   45 continue
   55 continue

 200  FORMAT (2I5,1PD20.12)

    
 

 50     DO I=1,JEDN
c        WRITE(IIZLAZ,*) SILE(I)
        WRITE(mufile2) SILE(I)
	ENDDO

c      rewind (mufile2)        
c	WRITE(mufile2) JEDN
c	WRITE(mufile2) iKc

c      REWIND (UNIT=iizlaz, IOSTAT=IOS)
c     REWIND (IIZLAZ)
c	WRITE(IIZLAZ,*) JEDN
c	WRITE(IIZLAZ,*) iKc


      CLOSE (mufile2)
       
      END
C=======================================================================
