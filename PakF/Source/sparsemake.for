C=======================================================================
      SUBROUTINE SPARSEmake(nel,id,net,npt,iks,ndim,neq,netip)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION nel(ndim+1,*),id(6,*)
	DIMENSION iks(2,*)
	dimension nn(24)

      


      do 5 nq=1,neq

	do 10 nbrel=1,net
	 k=0
	 do 20 i=1,ndim
        n=nel(i,nbrel)
	    do 30 j=1,netip
	      k=k+1
	      nn(k)=id(j,n)
30        continue   
20      continue
        do i=1,k
	    iks(1,nn(i))=1
	  enddo
10    continue
5     continue






	DO I=1,JEDN
       DO J=1,JEDN
           VREDN=0.D0
            IF(J.GE.I) THEN
            iKss=MaxKs(j)-j+i
            if (iKss.gt.MaxKs(j-1)) then
               if (DABS(ALEVO(iKs)).GT.1.D-12) then
                  iKc=iKc+1
	            vredn=alevo(iks)
c                  WRITE(IIZLAZ,*)I,J,ALEVO(IKS)
C                  Kc(iKc)=ALEVO(iKs)
               endif
            endif
            ELSE
             IKSs=MAXKS(I)-I+J
             if (iKss.gt.MaxKs(I-1)) then
               if (DABS(DESNO(iKs)).GT.1.D-12) then
                  iKc=iKc+1
                   vredn=desno(iks)
c                  WRITE(IIZLAZ,*)I,J,DESNO(IKS)
C                  Kc(iKc)=DESNO(iKs)
               endif
            endif
           ENDIF
C         write(iizlaz,200)i,j,vredn
       ENDDO
      ENDDO


	WRITE(IIZLAZ,*) JEDN
	WRITE(IIZLAZ,*) iKc
       
      iKc=0

C
CE Subroutine FULLMATR IS USED FOR PRINTING FULL MATRIX
C
 
 
 
 
 
      
	DO I=1,JEDN
       DO J=1,JEDN
           VREDN=0.D0
            IF(J.GE.I) THEN
            iKss=MaxKs(j)-j+i
            if (iKss.gt.MaxKs(j-1)) then
               if (DABS(ALEVO(iKs)).GT.1.D-12) then
                  iKc=iKc+1
	            vredn=alevo(iks)
                  WRITE(IIZLAZ,*)I,J,ALEVO(IKS)
C                  Kc(iKc)=ALEVO(iKs)
               endif
            endif
            ELSE
             IKSs=MAXKS(I)-I+J
             if (iKss.gt.MaxKs(I-1)) then
               if (DABS(DESNO(iKs)).GT.1.D-12) then
                  iKc=iKc+1
                   vredn=desno(ikss)
                  WRITE(IIZLAZ,*)I,J,DESNO(IKSs)
C                  Kc(iKc)=DESNO(iKs)
               endif
            endif
           ENDIF
C         write(iizlaz,200)i,j,vredn
       ENDDO
      ENDDO

 200  FORMAT (2I5,1PD20.12)

      DO I=1,JEDN
        WRITE(IIZLAZ,*) SILE(I)
	ENDDO

c      REWIND (UNIT=iizlaz, IOSTAT=IOS)
c     REWIND (IIZLAZ)
c	WRITE(IIZLAZ,*) JEDN
c	WRITE(IIZLAZ,*) iKc


      CLOSE (IIZLAZ)
       
      END
C=======================================================================
C==========================================================================
