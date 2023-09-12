C=========================================================================
      SUBROUTINE TetGen4_8(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(3,*)
	dimension nizf(*)
	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')
      
      k=1
	 write(iZ,*)'c tethra4 elements'
      do i=1,net
	 write(iz,'(11i5)') k,nel(1,i),nel(3,i),nel(8,i),nel(8,i),
     & nel(4,i),nel(4,i),nel(4,i),nel(4,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(2,i),nel(6,i),nel(6,i),
     & nel(3,i),nel(3,i),nel(3,i),nel(3,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(5,i),nel(8,i),nel(8,i),
     & nel(6,i),nel(6,i),nel(6,i),nel(6,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(3,i),nel(6,i),nel(6,i),
     & nel(8,i),nel(8,i),nel(8,i),nel(8,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(3,i),nel(6,i),nel(8,i),nel(8,i),
     & nel(7,i),nel(7,i),nel(7,i),nel(7,i)
	 k=k+1
	enddo 

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen4_6(CORD,NEL,NDIM,NPT,NET,ID,
     &niznew,nel1,ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(3,*)
	dimension nizf(*)
	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')
      
      k=1
	 write(iZ,*)'c tethra4 elements'
      do i=1,net
	 write(iz,'(11i5)') k,nel(1,i),nel(2,i),nel(4,i),nel(5,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(2,i),nel(3,i),nel(4,i),nel(5,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(2,i),nel(5,i),nel(6,i),nel(7,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(4,i),nel(5,i),nel(7,i),nel(8,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(2,i),nel(3,i),nel(5,i),nel(7,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(3,i),nel(4,i),nel(5,i),nel(7,i)
	 k=k+1
	enddo 


	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i5)') 3+(i-1)*6,nel(6,i),nel(5,i),nel(7,i)
	 write(iz,'(5i5)') 4+(i-1)*6,nel(7,i),nel(5,i),nel(8,i)
      enddo
      
c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 


      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen41_6(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(3,*)
	dimension nizf(*)
	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')
      

	 write(iZ,*)'c nodes'
      do i=1,net*6 
	  write(iz,'(6i5,3f10.4)') i+npt,1,1,1,0,1,0.d0,0.d0,0.d0
	enddo 


      k=1
	 write(iZ,*)'c tethra41 elements'
      do i=1,net
	 write(iz,'(11i5)') k,nel(1,i),nel(2,i),nel(4,i),nel(5,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(2,i),nel(3,i),nel(4,i),nel(5,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(2,i),nel(5,i),nel(6,i),nel(7,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(4,i),nel(5,i),nel(7,i),nel(8,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(2,i),nel(3,i),nel(5,i),nel(7,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(3,i),nel(4,i),nel(5,i),nel(7,i),npt+k
	 k=k+1
	enddo 
      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen41(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,ccord,IID
     &,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(3,*)
	dimension nizf(*)
	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')
      

	 write(iZ,*)'c nodes'
      do i=1,net*5 
	  write(iz,'(6i5,3f10.4)') i+npt,1,1,1,0,1,0.d0,0.d0,0.d0
	enddo 


      k=1
	 write(iZ,*)'c tethra41 elements'
      do i=1,net
	 write(iz,'(11i5)') k,nel(1,i),nel(3,i),nel(8,i),nel(4,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(2,i),nel(6,i),nel(3,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(5,i),nel(8,i),nel(6,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(3,i),nel(6,i),nel(8,i),npt+k
	 k=k+1
	 write(iz,'(11i5)') k,nel(3,i),nel(6,i),nel(8,i),nel(7,i),npt+k
	 k=k+1
	enddo 

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen4(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,ccord,IID,
     &nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(3,*)
	dimension nizf(*)
	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')
      
      k=1
	 write(iZ,*)'c tethra4 elements'
      do i=1,net
	 write(iz,'(11i5)') k,nel(1,i),nel(3,i),nel(8,i),nel(4,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(2,i),nel(6,i),nel(3,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(5,i),nel(8,i),nel(6,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(1,i),nel(3,i),nel(6,i),nel(8,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel(3,i),nel(6,i),nel(8,i),nel(7,i)
	 k=k+1
	enddo 

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen10(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,ccord,IID
     &,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo



	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=0
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo




      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
		 N9=NPT
		 CCORD(1,N9)=0.5D0*(CORD(1,N1)+CORD(1,N2))
		 CCORD(2,N9)=0.5D0*(CORD(2,N1)+CORD(2,N2))
		 CCORD(3,N9)=0.5D0*(CORD(3,N1)+CORD(3,N2))
	     IID(1,N9)=IID(1,N1).AND.IID(1,N2)
	     IID(2,N9)=IID(2,N1).AND.IID(2,N2)
	     IID(3,N9)=IID(3,N1).AND.IID(3,N2)
	     IID(4,N9)=1
         NPT=NPT+1
		 N10=NPT
		 CCORD(1,N10)=0.5D0*(CORD(1,N2)+CORD(1,N3))
		 CCORD(2,N10)=0.5D0*(CORD(2,N2)+CORD(2,N3))
		 CCORD(3,N10)=0.5D0*(CORD(3,N2)+CORD(3,N3))
	     IID(1,N10)=IID(1,N2).AND.IID(1,N3)
	     IID(2,N10)=IID(2,N2).AND.IID(2,N3)
	     IID(3,N10)=IID(3,N2).AND.IID(3,N3)
	     IID(4,N10)=1
         NPT=NPT+1
		 N11=NPT
		 CCORD(1,N11)=0.5D0*(CORD(1,N3)+CORD(1,N4))
		 CCORD(2,N11)=0.5D0*(CORD(2,N3)+CORD(2,N4))
		 CCORD(3,N11)=0.5D0*(CORD(3,N3)+CORD(3,N4))
	     IID(1,N11)=IID(1,N3).AND.IID(1,N4)
	     IID(2,N11)=IID(2,N3).AND.IID(2,N4)
	     IID(3,N11)=IID(3,N3).AND.IID(3,N4)
	     IID(4,N11)=1
         NPT=NPT+1
		 N12=NPT
		 CCORD(1,N12)=0.5D0*(CORD(1,N4)+CORD(1,N1))
		 CCORD(2,N12)=0.5D0*(CORD(2,N4)+CORD(2,N1))
		 CCORD(3,N12)=0.5D0*(CORD(3,N4)+CORD(3,N1))
	     IID(1,N12)=IID(1,N4).AND.IID(1,N1)
	     IID(2,N12)=IID(2,N4).AND.IID(2,N1)
	     IID(3,N12)=IID(3,N4).AND.IID(3,N1)
	     IID(4,N12)=1
         NPT=NPT+1
		 N13=NPT
		 CCORD(1,N13)=0.5D0*(CORD(1,N5)+CORD(1,N6))
		 CCORD(2,N13)=0.5D0*(CORD(2,N5)+CORD(2,N6))
		 CCORD(3,N13)=0.5D0*(CORD(3,N5)+CORD(3,N6))
	     IID(1,N13)=IID(1,N5).AND.IID(1,N6)
	     IID(2,N13)=IID(2,N5).AND.IID(2,N6)
	     IID(3,N13)=IID(3,N5).AND.IID(3,N6)
	     IID(4,N13)=1
         NPT=NPT+1
		 N14=NPT
		 CCORD(1,N14)=0.5D0*(CORD(1,N6)+CORD(1,N7))
		 CCORD(2,N14)=0.5D0*(CORD(2,N6)+CORD(2,N7))
		 CCORD(3,N14)=0.5D0*(CORD(3,N6)+CORD(3,N7))
	     IID(1,N14)=IID(1,N6).AND.IID(1,N7)
	     IID(2,N14)=IID(2,N6).AND.IID(2,N7)
	     IID(3,N14)=IID(3,N6).AND.IID(3,N7)
	     IID(4,N14)=1
         NPT=NPT+1
		 N15=NPT
		 CCORD(1,N15)=0.5D0*(CORD(1,N7)+CORD(1,N8))
		 CCORD(2,N15)=0.5D0*(CORD(2,N7)+CORD(2,N8))
		 CCORD(3,N15)=0.5D0*(CORD(3,N7)+CORD(3,N8))
	     IID(1,N15)=IID(1,N7).AND.IID(1,N8)
	     IID(2,N15)=IID(2,N7).AND.IID(2,N8)
	     IID(3,N15)=IID(3,N7).AND.IID(3,N8)
	     IID(4,N15)=1
         NPT=NPT+1
		 N16=NPT
		 CCORD(1,N16)=0.5D0*(CORD(1,N8)+CORD(1,N5))
		 CCORD(2,N16)=0.5D0*(CORD(2,N8)+CORD(2,N5))
		 CCORD(3,N16)=0.5D0*(CORD(3,N8)+CORD(3,N5))
	     IID(1,N16)=IID(1,N8).AND.IID(1,N5)
	     IID(2,N16)=IID(2,N8).AND.IID(2,N5)
	     IID(3,N16)=IID(3,N8).AND.IID(3,N5)
	     IID(4,N16)=1
         NPT=NPT+1
		 N17=NPT
		 CCORD(1,N17)=0.5D0*(CORD(1,N5)+CORD(1,N1))
		 CCORD(2,N17)=0.5D0*(CORD(2,N5)+CORD(2,N1))
		 CCORD(3,N17)=0.5D0*(CORD(3,N5)+CORD(3,N1))
	     IID(1,N17)=IID(1,N1).AND.IID(1,N5)
	     IID(2,N17)=IID(2,N1).AND.IID(2,N5)
	     IID(3,N17)=IID(3,N1).AND.IID(3,N5)
	     IID(4,N17)=1
         NPT=NPT+1
		 N18=NPT
		 CCORD(1,N18)=0.5D0*(CORD(1,N6)+CORD(1,N2))
		 CCORD(2,N18)=0.5D0*(CORD(2,N6)+CORD(2,N2))
		 CCORD(3,N18)=0.5D0*(CORD(3,N6)+CORD(3,N2))
	     IID(1,N18)=IID(1,N6).AND.IID(1,N2)
	     IID(2,N18)=IID(2,N6).AND.IID(2,N2)
	     IID(3,N18)=IID(3,N6).AND.IID(3,N2)
	     IID(4,N18)=1
         NPT=NPT+1
		 N19=NPT
		 CCORD(1,N19)=0.5D0*(CORD(1,N7)+CORD(1,N3))
		 CCORD(2,N19)=0.5D0*(CORD(2,N7)+CORD(2,N3))
		 CCORD(3,N19)=0.5D0*(CORD(3,N7)+CORD(3,N3))
	     IID(1,N19)=IID(1,N7).AND.IID(1,N3)
	     IID(2,N19)=IID(2,N7).AND.IID(2,N3)
	     IID(3,N19)=IID(3,N7).AND.IID(3,N3)
	     IID(4,N19)=1
         NPT=NPT+1
		 N20=NPT
		 CCORD(1,N20)=0.5D0*(CORD(1,N8)+CORD(1,N4))
		 CCORD(2,N20)=0.5D0*(CORD(2,N8)+CORD(2,N4))
		 CCORD(3,N20)=0.5D0*(CORD(3,N8)+CORD(3,N4))
	     IID(1,N20)=IID(1,N8).AND.IID(1,N4)
	     IID(2,N20)=IID(2,N8).AND.IID(2,N4)
	     IID(3,N20)=IID(3,N8).AND.IID(3,N4)
	     IID(4,N20)=1
         NPT=NPT+1
		 N21=NPT
       CCORD(1,N21)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N4))
       CCORD(2,N21)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N4))
       CCORD(3,N21)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4))
       IID(1,N21)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N3).AND.IID(1,N4)
       IID(2,N21)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N3).AND.IID(2,N4)
       IID(3,N21)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N3).AND.IID(3,N4)
       IID(4,N21)=1

         NPT=NPT+1
		 N22=NPT
       CCORD(1,N22)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N6)+CORD(1,N7))
       CCORD(2,N22)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N6)+CORD(2,N7))
       CCORD(3,N22)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N6)+CORD(3,N7))
       IID(1,N22)=IID(1,N2).AND.IID(1,N3).AND.IID(1,N6).AND.IID(1,N7)
       IID(2,N22)=IID(2,N2).AND.IID(2,N3).AND.IID(2,N6).AND.IID(2,N7)
       IID(3,N22)=IID(3,N2).AND.IID(3,N3).AND.IID(3,N6).AND.IID(3,N7)
       IID(4,N22)=1

         NPT=NPT+1
		 N23=NPT
       CCORD(1,N23)=0.25D0*(CORD(1,N5)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
       CCORD(2,N23)=0.25D0*(CORD(2,N5)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
       CCORD(3,N23)=0.25D0*(CORD(3,N5)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
       IID(1,N23)=IID(1,N5).AND.IID(1,N6).AND.IID(1,N7).AND.IID(1,N8)
       IID(2,N23)=IID(2,N5).AND.IID(2,N6).AND.IID(2,N7).AND.IID(2,N8)
       IID(3,N23)=IID(3,N5).AND.IID(3,N6).AND.IID(3,N7).AND.IID(3,N8)
       IID(4,N23)=1

         NPT=NPT+1
		 N24=NPT
       CCORD(1,N24)=0.25D0*(CORD(1,N1)+CORD(1,N4)+CORD(1,N5)+CORD(1,N8))
       CCORD(2,N24)=0.25D0*(CORD(2,N1)+CORD(2,N4)+CORD(2,N5)+CORD(2,N8))
       CCORD(3,N24)=0.25D0*(CORD(3,N1)+CORD(3,N4)+CORD(3,N5)+CORD(3,N8))
       IID(1,N24)=IID(1,N1).AND.IID(1,N4).AND.IID(1,N5).AND.IID(1,N8)
       IID(2,N24)=IID(2,N1).AND.IID(2,N4).AND.IID(2,N5).AND.IID(2,N8)
       IID(3,N24)=IID(3,N1).AND.IID(3,N4).AND.IID(3,N5).AND.IID(3,N8)
       IID(4,N24)=1

         NPT=NPT+1
		 N25=NPT
       CCORD(1,N25)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N5)+CORD(1,N6))
       CCORD(2,N25)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N5)+CORD(2,N6))
       CCORD(3,N25)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N5)+CORD(3,N6))
       IID(1,N25)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N5).AND.IID(1,N6)
       IID(2,N25)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N5).AND.IID(2,N6)
       IID(3,N25)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N5).AND.IID(3,N6)
       IID(4,N25)=1

         NPT=NPT+1
		 N26=NPT
       CCORD(1,N26)=0.25D0*(CORD(1,N3)+CORD(1,N4)+CORD(1,N7)+CORD(1,N8))
       CCORD(2,N26)=0.25D0*(CORD(2,N3)+CORD(2,N4)+CORD(2,N7)+CORD(2,N8))
       CCORD(3,N26)=0.25D0*(CORD(3,N3)+CORD(3,N4)+CORD(3,N7)+CORD(3,N8))
       IID(1,N26)=IID(1,N3).AND.IID(1,N4).AND.IID(1,N7).AND.IID(1,N8)
       IID(2,N26)=IID(2,N3).AND.IID(2,N4).AND.IID(2,N7).AND.IID(2,N8)
       IID(3,N26)=IID(3,N3).AND.IID(3,N4).AND.IID(3,N7).AND.IID(3,N8)
       IID(4,N26)=1


         
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
          NEL1 (14,NBREL)=N14
          NEL1 (15,NBREL)=N15
          NEL1 (16,NBREL)=N16
          NEL1 (17,NBREL)=N17
          NEL1 (18,NBREL)=N18
          NEL1 (19,NBREL)=N19
          NEL1 (20,NBREL)=N20
          NEL1 (21,NBREL)=N21
          NEL1 (22,NBREL)=N22
          NEL1 (23,NBREL)=N23
          NEL1 (24,NBREL)=N24
          NEL1 (25,NBREL)=N25
          NEL1 (26,NBREL)=N26
10    CONTINUE

      TOL=2.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=26

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (i10,5i5,3(x,1pe10.3))
      
      k=1
	 write(iZ,*)'c tethra elements'
      do i=1,net
	 write(iz,'(11i10)') k,nel1(1,i),nel1(3,i),nel1(8,i),nel1(4,i)
     &,nel1(21,i),nel1(26,i),nel1(24,i),nel1(12,i),nel1(11,i),nel1(20,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(1,i),nel1(2,i),nel1(6,i),nel1(3,i)
     &,nel1(9,i),nel1(18,i),nel1(25,i),nel1(21,i),nel1(10,i),nel1(22,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(1,i),nel1(5,i),nel1(8,i),nel1(6,i)
     &,nel1(17,i),nel1(16,i),nel1(24,i),nel1(25,i),nel1(13,i),nel1(23,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(1,i),nel1(3,i),nel1(6,i),nel1(8,i)
     &,nel1(21,i),nel1(22,i),nel1(25,i),nel1(24,i),nel1(26,i),nel1(23,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(3,i),nel1(6,i),nel1(8,i),nel1(7,i)
     &,nel1(22,i),nel1(23,i),nel1(26,i),nel1(19,i),nel1(14,i),nel1(15,i)
	 k=k+1
	enddo 



	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i10)') 3+(i-1)*5,nel(5,i),nel(8,i),nel(6,i)
	 write(iz,'(5i10)') 5+(i-1)*5,nel(6,i),nel(8,i),nel(7,i)
      enddo
      
c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 

      END
C=========================================================================
c1         3         8         4        21        26        24        12        11        20
C1         2         6         3         9        18        25        21        10        22
C1         5         8         6        17        16        24        25        13        23
C1         3         6         8        21        22        25        24        26        23
C3         6         8         7        22        23        26        19        14        15
C=========================================================================
C=========================================================================
      SUBROUTINE OctGen21(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,ccord,IID
     &,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(26,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo

	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=0
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo


      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
		 N9=NPT
		 CCORD(1,N9)=0.5D0*(CORD(1,N1)+CORD(1,N2))
		 CCORD(2,N9)=0.5D0*(CORD(2,N1)+CORD(2,N2))
		 CCORD(3,N9)=0.5D0*(CORD(3,N1)+CORD(3,N2))
	     IID(1,N9)=IID(1,N1).AND.IID(1,N2)
	     IID(2,N9)=IID(2,N1).AND.IID(2,N2)
	     IID(3,N9)=IID(3,N1).AND.IID(3,N2)
	     IID(4,N9)=1
         NPT=NPT+1
		 N10=NPT
		 CCORD(1,N10)=0.5D0*(CORD(1,N2)+CORD(1,N3))
		 CCORD(2,N10)=0.5D0*(CORD(2,N2)+CORD(2,N3))
		 CCORD(3,N10)=0.5D0*(CORD(3,N2)+CORD(3,N3))
	     IID(1,N10)=IID(1,N2).AND.IID(1,N3)
	     IID(2,N10)=IID(2,N2).AND.IID(2,N3)
	     IID(3,N10)=IID(3,N2).AND.IID(3,N3)
	     IID(4,N10)=1
         NPT=NPT+1
		 N11=NPT
		 CCORD(1,N11)=0.5D0*(CORD(1,N3)+CORD(1,N4))
		 CCORD(2,N11)=0.5D0*(CORD(2,N3)+CORD(2,N4))
		 CCORD(3,N11)=0.5D0*(CORD(3,N3)+CORD(3,N4))
	     IID(1,N11)=IID(1,N3).AND.IID(1,N4)
	     IID(2,N11)=IID(2,N3).AND.IID(2,N4)
	     IID(3,N11)=IID(3,N3).AND.IID(3,N4)
	     IID(4,N11)=1
         NPT=NPT+1
		 N12=NPT
		 CCORD(1,N12)=0.5D0*(CORD(1,N4)+CORD(1,N1))
		 CCORD(2,N12)=0.5D0*(CORD(2,N4)+CORD(2,N1))
		 CCORD(3,N12)=0.5D0*(CORD(3,N4)+CORD(3,N1))
	     IID(1,N12)=IID(1,N4).AND.IID(1,N1)
	     IID(2,N12)=IID(2,N4).AND.IID(2,N1)
	     IID(3,N12)=IID(3,N4).AND.IID(3,N1)
	     IID(4,N12)=1
         NPT=NPT+1
		 N13=NPT
		 CCORD(1,N13)=0.5D0*(CORD(1,N5)+CORD(1,N6))
		 CCORD(2,N13)=0.5D0*(CORD(2,N5)+CORD(2,N6))
		 CCORD(3,N13)=0.5D0*(CORD(3,N5)+CORD(3,N6))
	     IID(1,N13)=IID(1,N5).AND.IID(1,N6)
	     IID(2,N13)=IID(2,N5).AND.IID(2,N6)
	     IID(3,N13)=IID(3,N5).AND.IID(3,N6)
	     IID(4,N13)=1
         NPT=NPT+1
		 N14=NPT
		 CCORD(1,N14)=0.5D0*(CORD(1,N6)+CORD(1,N7))
		 CCORD(2,N14)=0.5D0*(CORD(2,N6)+CORD(2,N7))
		 CCORD(3,N14)=0.5D0*(CORD(3,N6)+CORD(3,N7))
	     IID(1,N14)=IID(1,N6).AND.IID(1,N7)
	     IID(2,N14)=IID(2,N6).AND.IID(2,N7)
	     IID(3,N14)=IID(3,N6).AND.IID(3,N7)
	     IID(4,N14)=1
         NPT=NPT+1
		 N15=NPT
		 CCORD(1,N15)=0.5D0*(CORD(1,N7)+CORD(1,N8))
		 CCORD(2,N15)=0.5D0*(CORD(2,N7)+CORD(2,N8))
		 CCORD(3,N15)=0.5D0*(CORD(3,N7)+CORD(3,N8))
	     IID(1,N15)=IID(1,N7).AND.IID(1,N8)
	     IID(2,N15)=IID(2,N7).AND.IID(2,N8)
	     IID(3,N15)=IID(3,N7).AND.IID(3,N8)
	     IID(4,N15)=1
         NPT=NPT+1
		 N16=NPT
		 CCORD(1,N16)=0.5D0*(CORD(1,N8)+CORD(1,N5))
		 CCORD(2,N16)=0.5D0*(CORD(2,N8)+CORD(2,N5))
		 CCORD(3,N16)=0.5D0*(CORD(3,N8)+CORD(3,N5))
	     IID(1,N16)=IID(1,N8).AND.IID(1,N5)
	     IID(2,N16)=IID(2,N8).AND.IID(2,N5)
	     IID(3,N16)=IID(3,N8).AND.IID(3,N5)
	     IID(4,N16)=1
         NPT=NPT+1
		 N17=NPT
		 CCORD(1,N17)=0.5D0*(CORD(1,N5)+CORD(1,N1))
		 CCORD(2,N17)=0.5D0*(CORD(2,N5)+CORD(2,N1))
		 CCORD(3,N17)=0.5D0*(CORD(3,N5)+CORD(3,N1))
	     IID(1,N17)=IID(1,N1).AND.IID(1,N5)
	     IID(2,N17)=IID(2,N1).AND.IID(2,N5)
	     IID(3,N17)=IID(3,N1).AND.IID(3,N5)
	     IID(4,N17)=1
         NPT=NPT+1
		 N18=NPT
		 CCORD(1,N18)=0.5D0*(CORD(1,N6)+CORD(1,N2))
		 CCORD(2,N18)=0.5D0*(CORD(2,N6)+CORD(2,N2))
		 CCORD(3,N18)=0.5D0*(CORD(3,N6)+CORD(3,N2))
	     IID(1,N18)=IID(1,N6).AND.IID(1,N2)
	     IID(2,N18)=IID(2,N6).AND.IID(2,N2)
	     IID(3,N18)=IID(3,N6).AND.IID(3,N2)
	     IID(4,N18)=1
         NPT=NPT+1
		 N19=NPT
		 CCORD(1,N19)=0.5D0*(CORD(1,N7)+CORD(1,N3))
		 CCORD(2,N19)=0.5D0*(CORD(2,N7)+CORD(2,N3))
		 CCORD(3,N19)=0.5D0*(CORD(3,N7)+CORD(3,N3))
	     IID(1,N19)=IID(1,N7).AND.IID(1,N3)
	     IID(2,N19)=IID(2,N7).AND.IID(2,N3)
	     IID(3,N19)=IID(3,N7).AND.IID(3,N3)
	     IID(4,N19)=1
         NPT=NPT+1
		 N20=NPT
		 CCORD(1,N20)=0.5D0*(CORD(1,N8)+CORD(1,N4))
		 CCORD(2,N20)=0.5D0*(CORD(2,N8)+CORD(2,N4))
		 CCORD(3,N20)=0.5D0*(CORD(3,N8)+CORD(3,N4))
	     IID(1,N20)=IID(1,N8).AND.IID(1,N4)
	     IID(2,N20)=IID(2,N8).AND.IID(2,N4)
	     IID(3,N20)=IID(3,N8).AND.IID(3,N4)
	     IID(4,N20)=1
         NPT=NPT+1
		 N21=NPT
      CCORD(1,N21)=0.125D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N4)
     &+CORD(1,N5)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
      CCORD(2,N21)=0.125D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N4)
     &+CORD(2,N5)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
      CCORD(3,N21)=0.125D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4)
     &+CORD(3,N5)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
      IID(1,N21)=0
      IID(2,N21)=0
      IID(3,N21)=0
      IID(4,N21)=1

         
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
          NEL1 (14,NBREL)=N14
          NEL1 (15,NBREL)=N15
          NEL1 (16,NBREL)=N16
          NEL1 (17,NBREL)=N17
          NEL1 (18,NBREL)=N18
          NEL1 (19,NBREL)=N19
          NEL1 (20,NBREL)=N20
          NEL1 (21,NBREL)=N21
10    CONTINUE

      TOL=1.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=21

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (6i5,3(1pe10.3))
      
      k=1
	 write(iZ,*)'c 21-node elements'
      do i=1,net
	 write(iz,'(13i5)') i,nel1(1,i),nel1(2,i),nel1(3,i),nel1(4,i)
     &,nel1(5,i),nel1(6,i),nel1(7,i),nel1(8,i)
       write(iz,'(13i5)')nel1(9,i),nel1(10,i),nel1(11,i),nel1(12,i)
     &,nel1(13,i),nel1(14,i),nel1(15,i),nel1(16,i),nel1(17,i)
     &,nel1(18,i),nel1(19,i),nel1(20,i),nel1(21,i)
	enddo 

	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
	do i=1,116
	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i5)') i,nel(5,i),nel(6,i),nel(7,i),nel(8,i)
      enddo
      
	do i=net-100+1,net
	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
      enddo

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen10_6(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(27,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo



	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=0
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo




      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
		 N9=NPT
		 CCORD(1,N9)=0.5D0*(CORD(1,N1)+CORD(1,N2))
		 CCORD(2,N9)=0.5D0*(CORD(2,N1)+CORD(2,N2))
		 CCORD(3,N9)=0.5D0*(CORD(3,N1)+CORD(3,N2))
	     IID(1,N9)=IID(1,N1).AND.IID(1,N2)
	     IID(2,N9)=IID(2,N1).AND.IID(2,N2)
	     IID(3,N9)=IID(3,N1).AND.IID(3,N2)
	     IID(4,N9)=1
         NPT=NPT+1
		 N10=NPT
		 CCORD(1,N10)=0.5D0*(CORD(1,N2)+CORD(1,N3))
		 CCORD(2,N10)=0.5D0*(CORD(2,N2)+CORD(2,N3))
		 CCORD(3,N10)=0.5D0*(CORD(3,N2)+CORD(3,N3))
	     IID(1,N10)=IID(1,N2).AND.IID(1,N3)
	     IID(2,N10)=IID(2,N2).AND.IID(2,N3)
	     IID(3,N10)=IID(3,N2).AND.IID(3,N3)
	     IID(4,N10)=1
         NPT=NPT+1
		 N11=NPT
		 CCORD(1,N11)=0.5D0*(CORD(1,N3)+CORD(1,N4))
		 CCORD(2,N11)=0.5D0*(CORD(2,N3)+CORD(2,N4))
		 CCORD(3,N11)=0.5D0*(CORD(3,N3)+CORD(3,N4))
	     IID(1,N11)=IID(1,N3).AND.IID(1,N4)
	     IID(2,N11)=IID(2,N3).AND.IID(2,N4)
	     IID(3,N11)=IID(3,N3).AND.IID(3,N4)
	     IID(4,N11)=1
         NPT=NPT+1
		 N12=NPT
		 CCORD(1,N12)=0.5D0*(CORD(1,N4)+CORD(1,N1))
		 CCORD(2,N12)=0.5D0*(CORD(2,N4)+CORD(2,N1))
		 CCORD(3,N12)=0.5D0*(CORD(3,N4)+CORD(3,N1))
	     IID(1,N12)=IID(1,N4).AND.IID(1,N1)
	     IID(2,N12)=IID(2,N4).AND.IID(2,N1)
	     IID(3,N12)=IID(3,N4).AND.IID(3,N1)
	     IID(4,N12)=1
         NPT=NPT+1
		 N13=NPT
		 CCORD(1,N13)=0.5D0*(CORD(1,N5)+CORD(1,N6))
		 CCORD(2,N13)=0.5D0*(CORD(2,N5)+CORD(2,N6))
		 CCORD(3,N13)=0.5D0*(CORD(3,N5)+CORD(3,N6))
	     IID(1,N13)=IID(1,N5).AND.IID(1,N6)
	     IID(2,N13)=IID(2,N5).AND.IID(2,N6)
	     IID(3,N13)=IID(3,N5).AND.IID(3,N6)
	     IID(4,N13)=1
         NPT=NPT+1
		 N14=NPT
		 CCORD(1,N14)=0.5D0*(CORD(1,N6)+CORD(1,N7))
		 CCORD(2,N14)=0.5D0*(CORD(2,N6)+CORD(2,N7))
		 CCORD(3,N14)=0.5D0*(CORD(3,N6)+CORD(3,N7))
	     IID(1,N14)=IID(1,N6).AND.IID(1,N7)
	     IID(2,N14)=IID(2,N6).AND.IID(2,N7)
	     IID(3,N14)=IID(3,N6).AND.IID(3,N7)
	     IID(4,N14)=1
         NPT=NPT+1
		 N15=NPT
		 CCORD(1,N15)=0.5D0*(CORD(1,N7)+CORD(1,N8))
		 CCORD(2,N15)=0.5D0*(CORD(2,N7)+CORD(2,N8))
		 CCORD(3,N15)=0.5D0*(CORD(3,N7)+CORD(3,N8))
	     IID(1,N15)=IID(1,N7).AND.IID(1,N8)
	     IID(2,N15)=IID(2,N7).AND.IID(2,N8)
	     IID(3,N15)=IID(3,N7).AND.IID(3,N8)
	     IID(4,N15)=1
         NPT=NPT+1
		 N16=NPT
		 CCORD(1,N16)=0.5D0*(CORD(1,N8)+CORD(1,N5))
		 CCORD(2,N16)=0.5D0*(CORD(2,N8)+CORD(2,N5))
		 CCORD(3,N16)=0.5D0*(CORD(3,N8)+CORD(3,N5))
	     IID(1,N16)=IID(1,N8).AND.IID(1,N5)
	     IID(2,N16)=IID(2,N8).AND.IID(2,N5)
	     IID(3,N16)=IID(3,N8).AND.IID(3,N5)
	     IID(4,N16)=1
         NPT=NPT+1
		 N17=NPT
		 CCORD(1,N17)=0.5D0*(CORD(1,N5)+CORD(1,N1))
		 CCORD(2,N17)=0.5D0*(CORD(2,N5)+CORD(2,N1))
		 CCORD(3,N17)=0.5D0*(CORD(3,N5)+CORD(3,N1))
	     IID(1,N17)=IID(1,N1).AND.IID(1,N5)
	     IID(2,N17)=IID(2,N1).AND.IID(2,N5)
	     IID(3,N17)=IID(3,N1).AND.IID(3,N5)
	     IID(4,N17)=1
         NPT=NPT+1
		 N18=NPT
		 CCORD(1,N18)=0.5D0*(CORD(1,N6)+CORD(1,N2))
		 CCORD(2,N18)=0.5D0*(CORD(2,N6)+CORD(2,N2))
		 CCORD(3,N18)=0.5D0*(CORD(3,N6)+CORD(3,N2))
	     IID(1,N18)=IID(1,N6).AND.IID(1,N2)
	     IID(2,N18)=IID(2,N6).AND.IID(2,N2)
	     IID(3,N18)=IID(3,N6).AND.IID(3,N2)
	     IID(4,N18)=1
         NPT=NPT+1
		 N19=NPT
		 CCORD(1,N19)=0.5D0*(CORD(1,N7)+CORD(1,N3))
		 CCORD(2,N19)=0.5D0*(CORD(2,N7)+CORD(2,N3))
		 CCORD(3,N19)=0.5D0*(CORD(3,N7)+CORD(3,N3))
	     IID(1,N19)=IID(1,N7).AND.IID(1,N3)
	     IID(2,N19)=IID(2,N7).AND.IID(2,N3)
	     IID(3,N19)=IID(3,N7).AND.IID(3,N3)
	     IID(4,N19)=1
         NPT=NPT+1
		 N20=NPT
		 CCORD(1,N20)=0.5D0*(CORD(1,N8)+CORD(1,N4))
		 CCORD(2,N20)=0.5D0*(CORD(2,N8)+CORD(2,N4))
		 CCORD(3,N20)=0.5D0*(CORD(3,N8)+CORD(3,N4))
	     IID(1,N20)=IID(1,N8).AND.IID(1,N4)
	     IID(2,N20)=IID(2,N8).AND.IID(2,N4)
	     IID(3,N20)=IID(3,N8).AND.IID(3,N4)
	     IID(4,N20)=1
         NPT=NPT+1
		 N21=NPT
       CCORD(1,N21)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N4))
       CCORD(2,N21)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N4))
       CCORD(3,N21)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4))
       IID(1,N21)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N3).AND.IID(1,N4)
       IID(2,N21)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N3).AND.IID(2,N4)
       IID(3,N21)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N3).AND.IID(3,N4)
       IID(4,N21)=1

         NPT=NPT+1
		 N22=NPT
       CCORD(1,N22)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N6)+CORD(1,N7))
       CCORD(2,N22)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N6)+CORD(2,N7))
       CCORD(3,N22)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N6)+CORD(3,N7))
       IID(1,N22)=IID(1,N2).AND.IID(1,N3).AND.IID(1,N6).AND.IID(1,N7)
       IID(2,N22)=IID(2,N2).AND.IID(2,N3).AND.IID(2,N6).AND.IID(2,N7)
       IID(3,N22)=IID(3,N2).AND.IID(3,N3).AND.IID(3,N6).AND.IID(3,N7)
       IID(4,N22)=1

         NPT=NPT+1
		 N23=NPT
       CCORD(1,N23)=0.25D0*(CORD(1,N5)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
       CCORD(2,N23)=0.25D0*(CORD(2,N5)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
       CCORD(3,N23)=0.25D0*(CORD(3,N5)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
       IID(1,N23)=IID(1,N5).AND.IID(1,N6).AND.IID(1,N7).AND.IID(1,N8)
       IID(2,N23)=IID(2,N5).AND.IID(2,N6).AND.IID(2,N7).AND.IID(2,N8)
       IID(3,N23)=IID(3,N5).AND.IID(3,N6).AND.IID(3,N7).AND.IID(3,N8)
       IID(4,N23)=1

         NPT=NPT+1
		 N24=NPT
       CCORD(1,N24)=0.25D0*(CORD(1,N1)+CORD(1,N4)+CORD(1,N5)+CORD(1,N8))
       CCORD(2,N24)=0.25D0*(CORD(2,N1)+CORD(2,N4)+CORD(2,N5)+CORD(2,N8))
       CCORD(3,N24)=0.25D0*(CORD(3,N1)+CORD(3,N4)+CORD(3,N5)+CORD(3,N8))
       IID(1,N24)=IID(1,N1).AND.IID(1,N4).AND.IID(1,N5).AND.IID(1,N8)
       IID(2,N24)=IID(2,N1).AND.IID(2,N4).AND.IID(2,N5).AND.IID(2,N8)
       IID(3,N24)=IID(3,N1).AND.IID(3,N4).AND.IID(3,N5).AND.IID(3,N8)
       IID(4,N24)=1

         NPT=NPT+1
		 N25=NPT
       CCORD(1,N25)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N5)+CORD(1,N6))
       CCORD(2,N25)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N5)+CORD(2,N6))
       CCORD(3,N25)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N5)+CORD(3,N6))
       IID(1,N25)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N5).AND.IID(1,N6)
       IID(2,N25)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N5).AND.IID(2,N6)
       IID(3,N25)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N5).AND.IID(3,N6)
       IID(4,N25)=1

         NPT=NPT+1
		 N26=NPT
       CCORD(1,N26)=0.25D0*(CORD(1,N3)+CORD(1,N4)+CORD(1,N7)+CORD(1,N8))
       CCORD(2,N26)=0.25D0*(CORD(2,N3)+CORD(2,N4)+CORD(2,N7)+CORD(2,N8))
       CCORD(3,N26)=0.25D0*(CORD(3,N3)+CORD(3,N4)+CORD(3,N7)+CORD(3,N8))
       IID(1,N26)=IID(1,N3).AND.IID(1,N4).AND.IID(1,N7).AND.IID(1,N8)
       IID(2,N26)=IID(2,N3).AND.IID(2,N4).AND.IID(2,N7).AND.IID(2,N8)
       IID(3,N26)=IID(3,N3).AND.IID(3,N4).AND.IID(3,N7).AND.IID(3,N8)
       IID(4,N26)=1

         NPT=NPT+1
		 N27=NPT
      CCORD(1,N27)=0.125D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N4)+
     &CORD(1,N5)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
      CCORD(2,N27)=0.125D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N4)+
     &CORD(2,N5)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
      CCORD(3,N27)=0.125D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4)+
     &CORD(3,N5)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
       IID(1,N27)=0
       IID(2,N27)=0
       IID(3,N27)=0
       IID(4,N27)=1

         
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
          NEL1 (14,NBREL)=N14
          NEL1 (15,NBREL)=N15
          NEL1 (16,NBREL)=N16
          NEL1 (17,NBREL)=N17
          NEL1 (18,NBREL)=N18
          NEL1 (19,NBREL)=N19
          NEL1 (20,NBREL)=N20
          NEL1 (21,NBREL)=N21
          NEL1 (22,NBREL)=N22
          NEL1 (23,NBREL)=N23
          NEL1 (24,NBREL)=N24
          NEL1 (25,NBREL)=N25
          NEL1 (26,NBREL)=N26
          NEL1 (27,NBREL)=N27
10    CONTINUE

      TOL=2.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=27

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (6i5,3(1pe10.3))
      
      k=1
	 write(iZ,*)'c tethra elements'


      do i=1,net
	 write(iz,'(11i5)') k,nel1(1,i),nel1(4,i),nel1(2,i),nel1(5,i)
     &,nel1(12,i),nel1(21,i),nel1(9,i),nel1(17,i),nel1(24,i),nel1(25,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(2,i),nel1(4,i),nel1(3,i),nel1(5,i)
     &,nel1(21,i),nel1(11,i),nel1(10,i),nel1(25,i),nel1(24,i),nel1(27,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(2,i),nel1(6,i),nel1(5,i),nel1(7,i)
     &,nel1(18,i),nel1(13,i),nel1(25,i),nel1(22,i),nel1(14,i),nel1(23,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(4,i),nel1(7,i),nel1(5,i),nel1(8,i)
     &,nel1(26,i),nel1(23,i),nel1(24,i),nel1(20,i),nel1(15,i),nel1(16,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(2,i),nel1(5,i),nel1(3,i),nel1(7,i)
     &,nel1(25,i),nel1(27,i),nel1(10,i),nel1(22,i),nel1(23,i),nel1(19,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(3,i),nel1(5,i),nel1(4,i),nel1(7,i)
     &,nel1(27,i),nel1(24,i),nel1(11,i),nel1(19,i),nel1(23,i),nel1(26,i)
	 k=k+1
	enddo 

c      do i=1,net
c	 write(iz,'(11i5)') k,nel1(1,i),nel1(2,i),nel1(4,i),nel1(5,i)
c    &,nel1(9,i),nel1(21,i),nel1(12,i),nel1(17,i),nel1(25,i),nel1(24,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(3,i),nel1(4,i),nel1(5,i)
c    &,nel1(10,i),nel1(11,i),nel1(21,i),nel1(25,i),nel1(27,i),nel1(24,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(5,i),nel1(6,i),nel1(7,i)
c     &,nel1(25,i),nel1(13,i),nel1(18,i),nel1(22,i),nel1(23,i),nel1(14,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(4,i),nel1(5,i),nel1(7,i),nel1(8,i)
c     &,nel1(24,i),nel1(23,i),nel1(26,i),nel1(20,i),nel1(16,i),nel1(15,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(3,i),nel1(5,i),nel1(7,i)
c     &,nel1(10,i),nel1(27,i),nel1(25,i),nel1(22,i),nel1(19,i),nel1(23,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(3,i),nel1(4,i),nel1(5,i),nel1(7,i)
c     &,nel1(11,i),nel1(24,i),nel1(27,i),nel1(19,i),nel1(26,i),nel1(23,i)
c	 k=k+1
c	enddo 


	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i5)') 3+(i-1)*6,nel(6,i),nel(5,i),nel(7,i)
	 write(iz,'(5i5)') 4+(i-1)*6,nel(7,i),nel(5,i),nel(8,i)
      enddo
      
c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen10_B5(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(14,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo



	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=0
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo




      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
		 N9=NPT
	 CCORD(1,N9)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N4)+CORD(1,N5))
	 CCORD(2,N9)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N4)+CORD(2,N5))
	 CCORD(3,N9)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N4)+CORD(3,N5))
	     IID(1,N9)=0
	     IID(2,N9)=0
	     IID(3,N9)=0
	     IID(4,N9)=1
         NPT=NPT+1
		 N10=NPT
	 CCORD(1,N10)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N4)+CORD(1,N5))
	 CCORD(2,N10)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N4)+CORD(2,N5))
	 CCORD(3,N10)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N4)+CORD(3,N5))
	     IID(1,N10)=0
	     IID(2,N11)=0
	     IID(3,N12)=0
	     IID(4,N13)=1
         NPT=NPT+1
		 N11=NPT
	 CCORD(1,N11)=0.25D0*(CORD(1,N2)+CORD(1,N5)+CORD(1,N6)+CORD(1,N7))
	 CCORD(2,N11)=0.25D0*(CORD(2,N2)+CORD(2,N5)+CORD(2,N6)+CORD(2,N7))
	 CCORD(3,N11)=0.25D0*(CORD(3,N2)+CORD(3,N5)+CORD(3,N6)+CORD(3,N7))
	     IID(1,N11)=0
	     IID(2,N11)=0
	     IID(3,N11)=0
	     IID(4,N11)=1
         NPT=NPT+1
		 N12=NPT
	 CCORD(1,N12)=0.25D0*(CORD(1,N4)+CORD(1,N5)+CORD(1,N7)+CORD(1,N8))
	 CCORD(2,N12)=0.25D0*(CORD(2,N4)+CORD(2,N5)+CORD(2,N7)+CORD(2,N8))
	 CCORD(3,N12)=0.25D0*(CORD(3,N4)+CORD(3,N5)+CORD(3,N7)+CORD(3,N8))
	     IID(1,N12)=0
	     IID(2,N12)=0
	     IID(3,N12)=0
	     IID(4,N12)=1
         NPT=NPT+1
		 N13=NPT
	 CCORD(1,N13)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N5)+CORD(1,N7))
	 CCORD(2,N13)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N5)+CORD(2,N7))
	 CCORD(3,N13)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N13)=0
	     IID(2,N13)=0
	     IID(3,N13)=0
	     IID(4,N13)=1
         NPT=NPT+1
		 N14=NPT
	 CCORD(1,N14)=0.25D0*(CORD(1,N3)+CORD(1,N4)+CORD(1,N5)+CORD(1,N7))
	 CCORD(2,N14)=0.25D0*(CORD(2,N3)+CORD(2,N4)+CORD(2,N5)+CORD(2,N7))
	 CCORD(3,N14)=0.25D0*(CORD(3,N3)+CORD(3,N4)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N14)=0
	     IID(2,N14)=0
	     IID(3,N14)=0
	     IID(4,N14)=1
         
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
          NEL1 (14,NBREL)=N14
10    CONTINUE

      TOL=2.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=14

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (i10,5i5,3(x,1pe13.5))
      
      k=1
	 write(iZ,*)'c tethra elements'


      do i=1,net
	 write(iz,'(11i10)') k,nel1(1,i),nel1(4,i),nel1(2,i),nel1(5,i)
     &,nel1(9,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(2,i),nel1(4,i),nel1(3,i),nel1(5,i)
     &,nel1(10,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(2,i),nel1(6,i),nel1(5,i),nel1(7,i)
     &,nel1(11,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(4,i),nel1(7,i),nel1(5,i),nel1(8,i)
     &,nel1(12,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(2,i),nel1(5,i),nel1(3,i),nel1(7,i)
     &,nel1(13,i)
	 k=k+1
	 write(iz,'(11i10)') k,nel1(3,i),nel1(5,i),nel1(4,i),nel1(7,i)
     &,nel1(14,i)
	 k=k+1
	enddo 

c      do i=1,net
c	 write(iz,'(11i5)') k,nel1(1,i),nel1(2,i),nel1(4,i),nel1(5,i)
c    &,nel1(9,i),nel1(21,i),nel1(12,i),nel1(17,i),nel1(25,i),nel1(24,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(3,i),nel1(4,i),nel1(5,i)
c    &,nel1(10,i),nel1(11,i),nel1(21,i),nel1(25,i),nel1(27,i),nel1(24,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(5,i),nel1(6,i),nel1(7,i)
c     &,nel1(25,i),nel1(13,i),nel1(18,i),nel1(22,i),nel1(23,i),nel1(14,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(4,i),nel1(5,i),nel1(7,i),nel1(8,i)
c     &,nel1(24,i),nel1(23,i),nel1(26,i),nel1(20,i),nel1(16,i),nel1(15,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(3,i),nel1(5,i),nel1(7,i)
c     &,nel1(10,i),nel1(27,i),nel1(25,i),nel1(22,i),nel1(19,i),nel1(23,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(3,i),nel1(4,i),nel1(5,i),nel1(7,i)
c     &,nel1(11,i),nel1(24,i),nel1(27,i),nel1(19,i),nel1(26,i),nel1(23,i)
c	 k=k+1
c	enddo 


	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure inlet' 
      do i=1,100
	 write(iz,'(5i10)') 3+(i-1)*6,nel(6,i),nel(5,i),nel(7,i)
	 write(iz,'(5i10)') 4+(i-1)*6,nel(7,i),nel(5,i),nel(8,i)
      enddo

	WRITE(IZ,*)'C boundary conditions pressure outlet' 
      do i=net-100+1,net
	 write(iz,'(5i10)') 1+(i-1)*6,nel(1,i),nel(4,i),nel(2,i)
	 write(iz,'(5i10)') 2+(i-1)*6,nel(2,i),nel(4,i),nel(3,i)
      enddo
      
c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen10_B55(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(13,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo



	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=0
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo




      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
		 N9=NPT
	 CCORD(1,N9)=0.25D0*(CORD(1,N1)+CORD(1,N3)+CORD(1,N4)+CORD(1,N8))
	 CCORD(2,N9)=0.25D0*(CORD(2,N1)+CORD(2,N3)+CORD(2,N4)+CORD(2,N8))
	 CCORD(3,N9)=0.25D0*(CORD(3,N1)+CORD(3,N3)+CORD(3,N4)+CORD(3,N8))
	     IID(1,N9)=0
	     IID(2,N9)=0
	     IID(3,N9)=0
	     IID(4,N9)=1
         NPT=NPT+1
		 N10=NPT
	 CCORD(1,N10)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N6))
	 CCORD(2,N10)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N6))
	 CCORD(3,N10)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N6))
	     IID(1,N10)=0
	     IID(2,N11)=0
	     IID(3,N12)=0
	     IID(4,N13)=1
         NPT=NPT+1
		 N11=NPT
	 CCORD(1,N11)=0.25D0*(CORD(1,N1)+CORD(1,N5)+CORD(1,N6)+CORD(1,N8))
	 CCORD(2,N11)=0.25D0*(CORD(2,N1)+CORD(2,N5)+CORD(2,N6)+CORD(2,N8))
	 CCORD(3,N11)=0.25D0*(CORD(3,N1)+CORD(3,N5)+CORD(3,N6)+CORD(3,N8))
	     IID(1,N11)=0
	     IID(2,N11)=0
	     IID(3,N11)=0
	     IID(4,N11)=1
         NPT=NPT+1
		 N12=NPT
	 CCORD(1,N12)=0.25D0*(CORD(1,N1)+CORD(1,N3)+CORD(1,N6)+CORD(1,N8))
	 CCORD(2,N12)=0.25D0*(CORD(2,N1)+CORD(2,N3)+CORD(2,N6)+CORD(2,N8))
	 CCORD(3,N12)=0.25D0*(CORD(3,N1)+CORD(3,N3)+CORD(3,N6)+CORD(3,N8))
	     IID(1,N12)=0
	     IID(2,N12)=0
	     IID(3,N12)=0
	     IID(4,N12)=1
         NPT=NPT+1
		 N13=NPT
	 CCORD(1,N13)=0.25D0*(CORD(1,N3)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
	 CCORD(2,N13)=0.25D0*(CORD(2,N3)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
	 CCORD(3,N13)=0.25D0*(CORD(3,N3)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
	     IID(1,N13)=0
	     IID(2,N13)=0
	     IID(3,N13)=0
	     IID(4,N13)=1
         
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
10    CONTINUE

      TOL=2.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=13

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (6i5,3(1pe10.3))
      
      k=1
	 write(iZ,*)'c tethra elements'
      do i=1,net
	 write(iz,'(11i5)') k,nel1(1,i),nel1(3,i),nel1(8,i),nel1(4,i)
     &,nel1(9,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(1,i),nel1(2,i),nel1(6,i),nel1(3,i)
     &,nel1(10,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(1,i),nel1(5,i),nel1(8,i),nel1(6,i)
     &,nel1(11,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(1,i),nel1(3,i),nel1(6,i),nel1(8,i)
     &,nel1(12,i)
	 k=k+1
	 write(iz,'(11i5)') k,nel1(3,i),nel1(6,i),nel1(8,i),nel1(7,i)
     &,nel1(13,i)
	 k=k+1
	enddo 



	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i5)') 3+(i-1)*5,nel(5,i),nel(8,i),nel(6,i)
	 write(iz,'(5i5)') 5+(i-1)*5,nel(6,i),nel(8,i),nel(7,i)
      enddo
      
c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen8FluxFace_6(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(32,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo



	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=1
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo


      t=1.d0/3.d0


      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
C===========================================================================
C FIRST TETRAHEDRAL ELEMENT 1,4,2,5
		 N9=NPT
		 CCORD(1,N9)=T*(CORD(1,N1)+CORD(1,N4)+CORD(1,N2))
		 CCORD(2,N9)=T*(CORD(2,N1)+CORD(2,N4)+CORD(2,N2))
		 CCORD(3,N9)=T*(CORD(3,N1)+CORD(3,N4)+CORD(3,N2))
	     IID(1,N9)=IID(1,N1).AND.IID(1,N4).AND.IID(1,N2)
	     IID(2,N9)=IID(2,N1).AND.IID(2,N4).AND.IID(2,N2)
	     IID(3,N9)=IID(3,N1).AND.IID(3,N4).AND.IID(3,N2)
	     IID(4,N9)=1
         NPT=NPT+1

		 N10=NPT
		 CCORD(1,N10)=T*(CORD(1,N1)+CORD(1,N4)+CORD(1,N5))
		 CCORD(2,N10)=T*(CORD(2,N1)+CORD(2,N4)+CORD(2,N5))
		 CCORD(3,N10)=T*(CORD(3,N1)+CORD(3,N4)+CORD(3,N5))
	     IID(1,N10)=IID(1,N1).AND.IID(1,N4).AND.IID(1,N5)
	     IID(2,N10)=IID(2,N1).AND.IID(2,N4).AND.IID(2,N5)
	     IID(3,N10)=IID(3,N1).AND.IID(3,N4).AND.IID(3,N5)
	     IID(4,N10)=1
         NPT=NPT+1

		 N11=NPT
		 CCORD(1,N11)=T*(CORD(1,N1)+CORD(1,N2)+CORD(1,N5))
		 CCORD(2,N11)=T*(CORD(2,N1)+CORD(2,N2)+CORD(2,N5))
		 CCORD(3,N11)=T*(CORD(3,N1)+CORD(3,N2)+CORD(3,N5))
	     IID(1,N11)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N5)
	     IID(2,N11)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N5)
	     IID(3,N11)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N5)
	     IID(4,N11)=1
         NPT=NPT+1

		 N12=NPT
		 CCORD(1,N12)=T*(CORD(1,N4)+CORD(1,N2)+CORD(1,N5))
		 CCORD(2,N12)=T*(CORD(2,N4)+CORD(2,N2)+CORD(2,N5))
		 CCORD(3,N12)=T*(CORD(3,N4)+CORD(3,N2)+CORD(3,N5))
	     IID(1,N12)=IID(1,N4).AND.IID(1,N2).AND.IID(1,N5)
	     IID(2,N12)=IID(2,N4).AND.IID(2,N2).AND.IID(2,N5)
	     IID(3,N12)=IID(3,N4).AND.IID(3,N2).AND.IID(3,N5)
	     IID(4,N12)=1
         NPT=NPT+1

C===========================================================================
C===========================================================================
C SECOND TETRAHEDRAL ELEMENT 2,4,3,5
		 N13=NPT
		 CCORD(1,N13)=T*(CORD(1,N2)+CORD(1,N4)+CORD(1,N3))
		 CCORD(2,N13)=T*(CORD(2,N2)+CORD(2,N4)+CORD(2,N3))
		 CCORD(3,N13)=T*(CORD(3,N2)+CORD(3,N4)+CORD(3,N3))
	     IID(1,N13)=IID(1,N2).AND.IID(1,N4).AND.IID(1,N3)
	     IID(2,N13)=IID(2,N2).AND.IID(2,N4).AND.IID(2,N3)
	     IID(3,N13)=IID(3,N2).AND.IID(3,N4).AND.IID(3,N3)
	     IID(4,N13)=1
         NPT=NPT+1

		 N14=NPT
		 CCORD(1,N14)=T*(CORD(1,N2)+CORD(1,N4)+CORD(1,N5))
		 CCORD(2,N14)=T*(CORD(2,N2)+CORD(2,N4)+CORD(2,N5))
		 CCORD(3,N14)=T*(CORD(3,N2)+CORD(3,N4)+CORD(3,N5))
	     IID(1,N14)=IID(1,N2).AND.IID(1,N4).AND.IID(1,N5)
	     IID(2,N14)=IID(2,N2).AND.IID(2,N4).AND.IID(2,N5)
	     IID(3,N14)=IID(3,N2).AND.IID(3,N4).AND.IID(3,N5)
	     IID(4,N14)=1
         NPT=NPT+1

		 N15=NPT
		 CCORD(1,N15)=T*(CORD(1,N2)+CORD(1,N3)+CORD(1,N5))
		 CCORD(2,N15)=T*(CORD(2,N2)+CORD(2,N3)+CORD(2,N5))
		 CCORD(3,N15)=T*(CORD(3,N2)+CORD(3,N3)+CORD(3,N5))
	     IID(1,N15)=IID(1,N2).AND.IID(1,N3).AND.IID(1,N5)
	     IID(2,N15)=IID(2,N2).AND.IID(2,N3).AND.IID(2,N5)
	     IID(3,N15)=IID(3,N2).AND.IID(3,N3).AND.IID(3,N5)
	     IID(4,N15)=1
         NPT=NPT+1

		 N16=NPT
		 CCORD(1,N16)=T*(CORD(1,N4)+CORD(1,N3)+CORD(1,N5))
		 CCORD(2,N16)=T*(CORD(2,N4)+CORD(2,N3)+CORD(2,N5))
		 CCORD(3,N16)=T*(CORD(3,N4)+CORD(3,N3)+CORD(3,N5))
	     IID(1,N16)=IID(1,N4).AND.IID(1,N3).AND.IID(1,N5)
	     IID(2,N16)=IID(2,N4).AND.IID(2,N3).AND.IID(2,N5)
	     IID(3,N16)=IID(3,N4).AND.IID(3,N3).AND.IID(3,N5)
	     IID(4,N16)=1
         NPT=NPT+1

C===========================================================================
C===========================================================================
C THIRD TETRAHEDRAL ELEMENT 2,6,5,7
		 N17=NPT
		 CCORD(1,N17)=T*(CORD(1,N2)+CORD(1,N6)+CORD(1,N5))
		 CCORD(2,N17)=T*(CORD(2,N2)+CORD(2,N6)+CORD(2,N5))
		 CCORD(3,N17)=T*(CORD(3,N2)+CORD(3,N6)+CORD(3,N5))
	     IID(1,N17)=IID(1,N2).AND.IID(1,N6).AND.IID(1,N5)
	     IID(2,N17)=IID(2,N2).AND.IID(2,N6).AND.IID(2,N5)
	     IID(3,N17)=IID(3,N2).AND.IID(3,N6).AND.IID(3,N5)
	     IID(4,N17)=1
         NPT=NPT+1

		 N18=NPT
		 CCORD(1,N18)=T*(CORD(1,N2)+CORD(1,N6)+CORD(1,N7))
		 CCORD(2,N18)=T*(CORD(2,N2)+CORD(2,N6)+CORD(2,N7))
		 CCORD(3,N18)=T*(CORD(3,N2)+CORD(3,N6)+CORD(3,N7))
	     IID(1,N18)=IID(1,N2).AND.IID(1,N6).AND.IID(1,N7)
	     IID(2,N18)=IID(2,N2).AND.IID(2,N6).AND.IID(2,N7)
	     IID(3,N18)=IID(3,N2).AND.IID(3,N6).AND.IID(3,N7)
	     IID(4,N18)=1
         NPT=NPT+1

		 N19=NPT
		 CCORD(1,N19)=T*(CORD(1,N2)+CORD(1,N5)+CORD(1,N7))
		 CCORD(2,N19)=T*(CORD(2,N2)+CORD(2,N5)+CORD(2,N7))
		 CCORD(3,N19)=T*(CORD(3,N2)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N19)=IID(1,N2).AND.IID(1,N5).AND.IID(1,N7)
	     IID(2,N19)=IID(2,N2).AND.IID(2,N5).AND.IID(2,N7)
	     IID(3,N19)=IID(3,N2).AND.IID(3,N5).AND.IID(3,N7)
	     IID(4,N19)=1
         NPT=NPT+1

		 N20=NPT
		 CCORD(1,N20)=T*(CORD(1,N6)+CORD(1,N5)+CORD(1,N7))
		 CCORD(2,N20)=T*(CORD(2,N6)+CORD(2,N5)+CORD(2,N7))
		 CCORD(3,N20)=T*(CORD(3,N6)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N20)=IID(1,N6).AND.IID(1,N5).AND.IID(1,N7)
	     IID(2,N20)=IID(2,N6).AND.IID(2,N5).AND.IID(2,N7)
	     IID(3,N20)=IID(3,N6).AND.IID(3,N5).AND.IID(3,N7)
	     IID(4,N20)=1
         NPT=NPT+1

C===========================================================================
C===========================================================================
C FOURTH TETRAHEDRAL ELEMENT 4,7,5,8

		 N21=NPT
		 CCORD(1,N21)=T*(CORD(1,N4)+CORD(1,N7)+CORD(1,N5))
		 CCORD(2,N21)=T*(CORD(2,N4)+CORD(2,N7)+CORD(2,N5))
		 CCORD(3,N21)=T*(CORD(3,N4)+CORD(3,N7)+CORD(3,N5))
	     IID(1,N21)=IID(1,N4).AND.IID(1,N7).AND.IID(1,N5)
	     IID(2,N21)=IID(2,N4).AND.IID(2,N7).AND.IID(2,N5)
	     IID(3,N21)=IID(3,N4).AND.IID(3,N7).AND.IID(3,N5)
	     IID(4,N21)=1
         NPT=NPT+1

		 N22=NPT
		 CCORD(1,N22)=T*(CORD(1,N4)+CORD(1,N7)+CORD(1,N8))
		 CCORD(2,N22)=T*(CORD(2,N4)+CORD(2,N7)+CORD(2,N8))
		 CCORD(3,N22)=T*(CORD(3,N4)+CORD(3,N7)+CORD(3,N8))
	     IID(1,N22)=IID(1,N4).AND.IID(1,N7).AND.IID(1,N8)
	     IID(2,N22)=IID(2,N4).AND.IID(2,N7).AND.IID(2,N8)
	     IID(3,N22)=IID(3,N4).AND.IID(3,N7).AND.IID(3,N8)
	     IID(4,N22)=1
         NPT=NPT+1

		 N23=NPT
		 CCORD(1,N23)=T*(CORD(1,N4)+CORD(1,N5)+CORD(1,N8))
		 CCORD(2,N23)=T*(CORD(2,N4)+CORD(2,N5)+CORD(2,N8))
		 CCORD(3,N23)=T*(CORD(3,N4)+CORD(3,N5)+CORD(3,N8))
	     IID(1,N23)=IID(1,N4).AND.IID(1,N5).AND.IID(1,N8)
	     IID(2,N23)=IID(2,N4).AND.IID(2,N5).AND.IID(2,N8)
	     IID(3,N23)=IID(3,N4).AND.IID(3,N5).AND.IID(3,N8)
	     IID(4,N23)=1
         NPT=NPT+1

		 N24=NPT
		 CCORD(1,N24)=T*(CORD(1,N7)+CORD(1,N5)+CORD(1,N8))
		 CCORD(2,N24)=T*(CORD(2,N7)+CORD(2,N5)+CORD(2,N8))
		 CCORD(3,N24)=T*(CORD(3,N7)+CORD(3,N5)+CORD(3,N8))
	     IID(1,N24)=IID(1,N7).AND.IID(1,N5).AND.IID(1,N8)
	     IID(2,N24)=IID(2,N7).AND.IID(2,N5).AND.IID(2,N8)
	     IID(3,N24)=IID(3,N7).AND.IID(3,N5).AND.IID(3,N8)
	     IID(4,N24)=1
         NPT=NPT+1

C===========================================================================
C===========================================================================
C FIFTH TETRAHEDRAL ELEMENT 2,5,3,7

		 N25=NPT
		 CCORD(1,N25)=T*(CORD(1,N2)+CORD(1,N5)+CORD(1,N3))
		 CCORD(2,N25)=T*(CORD(2,N2)+CORD(2,N5)+CORD(2,N3))
		 CCORD(3,N25)=T*(CORD(3,N2)+CORD(3,N5)+CORD(3,N3))
	     IID(1,N25)=IID(1,N2).AND.IID(1,N5).AND.IID(1,N3)
	     IID(2,N25)=IID(2,N2).AND.IID(2,N5).AND.IID(2,N3)
	     IID(3,N25)=IID(3,N2).AND.IID(3,N5).AND.IID(3,N3)
	     IID(4,N25)=1
         NPT=NPT+1

		 N26=NPT
		 CCORD(1,N26)=T*(CORD(1,N2)+CORD(1,N5)+CORD(1,N7))
		 CCORD(2,N26)=T*(CORD(2,N2)+CORD(2,N5)+CORD(2,N7))
		 CCORD(3,N26)=T*(CORD(3,N2)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N26)=IID(1,N2).AND.IID(1,N5).AND.IID(1,N7)
	     IID(2,N26)=IID(2,N2).AND.IID(2,N5).AND.IID(2,N7)
	     IID(3,N26)=IID(3,N2).AND.IID(3,N5).AND.IID(3,N7)
	     IID(4,N26)=1
         NPT=NPT+1

		 N27=NPT
		 CCORD(1,N27)=T*(CORD(1,N2)+CORD(1,N3)+CORD(1,N7))
		 CCORD(2,N27)=T*(CORD(2,N2)+CORD(2,N3)+CORD(2,N7))
		 CCORD(3,N27)=T*(CORD(3,N2)+CORD(3,N3)+CORD(3,N7))
	     IID(1,N27)=IID(1,N2).AND.IID(1,N3).AND.IID(1,N7)
	     IID(2,N27)=IID(2,N2).AND.IID(2,N3).AND.IID(2,N7)
	     IID(3,N27)=IID(3,N2).AND.IID(3,N3).AND.IID(3,N7)
	     IID(4,N27)=1
         NPT=NPT+1

		 N28=NPT
		 CCORD(1,N28)=T*(CORD(1,N5)+CORD(1,N3)+CORD(1,N7))
		 CCORD(2,N28)=T*(CORD(2,N5)+CORD(2,N3)+CORD(2,N7))
		 CCORD(3,N28)=T*(CORD(3,N5)+CORD(3,N3)+CORD(3,N7))
	     IID(1,N28)=IID(1,N5).AND.IID(1,N3).AND.IID(1,N7)
	     IID(2,N28)=IID(2,N5).AND.IID(2,N3).AND.IID(2,N7)
	     IID(3,N28)=IID(3,N5).AND.IID(3,N3).AND.IID(3,N7)
	     IID(4,N28)=1
         NPT=NPT+1

C===========================================================================
C===========================================================================
C SIXTH TETRAHEDRAL ELEMENT 3,5,4,7

		 N29=NPT
		 CCORD(1,N29)=T*(CORD(1,N3)+CORD(1,N5)+CORD(1,N4))
		 CCORD(2,N29)=T*(CORD(2,N3)+CORD(2,N5)+CORD(2,N4))
		 CCORD(3,N29)=T*(CORD(3,N3)+CORD(3,N5)+CORD(3,N4))
	     IID(1,N29)=IID(1,N3).AND.IID(1,N5).AND.IID(1,N4)
	     IID(2,N29)=IID(2,N3).AND.IID(2,N5).AND.IID(2,N4)
	     IID(3,N29)=IID(3,N3).AND.IID(3,N5).AND.IID(3,N4)
	     IID(4,N29)=1
         NPT=NPT+1

		 N30=NPT
		 CCORD(1,N30)=T*(CORD(1,N3)+CORD(1,N5)+CORD(1,N7))
		 CCORD(2,N30)=T*(CORD(2,N3)+CORD(2,N5)+CORD(2,N7))
		 CCORD(3,N30)=T*(CORD(3,N3)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N30)=IID(1,N3).AND.IID(1,N5).AND.IID(1,N7)
	     IID(2,N30)=IID(2,N3).AND.IID(2,N5).AND.IID(2,N7)
	     IID(3,N30)=IID(3,N3).AND.IID(3,N5).AND.IID(3,N7)
	     IID(4,N30)=1
         NPT=NPT+1

		 N31=NPT
		 CCORD(1,N31)=T*(CORD(1,N3)+CORD(1,N4)+CORD(1,N7))
		 CCORD(2,N31)=T*(CORD(2,N3)+CORD(2,N4)+CORD(2,N7))
		 CCORD(3,N31)=T*(CORD(3,N3)+CORD(3,N4)+CORD(3,N7))
	     IID(1,N31)=IID(1,N3).AND.IID(1,N4).AND.IID(1,N7)
	     IID(2,N31)=IID(2,N3).AND.IID(2,N4).AND.IID(2,N7)
	     IID(3,N31)=IID(3,N3).AND.IID(3,N4).AND.IID(3,N7)
	     IID(4,N31)=1
         NPT=NPT+1

		 N32=NPT
		 CCORD(1,N32)=T*(CORD(1,N5)+CORD(1,N4)+CORD(1,N7))
		 CCORD(2,N32)=T*(CORD(2,N5)+CORD(2,N4)+CORD(2,N7))
		 CCORD(3,N32)=T*(CORD(3,N5)+CORD(3,N4)+CORD(3,N7))
	     IID(1,N32)=IID(1,N5).AND.IID(1,N4).AND.IID(1,N7)
	     IID(2,N32)=IID(2,N5).AND.IID(2,N4).AND.IID(2,N7)
	     IID(3,N32)=IID(3,N5).AND.IID(3,N4).AND.IID(3,N7)
	     IID(4,N32)=1
c         NPT=NPT+1

C===========================================================================

         
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
          NEL1 (14,NBREL)=N14
          NEL1 (15,NBREL)=N15
          NEL1 (16,NBREL)=N16
          NEL1 (17,NBREL)=N17
          NEL1 (18,NBREL)=N18
          NEL1 (19,NBREL)=N19
          NEL1 (20,NBREL)=N20
          NEL1 (21,NBREL)=N21
          NEL1 (22,NBREL)=N22
          NEL1 (23,NBREL)=N23
          NEL1 (24,NBREL)=N24
          NEL1 (25,NBREL)=N25
          NEL1 (26,NBREL)=N26
          NEL1 (27,NBREL)=N27
          NEL1 (28,NBREL)=N28
          NEL1 (29,NBREL)=N29
          NEL1 (30,NBREL)=N30
          NEL1 (31,NBREL)=N31
          NEL1 (32,NBREL)=N32
10    CONTINUE

      TOL=2.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=32

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (I10,5i5,3(X,1pe13.5))
      
      k=1
	 write(iZ,*)'c tethra elements'


      do i=1,net
	 write(iz,'(9i10)') k,nel1(1,i),nel1(4,i),nel1(2,i),nel1(5,i)
     &,nel1(9,i),nel1(10,i),nel1(11,i),nel1(12,i)
	 k=k+1
	 write(iz,'(9i10)') k,nel1(2,i),nel1(4,i),nel1(3,i),nel1(5,i)
     &,nel1(13,i),nel1(14,i),nel1(15,i),nel1(16,i)
	 k=k+1
	 write(iz,'(9i10)') k,nel1(2,i),nel1(6,i),nel1(5,i),nel1(7,i)
     &,nel1(17,i),nel1(18,i),nel1(19,i),nel1(20,i)
     	 k=k+1
	 write(iz,'(9i10)') k,nel1(4,i),nel1(7,i),nel1(5,i),nel1(8,i)
     &,nel1(21,i),nel1(22,i),nel1(23,i),nel1(24,i)
	 k=k+1
	 write(iz,'(9i10)') k,nel1(2,i),nel1(5,i),nel1(3,i),nel1(7,i)
     &,nel1(25,i),nel1(26,i),nel1(27,i),nel1(28,i)
	 k=k+1
	 write(iz,'(9i10)') k,nel1(3,i),nel1(5,i),nel1(4,i),nel1(7,i)
     &,nel1(29,i),nel1(30,i),nel1(31,i),nel1(32,i)
	 k=k+1
	enddo 

      STOP

	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i5)') 3+(i-1)*6,nel(6,i),nel(5,i),nel(7,i)
	 write(iz,'(5i5)') 4+(i-1)*6,nel(7,i),nel(5,i),nel(8,i)
      enddo
      
c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
      SUBROUTINE TetGen11_6(CORD,NEL,NDIM,NPT,NET,ID,niznew,nel1,
     &ccord,IID,nizf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*),ccord(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NIZNEW(*),nel1(33,*),IID(4,*)
	dimension nizf(*)


      do nbrel=1,net
	   do ii=1,ndim
	      nel1(ii,nbrel)=nel(ii,nbrel)
	    enddo
	enddo



	do i=1,npt
	  ccord(1,i)=cord(1,i)
	  ccord(2,i)=cord(2,i)
	  ccord(3,i)=cord(3,i)
	  IID(1,I)=1
	  IF (ID(1,I).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I).NE.0) IID(3,I)=0
	  IID(4,I)=0
	enddo



	do i=1,116
	  IID(4,I)=1
	  IID(4,npt-I+1)=1
	  IID(1,I)=1
	  IF (ID(1,I+116).NE.0) IID(1,I)=0
	  IID(2,I)=1
	  IF (ID(2,I+116).NE.0) IID(2,I)=0
	  IID(3,I)=1
	  IF (ID(3,I+116).NE.0) IID(3,I)=0
	enddo




      DO 10 NBREL=1,NET
         N1=NEL1(1,NBREL)
         N2=NEL1(2,NBREL)
         N3=NEL1(3,NBREL)
         N4=NEL1(4,NBREL)
         N5=NEL1(5,NBREL)
         N6=NEL1(6,NBREL)
         N7=NEL1(7,NBREL)
         N8=NEL1(8,NBREL)
         NPT=NPT+1
		 N9=NPT
		 CCORD(1,N9)=0.5D0*(CORD(1,N1)+CORD(1,N2))
		 CCORD(2,N9)=0.5D0*(CORD(2,N1)+CORD(2,N2))
		 CCORD(3,N9)=0.5D0*(CORD(3,N1)+CORD(3,N2))
	     IID(1,N9)=IID(1,N1).AND.IID(1,N2)
	     IID(2,N9)=IID(2,N1).AND.IID(2,N2)
	     IID(3,N9)=IID(3,N1).AND.IID(3,N2)
	     IID(4,N9)=1
         NPT=NPT+1
		 N10=NPT
		 CCORD(1,N10)=0.5D0*(CORD(1,N2)+CORD(1,N3))
		 CCORD(2,N10)=0.5D0*(CORD(2,N2)+CORD(2,N3))
		 CCORD(3,N10)=0.5D0*(CORD(3,N2)+CORD(3,N3))
	     IID(1,N10)=IID(1,N2).AND.IID(1,N3)
	     IID(2,N10)=IID(2,N2).AND.IID(2,N3)
	     IID(3,N10)=IID(3,N2).AND.IID(3,N3)
	     IID(4,N10)=1
         NPT=NPT+1
		 N11=NPT
		 CCORD(1,N11)=0.5D0*(CORD(1,N3)+CORD(1,N4))
		 CCORD(2,N11)=0.5D0*(CORD(2,N3)+CORD(2,N4))
		 CCORD(3,N11)=0.5D0*(CORD(3,N3)+CORD(3,N4))
	     IID(1,N11)=IID(1,N3).AND.IID(1,N4)
	     IID(2,N11)=IID(2,N3).AND.IID(2,N4)
	     IID(3,N11)=IID(3,N3).AND.IID(3,N4)
	     IID(4,N11)=1
         NPT=NPT+1
		 N12=NPT
		 CCORD(1,N12)=0.5D0*(CORD(1,N4)+CORD(1,N1))
		 CCORD(2,N12)=0.5D0*(CORD(2,N4)+CORD(2,N1))
		 CCORD(3,N12)=0.5D0*(CORD(3,N4)+CORD(3,N1))
	     IID(1,N12)=IID(1,N4).AND.IID(1,N1)
	     IID(2,N12)=IID(2,N4).AND.IID(2,N1)
	     IID(3,N12)=IID(3,N4).AND.IID(3,N1)
	     IID(4,N12)=1
         NPT=NPT+1
		 N13=NPT
		 CCORD(1,N13)=0.5D0*(CORD(1,N5)+CORD(1,N6))
		 CCORD(2,N13)=0.5D0*(CORD(2,N5)+CORD(2,N6))
		 CCORD(3,N13)=0.5D0*(CORD(3,N5)+CORD(3,N6))
	     IID(1,N13)=IID(1,N5).AND.IID(1,N6)
	     IID(2,N13)=IID(2,N5).AND.IID(2,N6)
	     IID(3,N13)=IID(3,N5).AND.IID(3,N6)
	     IID(4,N13)=1
         NPT=NPT+1
		 N14=NPT
		 CCORD(1,N14)=0.5D0*(CORD(1,N6)+CORD(1,N7))
		 CCORD(2,N14)=0.5D0*(CORD(2,N6)+CORD(2,N7))
		 CCORD(3,N14)=0.5D0*(CORD(3,N6)+CORD(3,N7))
	     IID(1,N14)=IID(1,N6).AND.IID(1,N7)
	     IID(2,N14)=IID(2,N6).AND.IID(2,N7)
	     IID(3,N14)=IID(3,N6).AND.IID(3,N7)
	     IID(4,N14)=1
         NPT=NPT+1
		 N15=NPT
		 CCORD(1,N15)=0.5D0*(CORD(1,N7)+CORD(1,N8))
		 CCORD(2,N15)=0.5D0*(CORD(2,N7)+CORD(2,N8))
		 CCORD(3,N15)=0.5D0*(CORD(3,N7)+CORD(3,N8))
	     IID(1,N15)=IID(1,N7).AND.IID(1,N8)
	     IID(2,N15)=IID(2,N7).AND.IID(2,N8)
	     IID(3,N15)=IID(3,N7).AND.IID(3,N8)
	     IID(4,N15)=1
         NPT=NPT+1
		 N16=NPT
		 CCORD(1,N16)=0.5D0*(CORD(1,N8)+CORD(1,N5))
		 CCORD(2,N16)=0.5D0*(CORD(2,N8)+CORD(2,N5))
		 CCORD(3,N16)=0.5D0*(CORD(3,N8)+CORD(3,N5))
	     IID(1,N16)=IID(1,N8).AND.IID(1,N5)
	     IID(2,N16)=IID(2,N8).AND.IID(2,N5)
	     IID(3,N16)=IID(3,N8).AND.IID(3,N5)
	     IID(4,N16)=1
         NPT=NPT+1
		 N17=NPT
		 CCORD(1,N17)=0.5D0*(CORD(1,N5)+CORD(1,N1))
		 CCORD(2,N17)=0.5D0*(CORD(2,N5)+CORD(2,N1))
		 CCORD(3,N17)=0.5D0*(CORD(3,N5)+CORD(3,N1))
	     IID(1,N17)=IID(1,N1).AND.IID(1,N5)
	     IID(2,N17)=IID(2,N1).AND.IID(2,N5)
	     IID(3,N17)=IID(3,N1).AND.IID(3,N5)
	     IID(4,N17)=1
         NPT=NPT+1
		 N18=NPT
		 CCORD(1,N18)=0.5D0*(CORD(1,N6)+CORD(1,N2))
		 CCORD(2,N18)=0.5D0*(CORD(2,N6)+CORD(2,N2))
		 CCORD(3,N18)=0.5D0*(CORD(3,N6)+CORD(3,N2))
	     IID(1,N18)=IID(1,N6).AND.IID(1,N2)
	     IID(2,N18)=IID(2,N6).AND.IID(2,N2)
	     IID(3,N18)=IID(3,N6).AND.IID(3,N2)
	     IID(4,N18)=1
         NPT=NPT+1
		 N19=NPT
		 CCORD(1,N19)=0.5D0*(CORD(1,N7)+CORD(1,N3))
		 CCORD(2,N19)=0.5D0*(CORD(2,N7)+CORD(2,N3))
		 CCORD(3,N19)=0.5D0*(CORD(3,N7)+CORD(3,N3))
	     IID(1,N19)=IID(1,N7).AND.IID(1,N3)
	     IID(2,N19)=IID(2,N7).AND.IID(2,N3)
	     IID(3,N19)=IID(3,N7).AND.IID(3,N3)
	     IID(4,N19)=1
         NPT=NPT+1
		 N20=NPT
		 CCORD(1,N20)=0.5D0*(CORD(1,N8)+CORD(1,N4))
		 CCORD(2,N20)=0.5D0*(CORD(2,N8)+CORD(2,N4))
		 CCORD(3,N20)=0.5D0*(CORD(3,N8)+CORD(3,N4))
	     IID(1,N20)=IID(1,N8).AND.IID(1,N4)
	     IID(2,N20)=IID(2,N8).AND.IID(2,N4)
	     IID(3,N20)=IID(3,N8).AND.IID(3,N4)
	     IID(4,N20)=1
         NPT=NPT+1
		 N21=NPT
       CCORD(1,N21)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N4))
       CCORD(2,N21)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N4))
       CCORD(3,N21)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4))
       IID(1,N21)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N3).AND.IID(1,N4)
       IID(2,N21)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N3).AND.IID(2,N4)
       IID(3,N21)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N3).AND.IID(3,N4)
       IID(4,N21)=1

         NPT=NPT+1
		 N22=NPT
       CCORD(1,N22)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N6)+CORD(1,N7))
       CCORD(2,N22)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N6)+CORD(2,N7))
       CCORD(3,N22)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N6)+CORD(3,N7))
       IID(1,N22)=IID(1,N2).AND.IID(1,N3).AND.IID(1,N6).AND.IID(1,N7)
       IID(2,N22)=IID(2,N2).AND.IID(2,N3).AND.IID(2,N6).AND.IID(2,N7)
       IID(3,N22)=IID(3,N2).AND.IID(3,N3).AND.IID(3,N6).AND.IID(3,N7)
       IID(4,N22)=1

         NPT=NPT+1
		 N23=NPT
       CCORD(1,N23)=0.25D0*(CORD(1,N5)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
       CCORD(2,N23)=0.25D0*(CORD(2,N5)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
       CCORD(3,N23)=0.25D0*(CORD(3,N5)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
       IID(1,N23)=IID(1,N5).AND.IID(1,N6).AND.IID(1,N7).AND.IID(1,N8)
       IID(2,N23)=IID(2,N5).AND.IID(2,N6).AND.IID(2,N7).AND.IID(2,N8)
       IID(3,N23)=IID(3,N5).AND.IID(3,N6).AND.IID(3,N7).AND.IID(3,N8)
       IID(4,N23)=1

         NPT=NPT+1
		 N24=NPT
       CCORD(1,N24)=0.25D0*(CORD(1,N1)+CORD(1,N4)+CORD(1,N5)+CORD(1,N8))
       CCORD(2,N24)=0.25D0*(CORD(2,N1)+CORD(2,N4)+CORD(2,N5)+CORD(2,N8))
       CCORD(3,N24)=0.25D0*(CORD(3,N1)+CORD(3,N4)+CORD(3,N5)+CORD(3,N8))
       IID(1,N24)=IID(1,N1).AND.IID(1,N4).AND.IID(1,N5).AND.IID(1,N8)
       IID(2,N24)=IID(2,N1).AND.IID(2,N4).AND.IID(2,N5).AND.IID(2,N8)
       IID(3,N24)=IID(3,N1).AND.IID(3,N4).AND.IID(3,N5).AND.IID(3,N8)
       IID(4,N24)=1

         NPT=NPT+1
		 N25=NPT
       CCORD(1,N25)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N5)+CORD(1,N6))
       CCORD(2,N25)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N5)+CORD(2,N6))
       CCORD(3,N25)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N5)+CORD(3,N6))
       IID(1,N25)=IID(1,N1).AND.IID(1,N2).AND.IID(1,N5).AND.IID(1,N6)
       IID(2,N25)=IID(2,N1).AND.IID(2,N2).AND.IID(2,N5).AND.IID(2,N6)
       IID(3,N25)=IID(3,N1).AND.IID(3,N2).AND.IID(3,N5).AND.IID(3,N6)
       IID(4,N25)=1

         NPT=NPT+1
		 N26=NPT
       CCORD(1,N26)=0.25D0*(CORD(1,N3)+CORD(1,N4)+CORD(1,N7)+CORD(1,N8))
       CCORD(2,N26)=0.25D0*(CORD(2,N3)+CORD(2,N4)+CORD(2,N7)+CORD(2,N8))
       CCORD(3,N26)=0.25D0*(CORD(3,N3)+CORD(3,N4)+CORD(3,N7)+CORD(3,N8))
       IID(1,N26)=IID(1,N3).AND.IID(1,N4).AND.IID(1,N7).AND.IID(1,N8)
       IID(2,N26)=IID(2,N3).AND.IID(2,N4).AND.IID(2,N7).AND.IID(2,N8)
       IID(3,N26)=IID(3,N3).AND.IID(3,N4).AND.IID(3,N7).AND.IID(3,N8)
       IID(4,N26)=1

         NPT=NPT+1
		 N27=NPT
      CCORD(1,N27)=0.125D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N3)+CORD(1,N4)+
     &CORD(1,N5)+CORD(1,N6)+CORD(1,N7)+CORD(1,N8))
      CCORD(2,N27)=0.125D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N3)+CORD(2,N4)+
     &CORD(2,N5)+CORD(2,N6)+CORD(2,N7)+CORD(2,N8))
      CCORD(3,N27)=0.125D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N3)+CORD(3,N4)+
     &CORD(3,N5)+CORD(3,N6)+CORD(3,N7)+CORD(3,N8))
       IID(1,N27)=0
       IID(2,N27)=0
       IID(3,N27)=0
       IID(4,N27)=1

         NPT=NPT+1
		 N28=NPT
	 CCORD(1,N28)=0.25D0*(CORD(1,N1)+CORD(1,N2)+CORD(1,N4)+CORD(1,N5))
	 CCORD(2,N28)=0.25D0*(CORD(2,N1)+CORD(2,N2)+CORD(2,N4)+CORD(2,N5))
	 CCORD(3,N28)=0.25D0*(CORD(3,N1)+CORD(3,N2)+CORD(3,N4)+CORD(3,N5))
	     IID(1,N28)=0
	     IID(2,N28)=0
	     IID(3,N28)=0
	     IID(4,N28)=1
         NPT=NPT+1
		 N29=NPT
	 CCORD(1,N29)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N4)+CORD(1,N5))
	 CCORD(2,N29)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N4)+CORD(2,N5))
	 CCORD(3,N29)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N4)+CORD(3,N5))
	     IID(1,N29)=0
	     IID(2,N29)=0
	     IID(3,N29)=0
	     IID(4,N29)=1
         NPT=NPT+1
		 N30=NPT
	 CCORD(1,N30)=0.25D0*(CORD(1,N2)+CORD(1,N5)+CORD(1,N6)+CORD(1,N7))
	 CCORD(2,N30)=0.25D0*(CORD(2,N2)+CORD(2,N5)+CORD(2,N6)+CORD(2,N7))
	 CCORD(3,N30)=0.25D0*(CORD(3,N2)+CORD(3,N5)+CORD(3,N6)+CORD(3,N7))
	     IID(1,N30)=0
	     IID(2,N30)=0
	     IID(3,N30)=0
	     IID(4,N30)=1
         NPT=NPT+1
		 N31=NPT
	 CCORD(1,N31)=0.25D0*(CORD(1,N4)+CORD(1,N5)+CORD(1,N7)+CORD(1,N8))
	 CCORD(2,N31)=0.25D0*(CORD(2,N4)+CORD(2,N5)+CORD(2,N7)+CORD(2,N8))
	 CCORD(3,N31)=0.25D0*(CORD(3,N4)+CORD(3,N5)+CORD(3,N7)+CORD(3,N8))
	     IID(1,N31)=0
	     IID(2,N31)=0
	     IID(3,N31)=0
	     IID(4,N31)=1
         NPT=NPT+1
		 N32=NPT
	 CCORD(1,N32)=0.25D0*(CORD(1,N2)+CORD(1,N3)+CORD(1,N5)+CORD(1,N7))
	 CCORD(2,N32)=0.25D0*(CORD(2,N2)+CORD(2,N3)+CORD(2,N5)+CORD(2,N7))
	 CCORD(3,N32)=0.25D0*(CORD(3,N2)+CORD(3,N3)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N32)=0
	     IID(2,N32)=0
	     IID(3,N32)=0
	     IID(4,N32)=1
         NPT=NPT+1
		 N33=NPT
	 CCORD(1,N33)=0.25D0*(CORD(1,N3)+CORD(1,N4)+CORD(1,N5)+CORD(1,N7))
	 CCORD(2,N33)=0.25D0*(CORD(2,N3)+CORD(2,N4)+CORD(2,N5)+CORD(2,N7))
	 CCORD(3,N33)=0.25D0*(CORD(3,N3)+CORD(3,N4)+CORD(3,N5)+CORD(3,N7))
	     IID(1,N33)=0
	     IID(2,N33)=0
	     IID(3,N33)=0
	     IID(4,N33)=1
         
      
          NEL1 (9,NBREL)=N9
          NEL1 (10,NBREL)=N10
          NEL1 (11,NBREL)=N11
          NEL1 (12,NBREL)=N12
          NEL1 (13,NBREL)=N13
          NEL1 (14,NBREL)=N14
          NEL1 (15,NBREL)=N15
          NEL1 (16,NBREL)=N16
          NEL1 (17,NBREL)=N17
          NEL1 (18,NBREL)=N18
          NEL1 (19,NBREL)=N19
          NEL1 (20,NBREL)=N20
          NEL1 (21,NBREL)=N21
          NEL1 (22,NBREL)=N22
          NEL1 (23,NBREL)=N23
          NEL1 (24,NBREL)=N24
          NEL1 (25,NBREL)=N25
          NEL1 (26,NBREL)=N26
          NEL1 (27,NBREL)=N27
          NEL1 (28,NBREL)=N28
          NEL1 (29,NBREL)=N29
          NEL1 (30,NBREL)=N30
          NEL1 (31,NBREL)=N31
          NEL1 (32,NBREL)=N32
          NEL1 (33,NBREL)=N33
10    CONTINUE

      TOL=2.5D-5

	IZ=34
	OPEN (IZ,FILE='MERGE.DAT')

      DO I=1,NPT
	  NIZNEW(I)=I
	ENDDO 

      K=0

      DO I=1,NPT-1
	IF (NIZNEW(I).GT.0) THEN
	K=K+1
	NIZNEW(I)=K
	ENDIF
	 DO J=I+1,NPT
	   IF (DABS(CCORD(1,I)-CCORD(1,J)).LT.TOL .AND.
     & DABS(CCORD(2,I)-CCORD(2,J)).LT.TOL.AND.
     & DABS(CCORD(3,I)-CCORD(3,J)).LT.TOL.AND.NIZNEW(I).GT.0)
     & NIZNEW(J)=-NIZNEW(I)
	 ENDDO  
	
	ENDDO  

	IF (NIZNEW(NPT).GT.0) THEN
	K=K+1
	NIZNEW(NPT)=K
	ENDIF
      
      ndim1=33

      DO I=1,NET
	 DO J=1,NDIM1
	   NEL1(J,I)=ABS(NIZNEW(NEL1(J,I)))
	 ENDDO
	ENDDO
      

      do i=1,npt
	   nizf(i)=i
	enddo

c      goto 90

      do 200 i=1,npt-1
	if (nizf(i).lt.0) goto 200
	  do 100 j=i+1,npt
	if (nizf(j).lt.0) goto 100
	   z1=ccord(3,nizf(i))
	   z2=ccord(3,nizf(j))
	    if (z1.gt.z2) then
	      nn=nizf(i)
	      nizf(i)=nizf(j)
	      nizf(j)=nn
	    endif
 100      continue
 200      continue


  90  continue  




	WRITE(IZ,*)'C NODES' 
	DO 30 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 30
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	WRITE(IZ,1000) ABS(NIZNEW(II)),
     &IDX,IDY,IDZ,IDP,1,(CCORD(J,II),J=1,3)
30    continue



c      NPT=K
      
1000  format (i10,5i5,3(x,1pe10.3))
      
      k=1
	 write(iZ,*)'c tethra elements'


      do i=1,net
	 write(iz,'(12i10)') k,nel1(1,i),nel1(4,i),nel1(2,i),nel1(5,i)
     &,nel1(12,i),nel1(21,i),nel1(9,i),nel1(17,i),nel1(24,i),nel1(25,i)
     &,nel1(28,i)
	 k=k+1
	 write(iz,'(12i10)') k,nel1(2,i),nel1(4,i),nel1(3,i),nel1(5,i)
     &,nel1(21,i),nel1(11,i),nel1(10,i),nel1(25,i),nel1(24,i),nel1(27,i)
     &,nel1(29,i)
	 k=k+1
	 write(iz,'(12i10)') k,nel1(2,i),nel1(6,i),nel1(5,i),nel1(7,i)
     &,nel1(18,i),nel1(13,i),nel1(25,i),nel1(22,i),nel1(14,i),nel1(23,i)
     &,nel1(30,i)
	 k=k+1
	 write(iz,'(12i10)') k,nel1(4,i),nel1(7,i),nel1(5,i),nel1(8,i)
     &,nel1(26,i),nel1(23,i),nel1(24,i),nel1(20,i),nel1(15,i),nel1(16,i)
     &,nel1(31,i)
	 k=k+1
	 write(iz,'(12i10)') k,nel1(2,i),nel1(5,i),nel1(3,i),nel1(7,i)
     &,nel1(25,i),nel1(27,i),nel1(10,i),nel1(22,i),nel1(23,i),nel1(19,i)
     &,nel1(32,i)
	 k=k+1
	 write(iz,'(12i10)') k,nel1(3,i),nel1(5,i),nel1(4,i),nel1(7,i)
     &,nel1(27,i),nel1(24,i),nel1(11,i),nel1(19,i),nel1(23,i),nel1(26,i)
     &,nel1(33,i)
	 k=k+1
	enddo 

c      do i=1,net
c	 write(iz,'(11i5)') k,nel1(1,i),nel1(2,i),nel1(4,i),nel1(5,i)
c    &,nel1(9,i),nel1(21,i),nel1(12,i),nel1(17,i),nel1(25,i),nel1(24,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(3,i),nel1(4,i),nel1(5,i)
c    &,nel1(10,i),nel1(11,i),nel1(21,i),nel1(25,i),nel1(27,i),nel1(24,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(5,i),nel1(6,i),nel1(7,i)
c     &,nel1(25,i),nel1(13,i),nel1(18,i),nel1(22,i),nel1(23,i),nel1(14,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(4,i),nel1(5,i),nel1(7,i),nel1(8,i)
c     &,nel1(24,i),nel1(23,i),nel1(26,i),nel1(20,i),nel1(16,i),nel1(15,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(2,i),nel1(3,i),nel1(5,i),nel1(7,i)
c     &,nel1(10,i),nel1(27,i),nel1(25,i),nel1(22,i),nel1(19,i),nel1(23,i)
c	 k=k+1
c	 write(iz,'(11i5)') k,nel1(3,i),nel1(4,i),nel1(5,i),nel1(7,i)
c     &,nel1(11,i),nel1(24,i),nel1(27,i),nel1(19,i),nel1(26,i),nel1(23,i)
c	 k=k+1
c	enddo 


	WRITE(IZ,*)'C numzad pressure' 
      do 31 i=1,116
	 write(iz,'(3i5,f10.4)') i,4,1,1.d0
   31 continue
      
c	do i=1,116
c	 write(iz,'(3i5,f10.4)') npt-i+1,4,1,0.d0
c      enddo
	
	WRITE(IZ,*)'C boundary conditions pressure' 
      do i=1,100
	 write(iz,'(5i5)') 3+(i-1)*6,nel(6,i),nel(5,i),nel(7,i)
	 write(iz,'(5i5)') 4+(i-1)*6,nel(7,i),nel(5,i),nel(8,i)
      enddo
      
	WRITE(IZ,*)'C boundary conditions velocity' 
	DO 87 I=1,npt
	  
	  II=nizf(I)
c	  II=I
	  if (niznew(ii).lt.0) goto 87
	  if (CCORD(3,II).ne.0.d0) goto 87
	  IDX=IID(1,II)
	  IDY=IID(2,II)
	  IDZ=IID(3,II)
	  IDP=IID(4,II)
	  if (idx.eq.0.and.idy.eq.0.and.idz.eq.0.) then
	 write(iz,'(3i5,f10.4)') ABS(NIZNEW(II)),3,1,1.d0
	 endif
87    continue

c	do i=net-100+1,net
c	 write(iz,'(5i5)') i,nel(1,i),nel(2,i),nel(3,i),nel(4,i)
c      enddo

      close(iz)
      stop 

      END
C=========================================================================
C=========================================================================
C=========================================================================
