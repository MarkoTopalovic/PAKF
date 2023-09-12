C==========================================================================
      SUBROUTINE genercurzio(II)
c======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

c======================================================================
C WRITING NODES	
      
      NR=50
      NL=10

      WRITE(II,*)'C NODES'
      
      D=3.5D0
      H=0.26D0
      
c      BIAS=0.9999D0
c      A1=R*(1.D0-BIAS)/(1.D0-BIAS**(NR))

      DX=D/NR
      DY=H/NL
      K=0
      DO I=0,NR
        DO J=0,NL 
	   Y=J*DY 
	   X=I*DX
          IDX=0
          IDY=0
         IF (I.EQ.0.OR.I.EQ.NR.OR.J.EQ.0) THEN
           IDX=1
           IDY=1
         ENDIF
	   K=K+1
	   WRITE(II,1000)K,IDX,IDY,1,1,1,X,Y,0.D0
	  ENDDO
      ENDDO
	  
c======================================================================
C WRITING ELEMENTS
      K=0
      WRITE(II,*)'C ELEMENTS'
	 DO I=0,NR-1
        DO J=1,NL 
	   K=K+1
	   I1=J+I*(NL+1)
	   WRITE(II,'(5I10)')K,I1,I1+NL+1,I1+NL+2,I1+1
	  ENDDO
      ENDDO
	  
c======================================================================
	
1000  FORMAT (I10,5I5,3(X,E13.5))

	
	close(ii)
	STOP
c======================================================================
      END
C==========================================================================
C==========================================================================
      SUBROUTINE MovingSurfaceCurzio2D(npt,kkorak,cord,GNODE,time)
c======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*),GNODE(2,6,*)
c======================================================================
C WRITING NODES	
      

      NR=50
      NL=10
      v=3.d-0
      sign=1.d0
      if (mod(kkorak/10,2).eq.0) sign=-1.d0
      v=v*sign
c      write(32,*) 'v= ',v
      
      DX=R/NR
      DY=AL/NL
      K=0
      DO I=1,NR+1
         gnode(2,1,i)=v
      ENDDO
      DO I=npt-nl-1,npt
         gnode(2,1,i)=v
      ENDDO
      DO I=1,npt,nl+1
         gnode(2,1,i)=v
      ENDDO

      DO I=nl+1,npt,nl+1
c         write(32,*) 'v pre= ',cord(2,i),gnode(2,2,i)*time
         cord(2,i)=cord(2,i)+gnode(2,2,i)*time
c         write(32,*) 'v posle= ',cord(2,i)
      ENDDO
      
      do j=0,nr
       do I=0,nl-1
       y=cord(2,(j+1)*(nl+1))
       dy=y/nl
       k=i+1+(nl+1)*j
         cord(2,k)=i*dy
       ENDDO
      ENDDO
    	  
      END
C==========================================================================
      SUBROUTINE PrintCurzio(kkorak,pres,time,cord,npt,nz)
c======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION pres(3,*),cord(3,*)
c======================================================================
C WRITING NODES	
      
      
      
      do i=1,npt,nz
        s=sqrt(pres(1,i)**2+pres(2,i)**2+pres(3,i)**2)
        write(32,'(i10,3e13.5)') kkorak,cord(1,i),cord(2,i),s
      enddo 
      
      return

      do i=1,npt,nz
       if(dabs(cord(1,i)).lt.1.d-5) then
        s=sqrt(pres(1,i)**2+pres(2,i)**2+pres(3,i)**2)
        write(32,'(i10,3e13.5)') kkorak,cord(1,i),cord(2,i),s
       endif
       if(dabs(cord(2,i)).lt.1.d-5) then
        s=sqrt(pres(1,i)**2+pres(2,i)**2+pres(3,i)**2)
        write(32,'(i10,3e13.5)') kkorak,cord(1,i),cord(2,i),s
       endif
      enddo 



      NR=50
      NL=10
      lleft=1
      lcenter=1+(nl+1)*(nr)/2
      lright=1+(nl+1)*nr
      
      s1=sqrt(pres(1,lleft)**2+pres(2,lleft)**2+pres(3,lleft)**2)
      s2=sqrt(pres(1,lcenter)**2+pres(2,lcenter)**2+pres(3,lcenter)**2)
      s3=sqrt(pres(1,lright)**2+pres(2,lright)**2+pres(3,lright)**2)
     
      write(32,'(i10,3e13.5)') kkorak,s1,s2,s3
      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE MovSurfaceCurzio3D(npt,kkorak,ccord,GNODE,time,ID,nz)
c======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common /orbital/ ORBRAD,ORBFRE

      DIMENSION CCORD(3,*),GNODE(2,6,*)
      DIMENSION ID(6,*)
c======================================================================
C WRITING NODES	
      

      t=time*kkorak
      PI=4.D0*DATAN(1.D0)
      f=ORBFRE/60.d0
      w=2.d0*PI*f
      fi=w*t
      R=ORBRAD

      vx=R*w*sin(fi)
      vy=-R*w*cos(fi)
c      write(32,*) 'v= ',v
      do i=1,npt
       if(id(1,i).eq.0.and.id(2,i).eq.0) then
        gnode(2,1,i)=vx
        gnode(2,2,i)=vy
       endif
      enddo

      do i=nz,npt,nz
        vz=gnode(2,3,i)
        ccord(3,i)=ccord(3,i)+vz*time
        if (ccord(3,i).lt.0.d0) ccord(3,i)=0.26d0
c       if(id(1,i).ne.0.and.id(2,i).ne.0) then
c        vxx=gnode(2,1,i)
c        vyy=gnode(2,2,i)
c        ccord(1,i)=ccord(1,i)+vxx*time
c        ccord(2,i)=ccord(2,i)+vyy*time
c       endif
      enddo

      do j=0,npt/nz-1
       n=nz+nz*j
       z=ccord(3,n)
       dz=z/(nz-1)
       do i=2,nz-1
        n=i+nz*j
        ccord(3,n)=dz*(i-1)
       enddo
      enddo
    	  
      END
C==========================================================================
