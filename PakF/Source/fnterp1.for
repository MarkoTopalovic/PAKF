C======================================================================
      SUBROUTINE FNTERP1(R,S,NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,
     &ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,NBREL,IIZLAZ,HFI2,HFI3,NEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    Subroutine FNTERP is used for integration 2D finite element
C          
C      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
C     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /VISKOZ/ AMI,INDAMI
C      COMMON /PENALL/ PENALT,PRESS
       COMMON /SHEAR2/ SHEAR(2)

      DIMENSION V2(9),V3(9),T(9),ZVHR(9),ZVHS(9),FI(9)
      DIMENSION H(9),ZVHX(9),ZVHY(9)
      DIMENSION HP(4),X(9),Y(9)
      DIMENSION TT21(36+9)

      DIMENSION AJ(2,2)
      DIMENSION VTANG(2,9),AN(2),VSHX(2),VSHY(2)
	DIMENSION NEL(NDIM+1,*)


  5   DO 10 I=1,NDIM
      V2(I)=TT21(I)
      V3(I)=TT21(I+NDIM)
      T(I)=TT21(I+2*NDIM+4)
      FI(I)=TT21(I+3*NDIM+4)
  10  CONTINUE
      

      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIM.GT.4) THEN
      H(9)=RR*SS
      H(8)=0.5*RP*SS-0.5*H(9)
      H(7)=0.5*RR*SM-0.5*H(9)
      H(6)=0.5*RM*SS-0.5*H(9)
      H(5)=0.5*RR*SP-0.5*H(9)
      H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
      H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
      H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
      H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSE
      H(4)=0.25*RP*SM
      H(3)=0.25*RM*SM
      H(2)=0.25*RM*SP
      H(1)=0.25*RP*SP

      ENDIF


      HP(4)=0.25*RP*SM
      HP(3)=0.25*RM*SM
      HP(2)=0.25*RM*SP
      HP(1)=0.25*RP*SP

      IF (NDIM.EQ.4) HP(1)=1.D0 


      IF (NDIM.GT.4) THEN
      ZVHR(9)=-2.D0*R*SS
      ZVHS(9)=-2.D0*S*RR
      ZVHR(8)=0.5*SS-0.5*ZVHR(9)
      ZVHS(8)=-RP*S-0.5*ZVHS(9)
      ZVHR(7)=-R*SM-0.5*ZVHR(9)
      ZVHS(7)=-0.5*RR-0.5*ZVHS(9)
      ZVHR(6)=-0.5*SS-0.5*ZVHR(9)
      ZVHS(6)=-RM*S-0.5*ZVHS(9)
      ZVHR(5)=-R*SP-0.5*ZVHR(9)
      ZVHS(5)=0.5*RR-0.5*ZVHS(9)
      ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))-0.25*ZVHS(9)
      ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))-0.25*ZVHR(9)
      ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))-0.25*ZVHS(9)
      ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))-0.25*ZVHR(9)
      ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))-0.25*ZVHS(9)
      ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))-0.25*ZVHS(9)
      ELSE

      ZVHR(4)=0.25*SM
      ZVHS(4)=-0.25*RP
      ZVHR(3)=-0.25*SM
      ZVHS(3)=-0.25*RM
      ZVHR(2)=-0.25*SP
      ZVHS(2)=0.25*RM
      ZVHR(1)=0.25*SP
      ZVHS(1)=0.25*RP

      ENDIF
      AJ(1,1)=DOT(ZVHR,X,NDIM)
      AJ(1,2)=DOT(ZVHR,Y,NDIM)
      AJ(2,1)=DOT(ZVHS,X,NDIM)
      AJ(2,2)=DOT(ZVHS,Y,NDIM)

      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
C      WRITE(*,*)'DETERMINANTA MANJA OD NULE!!!'
      WRITE(*,*)'DETERMINANTE = ',DETJ
      WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT NUMBER  ',NBREL
      WRITE(IIZLAZ,*)'DETERMINANTE = ',DETJ
      WRITE(IIZLAZ,*)
     &'DETERMINANTE LESS THEN ZERO FOR ELEMENT NUMBER  ',NBREL
      WRITE(IIZLAZ,*)'      NODES:       COORDINATES (X,Y):  '
      DO I=1,NDIM
       WRITE(IIZLAZ,1000) I,X(I),Y(I)
      ENDDO
C      NN=NEL(2,NBREL)
C	NEL(2,NBREL)=NEL(4,NBREL)
C	NEL(4,NBREL)=NN

C	RETURN
C      STOP
      ENDIF

1000  FORMAT(I10,3X,2D13.5)

      DO 20 I=1,NDIM
      ZVHX(I)=(ZVHR(I)*AJ(2,2)-ZVHS(I)*AJ(1,2))/DETJ
      ZVHY(I)=(ZVHS(I)*AJ(1,1)-ZVHR(I)*AJ(2,1))/DETJ
  20  CONTINUE  

      IF (IUPWIN.EQ.1)
     &CALL INTEF2(X,Y,V2,V3,H,ZVHX,ZVHY,NDIM,AKT,GUSM,AMI)

      IF (NPARAM.NE.0) THEN
      VXX=0.D0
      VXY=0.D0
      VYX=0.D0
      VYY=0.D0
      PRESS=0.D0
      TX=0.D0
      TY=0.D0
          IF (DABS(R-1.).LT.1.E-5 .OR. DABS(R-(-1.)).LT.1.E-5) THEN
C          IF (R.EQ.1. .OR. R.EQ.-1.) THEN
           DETJS=DSQRT(AJ(2,1)**2+AJ(2,2)**2)
           ELX=R*AJ(2,2)/DETJS
           ELY=-R*AJ(2,1)/DETJS
           ETX=R*AJ(2,1)/DETJS
           ETY=R*AJ(2,2)/DETJS
          ELSE
           DETJS=DSQRT(AJ(1,2)**2+AJ(1,1)**2)
           ELX=-S*AJ(1,2)/DETJS
           ELY=S*AJ(1,1)/DETJS
           ETX=S*AJ(1,1)/DETJS
           ETY=S*AJ(1,2)/DETJS
          ENDIF

      DO 30 I=1,NDIM
      I1=I
      I2=I+NDIM
      I3=2*NDIM+I
      I4=2*NDIM+4+I
      VXX=VXX+ZVHX(I)*TT21(I1)
      VXY=VXY+ZVHY(I)*TT21(I1)
      VYX=VYX+ZVHX(I)*TT21(I2)
      VYY=VYY+ZVHY(I)*TT21(I2)
      TX=TX+ZVHX(I)*TT21(I4)
      TY=TY+ZVHY(I)*TT21(I4)
      IF (I.LE.4.AND. PENALT.LT.1.D0) PRESS=PRESS+HP(I)*TT21(I3)
  30  CONTINUE

      IF (PENALT.GT.1.D0) THEN
       IF(INDAX.EQ.1) THEN
        RR=0.D0
        DO I=1,4
         RR=RR+HP(I)*X(I)
        ENDDO
         HV2=DOT(H,V2,NDIM)
          PRESS=-PENALT*(VXX+VYY+HV2/RR)
       ELSE
          PRESS=-PENALT*(VXX+VYY)
       ENDIF
      ENDIF
     

       FS2=0.D0
       FS3=0.D0
       TAU=0.D0

      ENDIF

      IF (NPARAM.EQ.1) THEN
       FS2=AMI*(ELX*VXX+ELY*VXY)-ELX*PRESS
       FS3=AMI*(ELX*VYX+ELY*VYY)-ELY*PRESS
C FOR SHEAR STRESS CALCULATION
      CALL SHEARR(AN,TT21,VTANG,NDIM,ETX,ETY) 
      AN(1)=ELX
      AN(2)=ELY
      VSHX(1)=0.D0
      VSHX(2)=0.D0
      VSHY(1)=0.D0
      VSHY(2)=0.D0
      DO I=1,NDIM
       VSHX(1)=VSHX(1)+ZVHX(I)*VTANG(1,I)
       VSHX(2)=VSHX(2)+ZVHY(I)*VTANG(1,I)
       VSHY(1)=VSHY(1)+ZVHX(I)*VTANG(2,I)
       VSHY(2)=VSHY(2)+ZVHY(I)*VTANG(2,I)
      ENDDO
      SHEAR(1)=DOT(AN,VSHX,2)
      SHEAR(2)=DOT(AN,VSHY,2)
      ENDIF

      IF (NPARAM.EQ.2) THEN
C       TAU=-AMI*((VXX*ETX+VYX*ETY)*ELX+(VXY*ETX+VYY*ETY)*ELY)
C       FLUX=-AKT*(ELX*TX+ELY*TY)
       TAU=-AKT*(ELX*TX+ELY*TY)
      ENDIF
      



      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)

      HFI2=DOT(ZVHX,FI,NDIM)
      HFI3=DOT(ZVHY,FI,NDIM)


      ZVXT=DOT(ZVHX,T,NDIM)
      ZVYT=DOT(ZVHY,T,NDIM)
      ZVXV2=DOT(ZVHX,V2,NDIM)
      ZVYV3=DOT(ZVHY,V3,NDIM)
      ZVYV2=DOT(ZVHY,V2,NDIM)
      ZVXV3=DOT(ZVHX,V3,NDIM)

      END
C======================================================================
C==========================================================================
