C======================================================================
C MODIFIKOVAO ALEKSANDAR NIKOLIC, AVGUST 2016-FEBRUAR 2017.
C======================================================================
      SUBROUTINE FNTERPT(R,S,NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,
     &ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVXOM,ZVYOM,NBREL,
     &IIZLAZ,GNODE,NEL,CORD,ID,TT210)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    Subroutine FNTERPT is used for integration 2D finite element
C       
	COMMON /TURB/ ITURB   
	COMMON /SHEAR2/ SHEAR(2)
      DIMENSION V2(9),V3(9),T(9),ZVHR(9),ZVHS(9)
      DIMENSION H(9),ZVHX(9),ZVHY(9)
      DIMENSION HP(4),X(9),Y(9)
      DIMENSION TT21(54)
C--------------------------------------------------------------------
C     K-OMEGA TURB. MODEL
C--------------------------------------------------------------------
      DIMENSION AAKT(9)
      DIMENSION OMT(9)
      DIMENSION YP(9)
      DIMENSION YZV(9)
C--------------------------------------------------------------------
      DIMENSION AJ(2,2)
      DIMENSION GNODE(2,7,*)
      DIMENSION NEL(NDIM+1,*)
      DIMENSION TT210(54)
      DIMENSION CORD(3,*)
      DIMENSION ID(7,*)
      DIMENSION VTANG(2,9)
      DIMENSION AN(2),V(2),VSHX(2),VSHY(2),YPLUS(2),VPLUS(2)
      
      CMI=0.09
      GRANICA=11.225 
      GGRANICA=500.0 
      AKP=0.1
      E=9.793
      AKAPA=0.42
C    TT210(I+4*NDIM)**0.5  CORD(2,NEL(I,NBREL) GNODE(1,6,NEL(I,NBREL))

C LOKALNE KOORDINATE ELEMENTA 

    

      
      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S
      
      DO  8 I=1,NDIM
      H(4)=0.25*RP*SM
      H(3)=0.25*RM*SM
      H(2)=0.25*RM*SP
      H(1)=0.25*RP*SP
      T(I)=TT21(I+2*NDIM+4)
C DODATO KT I OMT (KIN. ENER. TURBULENCIJE I SPEC. DISIPACIJA KIN. ENER. TURBULENCIJE)      

  8   CONTINUE
      AAKT(I)=TT21(I+4*NDIM)
      OMT(I)=TT21(I+5*NDIM) 
C---------------------------------------------------
C     ODREDJIVANJE CVORA U ELEMENTU
C---------------------------------------------------
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)

C      DO 10 I=1,NDIM
       V2(1)=TT21(1)
       V2(2)=TT21(2)
       V2(3)=TT21(3)
       V2(4)=TT21(4)
       V3(1)=TT21(1+NDIM)
       V3(2)=TT21(2+NDIM)
       V3(3)=TT21(3+NDIM)
       V3(4)=TT21(4+NDIM)
      
C      write(iizlaz,*)'v2= ',v2(i), yzv(i), h(i),y(i),TT210(I+4*NDIM)
C  10  CONTINUE
C  11  continue
       
C INTERPOLACIONE FUNKCIJE

      IF (NDIM.EQ.9) THEN
      H(9)=RR*SS
      H(8)=0.5*RP*SS-0.5*H(9)
      H(7)=0.5*RR*SM-0.5*H(9)
      H(6)=0.5*RM*SS-0.5*H(9)
      H(5)=0.5*RR*SP-0.5*H(9)
      H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
      H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
      H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
      H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSEIF (NDIM.EQ.8) THEN

      H(8)=0.5*RP*SS
      H(7)=0.5*RR*SM
      H(6)=0.5*RM*SS
      H(5)=0.5*RR*SP
      H(4)=0.25*RP*SM-0.5*(H(7)+H(8))
      H(3)=0.25*RM*SM-0.5*(H(6)+H(7))
      H(2)=0.25*RM*SP-0.5*(H(5)+H(6))
      H(1)=0.25*RP*SP-0.5*(H(5)+H(8))
      
      ELSEIF (NDIM.EQ.4) THEN
      H(4)=0.25*RP*SM
      H(3)=0.25*RM*SM
      H(2)=0.25*RM*SP
      H(1)=0.25*RP*SP
      
      ELSE
      WRITE (*,*) 'BROJ COVOROVA NIJE 4,8 ILI 9'
      STOP

      ENDIF

      
C INTERPOLACIONE FUNKCIJE ZA PRITISAK, K I EPSILON

      HP(4)=0.25*RP*SM
      HP(3)=0.25*RM*SM
      HP(2)=0.25*RM*SP
      HP(1)=0.25*RP*SP
      
C IZVODI INTERPOLACIONIH FUNKCIJA PO LOKALNIM KOORDINATAMA

      IF (NDIM.EQ.9) THEN
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
      
      ELSEIF (NDIM.EQ.8) THEN
      
      ZVHR(8)=0.5*SS
      ZVHS(8)=-RP*S
      ZVHR(7)=-R*SM
      ZVHS(7)=-0.5*RR
      ZVHR(6)=-0.5*SS
      ZVHS(6)=-RM*S
      ZVHR(5)=-R*SP
      ZVHS(5)=0.5*RR
      ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))
      ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))
      ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))
      ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))
      ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))
      ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))
      ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))
      ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))

      ELSEIF (NDIM.EQ.4) THEN
      
      ZVHR(4)=0.25*SM
      ZVHS(4)=-0.25*RP
      ZVHR(3)=-0.25*SM
      ZVHS(3)=-0.25*RM
      ZVHR(2)=-0.25*SP
      ZVHS(2)=0.25*RM
      ZVHR(1)=0.25*SP
      ZVHS(1)=0.25*RP

      ENDIF
      
C     SKALARNI PROIZVOD IZVODA I X KOORDINATE      

      AJ(1,1)=DOT(ZVHR,X,NDIM)
      AJ(1,2)=DOT(ZVHR,Y,NDIM)
      AJ(2,1)=DOT(ZVHS,X,NDIM)
      AJ(2,2)=DOT(ZVHS,Y,NDIM)
      
      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
C      WRITE(*,*)'DETERMINANTA MANJA OD NULE!!!'
      WRITE(*,*)'DETERMINANTE LESS THEN ZERO!!!'
      STOP
      ENDIF

C IZVODI INTERPOLACIONIH FUNKCIJA PO X,Y KOORDINATAMA

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
      
      AKTX=0.D0
      AKTY=0.D0
      OMTX=0.D0
      OMTY=0.D0
      
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
      I3=I+2*NDIM
      I4=I+3*NDIM
           
C     IZVOD BRZINE I TEMPERATURE PO X I Y

      VXX=VXX+ZVHX(I)*TT21(I1)
      VXY=VXY+ZVHY(I)*TT21(I1)
      VYX=VYX+ZVHX(I)*TT21(I2)
      VYY=VYY+ZVHY(I)*TT21(I2)
      TX=TX+ZVHX(I)*TT21(I4)
      TY=TY+ZVHY(I)*TT21(I4)
C--------------------------------------------------------------------
C     K-OMEGA TURB. MODEL
C--------------------------------------------------------------------
      I5=I+4*NDIM
      I6=I+5*NDIM
      AKTX=AKTX+ZVHX(I)*TT21(I5)
      AKTY=AKTY+ZVHY(I)*TT21(I5)
      OMTX=OMTX+ZVHX(I)*TT21(I6)
      OMTY=OMTY+ZVHY(I)*TT21(I6)
C--------------------------------------------------------------------
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
C--------------------------------------------------------------------
C     K-OMEGA TURB. MODEL
C--------------------------------------------------------------------
      FSK1=(ELX*AKTX+ELY*AKTY)
      FSOMEGA1=(ELX*OMTX+ELY*OMTY)
C--------------------------------------------------------------------

       
      ENDIF

      IF (NPARAM.EQ.2) THEN
C       TAU=-AMI*((VXX*ETX+VYX*ETY)*ELX+(VXY*ETX+VYY*ETY)*ELY)
C       FLUX=-AKT*(ELX*TX+ELY*TY)
       TAU=-AKT*(ELX*TX+ELY*TY)
      ENDIF

      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)
C      HT=DOT(H,T,NDIM)
      ZVXT=DOT(ZVHX,T,NDIM)
      ZVYT=DOT(ZVHY,T,NDIM)
      ZVXV2=DOT(ZVHX,V2,NDIM)
      ZVYV3=DOT(ZVHY,V3,NDIM)
      ZVYV2=DOT(ZVHY,V2,NDIM)
      ZVXV3=DOT(ZVHX,V3,NDIM)
      
      ZVXK=DOT(ZVHX,AAKT,NDIM)
      ZVYK=DOT(ZVHY,AAKT,NDIM)
      
      ZVXOM=DOT(ZVHX,OMT,NDIM)
      ZVYOM=DOT(ZVHY,OMT,NDIM)
      
C PROIZVOD KT, OMT I H

      HKT=DOT(H,AAKT,NDIM)
      HOM=DOT(H,OMT,NDIM)
      
      
C-----------------------------------------------------          
CE THIS is used for calculation of shear stresses
C-----------------------------------------------------  
       DO I=1,NDIM
         V(1)=TT21(I)          
         V(2)=TT21(I+NDIM)          
         T(1)=ETX
         T(2)=ETY
         VV=DOT(T,V,2)
         VTANG(1,I)=VV*T(1)
         VTANG(2,I)=VV*T(2)
      ENDDO   
C----------------------------------------------------- 
C KOMBINACIJE KAD SU 2 CVORA = 0
C AKO SU N1, N2 = 0

      DO 40 I=1,NDIM    
      IF (ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(1,N2).EQ.0.AND.ID(2,N2)
     &.EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=V2(1)
      TT21(2)=V2(2)
      TT21(3)=VPLUS(1)
      TT21(4)=VPLUS(1)
      TT21(1+NDIM)=V3(1)
      TT21(2+NDIM)=V3(2)
      TT21(3+NDIM)=VPLUS(2)
      TT21(4+NDIM)=VPLUS(2)
      PRESS=0.D0
      
C AKO SU N3, N4 = 0
      
      ELSEIF (ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(1,N4).EQ.0.AND.
     &ID(2,N4).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=VPLUS(1)
      TT21(2)=VPLUS(1)
      TT21(3)=V2(3)
      TT21(4)=V2(4)
      TT21(1+NDIM)=VPLUS(2)
      TT21(2+NDIM)=VPLUS(2)
      TT21(3+NDIM)=V3(3)
      TT21(4+NDIM)=V3(4)
      PRESS=0.D0
      
C AKO SU N2, N3 = 0
      
      ELSEIF (ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(1,N3).EQ.0.AND.
     &ID(2,N3).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=VPLUS(1)
      TT21(2)=V2(2)
      TT21(3)=V2(3)
      TT21(4)=VPLUS(1)
      TT21(1+NDIM)=VPLUS(2)
      TT21(2+NDIM)=V3(2)
      TT21(3+NDIM)=V3(3)
      TT21(4+NDIM)=VPLUS(2)
      PRESS=0.D0
      
C AKO SU N1, N4 = 0

      ELSEIF (ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(1,N4).EQ.0.AND.
     &ID(2,N4).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=V2(1)
      TT21(2)=VPLUS(1)
      TT21(3)=VPLUS(1)
      TT21(4)=V2(4)
      TT21(1+NDIM)=V3(1)
      TT21(2+NDIM)=VPLUS(2)
      TT21(3+NDIM)=VPLUS(2)
      TT21(4+NDIM)=V3(4)
      PRESS=0.D0
      
C KOMBINACIJE KAD SU 3 CVORA = 0
C AKO SU N1, N2 I N3 = 0

      ELSEIF (ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(1,N2).EQ.0.AND.
     &ID(2,N2).EQ.0.AND.ID(1,N3).EQ.0.AND.ID(1,N3).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=V2(1)
      TT21(2)=V2(2)
      TT21(3)=V2(3)
      TT21(4)=VPLUS(1)
      TT21(1+NDIM)=V3(1)
      TT21(2+NDIM)=V3(2)
      TT21(3+NDIM)=V3(3)
      TT21(4+NDIM)=VPLUS(2)
      PRESS=0.D0
      
C AKO SU N2, N3 I N4 = 0

      ELSEIF (ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(1,N3).EQ.0.AND.
     &ID(2,N3).EQ.0.AND.ID(1,N4).EQ.0.AND.ID(1,N4).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=VPLUS(1)
      TT21(2)=V2(2)
      TT21(3)=V2(3)
      TT21(4)=V2(4)
      TT21(1+NDIM)=VPLUS(2)
      TT21(2+NDIM)=V3(2)
      TT21(3+NDIM)=V3(3)
      TT21(4+NDIM)=V3(4)
      PRESS=0.D0
      
C AKO SU N3, N4 I N1 = 0

      ELSEIF (ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(1,N4).EQ.0.AND.
     &ID(2,N4).EQ.0.AND.ID(1,N1).EQ.0.AND.ID(1,N1).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=V2(1)
      TT21(2)=VPLUS(1)
      TT21(3)=V2(3)
      TT21(4)=V2(4)
      TT21(1+NDIM)=V3(1)
      TT21(2+NDIM)=VPLUS(2)
      TT21(3+NDIM)=V3(3)
      TT21(4+NDIM)=V3(4)
      PRESS=0.D0
      
C AKO SU N3, N4 I N1 = 0

      ELSEIF (ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(1,N1).EQ.0.AND.
     &ID(2,N1).EQ.0.AND.ID(1,N2).EQ.0.AND.ID(1,N2).EQ.0) THEN
      YP(I)=H(I)*Y(I)
      YPLUS(1)=GUSM*YP(I)*VTANG(1,I)
      YPLUS(2)=GUSM*YP(I)*VTANG(2,I)
      VPLUS(1)=(YPLUS(1))/AKAPA+5.D0
      VPLUS(2)=(YPLUS(2))/AKAPA+5.D0
      TT21(1)=V2(1)
      TT21(2)=V2(2)
      TT21(3)=VPLUS(1)
      TT21(4)=V2(4)
      TT21(1+NDIM)=V3(1)
      TT21(2+NDIM)=V3(2)
      TT21(3+NDIM)=VPLUS(2)
      TT21(4+NDIM)=V3(4)
      PRESS=0.D0
      
      ENDIF
  40  CONTINUE 
      
      
      RETURN
      END
C=========================================================================