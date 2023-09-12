! Bogdan
! mart 2016        
         MODULE SERVIS
         IMPLICIT NONE
#define PHASE_MAIN_FULL 1
#define PHASE_SUBSTRUCT_FULL 2
#define PHASE_SUBSTRUCT_RHS 3
#define PHASE_MAIN_RHS 4  
#define PHASE_MAIN_VALS 5       	                
! trefs_dtmods se sastoji od 2 double precision promenljive          
         TYPE TREFS_DTMODS
         DOUBLE PRECISION TREF,DTMOD
         END TYPE       
! niz trefs_and_dtmods je pomocna struktura za reorganizovanu petlju
! cuvaju se promenljive tref i dtmod za svaku podstrukturu         
         TYPE(TREFS_DTMODS),ALLOCATABLE:: TREFS_AND_DTMODS(:)     
! niz cell_convergence_ind je niz indikatora 
! 0- sistem nije konvergirao 1-jeste               
         INTEGER, ALLOCATABLE :: CELL_CONVERGENCE_IND (:)  
! BROJ_ZAHTEVA broj sistema koji se salju servisu na resavanje     
         INTEGER BROJ_ZAHTEVA 
! ITERATION_SUB iteracija do koje su stigle podstukture koje nisu konvergirale         
         INTEGER ITERATION_SUB
! PHASE faza oznacava koji se podaci salju servisu 
! 1- salje se ceo glavni sistem, 2 - salju se celi sistemi podstruktura 
! 3- salju se samo desne strane podstruktura 4- salje se samo desna strana glavnog sistema          
         INTEGER PHASE  
         integer servis_save_lmax   
         real :: startG,endG   
         integer count_sis  
         CONTAINS
!===============================================================
! SOLVER_REMOTE je procedura za slanje glavnog sistema servisu, i dobavljanje resenja istog
! potpis procedure je isti kao SOLVER_INTERNAL uz dodatak argumenta PHASE 
! ova procedura poziva funkcije iz projekta communication : 
! send_conf_data sluzi sa slanje konfiguracionih podataka, 
! servisu se salje phase i broj sistema koji treba da resi 
! send_request- salje se ceo sistem servisu, a send_rhs samo desna strana
! recv_solution je funkcija za dobavljanje resenja od servisa
         SUBROUTINE SOLVER_REMOTE(SILE, JEDN,INFF)      
         IMPLICIT NONE
         DOUBLE PRECISION SILE (*)
         INTEGER JEDN, ISPARSE_N, ISPARSE_NZ,INFF,TREE_INDEX, IFILE, I
         DOUBLE PRECISION, ALLOCATABLE :: VALS (:)
         INTEGER, ALLOCATABLE :: ROWS (:)
         INTEGER, ALLOCATABLE :: COLS (:)

         ISPARSE_N = JEDN
         CALL SPARSEASSEMBLER_GETNONZERO(ISPARSE_NZ)
         ALLOCATE(COLS(ISPARSE_NZ))
         ALLOCATE(ROWS(ISPARSE_NZ))
         ALLOCATE(VALS(ISPARSE_NZ))
         
         CALL sparseassembler_savesparsefile()
         !CALL SPARSEASSEMBLER_GETSPARSE(ISPARSE_NZ, ROWS,COLS,VALS,INFF)
		 CALL SPARSEASSEMBLER_KILL(INFF)      
		 ! read sparse file 
            IFILE = 987
            OPEN (IFILE, FILE='sparse.bin', FORM='BINARY')
            DO I = 1, ISPARSE_NZ
                READ(IFILE) ROWS(I),COLS(I),VALS(I)
            END DO
            CLOSE (IFILE)	
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	    
         IF (INFF.EQ.1) THEN 
            PRINT *,' *MAIN* '
            PRINT *,'Sending data... '
         ENDIF
		 IF(ISPARSE_N.GT.0.AND.ISPARSE_NZ.GT.0)THEN
           IF (PHASE .EQ. PHASE_MAIN_FULL) THEN 
              CALL SEND_REQUEST(ISPARSE_N,ISPARSE_NZ,ROWS,COLS,VALS,SILE)
           ENDIF
           IF (PHASE.EQ. PHASE_MAIN_RHS) THEN
              CALL SEND_DOUBLE_ARRAY(ISPARSE_N,SILE)
           ENDIF
           IF (PHASE.EQ. PHASE_MAIN_VALS) THEN
              CALL SEND_DOUBLE_ARRAY(ISPARSE_NZ,VALS)
              CALL SEND_DOUBLE_ARRAY(ISPARSE_N,SILE)            
           ENDIF
           IF (INFF.EQ.1) THEN 
               PRINT *,'Data sent.Waiting for solution...'
           ENDIF  	   
           CALL RECV_SOLUTION(SILE,ISPARSE_N,INFF)
           IF (INFF.EQ.1) THEN 
               PRINT *,'Got solution!'
           ENDIF 	   
		 ELSE
			WRITE(*,*) 'ERROR IN SOLVER CALLING!'
			STOP
		 ENDIF
         !call print_system(ISPARSE_N,ISPARSE_NZ,rows,cols,vals,SILE)
		 DEALLOCATE(COLS)
		 DEALLOCATE(ROWS)
		 DEALLOCATE(VALS)

		 END SUBROUTINE
  !=============================================================  
  ! procedura za cuvanje sistema, na poziciji INDEX u nizu systems
  ! sile - rhs, jedn -broj jednacina, ostali podaci o sistemu se ucitavaju koriscenjem sparseassembler        
         SUBROUTINE SEND_SYSTEM(INDEX,SILE, JEDN)
         IMPLICIT NONE
         DOUBLE PRECISION SILE (*)
         INTEGER JEDN,I,INDEX,NFILE,TREE_INDEX, KOLIKO,INDEX_TO_SEND  
         DOUBLE PRECISION, ALLOCATABLE :: VALS (:)
         INTEGER, ALLOCATABLE :: IROWS (:)
         INTEGER, ALLOCATABLE :: JCOLS (:)         
		 CALL SPARSEASSEMBLER_GETNONZERO(KOLIKO)        
         IF (PHASE .EQ. PHASE_SUBSTRUCT_FULL) THEN           
		     ALLOCATE(IROWS(KOLIKO))
		     ALLOCATE(JCOLS(KOLIKO))
		     ALLOCATE(VALS(KOLIKO)) 
			 CALL SPARSEASSEMBLER_GETSPARSE(KOLIKO, &
	         & IROWS,JCOLS, VALS,0)		      
		 ENDIF	 
		 CALL SPARSEASSEMBLER_KILL(0) 
         IF (PHASE .EQ. PHASE_SUBSTRUCT_FULL) THEN
                CALL SEND_REQUEST(JEDN,KOLIKO, &
               & IROWS,JCOLS, VALS, SILE)
         ENDIF
        IF (PHASE .EQ. PHASE_SUBSTRUCT_RHS) THEN
                INDEX_TO_SEND = INDEX - 1
                CALL SEND_INT_ARRAY(1,INDEX_TO_SEND)
				CALL SEND_DOUBLE_ARRAY(JEDN,SILE)
         ENDIF                
         END SUBROUTINE            
!============================================================
 ! svi elementi niza convergence se postavljaju na nulu
        SUBROUTINE SET_NONE_CONVERGED(N)
        IMPLICIT NONE 
        INTEGER N,I
        DO I=1,N 
           CELL_CONVERGENCE_IND(I)=0     
        END DO
        END SUBROUTINE                               
!========================================================== 
!=============================================================   
! alokacija niza trefs_and_dtmods i postavljanje svih elemenata na nulu       
        SUBROUTINE ZERO_TREFS_AND_DTMODS(N)
        IMPLICIT NONE
        INTEGER N,I
        DO I=1,N
           TREFS_AND_DTMODS(I)%TREF=0.D0
           TREFS_AND_DTMODS(I)%DTMOD=0.D0       
        END DO
        END SUBROUTINE         
!============================================================= 
! procedura za cuvanje tref i dtmod         
        SUBROUTINE SAVE_TREFS_AND_DTMODS(INDEX,TREF,DTMOD)
        IMPLICIT NONE
        INTEGER INDEX
        DOUBLE PRECISION TREF,DTMOD
        TREFS_AND_DTMODS(INDEX)%TREF=TREF
        TREFS_AND_DTMODS(INDEX)%DTMOD=DTMOD
        END SUBROUTINE
 !============================================================ 
 ! procedura za ucitavanje vrednosti tref i dtmod 
        SUBROUTINE LOAD_TREFS_AND_DTMODS(INDEX,TREF,DTMOD)
        IMPLICIT NONE
        INTEGER INDEX
        DOUBLE PRECISION TREF,DTMOD
        TREF=TREFS_AND_DTMODS(INDEX)%TREF
        DTMOD=TREFS_AND_DTMODS(INDEX)%DTMOD
        END SUBROUTINE
!========================================================== 
! procedura za slanje signala servisu da je posao zavrsen       
		SUBROUTINE SERVIS_DONE()
		CALL SEND_INT_ARRAY(1,-1)
		END SUBROUTINE
!==========================================================		
! procedura za konektovanje na servis      
		SUBROUTINE SERVIS_START()
		integer port
		character(100) :: hostname
		open(unit = 100, file = 'CONF.txt',action = 'read')
		read (100,*), hostname
		read(100,*), port
		close(100)
		CALL CONNECT_TO_SERVIS(trim(hostname)//CHAR(0),port)
		END SUBROUTINE   
!=========================================================
        SUBROUTINE MEASURE_TIME()
        real :: elapsedtime,etime,exetime(2)
        elapsedtime=etime(exetime)
        print *,' TIME:'
        print *, 'elapsed:', elapsedtime, ' user:', exetime(1),' sys:', exetime(2)
#if defined (_WIN32) || defined (_WIN64)       
        print *,' Press any key to continue ...'
        call getchar()
#endif    
        END SUBROUTINE  
!=========================================================       
        SUBROUTINE MEASURE_TIME_START()
        real :: etime,exetime(2)
        startG=etime(exetime)   
        END SUBROUTINE  
!========================================================= 
        SUBROUTINE MEASURE_TIME_END()
        real :: elapsedtime,etime,exetime(2)
        endG=etime(exetime)
        print *,' TIME:'
        print *, endG-startG   
        END SUBROUTINE  
!=========================================================
        subroutine print_system(m,koliko,rows,cols,vals,sile)
        integer :: m, koliko,i
        DOUBLE PRECISION :: VALS (:)
        DOUBLE PRECISION SILE (*)
        INTEGER :: ROWS (:)
        INTEGER :: COLS (:)
        if(count_sis .eq. 58) return
        if(count_sis.eq.0) then 
        count_sis=49
        else 
        count_sis=count_sis+1
        endif
        open(unit=1001,file='SISTEM'//char(count_sis)//'.txt',action='write')
        write (1001,*) '%%MatrixMarket matrix coordinate real general '
        write(1001,*) m,m,koliko
        do i=1,koliko
         write(1001,*) rows(i),cols(i),vals(i)
        end do
        do i=1,m
        write (1001,*) sile(i)
        end do 
        close(1001)
        end subroutine  
!=========================================================
        END MODULE

        
