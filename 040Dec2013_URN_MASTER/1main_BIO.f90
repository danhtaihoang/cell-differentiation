!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!     Hogil Kim Memorial Building #501 POSTECH,
!     San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!     E-mail: hoangdanhtai@gmail.com    Personal site: hoangdanhtai.com 
!-----------------------------------------------------------------------------------------!
!     Replication and Transdifferentation by using Master Equation 
!!    21.02.2013
!!    03.03.2013: histogram PA, PB
!!    05.03.2013: Write <nA>/<nA+nB>; <FnA>/<nA+nB> ; .....  
!!    04.12.2013: Tinh Entropy S(t) va delS(t)     
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
      PROGRAM main_BIO
      IMPLICIT NONE

      CHARACTER (256)  :: Ligne21,Ligne22,Ligne23,Ligne31,Ligne32,Ligne33,Ligne41

      INTEGER (KIND=8) :: n_time,nA0,nB0,nAmax,nBmax,i,j,t,t0,iPAmin,iPAmax,iPBmin,iPBmax
      
      REAL    (KIND=8) :: Pa,Pab,Pb,Pba
      REAL    (KIND=8) :: tab_tmp

      INTEGER (KIND=8),DIMENSION(3) :: clock
     
      REAL (KIND=8),DIMENSION(:,:,:), ALLOCATABLE :: P
      
      REAL (KIND=8),DIMENSION(:), ALLOCATABLE :: nA_av,nB_av,nA2_av,nB2_av,F_nA,F_nB,S,delS
      REAL (KIND=8),DIMENSION(:), ALLOCATABLE :: nAB21_av,nAB22_av,F_nAB21,F_nAB22

      REAL (KIND=8),DIMENSION(:,:), ALLOCATABLE :: PA_av,PB_av
!!!=======================================================================================
!!!=======================================================================================
      CALL system('rm *.dat*')

      CALL ini_rdm_number()
      CALL read_input_parameter_file()
      
      
      nAmax=nA0+n_time ; nBmax=nB0+n_time

      !WRITE(*,*)Pa,Pab,Pb,Pba

      ALLOCATE(P(-1:nAmax+1,-1:nBmax+1,0:n_time+1))
      ALLOCATE(nA_av(0:n_time+1))
      ALLOCATE(nB_av(0:n_time+1))
      ALLOCATE(nA2_av(0:n_time+1))
      ALLOCATE(nB2_av(0:n_time+1))
      ALLOCATE(F_nA(0:n_time+1))
      ALLOCATE(F_nB(0:n_time+1))
      ALLOCATE(S(0:n_time+1))
      ALLOCATE(delS(0:n_time+1))
      
      ALLOCATE(nAB21_av(0:n_time+1))
      ALLOCATE(nAB22_av(0:n_time+1))
      ALLOCATE(F_nAB21(0:n_time+1))
      ALLOCATE(F_nAB22(0:n_time+1))
                     
      ALLOCATE(PA_av(-1:nAmax+1,0:n_time+1))
      ALLOCATE(PB_av(-1:nBmax+1,0:n_time+1))                
                     

      OPEN(unit=20,file='P.dat') 
      OPEN(unit=21,file='average_thermal.dat')
      OPEN(unit=22,file='PA_av.dat')
      OPEN(unit=23,file='PB_av.dat')
      
      OPEN(unit=31,file='Pt1.dat')
      OPEN(unit=32,file='Pt2.dat')
      OPEN(unit=33,file='Pt3.dat')
      
      OPEN(unit=41,file='entropy.dat')

      CALL value_thermal()
    

      CONTAINS

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ini_rdm_number()
      IMPLICIT NONE

      INTEGER (KIND=8) :: i_time,i

      CALL ITIME(clock)
      i_time=(clock(1)+1)*(clock(2)+1)*(clock(3)+1)
        
      DO i=1,i_time
            CALL random_number(tab_tmp)
      ENDDO 

      END SUBROUTINE ini_rdm_number

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter_file() 
!!! OPEN the parameter from file "1parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_input_parameter_file()
      IMPLICIT NONE

      CHARACTER (LEN=150) :: tamp
      OPEN(11,file='1parameter.in')
      
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I10))')   tamp, nA0
      READ(11, '(A30,(I10))')   tamp, nB0
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, Pa
      READ(11, '(A30,(F7.4))')  tamp, Pb
      READ(11, '(A30,(F7.4))')  tamp, Pab
      READ(11, '(A30,(F7.4))')  tamp, Pba
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I10))')   tamp, n_time
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I10))')   tamp, t0
      READ(11, '(A30,(I10))')   tamp, iPAmin
      READ(11, '(A30,(I10))')   tamp, iPAmax
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I10))')   tamp, iPBmin
      READ(11, '(A30,(I10))')   tamp, iPBmax
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      END SUBROUTINE read_input_parameter_file

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()
      IMPLICIT NONE

      P(:,:,:)=0. ; P(nA0,nB0,0)=1. 

      !!!------------------------------------
      !!! Calculate the probabilities
      
      DO t=0,n_time
            DO j=0,nBmax
            DO i=0,nAmax
      
            IF ((i+j)>1) THEN
                  
            P(i,j,t+1)=P(i,j,t)+Pa*(i-1)/(i+j-1)*P(i-1,j,t)+Pb*(j-1)/(i+j-1)*P(i,j-1,t)&
                      +Pba*(j+1)/(i+j)*P(i-1,j+1,t)+Pab*(i+1)/(i+j)*P(i+1,j-1,t)&
                      -((Pa+Pab)*i/(i+j)+(Pb+Pba)*j/(i+j))*P(i,j,t)
                      !  -(Pa+Pab)*P(i,j,t)
            END IF
            
            !!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            !!! 06.03.2013
            
            IF (t==50) THEN

                  WRITE(Ligne31,*)i,j,P(i,j,t)   
                  WRITE(31,'(a)') trim(Ligne31)
            
            END IF

            IF (t==100) THEN

                  WRITE(Ligne32,*)i,j,P(i,j,t)   
                  WRITE(32,'(a)') trim(Ligne32)
            
            END IF
            
            IF (t==200) THEN

                  WRITE(Ligne33,*)i,j,P(i,j,t)   
                  WRITE(33,'(a)') trim(Ligne33)
            
            END IF

            !!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            END DO
            END DO
      END DO

      !!!------------------------------------
      !!! Calculate the average values
      
      nA_av(:)=0. ; nB_av(:)=0. ; nA2_av(:)=0. ; nB2_av(:)=0. 
      nAB21_av(:)=0. ; nAB22_av(:)=0.
      
      DO t=0,n_time
      
            DO j=0,nBmax
            DO i=0,nAmax
      
            nA_av(t)=nA_av(t)+i*P(i,j,t)
            nB_av(t)=nB_av(t)+j*P(i,j,t)
            
            nA2_av(t)=nA2_av(t)+(i**2.)*P(i,j,t)
            nB2_av(t)=nB2_av(t)+(j**2.)*P(i,j,t)
            
            nAB21_av(t)=nAB21_av(t)+((i+j)**2.)*P(i,j,t) 
            nAB22_av(t)=nAB22_av(t)+((i-j)**2.)*P(i,j,t) 
             
            END DO
            END DO
            
            F_nA(t)=nA2_av(t)-nA_av(t)**2.
            F_nB(t)=nB2_av(t)-nB_av(t)**2.
            
            F_nAB21(t)=nAB21_av(t)-(nA_av(t)+nB_av(t))**2.
            F_nAB22(t)=nAB22_av(t)-(nA_av(t)-nB_av(t))**2.
            
            WRITE(Ligne21,*) t,nA_av(t)/(nA_av(t)+nB_av(t)),nB_av(t)/(nA_av(t)+nB_av(t)),&
            F_nA(t)/(nA_av(t)+nB_av(t)),F_nB(t)/(nA_av(t)+nB_av(t)),F_nAB21(t),&
            F_nAB22(t),nA_av(t)+nB_av(t),nA_av(t)-nB_av(t)
            WRITE(21,'(a)') trim(Ligne21)
            
      END DO

      !!!------------------------------------ 
      !!! Histogram
     
      PA_av(:,:)=0.

      DO i=iPAmin,iPAmax
      
            DO j=0,nBmax
                  PA_av(i,t0)=PA_av(i,t0)+P(i,j,t0)
            END DO
         
            WRITE(Ligne22,*) i,PA_av(i,t0)
            WRITE(22,'(a)') trim(Ligne22)
      END DO
      
      !!!-------------------------------
      PB_av(:,:)=0.

      DO j=iPBmin,iPBmax
      
            DO i=0,nAmax
                  PB_av(j,t0)=PB_av(j,t0)+P(i,j,t0)
            END DO
         
            WRITE(Ligne23,*) j,PB_av(j,t0)
            WRITE(23,'(a)') trim(Ligne23)
      END DO

      !!!=================================================================================
      !!! 04.12.2013: Entropy S(t)
      S(:)=0.
      DO t=0,n_time
      
            DO j=0,nBmax
            DO i=0,nAmax
      
            IF (P(i,j,t)/=0.) THEN
                  S(t)=S(t)-P(i,j,t)*log(P(i,j,t))
            END IF
             
            END DO
            END DO
         
            S(t)=S(t)/log(2.)
      END DO

      !!! delS(t)
      delS(:)=0.      
      DO t=0,n_time
            delS(t)=S(t+1)-S(t)
            
            WRITE(Ligne41,*)t,S(t),delS(t)
            WRITE(41,'(a)') trim(Ligne41)
            
      END DO

      !!!=================================================================================
      END SUBROUTINE value_thermal    


      END PROGRAM main_BIO

