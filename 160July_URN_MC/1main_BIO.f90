!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!     Hogil Kim Memorial Building #501 POSTECH,
!     San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!     E-mail: hoangdanhtai@gmail.com
!-----------------------------------------------------------------------------------------!
!!    05.03.2013: calcul to compare with Master Equation    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
      PROGRAM main_BIO
      IMPLICIT NONE

      CHARACTER (256)  :: Ligne21

      INTEGER (KIND=8) :: t,n_time,i_average,n_average,nA0,nB0
      
      REAL    (KIND=8) :: rdn_ab1,rdn_ab2,Pa,Pab,Pb,Pba,rdn_Pa,rdn_Pab,rdn_Pb,rdn_Pba
      REAL    (KIND=8) :: tab_tmp
      REAL    (KIND=8) :: nA,nA1,nB,nB1

      INTEGER (KIND=8),DIMENSION(3)           :: clock
     
      REAL (KIND=8),DIMENSION(:), ALLOCATABLE :: tab_nA,tab_nB,nA_av,nB_av,nAB_av
      REAL (KIND=8),DIMENSION(:), ALLOCATABLE :: nA2_av,nB2_av,F_nA,F_nB
      REAL (KIND=8),DIMENSION(:), ALLOCATABLE :: nAB21_av,nAB22_av,F_nAB21,F_nAB22

!!!=======================================================================================
!!!=======================================================================================
      CALL system('rm *.dat*')

      CALL ini_rdm_number()
      CALL read_input_parameter_file()

      ALLOCATE(tab_nA(0:n_time))
      ALLOCATE(tab_nB(0:n_time))
      ALLOCATE(nA_av(0:n_time))
      ALLOCATE(nB_av(0:n_time))
      ALLOCATE(nAB_av(0:n_time))

      ALLOCATE(nA2_av(0:n_time))
      ALLOCATE(nB2_av(0:n_time))
      ALLOCATE(F_nA(0:n_time))
      ALLOCATE(F_nB(0:n_time))
      
      ALLOCATE(nAB21_av(0:n_time))
      ALLOCATE(nAB22_av(0:n_time))
      ALLOCATE(F_nAB21(0:n_time))
      ALLOCATE(F_nAB22(0:n_time))

      OPEN(unit=21,file='average_thermal.dat')  

      !CALL value_thermal()

      !WRITE(*,*)Pa,Pab,Pb,Pba

      CALL average_thermal()

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
      READ(11, '(A30,(I5))')    tamp, nA0
      READ(11, '(A30,(I5))')    tamp, nB0
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, Pa
      READ(11, '(A30,(F7.4))')  tamp, Pb
      READ(11, '(A30,(F7.4))')  tamp, Pab
      READ(11, '(A30,(F7.4))')  tamp, Pba
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I10))')   tamp, n_time
      READ(11, '(A30,(I10))')   tamp, n_average
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

      nA=real(nA0) ; nB=real(nB0)

      DO t=1,n_time

      nA1=nA; nB1=nB

      !!!------------------------------------
      CALL random_number(rdn_ab1)
      
      IF (rdn_ab1<=nA/(nA+nB)) THEN
      
            CALL random_number(rdn_Pa)
            IF (rdn_Pa<=Pa) THEN
                  nA1=nA1+1.
            END IF
      ELSE
      
            CALL random_number(rdn_Pb)
            IF (rdn_Pb<=Pb) THEN
                  nB1=nB1+1.
            END IF
      END IF
      !!!------------------------------------      
      
      !!!------------------------------------
      CALL random_number(rdn_ab2)
      
      IF (rdn_ab2<=nA/(nA+nB)) THEN
            
            CALL random_number(rdn_Pab)
            IF (rdn_Pab<=Pab) THEN
                  nA1=nA1-1.
                  nB1=nB1+1.
            END IF

      ELSE
    
            CALL random_number(rdn_Pba)
            IF (rdn_Pba<=Pba) THEN
                  nB1=nB1-1.
                  nA1=nA1+1.
            END IF
            
      END IF
    
      nA=nA1 ; nB=nB1

      tab_nA(t)=nA ; tab_nB(t)=nB


      END DO

      END SUBROUTINE value_thermal    

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal()
      IMPLICIT NONE

      tab_nA(0)=nA0 ; tab_nB(0)=nB0
      nA_av(:)=0. ; nB_av(:)=0. ; nA2_av(:)=0. ; nB2_av(:)=0.
      nAB21_av(:)=0. ; nAB22_av(:)=0.

      !!!---------------------------------------
      DO i_average=1,n_average
            CALL value_thermal()

            DO t=0,n_time
                  nA_av(t)=nA_av(t)+tab_nA(t)
                  nB_av(t)=nB_av(t)+tab_nB(t)
                  
                  nA2_av(t)=nA2_av(t)+tab_nA(t)**2.
                  nB2_av(t)=nB2_av(t)+tab_nB(t)**2.
                  
                  nAB21_av(t)=nAB21_av(t)+(tab_nA(t)+tab_nB(t))**2.
                  nAB22_av(t)=nAB22_av(t)+(tab_nA(t)-tab_nB(t))**2.
                                    
            END DO
      END DO

      !!!---------------------------------------
      DO t=0,n_time
            nA_av(t)=nA_av(t)/real(n_average)
            nB_av(t)=nB_av(t)/real(n_average)
            

            nA2_av(t)=nA2_av(t)/real(n_average)
            nB2_av(t)=nB2_av(t)/real(n_average)
                  
            nAB21_av(t)=nAB21_av(t)/real(n_average)
            nAB22_av(t)=nAB22_av(t)/real(n_average)

            nAB_av(t)=nA_av(t)+nB_av(t)

            F_nA(t)=nA2_av(t)-nA_av(t)**2.
            F_nB(t)=nB2_av(t)-nB_av(t)**2.
            
            F_nAB21(t)=nAB21_av(t)-nAB_av(t)**2.
            F_nAB22(t)=nAB22_av(t)-(nA_av(t)-nB_av(t))**2.
                      
            WRITE(Ligne21,*) t,nA_av(t)/nAB_av(t),nB_av(t)/nAB_av(t),F_nA(t)/nAB_av(t),&
            F_nB(t)/nAB_av(t),nA_av(t)+nB_av(t),nA_av(t)-nB_av(t),F_nAB21(t),F_nAB22(t)
            WRITE(21,'(a)') trim(Ligne21)


      END DO     

      END SUBROUTINE average_thermal


      END PROGRAM main_BIO

