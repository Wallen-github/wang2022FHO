PROGRAM MAIN
	
    USE GlobalDef
    
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER(4) :: time_begin, time_end
    EXTERNAL :: YHCN
    REAL(8) :: EP(1),Coe(6)
    REAL(8) :: T0, TF, TP, HH, X0(6), X1(6), X2(6)
    INTEGER(4),PARAMETER :: NN=60
    REAL(8) :: Q0(6,NN+1),DT

103 format(1X, E25.18,1X,E25.18,1X,E25.18)
108 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,&
        & 1X,E25.18, 1X, E25.18)  
	
    CALL system_clock(time_begin)   ! time Keeper
    
    CALL GlobalPar()    ! import system parameters
    
    CALL NewRaf(1, (/1.6D0/), 1D-15, EP)    ! get EP location within error margin 'Err'
    rS = EP(1)     ! set global variable
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'SEP         =   ', rS
    WRITE(*,*) 'The trunction error o(P4)        =   ', Ratio_Mass*rS**4D0/OrbitN**5D0
    WRITE(*,*) 'The trunction error o(S2)        =   ', Beta_Sun*rS**2D0/R_Sun**3D0/2D0
    
    WRITE(*,*) '=================================Initial Values================================'
    
    Coe = 0D0
    CALL YHCA(Coe, 0D0, X0)
    X0 = X0+(/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
    WRITE(*,*) 'X0 = '
    WRITE(*,103) X0

    T0=0D0
	TP=PI*2D0
    WRITE(*,*) 'TP          =   ', TP
    
    ! Propagate Numerical and Analytical Method
    
    WRITE(*,*) '=================================Loading Points================================'
    
    Q0(1:6,1)=X0
	DT=(TP-T0)/3D0
    
    DO I =1,NN
        T0 = DBLE(I-1)*DT
        TF = DBLE(I-0)*DT
        CALL YHCA(Coe, TF, X1)
        X1 = X1 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
        Q0(:,I+1) = X1
    ENDDO
    
    WRITE(*,*) '=================================Loading Refine================================'
    CALL REFINE(YHCN, DT,Q0,NN+1)
    
    OPEN(11,file='ORBIT_Ref.DAT', status='REPLACE')
    DO I = 1,NN+1
        write(11, 108) Q0(:,I), (I-1D0)*DT, 0D0
    ENDDO
    CLOSE(11)
    WRITE(*,*) 'Q0 = '
    WRITE(*,103) Q0(:,NN+1)

    
    OPEN(11,file='ORBIT_Num.DAT', status='REPLACE')
    OPEN(12,file='ORBIT_Ana.DAT', status='REPLACE')
    T0 = 0D0
    T1 = 0D0

    X2 = Q0(:,1)
    write(11, 108) X2, T0, 0D0
    
    
    DO I =1,NN
        T0 = DBLE(I-1)*DT
        T1 = DBLE(I-1)*DT
        TF = DBLE(I-0)*DT
        HH = DT/1D4
        DD = 1D-14
	    DO WHILE(T0.LT.TF)
            CALL RKF78(YHCN,HH,T0,X2,EE,DD,6)
            write(11, 108) X2, T0, 0D0
        ENDDO
        CALL RKF78(YHCN,TF-T0,T0,X2,EE,1D0,6)
        write(11, 108) X2, T0, 0D0
        
        X1 = Q0(:,I)
        write(12, 108) X1, T1, 0D0
        DO WHILE(T1.LT.TF)
            CALL RKF78(YHCN,HH,T1,X1,EE,DD,6)
            write(12, 108) X1, T1, 0D0
        ENDDO
        CALL RKF78(YHCN,TF-T1,T1,X1,EE,1D0,6)
        write(12, 108) X1, T1, 0D0
    ENDDO
    CLOSE(11)
    CLOSE(12)
    
        
        
    WRITE(*,*) '=====================================End======================================='
    CALL system_clock(time_end)     ! time Keeper
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    
    END PROGRAM MAIN