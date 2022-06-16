    PROGRAM MAIN
	
    USE GlobalDef
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER(4) :: time_begin, time_end
    INTEGER(4),PARAMETER :: NN=15,ORDER=3, ORDERP=5
    EXTERNAL :: YHC
    REAL(8) :: EP(1)
    REAL(8) :: SolC(0:NN,-NN:NN,1:3,1:3),SolS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: X0(6), T0, Tp, T1, X1(6), X2(6), X3(6), HH, DD
    INTEGER(4),PARAMETER :: Node=600
    REAL(8) :: Q0(6,Node+1),DT
    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
108 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,&
        & 1X,E25.18, 1X, E25.18) 
    
    CALL system_clock(time_begin)   ! time Keeper
    CALL GlobalPar()    ! import system parameters
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Equilibrium Points'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    EE = 1E-15    ! set the error margin for all code
    CALL NewRaf(1, (/1.6D0/), EE, EP)    ! get EP location within error margin 'Err'
    rS = EP(1)     ! set global variable
    WRITE(*,*) 'SEP         =   ', rS
    WRITE(*,*) 'The trunction error o(P4)        =   ', Ratio_Mass*rS**5D0/OrbitN**5D0
    WRITE(*,*) 'The trunction error o(P5)        =   ', Ratio_Mass*rS**6D0/OrbitN**6D0
    WRITE(*,*) 'The trunction error o(P6)        =   ', Ratio_Mass*rS**7D0/OrbitN**7D0
    WRITE(*,*) 'The trunction error o(S2)        =   ', Beta_Sun*rS**3D0/R_Sun**3D0/2D0
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Initial Values'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    CALL Sol_Coe_Put(SolC,SolS,NN,ORDERP)
    CALL AnaSol_Boost(SolC,SolS,0D0,X0,NN,ORDER)
    X0 = X0 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
    WRITE(*,*) 'X0 = '
    WRITE(*,103) X0
    T0=0D0
	TP=PI*2D0
    WRITE(*,*) 'TP          =   ', TP
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Select Target Points'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    Q0(1:6,1)=X0
	DT=(TP-T0)/3D0
    WRITE(*,*) 'Total node number = ', Node
    WRITE(*,*) 'node numer per peroid = ', (TP-T0)/DT
    WRITE(*,*) 'Total propagation time (day) = ', (Node-1)*DT,(Node-1)*DT * Tunit/(36D2*24)
    WRITE(*,*) 'delta time = ', DT
    
    
    DO I =1,Node
        T0 = DBLE(I-1)*DT
        TF = DBLE(I-0)*DT
        CALL AnaSol_Boost(SolC,SolS,TF,X1,NN,ORDER)
        X1 = X1 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
        Q0(:,I+1) = X1
    ENDDO
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Refine'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    CALL REFINE(YHC, DT,Q0,Node+1)
    
    OPEN(11,file='ORBIT_Ref.DAT', status='REPLACE')
    DO I = 1,Node+1
        write(11, 108) Q0(:,I), (I-1D0)*DT, 0D0
    ENDDO
    CLOSE(11)
    WRITE(*,*) 'Q0 = '
    WRITE(*,103) Q0(:,Node+1)
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Propagation'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    OPEN(11,file='ORBIT_Num.DAT', status='REPLACE')
    OPEN(12,file='ORBIT_Ana.DAT', status='REPLACE')
    OPEN(13,file='ORBIT_AnaRef.DAT', status='REPLACE')
    T0 = 0D0
    T1 = 0D0

    X2 = Q0(:,1)
    write(11, 108) X2, T0, 0D0
    
    DO I =1,Node
        T0 = DBLE(I-1)*DT
        T1 = DBLE(I-1)*DT
        TF = DBLE(I-0)*DT
        HH = DT/1D4
        DD = 1D-14
	    DO WHILE(T0.LT.TF)
            CALL RKF78(YHC,HH,T0,X2,EE,DD,6)
            write(11, 108) X2, T0, 0D0
            CALL AnaSol_Boost(SolC,SolS,T0,X3,NN,ORDER)
            X3 = X3 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
            write(13, 108) X3, T0, 0D0
        ENDDO
        CALL RKF78(YHC,TF-T0,T0,X2,EE,1D0,6)
        write(11, 108) X2, T0, 0D0
        CALL AnaSol_Boost(SolC,SolS,T0,X3,NN,ORDER)
        X3 = X3 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
        write(13, 108) X3, T0, 0D0
        
        X1 = Q0(:,I)
        write(12, 108) X1, T1, 0D0
        DO WHILE(T1.LT.TF)
            CALL RKF78(YHC,HH,T1,X1,EE,DD,6)
            write(12, 108) X1, T1, 0D0
        ENDDO
        CALL RKF78(YHC,TF-T1,T1,X1,EE,1D0,6)
        write(12, 108) X1, T1, 0D0
    ENDDO
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    
    WRITE(*,*) '=====================================End======================================='
    CALL system_clock(time_end)     ! time Keeper
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    END PROGRAM
    
    
    SUBROUTINE Sol_Coe_Put(SolC,SolS,NN,ORDERP)
    
    INTEGER(4) :: NN,ORDERP
    REAL(8) :: Sol1C(0:NN,-NN:NN,1:3),Sol1S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol2C(0:NN,-NN:NN,1:3),Sol2S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol3C(0:NN,-NN:NN,1:3),Sol3S(0:NN,-NN:NN,1:3)
    REAL(8) :: SolC(0:NN,-NN:NN,1:3,1:3),SolS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: OmegaE2C(0:NN,-NN:NN,1:3),OmegaE2S(0:NN,-NN:NN,1:3)
    REAL(8) :: OmegaE3C(0:NN,-NN:NN,1:3),OmegaE3S(0:NN,-NN:NN,1:3)
    
    
    CALL Sol1_Coe(Sol1C,Sol1S,NN,ORDERP)
    CALL OmegaE2(OmegaE2C,OmegaE2S,NN,ORDERP)
    CALL Sol_Coe(OmegaE2C,OmegaE2S,Sol2C,Sol2S,NN)
    CALL OmegaE3(OmegaE3C,OmegaE3S,NN,ORDERP)
    CALL Sol_Coe(OmegaE3C,OmegaE3S,Sol3C,Sol3S,NN)
    SolC(:,:,:,1) =  Sol1C(0:NN,-NN:NN,1:3)
    SolS(:,:,:,1) =  Sol1S(0:NN,-NN:NN,1:3)
    SolC(:,:,:,2) =  Sol2C(0:NN,-NN:NN,1:3)
    SolS(:,:,:,2) =  Sol2S(0:NN,-NN:NN,1:3)
    SolC(:,:,:,3) =  Sol3C(0:NN,-NN:NN,1:3)
    SolS(:,:,:,3) =  Sol3S(0:NN,-NN:NN,1:3)
    SolS(0,0,:,:) =  0D0
    
    END SUBROUTINE
    