    PROGRAM MAIN
    ! ========================================================================================
    ! 比较数值积分，分析解，半分析解精度，搭配的绘图程序为Plot_orbit_num_Ana.m
    ! 参数说明：
    ! ORDER:    分析解阶数，目前可选择1,2,3阶
    ! ORDERP:   第三体摄动展开阶数，目前可选择保留P3,P4,P5项
    ! ---------------------------------------------------------------------------------------
    use CMD_Progress
    USE GlobalDef
    EXTERNAL :: YHCN,YHC
    INTEGER(4),PARAMETER :: NN=15,ORDER=3, ORDERP=5
    REAL(8) :: EE,Ep(1)
    REAL(8) :: Coe(6),PVN(6),DPV(6),DPV2(6)
    REAL(8) :: T0,TF,step,DD
    INTEGER(4) :: index
    INTEGER(4) :: time_begin, time_end
    REAL(8) :: Sol1C(0:NN,-NN:NN,1:3),Sol1S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol2C(0:NN,-NN:NN,1:3),Sol2S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol3C(0:NN,-NN:NN,1:3),Sol3S(0:NN,-NN:NN,1:3)
    REAL(8) :: SolC(0:NN,-NN:NN,1:3,1:3),SolS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: OmegaE2C(0:NN,-NN:NN,1:3),OmegaE2S(0:NN,-NN:NN,1:3)
    REAL(8) :: OmegaE3C(0:NN,-NN:NN,1:3),OmegaE3S(0:NN,-NN:NN,1:3)
    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
108 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    CALL GlobalPar()    ! import system parameters
    CALL system_clock(time_begin)   ! time Keeper
    
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
    
    Coe = 0D0
    !CALL YHCA(Coe, 0D0, PVN)
    !CALL YHCA_order(Coe, 0D0, PVN,ORDER)
    CALL AnaSol_Boost(SolC,SolS,0D0,PVN,NN,ORDER)
    PVN = PVN + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/)
    WRITE(*,*) 'PV0 = '
    WRITE(*,103) PVN
    
    OPEN(11,file='ORBIT_Numerical.DAT', status='REPLACE')
    OPEN(12,file='ORBIT_Analytical.DAT', status='REPLACE')
    OPEN(13,file='ORBIT_AnaNum.DAT', status='REPLACE')
    write(11, 108) (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), 0D0, 0D0
    write(12, 108) (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), 0D0, 0D0
    write(13, 108) (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), 0D0, 0D0
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Propagation'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    T0 = 0D0
    TF = 1D2!PI2/3D1
    step = 1D-2
    DD = 1E-15
    EE = 1E-15
    index = 0
    WRITE(*,'(A\)') 'Begin'
    DO WHILE(T0 < TF)
        ! Progress bar
        IF (INT(T0/TF * 7D1)==index+1) THEN
        write(*,'(A\)') '/'
        ENDIF
        index = INT(T0/TF * 7D1)
        
        ! Numerical Solution
        CALL RKF78(YHC,step,T0,PVN,EE,DD,6) ! 定步长与变步长，单摆积分不对
        write(11, 108) PVN, T0, 0D0
        
        ! Analytical Solution
        !CALL YHCA(Coe, T0, DPV)
        CALL YHCA_order(Coe, T0, DPV,ORDER)
        write(12, 108) DPV + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), t0, 0D0
        
        ! Analytical Solution by numerical coeff
        !CALL AnaSol(T0,DPV2)
        CALL AnaSol_Boost(SolC,SolS,T0,DPV2,NN,ORDER)
        write(13, 108) DPV2 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), t0, 0D0
    ENDDO
    CALL RKF78(YHC,TF-T0,T0,PVN,EE,1D0,6)
    write(11, 108) PVN, T0, 0D0
    CALL YHCA_order(Coe, T0, DPV,ORDER)
    !CALL YHCA(Coe, T0, DPV)
    write(12, 108) DPV + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), t0, 0D0
    CALL AnaSol_Boost(SolC,SolS,T0,DPV2,NN,ORDER)
    write(13, 108) DPV2 + (/0D0, rS, 0D0, 0D0, 0D0, 0D0/), t0, 0D0
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    WRITE(*,*) 'END'
    
    WRITE(*,*) '=====================================End======================================='
    CALL system_clock(time_end)     ! time Keeper
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    END PROGRAM