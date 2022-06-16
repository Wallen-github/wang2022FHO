    PROGRAM MAIN
    
    USE GlobalDef
    EXTERNAL :: YHC
    INTEGER(4) :: time_begin, time_end
    REAL(8) :: EE,EP(1)
    REAL(8) :: Coe(6),DPV(6),X0(6),Tp
    INTEGER(4) :: ORDER,flagPo,flagEig,Num
    REAL(8) :: Xp(6),WR(6),WI(6),DIS(6),MM(6,6)
    REAL(8) :: TF,DeltaY1
    CALL system_clock(time_begin)   ! time Keeper
    CALL GlobalPar()    ! import system parameters
    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Equilibrium Points'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    EE = 1E-15    ! set the error margin for all code
    CALL NewRaf(1, (/1.6D0/), EE, EP)    ! get EP location within error margin 'Err'
    rS = EP(1)     ! set global variable
    WRITE(*,*) 'SEP         =   ', rS
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Initial Values'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    
    !Tn = 5D1
    !OmegaN = PI2/Tn
    !OrbitN = ((Ratio_Mass+1D0) / OmegaN**2D0)**(1D0/3D0)
    
    ! 此处仅求解第三体摄动下的周期轨道，将太阳光压置为0
    ShadeFactor = 0D0
    
    Coe = 0D0
    ORDER = 1
    CALL YHCA_order(Coe, 0D0, DPV,ORDER)
    !CALL AnaSol(0D0,DPV) ! (3rd X P5)
    X0 = (/0D0,rS,0D0,0D0,0D0,0D0/) + DPV
    Tp = PI2 / ABS(OmegaN-OmegaA)
    WRITE(*,*) 'Tp (hour) = ', Tp * Tunit / 36D2,'Tn (hour) = ', Tn * Tunit / 36D2
    WRITE(*,*) 'SemiAxis = ',OrbitN
    WRITE(*,*) 'X0 = '
    WRITE(*,103) X0
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Correct State Vector'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    EE = 1D-14
    CALL PeriodOrbit_Tn(X0,Tp,EE, Xp,flagPo)
    CALL MonEig(Xp,Tp,MM,WR,WI,DIS,flagEig)
    WRITE(*,*) 'Xp = '
    WRITE(*,103) Xp
    WRITE(*,*) 'Eigen Values (REAL,IMAGINARY,VALUE) = '
    DO I = 1,6
        WRITE(*,103) WR(I),WI(I),DIS(I)
    ENDDO
    
    
    !WRITE(*,*) '==============================================================================='
    !WRITE(*,*) 'Find Boundary'
    !WRITE(*,*) '-------------------------------------------------------------------------------'
    !TF = 2D2*Tp
    !DeltaY1 = 3D-2
    !Num = 100
    !CALL DeltaY_Bound(YHC,Xp,TF,DeltaY1,Num)
    !TF = 2D2*Tp
    !CALL DeltaYup_max(YHC,Xp,TF,DeltaY1)
    
    TF = 4D0*Tp
    !Xp(2) = Xp(2)+0.02D0
    CALL WritePoData(YHC,Xp,TF)
    
    WRITE(*,*) '=====================================End======================================='
    CALL system_clock(time_end)     ! time Keeper
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    END PROGRAM
    
    SUBROUTINE WritePoData(YHC,Xp,TF)
    
    EXTERNAL :: YHC
    REAL(8) :: ERR,HH,T0,TF,Xp(6),EE
103 format(1X, E25.18,1X,E25.18,1X,E25.18)    
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    OPEN(14,file='ORBIT.DAT', status='REPLACE')
    err = 1D-14
    hh = 1D-5
    t0 = 0D0
    DO WHILE (t0<TF)
        write(14,107) Xp,t0
        Call RKF78(YHC,hh,t0,Xp,EE,err,6)
    ENDDO
    write(14,107) Xp,T0
    Call RKF78(YHC,TF-t0,t0,Xp,EE,1D0,6) ! 积分积满，恰好积分到tf
    CLOSE(14)
    WRITE(*,*) 'Xp = '
    WRITE(*,103) Xp
    
    END SUBROUTINE
    