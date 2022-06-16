    PROGRAM MAIN
    
    USE GlobalDef
    EXTERNAL :: YHC
    INTEGER(4) :: time_begin, time_end
    REAL(8) :: EE,EP(1),Fre(3)
    REAL(8) :: Coe(6),DPV(6),X0(6),Tp,Tn
    INTEGER(4) :: ORDER,flagPo,flagEig,Num
    REAL(8) :: Xp(6),WR(6),WI(6),DIS(6),MM(6,6),XBup(6),XBbo(6),DeltaY
    REAL(8) :: TNend,TNstart,TNstep
    CALL system_clock(time_begin)
    CALL GlobalPar()
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
113 format(1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,I5,1X,E25.18,1X,E25.18,1X,E25.18)  
115 format(1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,I5,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Equilibrium Points'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    EE = 1E-15
    CALL NewRaf(1, (/1.6D0/), EE, EP)
    rS = EP(1)
    WRITE(*,*) 'SEP         =   ', rS
    CALL GetFre(Fre)
    WRITE(*,*) 'Fre'
    WRITE(*,103) Fre
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Initial Values'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    
    Num = 50
    TNstart = 5D1
    TNend = 1D2
    TNstep = (TNend-TNstart)/Num
    Tn = TNstart
    ! 此处仅求解第三体摄动下的周期轨道，将太阳光压置为0
    ShadeFactor = 0D0

    OmegaN = PI2/Tn
    OrbitN = ((Ratio_Mass+1D0) / OmegaN**2D0)**(1D0/3D0)
    Tp = PI2 / ABS(OmegaN-OmegaA)
    
    Coe = 0D0
    ORDER = 1
    CALL YHCA_order(Coe, 0D0, DPV,ORDER)
    X0 = (/0D0,rS,0D0,0D0,0D0,0D0/) + DPV
    WRITE(*,*) 'TN (hour) = ', Tn
    WRITE(*,*) 'Tp (hour) = ', Tp
    WRITE(*,*) 'X0 = '
    WRITE(*,103) X0
    
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Solve Period Orbit Family'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    OPEN(14,file='POF.DAT', status='REPLACE')
    
    DO I = 1,Num
        CALL PeriodOrbit_Tn(X0,Tp,EE, Xp,flagPo)
        CALL MonEig(Xp,Tp,MM,WR,WI,DIS,flagEig)
        IF (flagPo==1) THEN
            WRITE(*,"(A,I4,A,I4)") 'I/Num=',I,'/',Num
            WRITE(*,"(A,E25.18,1X,A,E25.18,1X,A,I2)") 'Tn = ', Tn, 'SemiAxis(m) = ', OrbitN*Runit,'flagEig = ', flagEig
            !WRITE(14,113) Xp, Tp, Tn*36D2/Tunit,OrbitN, flagEig, Fre
            !CALL PO_boundary(YHC,200,1D-1,10*Tp,Xp,XBup,XBbo)
            !WRITE(*,*) Xp(2),XBup(2),XBbo(2)
            CALL DeltaYup_max(YHC,Xp,2D2*Tp,DeltaY)
            WRITE(14,115) Xp, Tp, Tn,OrbitN, flagEig, Fre,Xp(2)+DeltaY,Axis_b/Runit!XBup(2),XBbo(2)
            X0 = Xp
        ELSEIF (flagPo==0) THEN
            WRITE(*,*) 'Error in Main_POT_Tn: Failed'
        ENDIF
        
        Tn = TNstart + I*TNstep
        OmegaN = PI2/Tn
        Tp = PI2 / ABS(OmegaN-OmegaA)
        OrbitN = ((Ratio_Mass+1D0) / OmegaN**2D0)**(1D0/3D0)
    ENDDO
    CLOSE(14)

    
    WRITE(*,*) '=====================================End======================================='
    CALL system_clock(time_end)     ! time Keeper
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    END PROGRAM
    
    
    SUBROUTINE PO_boundary(YHC,Num,deltaY,TF,Xp,XBup,XBbo)
    EXTERNAL :: YHC
    REAL(8) :: ERR,HH,T0,TF,EE
    REAL(8) :: Xp(6),XBup(6),XBbo(6),deltaY
    INTEGER(4) :: FLAG,Num
    
    err = 1D-14
    hh = 1D-5
    step = deltaY/Num
    DO I = 1,Num
        FLAG = 1
        XBup = Xp
        T0 = 0D0
        XBup(2) = Xp(2)+I*step
        DO WHILE (T0<TF)
            Call RKF78(YHC,hh,T0,XBup,EE,err,6)
            CALL CheckImpact(XBup,FLAG)
            IF (XBup(2)<0D0 .OR. FLAG == 0) THEN
                XBup = Xp
                XBup = Xp(2)+(I-1)*step
                GOTO 201
            ENDIF
        ENDDO
        Call RKF78(YHC,TF-T0,T0,XBup,EE,1D0,6)
        CALL CheckImpact(XBup,FLAG)
        IF (XBup(2)<0D0 .OR. FLAG == 0) THEN
            XBup = Xp
            XBup = Xp(2)+(I-1)*step
            GOTO 201
        ENDIF
    ENDDO
    
201 DO I = 1,Num
        FLAG = 1
        XBbo = Xp
        T0 = 0D0
        XBbo(2) = Xp(2)-I*step
        DO WHILE (T0<TF)
            Call RKF78(YHC,hh,T0,XBbo,EE,err,6)
            CALL CheckImpact(XBbo,FLAG)
            IF (XBbo(2)<0D0 .OR. FLAG == 0) THEN
                XBbo = Xp
                XBbo = Xp(2)-(I-1)*step
                GOTO 202
            ENDIF
        ENDDO
        Call RKF78(YHC,TF-T0,T0,XBbo,EE,1D0,6)
        CALL CheckImpact(XBbo,FLAG)
        IF (XBbo(2)<0D0 .OR. FLAG == 0) THEN
            XBbo = Xp
            XBbo = Xp(2)-(I-1)*step
            GOTO 202
        ENDIF
    ENDDO
    
202 END SUBROUTINE
    