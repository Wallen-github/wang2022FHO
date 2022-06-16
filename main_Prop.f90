PROGRAM MAIN
    
    USE GlobalDef
	
    IMPLICIT REAL(8) (A-H,O-Z)
    EXTERNAL :: YHCN,YHC,YHCM
    REAL(8) :: EE, step
    REAL(8) :: DD, T0, TF, PV(6), EP(1),Coe(6),DPV(6),FLI
    INTEGER(4) :: FLAG
103 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    CALL GlobalPar()    ! import system parameters
    
    CALL NewRaf(1, (/1D0/), 1D-14, EP)    ! get EP location within error margin 'Err'
    rS = EP(1)     ! set global variable
    WRITE(*,*) 'SEP         =   ', rS
    
    T0 = 0D0
    TF = 100D0
    EE = 1D-13
    DD = 1D-13
    step = 1D-2
    
    !R_Sun =  0.34D0*AU/Runit 
    !OrbitN = 1.484D0
    !OmegaN = -DSQRT((Ratio_Mass + 1D0) / OrbitN ** 3D0)
    
    CALL YHCA_order(Coe, 0D0, DPV,1)
    PV = (/0D0,rS,0D0,0D0,0D0,0D0/) + DPV
    !PV=(/0D0,rS+7D-2,0D0,3D-1,0D0,0D0/)
    
    CALL LCE_FLI(YHCM,PV,FLI,TF,FLAG)
    WRITE(*,*) 'FLI = ',FLI
    !FLAG = 1
    
    OPEN(11,file='Orbit1.DAT', status='REPLACE')
    write(11, 103) (/0D0,rS,0D0,0D0,0D0,0D0/), T0 
    DO WHILE(T0 < TF)
        
        ! Numerical Solution
        CALL RKF78(YHCN,step,T0,PV,EE,DD,6)
        !CALL CheckTotal(T0,PV,Flag)
        !IF (Flag==0)THEN
        !    EXIT
        !ENDIF
        write(11, 103) PV, T0       
        
    ENDDO
    !CALL RKF78(YHCN,TF-T0,T0,PV,EE,1D0,6)
    !CALL CheckTotal(T0,PV,Flag)
    !write(11, 103) PV, T0 
    
    CLOSE(11)
    
    END PROGRAM MAIN