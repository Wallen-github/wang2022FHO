    PROGRAM MAIN
    USE GlobalDef
    EXTERNAL :: YHCM
    REAL(8) :: EP(1),EE
    REAL(8) :: Coe(6),DPV(6),X0(6)
    INTEGER(4) :: ORDER,KK
    REAL(8) :: FLI,TTOL
    
    CALL GlobalPar()
    WRITE(*,*) '==============================================================================='
    WRITE(*,*) 'Equilibrium Points'
    WRITE(*,*) '-------------------------------------------------------------------------------'
    EE = 1E-15    ! set the error margin for all code
    CALL NewRaf(1, (/1.6D0/), EE, EP)    ! get EP location within error margin 'Err'
    rS = EP(1)     ! set global variable
    WRITE(*,*) 'SEP         =   ', rS
    
    Coe = 0D0
    ORDER = 1
    CALL YHCA_order(Coe, 0D0, DPV,ORDER)
    X0 = (/0D0,rS,0D0,0D0,0D0,0D0/) + DPV
    WRITE(*,*) 'X0 = '
    WRITE(*,*) X0
    
    II=1
    JJ=2
    TTOL=PI2*10D0
    CALL LCE_FLI(YHCM,X0,FLI,TTOL,KK)
    IF(KK.EQ.0)THEN
        WRITE(*,'(I5,I5,F25.10)') II,JJ,FLI
    ELSE
        WRITE(*,'(I5,I5,F25.10)') II,JJ,0D0
    ENDIF
        
    END PROGRAM