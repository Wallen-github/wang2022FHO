    PROGRAM MAIN
    
    USE GlobalDef
    EXTERNAL :: YHCM
    REAL(8) :: EP(1),EE
    REAL(8) :: Coe(6),DPV(6),X0(6)
    INTEGER(4) :: ORDER,KK,Num
    REAL(8) :: Ista,Iend,Jsta,Jend,stepI,stepJ
    REAL(8) :: FLI,TTOL,V3rd,VSRP
    
    
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
    
    Ista = 2D-1
    Iend = 4D-1
    Jsta = rS+0.2
    Jend = 3.5D0
    Num = 50
    stepI = (Iend-Ista)/Num
    stepJ = (Jend-Jsta)/Num
    
    OPEN(1,FILE="RELATION.DAT")
    OPEN(2,FILE="RELATIONXY.DAT")
    DO I = 1,Num
        DO J =1,Num
            WRITE(*,*) Ista + (I-1)*stepI, Jsta + (J-1)*stepJ
            
            R_Sun =  (Ista + (I-1)*stepI)*AU/Runit
            VSRP = Beta_Sun*rS**2D0/R_Sun ** 2D0
            
            OrbitN = Jsta + (J-1)*stepJ
            OmegaN = -DSQRT((Ratio_Mass + 1D0) / OrbitN ** 3D0)
            V3rd = Ratio_Mass*rS**3D0/OrbitN**3D0

            CALL YHCA_order(Coe, 0D0, DPV,ORDER)
            X0 = (/0D0,rS,0D0,0D0,0D0,0D0/) + DPV
            TTOL=1D2!ABS(PI2/OmegaN)*1D1
            
            CALL LCE_FLI(YHCM,X0,FLI,TTOL,KK)
            WRITE(2,'(I5,I5,F25.10,F25.10)') I,J,Ista + (I-1)*stepI,OrbitN
	        IF(KK.EQ.0)THEN
	            WRITE(1,'(I5,I5,F25.10,F25.10,F25.10,F25.10)') I,J,FLI,TTOL,V3rd,VSRP
                WRITE(*,*) I,J,FLI
	        ELSEIF(KK.EQ.1)THEN
	            WRITE(1,'(I5,I5,F25.10,F25.10,F25.10,F25.10)') I,J,0D0,TTOL,V3rd,VSRP
                WRITE(*,*) I,J,0D0
            ELSEIF(KK.EQ.2)THEN
                WRITE(1,'(I5,I5,F25.10,F25.10,F25.10,F25.10)') I,J,2D3,TTOL,V3rd,VSRP
                WRITE(*,*) I,J,2D4
            ELSEIF(KK.EQ.3)THEN
                WRITE(1,'(I5,I5,F25.10,F25.10,F25.10,F25.10)') I,J,1D3,TTOL,V3rd,VSRP
                WRITE(*,*) I,J,1D3
            ENDIF
            
        ENDDO
    ENDDO
    CLOSE(1)
    CLOSE(2)
            
    
    
    END PROGRAM