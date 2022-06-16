    SUBROUTINE CheckImpact(X,FLAG)
    USE GlobalDef
    REAL(8) :: CRITERION,X(6)
    REAL(8) :: AL,BT,CA,SA2,SB2,SC2,R1,R2,R3,R5,X2,Y2,Z2
    INTEGER(4) :: FLAG
    FLAG = 1
    AL=Axis_a/Runit
	BT=Axis_b/Runit
	CA=Axis_c/Runit
	SA2=AL**2
	SB2=BT**2
	SC2=CA**2

    R1=DSQRT(X(1)**2+X(2)**2+X(3)**2)
    R2=R1*R1
    R3=R1*R2
    R5=R2*R3
    X2=X(1)**2
    Y2=X(2)**2
    Z2=X(3)**2
    CRITERION=X2/SA2+Y2/SB2+Z2/SC2-1D0
    IF(CRITERION.LT.0D0)THEN
        FLAG=0
    ENDIF
    
    END SUBROUTINE
    
    SUBROUTINE CheckImpact2(T0,X,FLAG)
    USE GlobalDef
    REAL(8) :: CRITERION,X(6),T0,X1(3),DX(3),phi,RR2
    REAL(8) :: AL,BT,CA,SA2,SB2,SC2,R1,R2,R3,R5,X2,Y2,Z2
    INTEGER(4) :: FLAG
    FLAG = 1
    AL=Axis_a/Runit
	BT=Axis_b/Runit
	CA=Axis_c/Runit
	SA2=AL**2
	SB2=BT**2
	SC2=CA**2

    R1=DSQRT(X(1)**2+X(2)**2+X(3)**2)
    R2=R1*R1
    R3=R1*R2
    R5=R2*R3
    X2=X(1)**2
    Y2=X(2)**2
    Z2=X(3)**2
    CRITERION=X2/SA2+Y2/SB2+Z2/SC2-1D0
    
    phi = (OmegaN - OmegaA) * T0 + phi0
	X1(1) = OrbitN * DCOS(phi)
	X1(2) = OrbitN * DSIN(phi)
	X1(3) = 0D0
    DX = X1-X(1:3)
    RR2 = DSQRT(DX(1)**2D0+DX(2)**2D0+DX(3)**2D0)
    
    IF(CRITERION.LT.0D0.OR.RR2.LT.RadiiB)THEN
        FLAG=0
    ENDIF
    
    END SUBROUTINE
    
    SUBROUTINE CheckEscape(X,FLAG)
    USE GlobalDef
    REAL(8) :: X(6)
    REAL(8) :: R1
    INTEGER(4) :: FLAG
    FLAG = 1

    R1=DSQRT(X(1)**2+X(2)**2+X(3)**2)

    IF(R1.GT.5D0*rS)THEN
        FLAG=0
    ENDIF
    
    END SUBROUTINE
    
    SUBROUTINE CheckOrbit(X,FLAG)
    USE GlobalDef
    REAL(8) :: X(6)
    REAL(8) :: R1,Theta
    INTEGER(4) :: FLAG
    FLAG = 1

    R1=DSQRT(X(1)**2+X(2)**2+X(3)**2)
    Theta = DATAN2(X(1),X(2))

    !IF(Theta.LT.-3D0)THEN
    !    FLAG=0
    !ENDIF
    
    END SUBROUTINE
    
    SUBROUTINE CheckTotal(T0,X,Flag)
    USE GlobalDef
    REAL(8) :: X(6),T0
    INTEGER(4) :: FLAG
    FLAG = 1
    
    CALL CheckImpact2(T0,X,FLAG)
    IF (FLAG==0) THEN
        !WRITE(*,*) 'Imapct'
        GOTO 201
    ENDIF
    CALL CheckEscape(X,FLAG)
    IF (FLAG==0) THEN
        WRITE(*,*) 'Escape'
        GOTO 201
    ENDIF
    CALL CheckOrbit(X,FLAG)
    IF (FLAG==0) THEN
        WRITE(*,*) 'Orbit'
        GOTO 201
    ENDIF
201 END