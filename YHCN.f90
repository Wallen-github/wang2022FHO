    !===================================================================
    !   v1.0:
    !       Change: 1.  A SUBROUTINE Force1B is added, which is used for 
    !               geting the 3-body force in body-fixed frame.
    !
    !               2.  the 'SUBROUTINE Force1B' is used in substituting
    !               origin code and bringing eccentricity
    !   v1.1:
    !       Change: 1.  The rotation matrix is corrected to the other direction
    !               in which theta -> -theta.
    !===================================================================
    
    SUBROUTINE YHC(T0, T1, X, Y)
    
    USE GlobalDef
        
    REAL(8) :: X(6), Y(6), T0, T1
    REAL(8) :: F(3), thetaA, F1(3), F1B(3), FS(3), SF
    
    ! 求解椭球引力
    !CALL Carlson(X(1:3), F)
    CALL Force0(X(1:3), F)
    
    ! 求第三体引力
    CALL Force1B(X(1:3), T1, F1B)
    !CALL Force1_P3(X(1:3), T1, F1B)
    !CALL Force1_P4(X(1:3), T1, F1B)
    !CALL Force1_P5(X(1:3), T1, F1B)
    
    ! 太阳光压
    thetaA = omegaA * T1
    CALL ForceS(thetaA, FS)
    
    Y(1) = X(4)
    Y(2) = X(5)
    Y(3) = X(6)
    Y(4) = 2D0 * OmegaA * X(5) + OmegaA**2D0 * X(1) + F(1) + F1B(1) + FS(1)
    Y(5) = - 2D0 * OmegaA * X(4) + OmegaA**2D0 * X(2) + F(2) + F1B(2) + FS(2)
    Y(6) = F(3) + F1B(3) + FS(3)

    END SUBROUTINE YHC
    
    SUBROUTINE RotationMatrix(thetaA,CC,CT)
    
    REAL(8) :: CC(3,3), CT(3,3),thetaA
    
    CT(1,1) = DCOS(thetaA)
    CT(1,2) = DSIN(thetaA)
    CT(1,3) = 0D0
    CT(2,1) = -DSIN(thetaA)
    CT(2,2) = DCOS(thetaA)
    CT(2,3) = 0D0
    CT(3,1) = 0D0
    CT(3,2) = 0D0
    CT(3,3) = 1D0
    CC(1:3,1)=CT(1,1:3)
	CC(1:3,2)=CT(2,1:3)
	CC(1:3,3)=CT(3,1:3) 
    
    END SUBROUTINE RotationMatrix
    
    SUBROUTINE Force1(XX, X1, F1)
    
    USE GlobalDef
    
    REAL(8) :: F1(3), XX(3), X1(3)
    REAL(8) :: RA1, RA2, RA3, XA(3)
    REAL(8) :: RB1, RB2, RB3, XB(3),RM
    
    RM = Ratio_Mass
    XA = X1 - XX
    XB = X1
    RA1=DSQRT(XA(1)**2+XA(2)**2+XA(3)**2)
	RA2=RA1*RA1
	RA3=RA1*RA2
	RB1=DSQRT(XB(1)**2+XB(2)**2+XB(3)**2)
	RB2=RB1*RB1
	RB3=RB1*RB2
	RB5=RB2*RB3
	F1(1)=RM*(XA(1)/RA3-XB(1)/RB3)
    F1(2)=RM*(XA(2)/RA3-XB(2)/RB3)
	F1(3)=RM*(XA(3)/RA3-XB(3)/RB3)
    
    END SUBROUTINE Force1
    
    SUBROUTINE Force1B(X, T1, F1B)
    
    USE GlobalDef
    
    REAL(8) :: F1B(3), X(3), T1
    REAL(8) :: phi, R, RR1, cosPhi, Delta, X1(3), RM

    RM = Ratio_Mass
    phi = (OmegaN - OmegaA) * T1 + phi0
    R = DSQRT(X(1) ** 2D0 + X(2) ** 2D0 + X(3) ** 2D0)
    
    X1(1) = OrbitN * COS(phi)
    X1(2) = OrbitN * SIN(phi)
    X1(3) = 0D0
    RR1 = (X(1) * X1(1) + X(2) * X1(2) + X(3) * X1(3))
    cosPhi = RR1 / (R * OrbitN)
	Delta = DSQRT((X1(1) - X(1)) ** 2D0 + (X1(2) - X(2)) ** 2D0 + (X1(3) - X(3)) ** 2D0)
    F1B(1) = RM * ((X1(1) - X(1)) / Delta ** 3D0 - X1(1) / OrbitN ** 3D0)
	F1B(2) = RM * ((X1(2) - X(2)) / Delta ** 3D0 - X1(2) / OrbitN ** 3D0)
	F1B(3) = RM * ((X1(3) - X(3)) / Delta ** 3D0 - X1(3) / OrbitN ** 3D0)
    
    END SUBROUTINE Force1B
    
    SUBROUTINE Force1_P3(X, T1, F1B)
    
    USE GlobalDef
    
    REAL(8) :: F1B(3), X(3), T1
    REAL(8) :: X1(3),phi,RM,r1
    
    RM = Ratio_Mass
    phi = (OmegaN - OmegaA) * T1 + phi0
    X1(1) = OrbitN * COS(phi)
    X1(2) = OrbitN * SIN(phi)
    X1(3) = 0D0
    
    r1 = OrbitN
    F1B(1) = -(1D0/2D0)*(2*r1**4*X(1)+9*r1**2*X(1)**2*X1(1)-6*r1**2*X(1)*X1(1)**2+6*r1**2*X(1)*X(2)*X1(2)+6*r1**2*X(1)*X(3)*X1(3)+3*r1**2*X1(1)*X(2)**2-6*r1**2*X1(1)*X(2)*X1(2)+3*r1**2*X1(1)*X(3)**2-6*r1**2*X1(1)*X(3)*X1(3)-15*X(1)**2*X1(1)**3-30*X(1)*X1(1)**2*X(2)*X1(2)-30*X(1)*X1(1)**2*X(3)*X1(3)-15*X1(1)*X(2)**2*X1(2)**2-30*X1(1)*X(2)*X1(2)*X(3)*X1(3)-15*X1(1)*X(3)**2*X1(3)**2)*RM/r1**7
    F1B(2) = -(1D0/2D0)*(2*r1**4*X(2)+3*r1**2*X(1)**2*X1(2)+6*r1**2*X(1)*X1(1)*X(2)-6*r1**2*X(1)*X1(1)*X1(2)+9*r1**2*X(2)**2*X1(2)-6*r1**2*X(2)*X1(2)**2+6*r1**2*X(2)*X(3)*X1(3)+3*r1**2*X1(2)*X(3)**2-6*r1**2*X1(2)*X(3)*X1(3)-15*X(1)**2*X1(1)**2*X1(2)-30*X(1)*X1(1)*X(2)*X1(2)**2-30*X(1)*X1(1)*X1(2)*X(3)*X1(3)-15*X(2)**2*X1(2)**3-30*X(2)*X1(2)**2*X(3)*X1(3)-15*X1(2)*X(3)**2*X1(3)**2)*RM/r1**7
    F1B(3) = -(1D0/2D0)*(2*r1**4*X(3)+3*r1**2*X(1)**2*X1(3)+6*r1**2*X(1)*X1(1)*X(3)-6*r1**2*X(1)*X1(1)*X1(3)+3*r1**2*X(2)**2*X1(3)+6*r1**2*X(2)*X1(2)*X(3)-6*r1**2*X(2)*X1(2)*X1(3)+9*r1**2*X(3)**2*X1(3)-6*r1**2*X(3)*X1(3)**2-15*X(1)**2*X1(1)**2*X1(3)-30*X(1)*X1(1)*X(2)*X1(2)*X1(3)-30*X(1)*X1(1)*X(3)*X1(3)**2-15*X(2)**2*X1(2)**2*X1(3)-30*X(2)*X1(2)*X(3)*X1(3)**2-15*X(3)**2*X1(3)**3)*RM/r1**7
    
    END SUBROUTINE
    
    SUBROUTINE Force1_P4(X, T1, F1B)
    
    USE GlobalDef
    
    REAL(8) :: F1B(3), X(3), T1
    REAL(8) :: X1(3),phi,RM,r1
    
    RM = Ratio_Mass
    phi = (OmegaN - OmegaA) * T1 + phi0
    X1(1) = OrbitN * COS(phi)
    X1(2) = OrbitN * SIN(phi)
    X1(3) = 0D0
    
    r1 = OrbitN
    F1B(1) = -(1D0/2D0)*(2*r1**6*X(1)-3*r1**4*X(1)**3+9*r1**4*X(1)**2*X1(1)-6*r1**4*X(1)*X1(1)**2-3*r1**4*X(1)*X(2)**2+6*r1**4*X(1)*X(2)*X1(2)-3*r1**4*X(1)*X(3)**2+6*r1**4*X(1)*X(3)*X1(3)+3*r1**4*X1(1)*X(2)**2-6*r1**4*X1(1)*X(2)*X1(2)+3*r1**4*X1(1)*X(3)**2-6*r1**4*X1(1)*X(3)*X1(3)+30*r1**2*X(1)**3*X1(1)**2-15*r1**2*X(1)**2*X1(1)**3+45*r1**2*X(1)**2*X1(1)*X(2)*X1(2)+45*r1**2*X(1)**2*X1(1)*X(3)*X1(3)+15*r1**2*X(1)*X1(1)**2*X(2)**2-30*r1**2*X(1)*X1(1)**2*X(2)*X1(2)+15*r1**2*X(1)*X1(1)**2*X(3)**2-30*r1**2*X(1)*X1(1)**2*X(3)*X1(3)+15*r1**2*X(1)*X(2)**2*X1(2)**2+30*r1**2*X(1)*X(2)*X1(2)*X(3)*X1(3)+15*r1**2*X(1)*X(3)**2*X1(3)**2+15*r1**2*X1(1)*X(2)**3*X1(2)-15*r1**2*X1(1)*X(2)**2*X1(2)**2+15*r1**2*X1(1)*X(2)**2*X(3)*X1(3)+15*r1**2*X1(1)*X(2)*X1(2)*X(3)**2-30*r1**2*X1(1)*X(2)*X1(2)*X(3)*X1(3)+15*r1**2*X1(1)*X(3)**3*X1(3)-15*r1**2*X1(1)*X(3)**2*X1(3)**2-35*X(1)**3*X1(1)**4-105*X(1)**2*X1(1)**3*X(2)*X1(2)-105*X(1)**2*X1(1)**3*X(3)*X1(3)-105*X(1)*X1(1)**2*X(2)**2*X1(2)**2-210*X(1)*X1(1)**2*X(2)*X1(2)*X(3)*X1(3)-105*X(1)*X1(1)**2*X(3)**2*X1(3)**2-35*X1(1)*X(2)**3*X1(2)**3-105*X1(1)*X(2)**2*X1(2)**2*X(3)*X1(3)-105*X1(1)*X(2)*X1(2)*X(3)**2*X1(3)**2-35*X1(1)*X(3)**3*X1(3)**3)*RM/r1**9

    F1B(2) = -(1D0/2D0)*(2*r1**6*X(2)-3*r1**4*X(1)**2*X(2)+3*r1**4*X(1)**2*X1(2)+6*r1**4*X(1)*X1(1)*X(2)-6*r1**4*X(1)*X1(1)*X1(2)-3*r1**4*X(2)**3+9*r1**4*X(2)**2*X1(2)-6*r1**4*X(2)*X1(2)**2-3*r1**4*X(2)*X(3)**2+6*r1**4*X(2)*X(3)*X1(3)+3*r1**4*X1(2)*X(3)**2-6*r1**4*X1(2)*X(3)*X1(3)+15*r1**2*X(1)**3*X1(1)*X1(2)+15*r1**2*X(1)**2*X1(1)**2*X(2)-15*r1**2*X(1)**2*X1(1)**2*X1(2)+15*r1**2*X(1)**2*X(2)*X1(2)**2+15*r1**2*X(1)**2*X1(2)*X(3)*X1(3)+45*r1**2*X(1)*X1(1)*X(2)**2*X1(2)-30*r1**2*X(1)*X1(1)*X(2)*X1(2)**2+30*r1**2*X(1)*X1(1)*X(2)*X(3)*X1(3)+15*r1**2*X(1)*X1(1)*X1(2)*X(3)**2-30*r1**2*X(1)*X1(1)*X1(2)*X(3)*X1(3)+30*r1**2*X(2)**3*X1(2)**2-15*r1**2*X(2)**2*X1(2)**3+45*r1**2*X(2)**2*X1(2)*X(3)*X1(3)+15*r1**2*X(2)*X1(2)**2*X(3)**2-30*r1**2*X(2)*X1(2)**2*X(3)*X1(3)+15*r1**2*X(2)*X(3)**2*X1(3)**2+15*r1**2*X1(2)*X(3)**3*X1(3)-15*r1**2*X1(2)*X(3)**2*X1(3)**2-35*X(1)**3*X1(1)**3*X1(2)-105*X(1)**2*X1(1)**2*X(2)*X1(2)**2-105*X(1)**2*X1(1)**2*X1(2)*X(3)*X1(3)-105*X(1)*X1(1)*X(2)**2*X1(2)**3-210*X(1)*X1(1)*X(2)*X1(2)**2*X(3)*X1(3)-105*X(1)*X1(1)*X1(2)*X(3)**2*X1(3)**2-35*X(2)**3*X1(2)**4-105*X(2)**2*X1(2)**3*X(3)*X1(3)-105*X(2)*X1(2)**2*X(3)**2*X1(3)**2-35*X1(2)*X(3)**3*X1(3)**3)*RM/r1**9

    F1B(3) = -(1D0/2D0)*(2*r1**6*X(3)-3*r1**4*X(1)**2*X(3)+3*r1**4*X(1)**2*X1(3)+6*r1**4*X(1)*X1(1)*X(3)-6*r1**4*X(1)*X1(1)*X1(3)-3*r1**4*X(2)**2*X(3)+3*r1**4*X(2)**2*X1(3)+6*r1**4*X(2)*X1(2)*X(3)-6*r1**4*X(2)*X1(2)*X1(3)-3*r1**4*X(3)**3+9*r1**4*X(3)**2*X1(3)-6*r1**4*X(3)*X1(3)**2+15*r1**2*X(1)**3*X1(1)*X1(3)+15*r1**2*X(1)**2*X1(1)**2*X(3)-15*r1**2*X(1)**2*X1(1)**2*X1(3)+15*r1**2*X(1)**2*X(2)*X1(2)*X1(3)+15*r1**2*X(1)**2*X(3)*X1(3)**2+15*r1**2*X(1)*X1(1)*X(2)**2*X1(3)+30*r1**2*X(1)*X1(1)*X(2)*X1(2)*X(3)-30*r1**2*X(1)*X1(1)*X(2)*X1(2)*X1(3)+45*r1**2*X(1)*X1(1)*X(3)**2*X1(3)-30*r1**2*X(1)*X1(1)*X(3)*X1(3)**2+15*r1**2*X(2)**3*X1(2)*X1(3)+15*r1**2*X(2)**2*X1(2)**2*X(3)-15*r1**2*X(2)**2*X1(2)**2*X1(3)+15*r1**2*X(2)**2*X(3)*X1(3)**2+45*r1**2*X(2)*X1(2)*X(3)**2*X1(3)-30*r1**2*X(2)*X1(2)*X(3)*X1(3)**2+30*r1**2*X(3)**3*X1(3)**2-15*r1**2*X(3)**2*X1(3)**3-35*X(1)**3*X1(1)**3*X1(3)-105*X(1)**2*X1(1)**2*X(2)*X1(2)*X1(3)-105*X(1)**2*X1(1)**2*X(3)*X1(3)**2-105*X(1)*X1(1)*X(2)**2*X1(2)**2*X1(3)-210*X(1)*X1(1)*X(2)*X1(2)*X(3)*X1(3)**2-105*X(1)*X1(1)*X(3)**2*X1(3)**3-35*X(2)**3*X1(2)**3*X1(3)-105*X(2)**2*X1(2)**2*X(3)*X1(3)**2-105*X(2)*X1(2)*X(3)**2*X1(3)**3-35*X(3)**3*X1(3)**4)*RM/r1**9
    
    END SUBROUTINE
    
    SUBROUTINE Force1_P5(X, T1, F1B)
    
    USE GlobalDef
    
    REAL(8) :: F1B(3), X(3), T1
    REAL(8) :: X1(3),phi,RM,r1,F1P4(3)
    
    RM = Ratio_Mass
    phi = (OmegaN - OmegaA) * T1 + phi0
    X1(1) = OrbitN * COS(phi)
    X1(2) = OrbitN * SIN(phi)
    X1(3) = 0D0
    
    CALL Force1_P4(X, T1, F1P4)
    
    r1 = OrbitN
    F1B(1) = F1P4(1) + (5D0/8D0)*RM*(15*r1**4*X(1)**4*X1(1)+12*r1**4*X(1)**3*X(2)*X1(2)+12*r1**4*X(1)**3*X(3)*X1(3)+18*r1**4*X(1)**2*X1(1)*X(2)**2+18*r1**4*X(1)**2*X1(1)*X(3)**2+12*r1**4*X(1)*X(2)**3*X1(2)+12*r1**4*X(1)*X(2)**2*X(3)*X1(3)+12*r1**4*X(1)*X(2)*X1(2)*X(3)**2+12*r1**4*X(1)*X(3)**3*X1(3)+3*r1**4*X1(1)*X(2)**4+6*r1**4*X1(1)*X(2)**2*X(3)**2+3*r1**4*X1(1)*X(3)**4-70*r1**2*X(1)**4*X1(1)**3-168*r1**2*X(1)**3*X1(1)**2*X(2)*X1(2)-168*r1**2*X(1)**3*X1(1)**2*X(3)*X1(3)-42*r1**2*X(1)**2*X1(1)**3*X(2)**2-42*r1**2*X(1)**2*X1(1)**3*X(3)**2-126*r1**2*X(1)**2*X1(1)*X(2)**2*X1(2)**2-252*r1**2*X(1)**2*X1(1)*X(2)*X1(2)*X(3)*X1(3)-126*r1**2*X(1)**2*X1(1)*X(3)**2*X1(3)**2-84*r1**2*X(1)*X1(1)**2*X(2)**3*X1(2)-84*r1**2*X(1)*X1(1)**2*X(2)**2*X(3)*X1(3)-84*r1**2*X(1)*X1(1)**2*X(2)*X1(2)*X(3)**2-84*r1**2*X(1)*X1(1)**2*X(3)**3*X1(3)-28*r1**2*X(1)*X(2)**3*X1(2)**3-84*r1**2*X(1)*X(2)**2*X1(2)**2*X(3)*X1(3)-84*r1**2*X(1)*X(2)*X1(2)*X(3)**2*X1(3)**2-28*r1**2*X(1)*X(3)**3*X1(3)**3-42*r1**2*X1(1)*X(2)**4*X1(2)**2-84*r1**2*X1(1)*X(2)**3*X1(2)*X(3)*X1(3)-42*r1**2*X1(1)*X(2)**2*X1(2)**2*X(3)**2-42*r1**2*X1(1)*X(2)**2*X(3)**2*X1(3)**2-84*r1**2*X1(1)*X(2)*X1(2)*X(3)**3*X1(3)-42*r1**2*X1(1)*X(3)**4*X1(3)**2+63*X(1)**4*X1(1)**5+252*X(1)**3*X1(1)**4*X(2)*X1(2)+252*X(1)**3*X1(1)**4*X(3)*X1(3)+378*X(1)**2*X1(1)**3*X(2)**2*X1(2)**2+756*X(1)**2*X1(1)**3*X(2)*X1(2)*X(3)*X1(3)+378*X(1)**2*X1(1)**3*X(3)**2*X1(3)**2+252*X(1)*X1(1)**2*X(2)**3*X1(2)**3+756*X(1)*X1(1)**2*X(2)**2*X1(2)**2*X(3)*X1(3)+756*X(1)*X1(1)**2*X(2)*X1(2)*X(3)**2*X1(3)**2+252*X(1)*X1(1)**2*X(3)**3*X1(3)**3+63*X1(1)*X(2)**4*X1(2)**4+252*X1(1)*X(2)**3*X1(2)**3*X(3)*X1(3)+378*X1(1)*X(2)**2*X1(2)**2*X(3)**2*X1(3)**2+252*X1(1)*X(2)*X1(2)*X(3)**3*X1(3)**3+63*X1(1)*X(3)**4*X1(3)**4)/r1**11

    F1B(2) = F1P4(2) + (5D0/8D0)*RM*(3*r1**4*X(1)**4*X1(2)+12*r1**4*X(1)**3*X1(1)*X(2)+18*r1**4*X(1)**2*X(2)**2*X1(2)+12*r1**4*X(1)**2*X(2)*X(3)*X1(3)+6*r1**4*X(1)**2*X1(2)*X(3)**2+12*r1**4*X(1)*X1(1)*X(2)**3+12*r1**4*X(1)*X1(1)*X(2)*X(3)**2+15*r1**4*X(2)**4*X1(2)+12*r1**4*X(2)**3*X(3)*X1(3)+18*r1**4*X(2)**2*X1(2)*X(3)**2+12*r1**4*X(2)*X(3)**3*X1(3)+3*r1**4*X1(2)*X(3)**4-42*r1**2*X(1)**4*X1(1)**2*X1(2)-28*r1**2*X(1)**3*X1(1)**3*X(2)-84*r1**2*X(1)**3*X1(1)*X(2)*X1(2)**2-84*r1**2*X(1)**3*X1(1)*X1(2)*X(3)*X1(3)-126*r1**2*X(1)**2*X1(1)**2*X(2)**2*X1(2)-84*r1**2*X(1)**2*X1(1)**2*X(2)*X(3)*X1(3)-42*r1**2*X(1)**2*X1(1)**2*X1(2)*X(3)**2-42*r1**2*X(1)**2*X(2)**2*X1(2)**3-84*r1**2*X(1)**2*X(2)*X1(2)**2*X(3)*X1(3)-42*r1**2*X(1)**2*X1(2)*X(3)**2*X1(3)**2-168*r1**2*X(1)*X1(1)*X(2)**3*X1(2)**2-252*r1**2*X(1)*X1(1)*X(2)**2*X1(2)*X(3)*X1(3)-84*r1**2*X(1)*X1(1)*X(2)*X1(2)**2*X(3)**2-84*r1**2*X(1)*X1(1)*X(2)*X(3)**2*X1(3)**2-84*r1**2*X(1)*X1(1)*X1(2)*X(3)**3*X1(3)-70*r1**2*X(2)**4*X1(2)**3-168*r1**2*X(2)**3*X1(2)**2*X(3)*X1(3)-42*r1**2*X(2)**2*X1(2)**3*X(3)**2-126*r1**2*X(2)**2*X1(2)*X(3)**2*X1(3)**2-84*r1**2*X(2)*X1(2)**2*X(3)**3*X1(3)-28*r1**2*X(2)*X(3)**3*X1(3)**3-42*r1**2*X1(2)*X(3)**4*X1(3)**2+63*X(1)**4*X1(1)**4*X1(2)+252*X(1)**3*X1(1)**3*X(2)*X1(2)**2+252*X(1)**3*X1(1)**3*X1(2)*X(3)*X1(3)+378*X(1)**2*X1(1)**2*X(2)**2*X1(2)**3+756*X(1)**2*X1(1)**2*X(2)*X1(2)**2*X(3)*X1(3)+378*X(1)**2*X1(1)**2*X1(2)*X(3)**2*X1(3)**2+252*X(1)*X1(1)*X(2)**3*X1(2)**4+756*X(1)*X1(1)*X(2)**2*X1(2)**3*X(3)*X1(3)+756*X(1)*X1(1)*X(2)*X1(2)**2*X(3)**2*X1(3)**2+252*X(1)*X1(1)*X1(2)*X(3)**3*X1(3)**3+63*X(2)**4*X1(2)**5+252*X(2)**3*X1(2)**4*X(3)*X1(3)+378*X(2)**2*X1(2)**3*X(3)**2*X1(3)**2+252*X(2)*X1(2)**2*X(3)**3*X1(3)**3+63*X1(2)*X(3)**4*X1(3)**4)/r1**11

    F1B(3) = F1P4(3) + (5D0/8D0)*RM*(3*r1**4*X(1)**4*X1(3)+12*r1**4*X(1)**3*X1(1)*X(3)+6*r1**4*X(1)**2*X(2)**2*X1(3)+12*r1**4*X(1)**2*X(2)*X1(2)*X(3)+18*r1**4*X(1)**2*X(3)**2*X1(3)+12*r1**4*X(1)*X1(1)*X(2)**2*X(3)+12*r1**4*X(1)*X1(1)*X(3)**3+3*r1**4*X(2)**4*X1(3)+12*r1**4*X(2)**3*X1(2)*X(3)+18*r1**4*X(2)**2*X(3)**2*X1(3)+12*r1**4*X(2)*X1(2)*X(3)**3+15*r1**4*X(3)**4*X1(3)-42*r1**2*X(1)**4*X1(1)**2*X1(3)-28*r1**2*X(1)**3*X1(1)**3*X(3)-84*r1**2*X(1)**3*X1(1)*X(2)*X1(2)*X1(3)-84*r1**2*X(1)**3*X1(1)*X(3)*X1(3)**2-42*r1**2*X(1)**2*X1(1)**2*X(2)**2*X1(3)-84*r1**2*X(1)**2*X1(1)**2*X(2)*X1(2)*X(3)-126*r1**2*X(1)**2*X1(1)**2*X(3)**2*X1(3)-42*r1**2*X(1)**2*X(2)**2*X1(2)**2*X1(3)-84*r1**2*X(1)**2*X(2)*X1(2)*X(3)*X1(3)**2-42*r1**2*X(1)**2*X(3)**2*X1(3)**3-84*r1**2*X(1)*X1(1)*X(2)**3*X1(2)*X1(3)-84*r1**2*X(1)*X1(1)*X(2)**2*X1(2)**2*X(3)-84*r1**2*X(1)*X1(1)*X(2)**2*X(3)*X1(3)**2-252*r1**2*X(1)*X1(1)*X(2)*X1(2)*X(3)**2*X1(3)-168*r1**2*X(1)*X1(1)*X(3)**3*X1(3)**2-42*r1**2*X(2)**4*X1(2)**2*X1(3)-28*r1**2*X(2)**3*X1(2)**3*X(3)-84*r1**2*X(2)**3*X1(2)*X(3)*X1(3)**2-126*r1**2*X(2)**2*X1(2)**2*X(3)**2*X1(3)-42*r1**2*X(2)**2*X(3)**2*X1(3)**3-168*r1**2*X(2)*X1(2)*X(3)**3*X1(3)**2-70*r1**2*X(3)**4*X1(3)**3+63*X(1)**4*X1(1)**4*X1(3)+252*X(1)**3*X1(1)**3*X(2)*X1(2)*X1(3)+252*X(1)**3*X1(1)**3*X(3)*X1(3)**2+378*X(1)**2*X1(1)**2*X(2)**2*X1(2)**2*X1(3)+756*X(1)**2*X1(1)**2*X(2)*X1(2)*X(3)*X1(3)**2+378*X(1)**2*X1(1)**2*X(3)**2*X1(3)**3+252*X(1)*X1(1)*X(2)**3*X1(2)**3*X1(3)+756*X(1)*X1(1)*X(2)**2*X1(2)**2*X(3)*X1(3)**2+756*X(1)*X1(1)*X(2)*X1(2)*X(3)**2*X1(3)**3+252*X(1)*X1(1)*X(3)**3*X1(3)**4+63*X(2)**4*X1(2)**4*X1(3)+252*X(2)**3*X1(2)**3*X(3)*X1(3)**2+378*X(2)**2*X1(2)**2*X(3)**2*X1(3)**3+252*X(2)*X1(2)*X(3)**3*X1(3)**4+63*X(3)**4*X1(3)**5)/r1**11
    
    END SUBROUTINE
    
    SUBROUTINE Force0(X, V0)
    
    USE GlobalDef
    
    REAL(8) :: X(3), V0(3), R, R2, R3, R5
    REAL(8) :: V0x, V0y, V0z, V1x, V1y, V1z
    REAL(8) :: F0(3), F2(3), TEM1, TEM2, X2, Y2, Z2
    
    R = DSQRT(X(1) ** 2D0 + X(2) ** 2D0 + X(3) ** 2D0)
    R2 = R ** 2D0
    R3 = R2 * R
    R5 = R3 * R2
    X2 = X(1) ** 2D0
    Y2 = X(2) ** 2D0
    Z2 = X(3) ** 2D0
    
    F0(1)=-X(1)/R3
	F0(2)=-X(2)/R3
	F0(3)=-X(3)/R3
	TEM1=5D0*Z2/R2-1D0
	TEM2=5D0*(X2-Y2)/R2
	F2(1)=-X(1)/R5*(3D0/2D0*C20*TEM1+3D0*C22*(TEM2-2D0))
	F2(2)=-X(2)/R5*(3D0/2D0*C20*TEM1+3D0*C22*(TEM2+2D0))
	F2(3)=-X(3)/R5*(3D0/2D0*C20*(TEM1-2D0)+3D0*C22*TEM2)
    
    V0 = F0 + F2
    
    END SUBROUTINE Force0
    
    SUBROUTINE ForceS(thetaA, FS)
    
    USE GlobalDef
    
    REAL(8) :: thetaA, phiS, FS(3),SF
    SF = ShadeFactor
    phiS = -thetaA + phiS0
    
	FS(1) = - SF * Beta_Sun / R_sun ** 2D0 * DCOS(theta_Sun) * DCOS(phiS)
	FS(2) = - SF * Beta_Sun / R_sun ** 2D0 * DCOS(theta_Sun) * DSIN(phiS)
	FS(3) = - SF * Beta_Sun / R_sun ** 2D0 * DSIN(theta_Sun)
    
    END SUBROUTINE ForceS