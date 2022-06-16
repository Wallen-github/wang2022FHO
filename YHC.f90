    SUBROUTINE YHCN(T0, T1, X, Y)
    
    USE GlobalDef   ! Define glocal variables
        
    REAL(8) :: X(6), Y(6), T0, T1
    REAL(8) :: R, V0x, V0y, V0z, V1x, V1y, V1z
    REAL(8) :: phi, X_Prime, Y_Prime, Z_Prime, R_Prime, Delta, V2x, V2y, V2z, RR1, cosPhi
    REAL(8) :: phi_Sun, VxSRP, VySRP, VzSRP, r1, RM, SF
        
    R = DSQRT(X(1) ** 2D0 + X(2) ** 2D0 + X(3) ** 2D0)
    
    ! Two-body gravitation term
    
    V0x = - X(1) / R ** 3D0
    V0y = - X(2) / R ** 3D0
    V0z = - X(3) / R ** 3D0
        
    ! Non-spherical gravitation term truncted at 2-order
    
    V1x = - 3D0 * X(1) * C20 * (3D0 / 2D0 * X(3) ** 2D0 / R ** 2D0 - 1D0 / 2D0) / R ** 5D0 &
        & - 3D0 * C20 * X(3) ** 2D0 * X(1) / R ** 7D0 - 15D0 * C22 * X(1) * (X(1) ** 2D0 - X(2) ** 2D0) / R ** 7D0 &
        & + 6D0 * C22 * X(1) / R ** 5D0
        
    V1y = - 3D0 * X(2) * C20 * (3D0 / 2D0 * X(3) ** 2D0 / R ** 2D0 - 1D0 / 2D0) / R ** 5D0 &
        & - 3D0 * C20 * X(3) ** 2D0 * X(2) / R ** 7D0 - 15D0 * C22 * X(2) * (X(1) ** 2D0 - X(2) ** 2D0) / R ** 7D0 &
        & - 6D0 * C22 * X(2) / R ** 5D0
        
    V1z = - 3D0 * X(3) * C20 * (3D0 / 2D0 * X(3) ** 2D0 / R ** 2D0 - 1D0 / 2D0) / R ** 5D0 &
        & - 3D0 * C20 * X(3) ** 3D0 / R ** 7D0 - 15D0 * C22 * X(3) * (X(1) ** 2D0 - X(2) ** 2D0) / R ** 7D0 &
        & + 3D0 * C20 * X(3) / R ** 5D0
    
    ! Third-body perturbation term
    
    RM = Ratio_Mass
    phi = (OmegaN - OmegaA) * T1 + phi0
	X_Prime = OrbitN * DCOS(phi)
	Y_Prime = OrbitN * DSIN(phi)
	Z_Prime = 0D0
    R_Prime = OrbitN
    RR1 = (X(1) * X_Prime + X(2) * Y_Prime + X(3) * Z_Prime)
    cosPhi = RR1 / (R * OrbitN)
	Delta = DSQRT((X_Prime - X(1)) ** 2D0 + (Y_Prime - X(2)) ** 2D0 + (Z_Prime - X(3)) ** 2D0)
	
    V2x = RM * ((X_Prime - X(1)) / Delta ** 3D0 - X_Prime / OrbitN ** 3D0)
	V2y = RM * ((Y_Prime - X(2)) / Delta ** 3D0 - Y_Prime / OrbitN ** 3D0)
	V2z = RM * ((Z_Prime - X(3)) / Delta ** 3D0 - Z_Prime / OrbitN ** 3D0)
    
    !r1 = OrbitN
    !V2x = -(1D0/2D0)*(2*r1**4*X(1)+9*r1**2*X(1)**2*X_Prime-6*r1**2*X(1)*X_Prime**2+6*r1**2*X(1)*X(2)*Y_Prime+6*r1**2*X(1)*X(3)*Z_Prime+3*r1**2*X_Prime*X(2)**2-6*r1**2*X_Prime*X(2)*Y_Prime+3*r1**2*X_Prime*X(3)**2-6*r1**2*X_Prime*X(3)*Z_Prime-15*X(1)**2*X_Prime**3-30*X(1)*X_Prime**2*X(2)*Y_Prime-30*X(1)*X_Prime**2*X(3)*Z_Prime-15*X_Prime*X(2)**2*Y_Prime**2-30*X_Prime*X(2)*Y_Prime*X(3)*Z_Prime-15*X_Prime*X(3)**2*Z_Prime**2)*RM/r1**7
    !V2y = -(1D0/2D0)*(2*r1**4*X(2)+3*r1**2*X(1)**2*Y_Prime+6*r1**2*X(1)*X_Prime*X(2)-6*r1**2*X(1)*X_Prime*Y_Prime+9*r1**2*X(2)**2*Y_Prime-6*r1**2*X(2)*Y_Prime**2+6*r1**2*X(2)*X(3)*Z_Prime+3*r1**2*Y_Prime*X(3)**2-6*r1**2*Y_Prime*X(3)*Z_Prime-15*X(1)**2*X_Prime**2*Y_Prime-30*X(1)*X_Prime*X(2)*Y_Prime**2-30*X(1)*X_Prime*Y_Prime*X(3)*Z_Prime-15*X(2)**2*Y_Prime**3-30*X(2)*Y_Prime**2*X(3)*Z_Prime-15*Y_Prime*X(3)**2*Z_Prime**2)*RM/r1**7
    !V2z = -(1D0/2D0)*(2*r1**4*X(3)+3*r1**2*X(1)**2*Z_Prime+6*r1**2*X(1)*X_Prime*X(3)-6*r1**2*X(1)*X_Prime*Z_Prime+3*r1**2*X(2)**2*Z_Prime+6*r1**2*X(2)*Y_Prime*X(3)-6*r1**2*X(2)*Y_Prime*Z_Prime+9*r1**2*X(3)**2*Z_Prime-6*r1**2*X(3)*Z_Prime**2-15*X(1)**2*X_Prime**2*Z_Prime-30*X(1)*X_Prime*X(2)*Y_Prime*Z_Prime-30*X(1)*X_Prime*X(3)*Z_Prime**2-15*X(2)**2*Y_Prime**2*Z_Prime-30*X(2)*Y_Prime*X(3)*Z_Prime**2-15*X(3)**2*Z_Prime**3)*RM/r1**7
    
    ! Solar radiation pressure linearical term
	
    SF = ShadeFactor
    phi_Sun = - OmegaA * T1
	VxSRP = - SF * Beta_Sun / R_sun ** 2D0 * DCOS(theta_Sun) * DCOS(phi_Sun)
	VySRP = - SF * Beta_Sun / R_sun ** 2D0 * DCOS(theta_Sun) * DSIN(phi_Sun)
	VzSRP = - SF * Beta_Sun / R_sun ** 2D0 * DSIN(theta_Sun)
        
    Y(1) = X(4)
    Y(2) = X(5)
    Y(3) = X(6)
    Y(4) = 2D0 * X(5) + X(1) + V0x + V1x + V2x + VxSRP
    Y(5) = - 2D0 * X(4) + X(2) + V0y + V1y + V2y + VySRP
    Y(6) = V0z + V1z + V2z + VzSRP
    
    END SUBROUTINE YHCN
    
    SUBROUTINE YHCA_order(Coe, T, DPV,ORDER)
        
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: ORDER
    REAL(8) :: Coe(6), T, DPV(6)
    REAL(8) :: DPV0(6), DPV1(6), DPV2(6), DPVS2(6), DPV3(6)
        
    CALL YHC0(Coe, T, DPV0)
    CALL YHC1(T, DPV1)
    CALL YHC2_v2(T, DPV2)
    CALL YHC3_v2(T, DPV3)
    IF (ORDER==0)THEN 
        DPV = DPV0
    ELSEIF (ORDER==1)THEN
        DPV =  DPV1
    ELSEIF (ORDER==2)THEN
        DPV =  DPV1 + DPV2
    ELSEIF (ORDER==3)THEN
        DPV =  DPV1 + DPV2 + DPV3
    ENDIF
    
    END  SUBROUTINE
    
    SUBROUTINE YHCA(Coe, T, DPV)
        
    USE GlobalDef   ! Define glocal variables
        
    REAL(8) :: Coe(6), T, DPV(6)
    REAL(8) :: DPV0(6), DPV1(6), DPV2(6), DPVS2(6), DPV3(6)
        
    CALL YHC0(Coe, T, DPV0)
    CALL YHC1(T, DPV1)
    !CALL YHC2(T, DPV2)
    !CALL YHC3(T, DPV3)
    CALL YHC2_v2(T, DPV2)
    CALL YHC3_v2(T, DPV3)
    
    DPV =  DPV1 !+ DPV2 !+ DPV3!+ DPV0
    
    END  SUBROUTINE YHCA
    
    SUBROUTINE YHC0(Coe, T, DPV)
        
    USE GlobalDef   ! Define glocal variables
        
    REAL(8) :: DPV(6), Coe(6), Fre(3), T
    REAL(8) :: V0xx, kl, ks
        
    CALL GetFre(Fre)    ! Get eigen frequencies of dynamical system
        
    V0xx = - 1D0 / rS ** 3D0 + 3D0 / 2D0 * (C20 + 14D0 * C22) / rS ** 5D0
        
    kl = - (1D0 + V0xx + Fre(1) ** 2D0) / (2D0 * Fre(1))
    ks = - (1D0 + V0xx + Fre(2) ** 2D0) / (2D0 * Fre(2))
        
    DPV(1) = Coe(1) * DCOS(Fre(1) * T + Coe(4))     &
        &   + Coe(2) * DCOS(Fre(2) * T + Coe(5))
        
    DPV(2) = kl * Coe(1) * DSIN(Fre(1) * T + Coe(4))    &
        &   + ks * Coe(2) * DSIN(Fre(2) * T + Coe(5))
        
    DPV(3) = Coe(3) * DCOS(Fre(3) * T + Coe(6))
        
    DPV(4) = - Coe(1) * Fre(1) * DSIN(Fre(1) * T + Coe(4))  &
        &   - Coe(2) * Fre(2) * DSIN(Fre(2) * T + Coe(5))
        
    DPV(5) = kl * Coe(1) * Fre(1) * DCOS(Fre(1) * T + Coe(4))   &
        &   + ks * Coe(2) * Fre(2) * DCOS(Fre(2) * T + Coe(5))
        
    DPV(6) = - Coe(3) * Fre(3) * DSIN(Fre(3) * T + Coe(6))
    
    END SUBROUTINE YHC0
    
    SUBROUTINE YHC1(T, Y1)
        
    USE GlobalDef   ! Define glocal variables
        
    REAL(8) :: T, Y1(6)
    REAL(8) :: phi, phiA, R1, Omega0xx, Omega0yy, Omega0zz, Omega0yyy, omegaSun, omega, Omegaez
    REAL(8) :: RM, SF
    REAL(8) :: AlphaSxi, BetaSeta, Czeta
    REAL(8) :: Beta2xi_P2, Alpha2eta_P2, Ceta_P2
    REAL(8) :: Alpha1xi_P3, Beta1eta_P3, Alpha3xi_P3, Beta3eta_P3
    REAL(8) :: Beta2xi_P4, Alpha2eta_P4, Beta4xi_P4, Alpha4eta_P4, Ceta_P4
    REAL(8) :: Alpha1xi_P5, Beta1eta_P5, Alpha3xi_P5, Beta3eta_P5, Alpha5xi_P5, Beta5eta_P5
    
    ! Parameters
    
    SF = ShadeFactor
    RM = Ratio_Mass
    R1 = OrbitN
    omega = OmegaN - OmegaA
    phi =  omega * T + phi0
    omegaSun = - OmegaA
    phiA = omegaSun * T
    
    ! Difference of poteintal function
    
    Omega0xx = 1D0 - 1D0 / rS ** 3D0 + 3D0 / 2D0 * (C20 + 14D0 * C22) / rS ** 5D0
    Omega0yy = 1D0 + 2D0 / rS ** 3D0 - 6D0 * (C20 + 6D0 * C22) / rS ** 5D0
    Omega0zz = - 1D0 / rS ** 3D0 + 3D0 / 2D0 * (3D0 * C20 + 10D0 * C22) / rS ** 5D0
    Omegaez = - SF * Beta_Sun / R_sun ** 2D0 * DSIN(theta_Sun)
    Omega0yyy = - 6D0 / rS ** 4D0 + 30D0 * (C20 + 6D0 * C22) / rS ** 6D0
	
    ! Coefficients
    
    Alpha1xi_P3 = - 3D0 / 8D0 * RM * rS ** 2D0 * (omega ** 2D0 - 6D0 * omega &
        & + Omega0yy) / (r1 ** 4D0 * (omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta1eta_P3 = - 3D0 / 8D0 * RM * rS ** 2D0 * (3D0 * omega ** 2D0 - 2D0 * omega &
        & + 3D0 * Omega0xx) / (r1 ** 4D0 * (omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Alpha1xi_P5 = - 15D0 / 64D0 * RM * rS ** 4D0 * (omega ** 2D0 - 1D1 * omega &
        & + Omega0yy) / (r1 ** 6D0 * (omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta1eta_P5 = - 15D0 / 64D0 * RM * rS ** 4D0 * (5D0 * omega ** 2D0 - 2D0 * omega &
        & + 5D0 * Omega0xx) / (r1 ** 6D0 * (omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta2xi_P2 = - 3D0 / 2D0 * RM * rS * (4D0 * omega ** 2D0 - 4D0 * omega &
        & + Omega0yy) / (r1 ** 3D0 * (16D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 4D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Alpha2eta_P2 = 3D0 / 2D0 * RM * rS * (4D0 * omega ** 2D0 - 4D0 * omega &
        & + Omega0xx) / (r1 ** 3D0 * (16D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 4D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta2xi_P4 = - 5D0 / 8D0 * RM * rS ** 3D0 * (4D0 * omega ** 2D0 - 8D0 * omega &
        & + Omega0yy) / (r1 ** 5D0 * (16D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 4D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Alpha2eta_P4 = (5D0/4D0)*rS**3D0*RM*(4D0*omega**2D0-2D0*omega+Omega0xx) &
        & /(r1**5D0*(16D0*omega**4D0+4D0*omega**2D0*Omega0xx+4D0*omega**2D0*Omega0yy&
        & -16D0*omega**2D0+Omega0xx*Omega0yy))
    
    Alpha3xi_P3 = 15D0 * RM * rS ** 2D0 * (9D0 * omega ** 2D0 - 6D0 * omega &
        & + Omega0yy) / (8D0 * r1 ** 4D0 * (81D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 9D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta3eta_P3 = 15D0 * RM * rS ** 2D0 * (9D0 * omega ** 2D0 - 6D0 * omega &
        & + Omega0xx) / (8D0 * r1 ** 4D0 * (81D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 9D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Alpha3xi_P5 = 105D0 * RM * rS ** 4D0 * (9D0 * omega ** 2D0 - 1D1 * omega &
        & + Omega0yy) / (128D0 * r1 ** 6D0 * (81D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 9D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta3eta_P5 = 35D0 * RM * rS ** 4D0 * (45D0 * omega ** 2D0 - 18D0 * omega &
        & + 5D0 * Omega0xx) / (128D0 * r1 ** 6D0 * (81D0 * omega ** 4D0 + (Omega0xx + Omega0yy - 4D0) &
        & * 9D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    AlphaSxi = SF * DCOS(theta_Sun) * Beta_Sun * (omegaSun ** 2D0 + Omega0yy - 2D0 * omegaSun) &
        & / R_Sun ** 2D0 / (omegaSun ** 4D0 + (Omega0xx + Omega0yy - 4D0) * omegaSun ** 2D0 & 
        & + Omega0xx * Omega0yy)
    
    BetaSeta = SF * DCOS(theta_Sun) * Beta_Sun * (omegaSun ** 2D0 + Omega0xx - 2D0 * omegaSun) &
        & / R_Sun ** 2D0 / (omegaSun ** 4D0 + (Omega0xx + Omega0yy - 4D0) * omegaSun ** 2D0 & 
        & + Omega0xx * Omega0yy)
    
    Beta4xi_P4 = 35D0 * RM * rS ** 3D0 * (16D0 * omega ** 2D0 - 8D0 * omega + Omega0yy) &
        & / (16D0 * r1 ** 5D0 * (256D0 * omega ** 4D0 + 16D0 * omega ** 2D0 * Omega0xx &
        & + 16D0 * omega ** 2D0 * Omega0yy - 64D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Alpha4eta_P4 = - 35D0 * RM * rS ** 3D0 * (16D0 * omega ** 2D0 - 8D0 * omega + Omega0xx) &
        & / (16D0 * r1 ** 5D0 * (256D0 * omega ** 4D0 + 16D0 * omega ** 2D0 * Omega0xx &
        & + 16D0 * omega ** 2D0 * Omega0yy - 64D0 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Alpha5xi_P5 = - 315D0 * RM * rS ** 4D0 * (25D0 * omega ** 2D0 - 1D1 * omega + Omega0yy) &
        & / (128D0 * r1 ** 6D0 * (625D0 * omega ** 4D0 + 25D0 * omega ** 2D0 * Omega0xx &
        & + 25D0 * omega ** 2D0 * Omega0yy - 1D2 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Beta5eta_P5 = - 315D0 * RM * rS ** 4D0 * (25D0 * omega ** 2D0 - 1D1 * omega + Omega0xx) &
        & / (128D0 * r1 ** 6D0 * (625D0 * omega ** 4D0 + 25D0 * omega ** 2D0 * Omega0xx &
        & + 25D0 * omega ** 2D0 * Omega0yy - 1D2 * omega ** 2D0 + Omega0xx * Omega0yy))
    
    Ceta_P2 = - RM * rS / (2D0 * r1 ** 3D0 * Omega0yy)
    
    Ceta_P4 = - 9D0 * RM * rS ** 3D0 / (16D0 * r1 ** 5D0 * Omega0yy)

    Czeta = - Omegaez / Omega0zz
    
    ! First-oder solution
    
    Y1(1) = Beta2xi_P2 * DSIN(2D0 * phi) &
        !& + Alpha1xi_P3 * DCOS(phi) + Alpha3xi_P3 * DCOS(3D0 * phi) &
        !& + Beta2xi_P4 * DSIN(2D0 * phi) + Beta4xi_P4 * DSIN(4D0 * phi) &
        !& + Alpha1xi_P5 * DCOS(phi) + Alpha3xi_P5 * DCOS(3D0 * phi) + Alpha5xi_P5 * DCOS(5D0 * phi) &
        & + AlphaSxi * DCOS(phiA)
    
    Y1(2) = Alpha2eta_P2 * DCOS(2D0 * phi) + Ceta_P2 &
        !& + Beta1eta_P3 * DSIN(phi) + Beta3eta_P3 * DSIN(3D0 * phi) &
        !& + Alpha2eta_P4 * DCOS(2D0 * phi) + Alpha4eta_P4 * DCOS(4D0 * phi) + Ceta_P4 &
        !& + Beta1eta_P5 * DSIN(phi) + Beta3eta_P5 * DSIN(3D0 * phi) + Beta5eta_P5 * DSIN(5D0 * phi)&
        & + BetaSeta * DSIN(phiA) 
    
    Y1(3) = Czeta
    
    Y1(4) = 2D0 * Beta2xi_P2 * omega * DCOS(2D0 * phi) &
        !& - Alpha1xi_P3 * omega * DSIN(phi) - 3D0 * Alpha3xi_P3 * omega * DSIN(3D0 * phi) &
        !& + 2D0 * Beta2xi_P4 * omega * DCOS(2D0 * phi) + 4D0 * Beta4xi_P4 * omega * DCOS(4D0 * phi) &
        !& - Alpha1xi_P5 * omega * DSIN(phi) - 3D0 * omega * Alpha3xi_P5 * DSIN(3D0 * phi) - 5D0 * omega * Alpha5xi_P5 * DSIN(5D0 * phi) &
        & - AlphaSxi * omegaSun * DSIN(phiA)
    
    Y1(5) = - 2D0 * Alpha2eta_P2 * omega * DSIN(2D0 * phi) &
        !& + Beta1eta_P3 * omega * DCOS(phi) + 3D0 * Beta3eta_P3 * omega * DCOS(3D0 * phi) &
        !& - 2D0 * Alpha2eta_P4 * omega * DSIN(2D0 * phi) - 4D0 * Alpha4eta_P4 * omega * DSIN(4D0 * phi) &
        !& + Beta1eta_P5 * omega * DCOS(phi) + 3D0 * omega * Beta3eta_P5 * DCOS(3D0 * phi) + 5D0 * omega * Beta5eta_P5 * DCOS(5D0 * phi) &
        & + BetaSeta * omegaSun * DCOS(phiA)
    
    Y1(6) = 0D0
    
    END SUBROUTINE YHC1