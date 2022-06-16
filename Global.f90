    MODULE GlobalDef
    
        REAL(8),PARAMETER :: Gravity = 6.67D-11, AU=1.4959787D011
        REAL(8),PARAMETER :: PI=3.14159265358979323846D0,PI2=PI*2D0
        REAL(8) :: C20, C22, OmegaA, rS, RadiiA
        REAL(8) :: Tunit, Munit, Runit
        REAL(8) :: RadiiB, Ratio_Mass, OmegaN, OrbitN, phi0
        REAL(8) :: R_Sun, Beta_Sun, ShadeFactor, theta_Sun, phiS0
        REAL(8) :: Axis_a, Axis_b, Axis_c
    
    END MODULE

    SUBROUTINE GlobalPar()
    
        USE GlobalDef
    
        REAL(8) :: rhoA, rhoB, Omega_Orbit, ratio_ab, ratio_bc
        REAL(8) :: Mass_Primary, Omega_Primary, Radii_Primary, C20unit, C22unit 
        REAL(8) :: Mass_Secondary, Omega_Secondary, Radii_Secondary, SemiAxis
        REAL(8) :: Area_Mass, kappa, PSRP, rSun
        
        !WRITE(*,*) '#01 (22) Kalliope and Linus       Unstable but temperary maintains'
        !Radii_Primary = 166.2D3/2D0
        !ratio_ab = 1.32D0
        !ratio_bc = 1.2D0
        !rhoA = 3380D0
        !rhoB = 3380D0
        !RP_Primary = 4.15D0
        !Radii_Secondary = 28D3/2D0
        !SemiAxis = 1098D3
        !rSun = 2.9D0*AU
        !theta_Sun = 13.70D0 * PI / 18D1
        
        !WRITE(*,*) '#02 (41) Daphne and Peneius        Unstable'
        !Radii_Primary = 174D3/2D0
        !ratio_ab = 1.37D0
        !ratio_bc = 1.0D0
        !rhoA = 1950D0
        !rhoB = 1950D0
        !RP_Primary = 5.98D0
        !Radii_Secondary = 2D3/2D0
        !SemiAxis = 443D3
        !rSun = 2.76D0*AU
        !theta_Sun = 15.79D0 * PI / 18D1
        
        !WRITE(*,*) '#03 (317) Roxane and Olympias      Unstable but temperary maintains'
        !Radii_Primary = 19.86D3/2D0
        !ratio_ab = 1.71D0
        !ratio_bc = 1.18D0
        !rhoA = 2160D0
        !rhoB = 2160D0
        !RP_Primary = 8.169D0
        !Radii_Secondary = 5.2D3/2D0
        !SemiAxis = 245D3
        !rSun = 2.28D0*AU
        !theta_Sun = 1.76D0 * PI / 18D1
        
        !WRITE(*,*) '#04 (939) Isberga                   Unstable but temperary maintains'
        !Radii_Primary = 12.4D3/2D0
        !ratio_ab = 1.23D0
        !ratio_bc = 1.06D0
        !rhoA = 2910D0
        !rhoB = 2910D0
        !RP_Primary = 2.92D0
        !Radii_Secondary = 3.6D3/2D0
        !SemiAxis = 33D3
        !rSun = 2.25D0*AU
        !theta_Sun = 2.59D0 * PI / 18D1
        
        !WRITE(*,*) '#05 (1862) Apollo                   Unstable but temperary maintains'
        !Radii_Primary = 1.55D3/2D0
        !ratio_ab = 1.2D0
        !ratio_bc = 1.15D0
        !rhoA = 2050D0
        !rhoB = 2050D0
        !RP_Primary = 3.065D0
        !Radii_Secondary = 0.08D3/2D0
        !SemiAxis = 3.75D3
        !rSun = 1.47D0*AU
        !theta_Sun = 6.35D0 * PI / 18D1
        
        !WRITE(*,*) '#06 (65803) Didymos and Dimorphos      Inside Primary'
        !Radii_Primary = 0.78D3/2D0
        !ratio_ab = 1.01D0
        !ratio_bc = 1.06D0
        !rhoA = 2170D0
        !rhoB = 2170D0
        !RP_Primary = 2.26D0
        !Radii_Secondary = 0.17D3/2D0
        !SemiAxis = 1.19D3
        !rSun = 1.644D0*AU
        !theta_Sun = 3.41D0 * PI / 18D1
        
        WRITE(*,*) '#07 (66391) Moshup and Squannit      Stable'
        Radii_Primary = 1.317D3/2D0
        ratio_ab = 1.02D0
        ratio_bc = 1.11D0
        rhoA = 1970D0
        rhoB = 1300D0
        RP_Primary = 3.7645D0
        Radii_Secondary = 0.59D3/2D0
        SemiAxis = 2.548D3!1.4D3!
        rSun = 0.6423D0*AU!0.4*AU!
        theta_Sun = 38.87D0 * PI / 18D1
        
        !WRITE(*,*) 'Test'
        !Radii_Primary = 389.91451117656855406D0!454.545454545455D0!389.91451117656855406D0
        !ratio_ab = 4D2/39D1
        !ratio_bc = 39D1/38D1
        !rhoA = 2D3
        !rhoB = 2D3
        !RP_Primary = 3D0
        !Radii_Secondary = 8D1
        !SemiAxis = 2D3
        !rSun = AU
        !theta_Sun = 10D0 * PI / 18D1
        
        !WRITE(*,*) '(17260) Kusnirak      UnStable but temperary maintains'
        !Radii_Primary = 4.62D3/2D0
        !ratio_ab = 1.05D0
        !ratio_bc = 1.1D0
        !rhoA = 1600D0
        !rhoB = 1600D0
        !RP_Primary = 3.1287D0
        !Radii_Secondary = 1.14D3/2D0
        !SemiAxis = 7.4D3
        !rSun = 2.204534*AU
        !theta_Sun = 5.283D0 * PI / 18D1
        
        !WRITE(*,*) '(1052) Belgica                  Stable'
        !Radii_Primary = 9.79D3/2D0
        !ratio_ab = 1.02D0
        !ratio_bc = 1.11D0
        !rhoA = 1600D0
        !rhoB = 1600D0
        !RP_Primary = 2.7097D0
        !Radii_Secondary = 3.53D3/2D0
        !SemiAxis = 34D3
        !rSun = 2.2358*AU
        !theta_Sun = 4.695D0 * PI / 18D1
        
        !WRITE(*,*) '(702) Alauda and Pichi unem      Stable  ignore the secondary gravity'
        !Radii_Primary = 201.96D3/2D0
        !ratio_ab = 1.02D0
        !ratio_bc = 1.11D0
        !rhoA = 1570D0
        !rhoB = 1570D0
        !RP_Primary = 8.354D0
        !Radii_Secondary = 3.51D3/2D0
        !SemiAxis = 1227D3
        !rSun = 3.19*AU
        !theta_Sun = 20.611D0 * PI / 18D1
        
        !WRITE(*,*) '(48639) 1995 TL8      Stable  ignore the SRP'
        !Radii_Primary = 176D3/2D0
        !ratio_ab = 1.02D0
        !ratio_bc = 1.11D0
        !rhoA = 1000D0
        !rhoB = 1000D0
        !RP_Primary = 4.354D0
        !Radii_Secondary = 80D3/2D0
        !SemiAxis = 420D3
        !rSun = 52.560*AU
        !theta_Sun = 0.24245D0 * PI / 18D1
        
        !WRITE(*,*) '#08 (88710) 2001 SL9                    Inside Primary'
        !Radii_Primary = 0.77D3/2D0
        !ratio_ab = 1.07D0
        !ratio_bc = 1.64D0
        !rhoA = 1800D0
        !rhoB = 1800D0
        !RP_Primary = 2.4D0
        !Radii_Secondary = 0.18D3/2D0
        !SemiAxis = 1.35D3
        !rSun = 1.06D0*AU
        !theta_Sun = 21.90D0 * PI / 18D1
        
        !WRITE(*,*) '#09 (175706) 1996 FG3                   Stable'
        !Radii_Primary = 1.69D3/2D0
        !ratio_ab = 1.06D0
        !ratio_bc = 1.23D0
        !rhoA = 1300D0
        !rhoB = 1300D0
        !RP_Primary = 3.595D0
        !Radii_Secondary = 0.49D3/2D0
        !SemiAxis = 3D3
        !rSun = 1.054D0*AU
        !theta_Sun = 1.99D0 * PI / 18D1
        
        !WRITE(*,*) '#10 (617) Patroclus and Menoetius         V3rd is too large'
        !Radii_Primary = 126.6D3/2D0
        !ratio_ab = 1.036D0
        !ratio_bc = 1.028D0
        !rhoA = 810D0
        !rhoB = 810D0
        !RP_Primary = 102.595D0
        !Radii_Secondary = 111.8D3/2D0
        !SemiAxis = 695D3
        !rSun = 22.06D0*AU
        
        !WRITE(*,*) '(18503) 1996 PY4         '
        !Radii_Primary = 3.43D3/2D0
        !ratio_ab = 1.02D0
        !ratio_bc = 1.11D0
        !rhoA = 1.6D3
        !rhoB = 1.6D3
        !RP_Primary = 3.4391D0
        !Radii_Secondary = 0.82D3/2D0
        !SemiAxis = 6.7D3
        !rSun = 2.3611574399D0*AU
        
        ! Solar Pressure Parameters
        
        ShadeFactor = 1D0
        !theta_Sun = 10D0 * PI / 18D1
        Area_Mass = 1D-2
        kappa = 1.44D0
        PSRP = 4.5605D-6
        phiS0 = 0D0
        phi0 = 0D0!PI/2D0
        
        ! Other Parameters
        
        !Radii_Primary = (Axis_c * Axis_b * Axis_a) ** (1D0 / 3D0)
        Omega_Primary = 2D0 * PI / (RP_Primary * 6D1 * 6D1)
        Axis_c = (Radii_Primary**3D0/(ratio_ab*ratio_bc**2D0))**(1D0/3D0)
        Axis_b = Axis_c * ratio_bc
        Axis_a = Axis_b * ratio_ab
        Mass_Primary = 4D0*PI/3D0 * Radii_Primary ** 3D0 * rhoA
        Mass_Secondary = 4D0*PI/3D0 * Radii_Secondary ** 3D0 * rhoB
        C20unit = (Axis_c ** 2D0 - (Axis_a ** 2D0 + Axis_b ** 2D0) / 2D0) / 5D0
        C22unit = 1D0 / 2D1 * (Axis_a ** 2D0 - Axis_b ** 2D0)
        Omega_Orbit = -DSQRT(Gravity * (Mass_Primary + Mass_Secondary) / SemiAxis ** 3D0)
        
        ! Unit Parameters
        
        Tunit = 1D0 / Omega_Primary
        Munit = Mass_Primary
        Runit = (Gravity * Mass_Primary * Tunit ** 2D0) ** (1D0 / 3D0)
		
        ! Normalisation
        
        C20 = C20unit / Runit ** 2D0
        C22 = C22unit / Runit ** 2D0
        RadiiA = Radii_Primary / Runit
        RadiiB = Radii_Secondary / Runit
        OmegaN = Omega_Orbit * Tunit
        OmegaA = Omega_Primary * Tunit
        Ratio_Mass = Mass_Secondary / Munit
        OrbitN = SemiAxis / Runit
        R_Sun = rSun / Runit
        Beta_Sun = kappa * Area_Mass * PSRP * AU ** 2D0 * Tunit ** 2D0 / Runit ** 3D0
        
        WRITE(*,*) '==============================================================================='
        WRITE(*,*) 'Axis_a(1&m) =    ', Axis_a/Runit,Axis_a
        WRITE(*,*) 'Axis_b(1&m) =    ', Axis_b/Runit,Axis_b
        WRITE(*,*) 'Axis_c(1&m) =    ', Axis_c/Runit,Axis_c
        WRITE(*,*) 'TA(hour)    =    ', 2D0 * PI/Omega_Primary/36D2
        WRITE(*,*) 'Tn(hour)    =    ', 2D0 * PI/Omega_Orbit/36D2
        WRITE(*,*) 'r1(m)       =    ', OrbitN
        WRITE(*,*) 'rSun(AU)    =    ', rSun/AU
        WRITE(*,*) ' '
        WRITE(*,*) 'C20         =    ', C20
        WRITE(*,*) 'C22         =    ', C22
        WRITE(*,*) ' '
        WRITE(*,*) 'RadiiA      =    ', RadiiA
        WRITE(*,*) 'RadiiB      =    ', RadiiB
        WRITE(*,*) 'OmegaN      =    ', OmegaN
        WRITE(*,*) 'OmegaA      =    ', OmegaA
        WRITE(*,*) 'Ratio_Mass  =    ', Ratio_Mass
        WRITE(*,*) 'Beta_Sun    =    ', Beta_Sun
        WRITE(*,*) 'kappa_Sun   =    ', Beta_Sun / R_sun ** 2D0 * 1D3
        WRITE(*,*) ' '
        WRITE(*,*) 'Tunit(hour) =    ', Tunit/36D2
        WRITE(*,*) 'Runit(m)    =    ', Runit
        WRITE(*,*) 'Munit(kg)    =    ', Munit
        
    END SUBROUTINE GlobalPar
    
