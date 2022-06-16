    SUBROUTINE AnaSol_Boost(SolC,SolS,T,YHC,NN,ORDER)
    USE GlobalDef   ! Define glocal variables
    IMPLICIT REAL*8(A-H,O-Z)
    REAL(8) :: YHC(6),T
    REAL(8) :: phi, phiA, omegaSun, omega
    INTEGER(4) :: ORDER,NN
    REAL(8) :: SolC(0:NN,-NN:NN,1:3,1:3),SolS(0:NN,-NN:NN,1:3,1:3)
    
    omega = OmegaN - OmegaA
    phi =  omega * T + phi0
    omegaSun = - OmegaA
    phiA = omegaSun * T
    
    YHC = 0D0
    
    DO K = 1,3
        DO I = 0,NN
            DO J = -NN,NN
                DO N = 1,ORDER
                    IF (SolC(I,J,K,N).NE.0D0 .OR. SolS(I,J,K,N).NE.0D0) THEN
                        YHC(K) = YHC(K) + SolC(I,J,K,N)*DCOS(I*phi + J*phiA) + SolS(I,J,K,N)*DSIN(I*phi + J*phiA)
                        YHC(K+3) = YHC(K+3) + (-(I*omega+J*omegaSun)*SolC(I,J,K,N)*DSIN(I*phi + J*phiA) + (I*omega+J*omegaSun)*SolS(I,J,K,N)*DCOS(I*phi + J*phiA))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    
    END SUBROUTINE
    
    SUBROUTINE AnaSol(T,YHC)
    USE GlobalDef   ! Define glocal variables
    IMPLICIT REAL*8(A-H,O-Z)
    INTEGER(4),PARAMETER :: NN=15,ORDERP=5
    REAL(8) :: Sol1C(0:NN,-NN:NN,1:3),Sol1S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol2C(0:NN,-NN:NN,1:3),Sol2S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol3C(0:NN,-NN:NN,1:3),Sol3S(0:NN,-NN:NN,1:3)
    REAL(8) :: OmegaE2C(0:NN,-NN:NN,1:3),OmegaE2S(0:NN,-NN:NN,1:3)
    REAL(8) :: OmegaE3C(0:NN,-NN:NN,1:3),OmegaE3S(0:NN,-NN:NN,1:3)
    REAL(8) :: SolC(0:NN,-NN:NN,1:3,1:3),SolS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: YHC(6),T
    REAL(8) :: phi, phiA, omegaSun, omega
    
    omega = OmegaN - OmegaA
    phi =  omega * T + phi0
    omegaSun = - OmegaA
    phiA = omegaSun * T
    
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
    
    SolS(0,0,:,:) = 0D0
    
    YHC = 0D0
    
    DO K = 1,3
        DO I = 0,NN
            DO J = -NN,NN
                YHC(K) = YHC(K) + Sol1C(I,J,K)*DCOS(I*phi + J*phiA) + Sol1S(I,J,K)*DSIN(I*phi + J*phiA)
                YHC(K+3) = YHC(K+3) + (-(I*omega+J*omegaSun)*Sol1C(I,J,K)*DSIN(I*phi + J*phiA) + (I*omega+J*omegaSun)*Sol1S(I,J,K)*DCOS(I*phi + J*phiA))
                
                YHC(K) = YHC(K) + Sol2C(I,J,K)*DCOS(I*phi + J*phiA) + Sol2S(I,J,K)*DSIN(I*phi + J*phiA)
                YHC(K+3) = YHC(K+3) + (-(I*omega+J*omegaSun)*Sol2C(I,J,K)*DSIN(I*phi + J*phiA) + (I*omega+J*omegaSun)*Sol2S(I,J,K)*DCOS(I*phi + J*phiA))
                
                YHC(K) = YHC(K) + Sol3C(I,J,K)*DCOS(I*phi + J*phiA) + Sol3S(I,J,K)*DSIN(I*phi + J*phiA)
                YHC(K+3) = YHC(K+3) + (-(I*omega+J*omegaSun)*Sol3C(I,J,K)*DSIN(I*phi + J*phiA) + (I*omega+J*omegaSun)*Sol3S(I,J,K)*DCOS(I*phi + J*phiA))
           
            ENDDO
        ENDDO
    ENDDO
    
    END SUBROUTINE
    
    SUBROUTINE Sol_Coe(OmegaEC,OmegaES,SolC,SolS,NN)
    
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: NN
    REAL(8) :: SolC(0:NN,-NN:NN,1:3),SolS(0:NN,-NN:NN,1:3)
    REAL(8) :: OmegaEC(0:NN,-NN:NN,1:3),OmegaES(0:NN,-NN:NN,1:3)
    REAL(8) :: CoeMC(3,3),CoeMS(3,3),RfVec(2),ReVec(2)
    
    SolC = 0D0
    SolS = 0D0
    
    DO I = 0,NN,1
        DO J = -NN,NN,1
            ! x,y方向求解
            IF (OmegaEC(I,J,1)==0 .AND. OmegaES(I,J,1)==0 .AND. OmegaEC(I,J,2)==0 .AND. OmegaES(I,J,2)==0) THEN
                CONTINUE
            ELSE
                CALL Coeff_Matrix(CoeMC,CoeMS,I,J)
                
                RfVec(1) = OmegaEC(I,J,1)
                RfVec(2) = OmegaES(I,J,2)
                CALL SolFun(CoeMC(1:2,1:2),RfVec,ReVec,2)
                SolC(I,J,1) = ReVec(1)
                SolS(I,J,2) = ReVec(2)
                
                RfVec(1) = OmegaES(I,J,1)
                RfVec(2) = OmegaEC(I,J,2)
                CALL SolFun(CoeMS(1:2,1:2),RfVec,ReVec,2)
                SolS(I,J,1) = ReVec(1)
                SolC(I,J,2) = ReVec(2)
            ENDIF
            ! z方向求解
            IF (OmegaEC(I,J,3)==0 .AND. OmegaES(I,J,3)==0) THEN
                CONTINUE
            ELSE
                CALL Coeff_Matrix(CoeMC,CoeMS,I,J)
                SolC(I,J,3) = OmegaEC(I,J,3)/CoeMC(3,3)
                SolS(I,J,3) = OmegaES(I,J,3)/CoeMS(3,3)
            ENDIF
        ENDDO
    ENDDO
    
    END SUBROUTINE
    
    SUBROUTINE OmegaE2(OmegaE2C,OmegaE2S,NN,ORDERP)
    
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: NN,ORDERP
    REAL(8) :: OmegaE2C(0:NN,-NN:NN,1:3),OmegaE2S(0:NN,-NN:NN,1:3)
    REAL(8) :: F0xi1C(0:NN,-NN:NN,1:3),F0xi1S(0:NN,-NN:NN,1:3)
    REAL(8) :: Fexi1C(0:NN,-NN:NN,1:3),Fexi1S(0:NN,-NN:NN,1:3)
	REAL(8) :: CC(0:NN,-NN:NN),CS(0:NN,-NN:NN)
    DIMENSION IDA(2),IDB(2),IDC(2)
    REAL(8) :: FeXC(0:NN,-NN:NN,1:3,1:3),FeXS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: Sol1C(0:NN,-NN:NN,1:3),Sol1S(0:NN,-NN:NN,1:3)
    REAL(8) :: F0X2(1:3,1:3,1:3)
    
    CALL Sol1_Coe(Sol1C,Sol1S,NN,ORDERP)
    !CALL WriteToScreen(Sol1C,Sol1S,NN)
    CALL FeX(FeXC,FeXS,NN,ORDERP)
    CALL F0X2Matri(F0X2)
    
    IDA=NN
	IDB=NN
	IDC=NN
    
    F0xi1C = 0D0
    F0xi1S = 0D0
    Fexi1C = 0D0
    Fexi1S = 0D0
    DO K = 1,3
        DO I = 1,3
            DO J = 1,3
                CALL MULTIPLY(Sol1C(:,:,I),Sol1S(:,:,I),Sol1C(:,:,J),Sol1S(:,:,J),CC,CS,IDA,IDB,IDC)
                F0xi1C(:,:,K) = F0xi1C(:,:,K) + CC*F0X2(K,I,J)/2D0
                F0xi1S(:,:,K) = F0xi1S(:,:,K) + CS*F0X2(K,I,J)/2D0
                !WRITE(*,*) 'lalalal'
            END DO
            CALL MULTIPLY(FeXC(:,:,K,I),FeXS(:,:,K,I),Sol1C(:,:,I),Sol1S(:,:,I),CC,CS,IDA,IDB,IDC)
            Fexi1C(:,:,K) = Fexi1C(:,:,K)+CC
            Fexi1S(:,:,K) = Fexi1S(:,:,K)+CS
            !WRITE(*,*) 'lalalal'
        END DO
    END DO
    
    OmegaE2C = F0xi1C + Fexi1C
    OmegaE2S = F0xi1S + Fexi1S
    
    END SUBROUTINE
    
    SUBROUTINE OmegaE3(OmegaE3C,OmegaE3S,NN,ORDERP)
    
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: NN,ORDERP
    DIMENSION IDA(2),IDB(2),IDC(2)
    REAL(8) :: OmegaE3C(0:NN,-NN:NN,1:3),OmegaE3S(0:NN,-NN:NN,1:3)
    REAL(8) :: OmegaE2C(0:NN,-NN:NN,1:3),OmegaE2S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol1C(0:NN,-NN:NN,1:3),Sol1S(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol2C(0:NN,-NN:NN,1:3),Sol2S(0:NN,-NN:NN,1:3)
    REAL(8) :: F0X2(1:3,1:3,1:3),F0X3(1:3,1:3,1:3,1:3)
    REAL(8) :: FeXC(0:NN,-NN:NN,1:3,1:3),FeXS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: FeXXC(0:NN,-NN:NN,1:3,1:3,1:3),FeXXS(0:NN,-NN:NN,1:3,1:3,1:3)
    
    REAL(8) :: F0XXXSol13SumC(0:NN,-NN:NN,1:3),F0XXXSol13SumS(0:NN,-NN:NN,1:3)
    REAL(8) :: FeXXSol12SumC(0:NN,-NN:NN,1:3),FeXXSol12SumS(0:NN,-NN:NN,1:3)
    REAL(8) :: F0XXSol12SumC(0:NN,-NN:NN,1:3),F0XXSol12SumS(0:NN,-NN:NN,1:3)
    REAL(8) :: FeXSol2SumC(0:NN,-NN:NN,1:3),FeXSol2SumS(0:NN,-NN:NN,1:3)
    REAL(8) :: Sol12C(0:NN,-NN:NN),Sol12S(0:NN,-NN:NN)
    REAL(8) :: FeXXSol12C(0:NN,-NN:NN),FeXXSol12S(0:NN,-NN:NN)
    REAL(8) :: FeXSol2C(0:NN,-NN:NN),FeXSol2S(0:NN,-NN:NN)
    REAL(8) :: Sol13C(0:NN,-NN:NN),Sol13S(0:NN,-NN:NN)
    
    
    ! 计算一阶解
    CALL Sol1_Coe(Sol1C,Sol1S,NN,ORDERP)
    ! 计算二阶解
    CALL OmegaE2(OmegaE2C,OmegaE2S,NN,ORDERP)
    CALL Sol_Coe(OmegaE2C,OmegaE2S,Sol2C,Sol2S,NN)
    ! 计算偏导系数
    CALL FeX(FeXC,FeXS,NN,ORDERP)
    CALL FeXX(FeXXC,FeXXS,NN,ORDERP)
    CALL F0X2Matri(F0X2)
    CALL F0X3Matri(F0X3)
    
    IDA=NN
	IDB=NN
	IDC=NN
    F0XXXSol13SumC = 0D0
    F0XXXSol13SumS = 0D0
    FeXXSol12SumC = 0D0
    FeXXSol12SumS = 0D0
    F0XXSol12SumC = 0D0
    F0XXSol12SumS = 0D0
    FeXSol2SumC = 0D0
    FeXSol2SumS = 0D0
    OmegaE3C = 0D0
    OmegaE3S = 0D0
    
    DO L = 1,3
        DO I = 1,3
            DO J = 1,3
                DO K = 1,3
                    CALL MULTIPLY(Sol1C(:,:,I),Sol1S(:,:,I),Sol1C(:,:,J),Sol1S(:,:,J),Sol12C,Sol12S,IDA,IDB,IDC)
                    CALL MULTIPLY(Sol12C,Sol12S,Sol1C(:,:,K),Sol1S(:,:,K),Sol13C,Sol13S,IDA,IDB,IDC)
                    F0XXXSol13SumC(:,:,L) = F0XXXSol13SumC(:,:,L) + F0X3(L,I,J,K)*Sol13C/6D0
                    F0XXXSol13SumS(:,:,L) = F0XXXSol13SumS(:,:,L) + F0X3(L,I,J,K)*Sol13S/6D0
                    
                ENDDO
                CALL MULTIPLY(Sol1C(:,:,I),Sol1S(:,:,I),Sol1C(:,:,J),Sol1S(:,:,J),Sol12C,Sol12S,IDA,IDB,IDC)
                CALL MULTIPLY(FeXXC(:,:,L,I,J),FeXXS(:,:,L,I,J),Sol12C,Sol12S,FeXXSol12C,FeXXSol12S,IDA,IDB,IDC)
                FeXXSol12SumC(:,:,L) = FeXXSol12SumC(:,:,L) + FeXXSol12C/2D0
                FeXXSol12SumS(:,:,L) = FeXXSol12SumS(:,:,L) + FeXXSol12S/2D0
                
                CALL MULTIPLY(Sol1C(:,:,I),Sol1S(:,:,I),Sol2C(:,:,J),Sol2S(:,:,J),Sol12C,Sol12S,IDA,IDB,IDC)
                F0XXSol12SumC(:,:,L) = F0XXSol12SumC(:,:,L) + F0X2(L,I,J)*Sol12C
                F0XXSol12SumS(:,:,L) = F0XXSol12SumS(:,:,L) + F0X2(L,I,J)*Sol12S
            ENDDO
            CALL MULTIPLY(FeXC(:,:,L,I),FeXS(:,:,L,I),Sol2C(:,:,I),Sol2S(:,:,I),FeXSol2C,FeXSol2S,IDA,IDB,IDC)
            FeXSol2SumC(:,:,L) = FeXSol2SumC(:,:,L) + FeXSol2C
            FeXSol2SumS(:,:,L) = FeXSol2SumS(:,:,L) + FeXSol2S
        ENDDO
    ENDDO
    
    OmegaE3C = F0XXSol12SumC + F0XXXSol13SumC + FeXSol2SumC + FeXXSol12SumC
    OmegaE3S = F0XXSol12SumS + F0XXXSol13SumS + FeXSol2SumS + FeXXSol12SumS
    
    END SUBROUTINE
    
    SUBROUTINE Coeff_Matrix(CoeMC,CoeMS,L,M)
    
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: L,M
    REAL(8) :: CoeMC(3,3),CoeMS(3,3)
    REAL(8) :: Omega0xx, Omega0yy, Omega0zz,omega,omegaSun
    
    ! Parameters
    omega = OmegaN - OmegaA
    omegaSun = - OmegaA
    ! Difference of poteintal function
    Omega0xx = 1D0 - 1D0 / rS ** 3D0 + 3D0 / 2D0 * (C20 + 14D0 * C22) / rS ** 5D0
    Omega0yy = 1D0 + 2D0 / rS ** 3D0 - 6D0 * (C20 + 6D0 * C22) / rS ** 5D0
    Omega0zz = - 1D0 / rS ** 3D0 + 3D0 / 2D0 * (3D0 * C20 + 10D0 * C22) / rS ** 5D0
    
    CoeMC(1,1) = -(L*omega+M*omegaSun)**2D0 - Omega0xx
    CoeMC(1,2) = -2D0*(L*omega+M*omegaSun)
    CoeMC(2,1) = -2D0*(L*omega+M*omegaSun)
    CoeMC(2,2) = -(L*omega+M*omegaSun)**2D0 - Omega0yy
    
    CoeMS = CoeMC
    CoeMS(1,2) = -CoeMC(1,2)
    CoeMS(2,1) = -CoeMC(2,1)
    
    CoeMS(3,3) = -(L*omega+M*omegaSun)**2D0 - Omega0zz
    CoeMC(3,3) = CoeMS(3,3)
    
    END SUBROUTINE
    
    SUBROUTINE F0X2Matri(F0X2)
    
    USE GlobalDef   ! Define glocal variables
    REAL(8) :: F0X2(1:3,1:3,1:3)
    
    F0X2 = 0D0
    
    F0X2(1,1,2) = -(3D0/2D0) * (-2D0 * rS ** 2D0 + 5D0*C20 + 7D1*C22)/rS ** 6D0
    F0X2(1,2,1) = -(3D0/2D0) * (-2D0 * rS ** 2D0 + 5D0*C20 + 7D1*C22)/rS ** 6D0
    F0X2(2,1,1) = -(3D0/2D0) * (-2D0 * rS ** 2D0 + 5D0*C20 + 7D1*C22)/rS ** 6D0
    F0X2(2,2,2) = (6D0 * (-rS ** 2D0 + 5D0*C20 + 3D1 * C22))/rS ** 6D0
    F0X2(2,3,3)= -(3D0/2D0) * (-2D0 *rS ** 2D0 + 15D0*C20 + 5D1*C22)/rS ** 6D0
    F0X2(3,2,3)= -(3D0/2D0) * (-2D0 *rS ** 2D0 + 15D0*C20 + 5D1*C22)/rS ** 6D0
    F0X2(3,3,2)= -(3D0/2D0) * (-2D0 *rS ** 2D0 + 15D0*C20 + 5D1*C22)/rS ** 6D0
    
    END SUBROUTINE
    
    SUBROUTINE F0X3Matri(F0X3)
    
    USE GlobalDef   ! Define glocal variables
    REAL(8) :: F0X3(1:3,1:3,1:3,1:3)
    
    F0X3 = 0D0
    
    F0X3(1,1,1,1) = -(9D0/2D0)*(-2*rS**2+5*C20+110*C22)/rS**7
    F0X3(1,1,2,2) = (3D0*(-4D0*rS**2+15*C20+210*C22))/rS**7
    F0X3(1,2,1,2) = (3D0*(-4D0*rS**2+15*C20+210*C22))/rS**7
    F0X3(1,2,2,1) = (3D0*(-4D0*rS**2+15*C20+210*C22))/rS**7
    F0X3(1,3,3,1) = -(3D0/2D0)*(-2*rS**2+15*C20+90*C22)/rS**7
    F0X3(1,1,3,3) = -(3D0/2D0)*(-2*rS**2+15*C20+90*C22)/rS**7
    F0X3(1,3,1,3) = -(3D0/2D0)*(-2*rS**2+15*C20+90*C22)/rS**7
    F0X3(2,2,2,2) = -(12D0*(-2*rS**2+15*C20+90*C22))/rS**7
    F0X3(2,1,1,2) = (3D0*(-4D0*rS**2+15*C20+210*C22))/rS**7
    F0X3(2,1,2,1) = (3D0*(-4D0*rS**2+15*C20+210*C22))/rS**7
    F0X3(2,2,1,1) = (3D0*(-4D0*rS**2+15*C20+210*C22))/rS**7
    F0X3(2,3,3,2) = (3D0*(-4*rS**2+45*C20+150*C22))/rS**7
    F0X3(2,2,3,3) = (3D0*(-4*rS**2+45*C20+150*C22))/rS**7
    F0X3(2,3,2,3) = (3D0*(-4*rS**2+45*C20+150*C22))/rS**7
    F0X3(3,3,3,3) = -(9D0/2D0)*(-2*rS**2+25*C20+70*C22)/rS**7
    F0X3(3,1,1,3) = -(3D0/2D0)*(-2*rS**2+15*C20+90*C22)/rS**7
    F0X3(3,3,1,1) = -(3D0/2D0)*(-2*rS**2+15*C20+90*C22)/rS**7
    F0X3(3,1,3,1) = -(3D0/2D0)*(-2*rS**2+15*C20+90*C22)/rS**7
    F0X3(3,2,2,3) = (3D0*(-4*rS**2+45*C20+150*C22))/rS**7
    F0X3(3,2,3,2) = (3D0*(-4*rS**2+45*C20+150*C22))/rS**7
    F0X3(3,3,2,2) = (3D0*(-4*rS**2+45*C20+150*C22))/rS**7
    
    END SUBROUTINE
    
    SUBROUTINE FeX(FeXC,FeXS,NN,ORDERP)
    !变量说明
    !FeXC(l,m,i,j),FeXS(l,m,i,j)
    !   l: 三角函数第一个幅角阶数
    !   m: 三角函数第二个幅角阶数
    !   i: 摄动力在i方向上的分量,i={x,y,z},分别对应方程组中第1,2,3个方程
    !   j: 对第j个分量求导,j={x,y,z}
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: NN,ORDERP
    REAL(8) :: FeXC(0:NN,-NN:NN,1:3,1:3),FeXS(0:NN,-NN:NN,1:3,1:3)
    REAL(8) :: RM, SF, R1
    
    SF = ShadeFactor
    RM = Ratio_Mass
    R1 = OrbitN
    
    FeXC = 0D0
    FeXS = 0D0
    
    ! P2+P3
    IF (ORDERP>=3)THEN
        FeXS(1,0,1,1) = (3D0/4D0)*RM*rS/R1**4D0
        FeXC(2,0,1,1) = 3D0/2D0*RM/R1**3D0
        FeXS(3,0,1,1) = (15D0/4D0)*RM*rS/R1**4D0
        FeXC(0,0,1,1) = 1D0/2D0*RM/R1**3D0
        FeXC(1,0,1,2) = (3D0/4D0)*RM*rS/R1**4D0
        FeXS(2,0,1,2) = (3D0/2D0)*RM/R1**3D0
        FeXC(3,0,1,2) = -(15D0/4D0)*RM*rS/R1**4D0
    
        FeXC(3,0,2,1) = -(15D0/4D0)*RM*rS/R1**4D0
        FeXC(1,0,2,1) = (3D0/4D0)*RM*rS/R1**4D0
        FeXS(2,0,2,1) = (3D0/2D0)*RM/r1**3D0
        FeXC(0,0,2,2) = 1D0/2D0*RM/R1**3D0
        FeXC(2,0,2,2) = -(3D0/2D0)*RM/R1**3D0
        FeXS(1,0,2,2) = (9D0/4D0)*RM*rS/R1**4D0
        FeXS(3,0,2,2) = -(15D0/4D0)*RM*rS/R1**4D0
    
        FeXS(1,0,3,3) = -3D0*RM*rS/R1**4D0
        FeXC(0,0,3,3) = -RM/R1**3D0
    ENDIF
    
    !P4
    IF(ORDERP>=4)THEN
        FeXC(4,0,1,1) = FeXC(4,0,1,1) - (105D0/16D0)*RM*rS**2D0/R1**5D0
        FeXC(0,0,1,1) = FeXC(0,0,1,1) + (9D0/16D0)*RM*rS**2D0/R1**5D0
        FeXS(4,0,1,2) = FeXS(4,0,1,2) - (105D0/16D0)*RM*rS**2D0/R1**5D0
        FeXS(2,0,1,2) = FeXS(2,0,1,2) + (15D0/8D0)*RM*rS**2D0/R1**5D0
    
        FeXS(2,0,2,1) = FeXS(2,0,2,1) + (15D0/8D0)*RM*rS**2D0/R1**5D0
        FeXS(4,0,2,1) = FeXS(4,0,2,1) - (105D0/16D0)*RM*rS**2D0/R1**5D0
        FeXC(0,0,2,2) = FeXC(0,0,2,2) + 27D0/16D0*RM*rS**2D0/R1**5D0
        FeXC(2,0,2,2) = FeXC(2,0,2,2) - (15D0/4D0)*RM*rS**2D0/R1**5D0
        FeXC(4,0,2,2) = FeXC(4,0,2,2) + (105D0/16D0)*RM*rS**2D0/R1**5D0
    
        FeXC(2,0,3,3) = FeXC(2,0,3,3) + (15D0/4D0)*RM*rS**2D0/R1**5D0
        FeXC(0,0,3,3) = FeXC(0,0,3,3) - 9D0/4D0*RM*rS**2D0/R1**5D0
    ENDIF
    
    !P5
    IF (ORDERP>=5)THEN
        FeXS(1,0,1,1) = FeXS(1,0,1,1) + (15D0/16D0)*RM*rS**3D0/R1**6D0
        FeXS(3,0,1,1) = FeXS(3,0,1,1) + (35D0/32D0)*RM*rS**3D0/R1**6D0
        FeXS(5,0,1,1) = FeXS(5,0,1,1) - (315D0/32D0)*RM*rS**3D0/R1**6D0
        FeXC(1,0,1,2) = FeXC(1,0,1,2) + (15D0/16D0)*RM*rS**3D0/R1**6D0
        FeXC(3,0,1,2) = FeXC(3,0,1,2) - (105D0/32D0)*RM*rS**3D0/R1**6D0
        FeXC(5,0,1,2) = FeXC(5,0,1,2) + (315D0/32D0)*RM*rS**3D0/R1**6D0
    
        FeXC(1,0,2,1) = FeXC(1,0,2,1) + (15D0/16D0)*RM*rS**3D0/R1**6D0
        FeXC(3,0,2,1) = FeXC(3,0,2,1) - (105D0/32D0)*RM*rS**3D0/R1**6D0
        FeXC(5,0,2,1) = FeXC(5,0,2,1) + (315D0/32D0)*RM*rS**3D0/R1**6D0
        FeXS(1,0,2,2) = FeXS(1,0,2,2) + (75D0/16D0)*RM*rS**3D0/R1**6D0
        FeXS(3,0,2,2) = FeXS(3,0,2,2) - (175D0/32D0)*RM*rS**3D0/R1**6D0
        FeXS(5,0,2,2) = FeXS(5,0,2,2) + (315D0/32D0)*RM*rS**3D0/R1**6D0
    
        FeXS(1,0,3,3) = FeXS(1,0,3,3) - (45D0/8D0)*RM*rS**3D0/R1**6D0
        FeXS(3,0,3,3) = FeXS(3,0,3,3) + (35D0/8D0)*RM*rS**3D0/R1**6D0
    ENDIF
    
    END
    
    SUBROUTINE FeXX(FeXXC,FeXXS,NN,ORDERP)
    !变量说明
    !FeXC(l,m,i,j,k),FeXS(l,m,i,j,k)
    !   l: 三角函数第一个幅角阶数
    !   m: 三角函数第二个幅角阶数
    !   i: 摄动力在i方向上的分量,i={x,y,z},分别对应方程组中第1,2,3个方程
    !   j: 对第j个分量求导,j={x,y,z}
    !   k: 对第k个分量求导,k={x,y,z}
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: NN,ORDERP
    REAL(8) :: FeXXC(0:NN,-NN:NN,1:3,1:3,1:3),FeXXS(0:NN,-NN:NN,1:3,1:3,1:3)
    REAL(8) :: RM, SF, R1
    
    SF = ShadeFactor
    RM = Ratio_Mass
    R1 = OrbitN
    
    FeXXC = 0D0
    FeXXS = 0D0
    
    ! P2+P3
    IF (ORDERP>=3)THEN
        FeXXC(1,0,1,1,1) = 9D0/4D0*RM/R1**4D0
        FeXXC(3,0,1,1,1) = (15D0/4D0)*RM/R1**4D0
        FeXXS(1,0,1,1,2) = 3D0/4D0*RM/R1**4D0
        FeXXS(3,0,1,1,2) = (15D0/4D0)*RM/R1**4D0
        FeXXS(1,0,1,2,1) = 3D0/4D0*RM/R1**4D0
        FeXXS(3,0,1,2,1) = (15D0/4D0)*RM/R1**4D0
        FeXXC(1,0,1,2,2) = 3D0/4D0*RM/R1**4D0
        FeXXC(3,0,1,2,2) = -(15D0/4D0)*RM/R1**4D0
        FeXXC(1,0,1,3,3) = -3D0*RM/R1**4D0

        FeXXS(1,0,2,1,1) = 3D0/4D0*RM/R1**4D0
        FeXXS(3,0,2,1,1) = (15D0/4D0)*RM/R1**4D0
        FeXXC(1,0,2,1,2) = 3D0/4D0*RM/R1**4D0
        FeXXC(3,0,2,1,2) = -(15D0/4D0)*RM/R1**4D0
        FeXXC(1,0,2,2,1) = 3D0/4D0*RM/R1**4D0
        FeXXC(3,0,2,2,1) = -(15D0/4D0)*RM/R1**4D0
        FeXXS(1,0,2,2,2) = 9D0/4D0*RM/R1**4D0
        FeXXS(3,0,2,2,2) = -(15D0/4D0)*RM/R1**4D0
        FeXXS(1,0,2,3,3) = -3D0*RM/R1**4D0

        FeXXC(1,0,3,1,3) = -3D0*RM/R1**4D0
        FeXXC(1,0,3,3,1) = -3D0*RM/R1**4D0
        FeXXS(1,0,3,2,3) = -3D0*RM/R1**4D0
        FeXXS(1,0,3,3,2) = -3D0*RM/R1**4D0
    ENDIF
    
    ! P4
    IF (ORDERP>=4)THEN
        FeXXS(2,0,1,1,1) = FeXXS(2,0,1,1,1) + (15D0/4D0)*RM*rS/R1**5D0
        FeXXS(4,0,1,1,1) = FeXXS(4,0,1,1,1) + (105D0/8D0)*RM*rS/R1**5D0
        FeXXC(2,0,1,1,2) = FeXXC(2,0,1,1,2) - (105D0/8D0)*RM*rS/R1**5D0 
        FeXXC(0,0,1,1,2) = FeXXC(0,0,1,1,2) - (9D0/8D0)*RM*rS/R1**5D0
        FeXXS(2,0,1,2,2) = FeXXS(2,0,1,2,2) + (15D0/4D0)*RM*rS/R1**5D0 
        FeXXS(4,0,1,2,2) = FeXXS(4,0,1,2,2) - (105D0/8D0)*RM*rS/R1**5D0
        FeXXS(2,0,1,3,3) = FeXXS(2,0,1,3,3) - (15D0/2D0)*RM*rS/R1**5D0
    
        FeXXC(4,0,2,1,1) = FeXXC(4,0,2,1,1) - (105D0/8D0)*RM*rS/R1**5D0
        FeXXC(0,0,2,1,1) = FeXXC(0,0,2,1,1) + (9D0/8D0)*RM*rS/R1**5D0
        FeXXS(4,0,2,1,2) = FeXXS(4,0,2,1,2) - (105D0/8D0)*RM*rS/R1**5D0
        FeXXS(2,0,2,1,2) = FeXXS(2,0,2,1,2) + (15D0/4D0)*RM*rS/R1**5D0
        FeXXC(2,0,2,2,2) = FeXXC(2,0,2,2,2) - (15D0/2D0)*RM*rS/R1**5D0
        FeXXC(4,0,2,2,2) = FeXXC(4,0,2,2,2) + (105D0/8D0)*RM*rS/R1**5D0
        FeXXC(0,0,2,2,2) = FeXXC(0,0,2,2,2) + (27D0/8D0)*RM*rS/R1**5D0
    
        FeXXC(0,0,2,3,3) = FeXXC(0,0,2,3,3) - (9D0/2D0)*RM*rS/R1**5D0
        FeXXC(2,0,2,3,3) = FeXXC(2,0,2,3,3) + (15D0/2D0)*RM*rS/R1**5D0
    
        FeXXS(2,0,3,1,3) = FeXXS(2,0,3,1,3) - (15D0/2D0)*RM*rS/R1**5D0
        FeXXC(2,0,3,2,3) = FeXXC(2,0,3,2,3) + (15D0/2D0)*RM*rS/R1**5D0
        FeXXC(0,0,3,2,3) = FeXXC(0,0,3,2,3) - 9D0/2D0*RM*rS/R1**5D0
    ENDIF
    
    !P5
    IF (ORDERP>=5)THEN
        FeXXC(1,0,1,1,1) = FeXXC(1,0,1,1,1) + (45D0/16D0)*RM*rS**2D0/R1**6D0
        FeXXC(3,0,1,1,1) = FeXXC(3,0,1,1,1) - (105D0/32D0)*RM*rS**2D0/R1**6D0
        FeXXC(5,0,1,1,1) = FeXXC(5,0,1,1,1) - (945D0/32D0)*RM*rS**2D0/R1**6D0
        FeXXS(1,0,1,1,2) = FeXXS(1,0,1,1,2) + (45D0/16D0)*RM*rS**2D0/R1**6D0
        FeXXS(3,0,1,1,2) = FeXXS(3,0,1,1,2) + (105D0/32D0)*RM*rS**2D0/R1**6D0 
        FeXXS(5,0,1,1,2) = FeXXS(5,0,1,1,2) - (945D0/32D0)*RM*rS**2D0/R1**6D0
        FeXXC(1,0,1,3,3) = FeXXC(1,0,1,3,3) - (45D0/8D0)*RM*rS**2D0/R1**6D0
        FeXXC(3,0,1,3,3) = FeXXC(3,0,1,3,3) + (105D0/8D0)*RM*rS**2D0/R1**6D0
    
        FeXXS(1,0,2,1,1) = FeXXS(1,0,2,1,1) + (45D0/16D0)*RM*rS**2D0/R1**6D0
        FeXXS(3,0,2,1,1) = FeXXS(3,0,2,1,1) + (105D0/32D0)*RM*rS**2D0/R1**6D0 
        FeXXS(5,0,2,1,1) = FeXXS(5,0,2,1,1) - (945D0/32D0)*RM*rS**2D0/R1**6D0
        FeXXC(1,0,2,1,2) = FeXXC(1,0,2,1,2) + (45D0/16D0)*RM*rS**2D0/R1**6D0
        FeXXC(3,0,2,1,2) = FeXXC(3,0,2,1,2) - (315D0/32D0)*RM*rS**2D0/R1**6D0
        FeXXC(5,0,2,1,2) = FeXXC(5,0,2,1,2) + (945D0/32D0)*RM*rS**2D0/R1**6D0
        FeXXS(1,0,2,3,3) = FeXXS(1,0,2,3,3) - (135D0/8D0)*RM*rS**2D0/R1**6D0
        FeXXS(3,0,2,3,3) = FeXXS(3,0,2,3,3) + (105D0/8D0)*RM*rS**2D0/R1**6D0
    
        FeXXC(1,0,3,1,3) = FeXXC(1,0,3,1,3) - (45D0/8D0)*RM*rS**2D0/R1**6D0
        FeXXC(3,0,3,1,3) = FeXXC(3,0,3,1,3) + (105D0/8D0)*RM*rS**2D0/R1**6D0
        FeXXS(1,0,3,2,3) = FeXXS(1,0,3,2,3) - (135D0/8D0)*RM*rS**2D0/R1**6D0
        FeXXS(3,0,3,2,3) = FeXXS(3,0,3,2,3) + (105D0/8D0)*RM*rS**2D0/R1**6D0
    ENDIF
    
    END
    
    SUBROUTINE Sol1_Coe(Sol1C,Sol1S,NN,ORDERP)
        
    USE GlobalDef   ! Define glocal variables
    INTEGER(4) :: NN,ORDERP
    REAL(8) :: Sol1C(0:NN,-NN:NN,1:3),Sol1S(0:NN,-NN:NN,1:3)
    REAL(8) :: R1, Omega0xx, Omega0yy, Omega0zz, Omega0yyy, omegaSun, omega, Omegaez
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
    omegaSun = - OmegaA
    
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
    
    ! First-oder coeffcients
    ! P3
    Sol1C = 0D0
    Sol1S = 0D0
    
    IF (ORDERP>=3)THEN
        Sol1C(1,0,1) = Alpha1xi_P3
        Sol1C(3,0,1) = Alpha3xi_P3
        Sol1C(0,1,1) = AlphaSxi
        Sol1S(2,0,1) = Beta2xi_P2
    
        Sol1S(1,0,2) = Beta1eta_P3
        Sol1S(3,0,2) = Beta3eta_P3
        Sol1S(0,1,2) = BetaSeta
        Sol1C(2,0,2) = Alpha2eta_P2
        Sol1C(0,0,2) = Ceta_P2
    
        Sol1C(0,0,3) = Czeta
    ENDIF
    
    ! P4
    IF (ORDERP>=4)THEN
        Sol1S(2,0,1) = Sol1S(2,0,1) + Beta2xi_P4
        Sol1S(4,0,1) = Beta4xi_P4
    
        Sol1C(2,0,2) = Sol1C(2,0,2) + Alpha2eta_P4
        Sol1C(4,0,2) = Alpha4eta_P4
        Sol1C(0,0,2) = Sol1C(0,0,2) + Ceta_P4
    ENDIF
    
    ! P5
    IF (ORDERP>=5)THEN
        Sol1C(1,0,1) = Sol1C(1,0,1) + Alpha1xi_P5
        Sol1C(3,0,1) = Sol1C(3,0,1) + Alpha3xi_P5
        Sol1C(5,0,1) = Alpha5xi_P5
    
        Sol1S(1,0,2) = Sol1S(1,0,2) + Beta1eta_P5
        Sol1S(3,0,2) = Sol1S(3,0,2) + Beta3eta_P5
        Sol1S(5,0,2) = Beta5eta_P5
    ENDIF
    
    END
    
    SUBROUTINE WriteToScreen(SolC,SolS,NN)
    
    IMPLICIT REAL*8(A-H,O-Z)
    INTEGER(4) :: NN
    REAL(8) :: SolC(0:NN,-NN:NN,1:3),SolS(0:NN,-NN:NN,1:3)
    
    WRITE(*,*) 'x axis'
    WRITE(*,*) '            COS terms                       SIN terms'
    DO I = 0,NN
        DO J = -NN,NN
            IF (SolC(I,J,1).NE.0D0 .OR. SolS(I,J,1).NE.0D0) THEN
                write(*,'(A,I2,A,I2,A,E25.18,A,1X,A,I2,A,I2,A,E25.18,A)') '[',I,',',J,',',SolC(I,J,1),']','[',I,',',J,',',SolS(I,J,1),']'
            ENDIF
        ENDDO
    ENDDO
    WRITE(*,*) 'y axis'
    WRITE(*,*) '            COS terms                       SIN terms'
    DO I = 0,NN
        DO J = -NN,NN
            IF (SolC(I,J,2).NE.0D0 .OR. SolS(I,J,2).NE.0D0) THEN
                write(*,'(A,I2,A,I2,A,E25.18,A,1X,A,I2,A,I2,A,E25.18,A)') '[',I,',',J,',',SolC(I,J,2),']','[',I,',',J,',',SolS(I,J,2),']'
            ENDIF
        ENDDO
    ENDDO
    WRITE(*,*) 'z axis'
    WRITE(*,*) '            COS terms                       SIN terms'
    DO I = 0,NN
        DO J = -NN,NN
            IF (SolC(I,J,3).NE.0D0 .OR. SolS(I,J,3).NE.0D0) THEN
                write(*,'(A,I2,A,I2,A,E25.18,A,1X,A,I2,A,I2,A,E25.18,A)') '[',I,',',J,',',SolC(I,J,3),']','[',I,',',J,',',SolS(I,J,3),']'
            ENDIF
        ENDDO
    ENDDO
    
    END SUBROUTINE
    
    