    SUBROUTINE PeriodOrbit_Tn(X0,Tp,TOF, Xp,flag)
    
    REAL(8) :: X0(6),Tp,TOF
    Integer(4) :: I,flag,Num
    REAL(8) :: Meye(6,6), Xp(6), XT(6), MM(6,6)
    REAL(8) :: Coe(2,2), Anti_Coe(2,2), dX(2),dF(2)
    REAL(8) :: dX2,dF2, TOF0,be
    
    ! 初始化Monodromy矩阵，单位矩阵
    Meye = 0D0
    DO I = 1, 6
        Meye(I, I) = 1D0
    ENDDO
    
    Xp = X0
    be = 1D0
    Num = 1000
    
    DO I = 1,Num
        
        CALL MonMatrix(XP,Tp,XT,MM)
        
        Coe(1,1) = MM(1,2) - Meye(1,2)
        Coe(1,2) = MM(1,4) - Meye(1,4)
        Coe(2,1) = MM(5,2) - Meye(5,2)
        Coe(2,2) = MM(5,4) - Meye(5,4)
        dF(1) = XT(1) - Xp(1)
        dF(2) = XT(5) - Xp(5)
        CALL ANTIMATRIX(Coe,Anti_Coe,2,flag)
        CALL MATRIXVECTOR(Anti_Coe,dF,dX,2,2)
        CALL VECTOR(dF,dF,dF2,2)
        TOF0 = DSQRT(dF2)
        
        Xp(2) = Xp(2) - dX(1)/(1D0+be*TOF0)
        Xp(4) = Xp(4) - dX(2)/(1D0+be*TOF0)
        
        !WRITE(*,*) I,TOF0
        IF(TOF0>1D3)THEN
            write(*,*) 'Error in Sub PeriodOrbitT_Tn: the correction is not convergent.'
            flag = 0
            EXIT
        ELSEIF (TOF0<TOF) THEN
            WRITE(*,*) 'Success !'
            flag = 1
            EXIT
        ENDIF
        
    ENDDO
    
    IF (I.GT.Num) THEN
        WRITE(*,*) 'Error in Sub PeriodOrbitT_Tn: Exceeds the maximum number of cycles'
        flag = 0
    ENDIF
    
    
    END SUBROUTINE
    
    SUBROUTINE MonMatrix(X0,tf,XT,MM)
    EXTERNAL :: YHCM
    REAL(8) :: M0(42),X0(6),XT(6),t0,tf,hh,err,EE,MM(6,6),Meye(6,6)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18)
    N = 6
    
    DO I = 1,N
        DO J = 1,N
            IF(I == J)THEN
                M0((I-1)*N+J) = 1.0D0
            ELSE
                M0((I-1)*N+J) = 0.0D0
            END IF
        END DO
    END DO
    
    M0(N*N+1:N*N+N) = X0
    
    !积分
    !OPEN(14,file='ORBIT.DAT', status='REPLACE')
    err = 1D-15
    hh = 1D-2
    t0 = 0D0
    EE = 1D-15
    DO WHILE (t0<tf)
        Call RKF78(YHCM,hh,t0,M0,EE,err,N*N+N)
        !write(14,107) M0(N*N+1:N*N+N),t0*Tunit
    ENDDO
    Call RKF78(YHCM,tf-t0,t0,M0,EE,1D0,N*N+N) ! 积分积满，恰好积分到tf
    !write(14,107) M0(N*N+1:N*N+N),t0*Tunit
    !CLOSE(14)
    

    DO I=1,N*N
	    J=INT((I+N-1)/N)
	    MM(J,I-N*(J-1)) = M0(I)
    ENDDO
    
    MM = MM !- Meye
    XT = M0(N*N+1:N*N+N)
    
    END SUBROUTINE
    
    SUBROUTINE MonEig(Xp,Tp,MM,WR,WI,DIS,flag)
    
    IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION WR(6),WI(6)
    REAL(8) :: Xp(6),Tp,MM(6,6),XT(6)
    Integer(4) :: I,flag
    REAL(8) :: DIS(6)
    
    CALL MonMatrix(XP,Tp,XT,MM)
    CALL EIGENVALUE(6,MM,WR,WI)
    
    DIS = DSQRT(WR**2D0+WI**2D0)
    
    flag = 1
    DO I = 1,6
        IF (ABS(DIS(I)-1D0)>1D-10) THEN
            flag = 0
        ENDIF
    ENDDO
    
    END SUBROUTINE MonEig
    
    SUBROUTINE YHCM(T0, T1, X, Y)
    
    USE GlobalDef
        
    REAL(8) :: X(42), Y(42), T0, T1
    REAL(8) :: F(3), F1B(3), FS(3),XX(6), RM
    REAL(8) :: AA(6,6), AMM(6,6), MM(6,6),DDU(3,3),DDK(3,3)
    
    RM = Ratio_Mass
    XX(1:6) = X(37:42)
    
    CALL DVRRC(XX(1:3),DDU)
    CALL DKRR(XX(1:3), DDK)
    
    AA = 0D0
    AA(1,4) = 1D0
    AA(2,5) = 1D0
    AA(3,6) = 1D0
    AA(4,1) = 1D0 + DDU(1,1) + RM*DDK(1,1)
    AA(4,2) = DDU(1,2) + RM*DDK(1,2)
    AA(4,3) = DDU(1,3) + RM*DDK(1,3)
    AA(5,1) = DDU(2,1) + RM*DDK(2,1)
    AA(5,2) = 1D0 + DDU(2,2) + RM*DDK(2,2)
    AA(5,3) = DDU(2,3) + RM*DDK(2,3)
    AA(6,1) = DDU(3,1) + RM*DDK(3,1)
    AA(6,2) = DDU(3,2) + RM*DDK(3,2)
    AA(6,3) = DDU(3,3) + RM*DDK(3,3)
    AA(4,5) = 2D0
    AA(5,4) = -2D0
    
    DO I = 1,6
        DO J = 1,6
            MM(I,J) = X((I-1)*6+J)
        ENDDO
    ENDDO
    
    CALL MATRIXMATRIX(AA,MM,AMM,6,6,6)
    
    DO I = 1,6
        DO J = 1,6
            Y((I-1)*6+J) = AMM(I,J)
        ENDDO
    ENDDO
    
    ! 求解椭球引力
    CALL Force0(XX(1:3), F)
    
    ! 求第三体引力
    CALL Force1B(XX(1:3), T1, F1B)
    
    ! 太阳光压
    CALL ForceS(omegaA * T1, FS)
    
    Y(37) = XX(4)
    Y(38) = XX(5)
    Y(39) = XX(6)
    Y(40) = 2D0 * OmegaA * XX(5) + OmegaA**2D0*XX(1) + F(1) + F1B(1) + FS(1)
    Y(41) = - 2D0 * OmegaA * XX(4) + OmegaA**2D0*XX(2) + F(2) + F1B(2) + FS(2)
    Y(42) = F(3) + F1B(3) + FS(3)

    END SUBROUTINE YHCM
    
    SUBROUTINE DVRRC(X, DV)
    !解析解给出偏导数
!-----------------------------------------------------------------------
    USE GlobalDef
    IMPLICIT REAL*8(A-H,O-Z)
    DIMENSION X(3),DV(3,3)
    
    R = DSQRT(X(1)**2D0+X(2)**2D0+X(3)**2D0)
    R5 = R**5D0
    R9 = R**9D0
    X2 = X(1)**2D0
    Y2 = X(2)**2D0
    Z2 = X(3)**2D0
    X4 = X2**2D0
    Y4 = Y2**2D0
    Z4 = Z2**2D0
    
    DV(1,1) = (2D0*X2-Y2-Z2)/R5 - 3D0/2D0*C20*(4D0*X4+3D0*X2*Y2-27D0*X2*Z2-Y4+3D0*Y2*Z2+4D0*Z4)/R9 + 3D0*C22*(12D0*X4-51D0*X2*Y2-21D0*X2*Z2+7D0*Y4+9D0*Y2*Z2+2D0*Z4)/R9
    DV(1,2) = 3D0*X(1)*X(2)/R5 - 15D0/2D0*C20*X(1)*X(2)*(X2+Y2-6D0*Z2)/R9 + 105D0*C22*(X2-Y2)*X(1)*X(2)/R9
    DV(1,3) = 3D0*X(1)*X(3)/R5 - 15D0/2D0*C20*X(1)*X(3)*(3D0*X2+3D0*Y2-4D0*Z2)/R9 + 15D0*C22*X(1)*X(3)*(5D0*X2-9D0*Y2-2D0*Z2)/R9
    DV(2,1) = DV(1,2)
    DV(2,2) = -(X2-2D0*Y2+Z2)/R5 + 3D0/2D0*C20*(X4-3D0*X2*Y2-3D0*X2*Z2-4D0*Y4+27D0*Y2*Z2-4D0*Z4)/R9 - 3D0*C22*(7D0*X4-51D0*X2*Y2+9D0*X2*Z2+12D0*Y4-21D0*Y2*Z2+2D0*Z4)/R9
    DV(2,3) = 3D0*X(2)*X(3)/R5 - 15D0/2D0*C20*X(2)*X(3)*(3D0*X2+3D0*Y2-4D0*Z2)/R9+15D0*C22*X(2)*X(3)*(9D0*X2-5D0*Y2+2D0*Z2)/R9
    DV(3,1) = DV(1,3)
    DV(3,2) = DV(2,3)
    DV(3,3) = - (X2+Y2-2D0*Z2)/R5 + 3D0/2D0*C20*(3D0*X4+6D0*X2*Y2-24D0*X2*Z2+3D0*Y4-24D0*Y2*Z2+8D0*Z4)/R9-15D0*C22*(X2-Y2)*(X2+Y2-6D0*Z2)/R9
    
    END SUBROUTINE DVRRC
    
    SUBROUTINE DKRR(X, DK)
    !解析解给出偏导数
!-----------------------------------------------------------------------
    USE GlobalDef
    IMPLICIT REAL*8(A-H,O-Z)
    DIMENSION X(3),DK(3,3)
    
    R = DSQRT(X(1)**2D0+X(2)**2D0+X(3)**2D0)
    R2 = R*R
    X2 = X(1)*X(1)
    Y2 = X(2)*X(2)
    Z2 = X(3)*X(3)
    RB = Axis
    RB2 = RB*RB
    
    thetaA = OmegaA * T1
    thetaN = OmegaN * T1 + thetaN0
    theta = thetaA+thetaN
    XB = RB*DCOS(theta)
    XB2 = XB*XB
    YB = RB*DSIN(theta)
    YB2 = YB*YB
    RBR5 = DSQRT(RB2+R2-2D0*X(1)*XB-2D0*X(2)*YB)**5D0
    
    DK(1,1) = -(RB2-2D0*X2+4D0*X(1)*XB-3D0*XB2+Y2-2D0*X(2)*YB2+Z2)/RBR5
    DK(1,2) = 3D0*(-YB+X(2))*(X(1)-XB)/RBR5
    DK(1,3) = 3D0*X(3)*(X(1)-XB)/RBR5
    DK(2,1) = DK(1,2)
    DK(2,2) = -(RB2+X2-2D0*X(1)*XB-2D0*Y2+4D0*X(2)*YB-3D0*YB2+Z2)/RBR5
    DK(2,3) = 3D0*X(3)*(X(2)-YB)/RBR5
    DK(3,1) = DK(1,3)
    DK(3,2) = DK(2,3)
    DK(3,3) = -(RB2+X2-2D0*X(1)*XB+Y2-2D0*X(2)*YB-2D0*Z2)/RBR5
    
    END SUBROUTINE DKRR