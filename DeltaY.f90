SUBROUTINE DeltaYup_max(YHC,Xp,TF,DeltaYu)
    ! 这个子函数可以找到在y方向改变初始状态量后轨迹在TF时间积分内的相较于平衡点最大y轴偏移量的最大值
    USE GlobalDef
    EXTERNAL :: YHC
    INTEGER(4) :: I,Num,Flag
    REAL(8) :: step,step0,DeltaYu, X0(6),Xp(6),XX(6),T0,TF,hh,EE,err,Xmid
    REAL(8) :: Xmid6(6)
    
    Num = 100
    step0 = 1D-2
    step = step0
    X0 = Xp
    
    
    DO I = 0,Num
        !step = step0
        Xmid6 = X0
        
201     X0(2) = X0(2) + step
        XX = X0
        Xmid = XX(1)
        
        
        ! 积分状态量XX，若撞击或逃逸则中断积分
        err = 1D-14
        hh = 1D-5
        t0 = 0D0
        DeltaYu = 0D0
        Flag = 1
        DO WHILE (t0<TF)
            Call RKF78(YHC,hh,t0,XX,EE,err,6)
            CALL CheckTotal(T0,XX,Flag)
            IF (Flag==0) THEN
                !WRITE(*,*) 'Flag = ',Flag
                EXIT
            ENDIF 
            IF (Xmid*XX(1)<0 .AND. (XX(2)-rS)>0 .AND. DeltaYu<(XX(2)-rS) )THEN
                DeltaYu = XX(2)-rS
            ENDIF
            Xmid = XX(1)
        ENDDO
        
        !WRITE(*,*) step,X0(2)+step-Xp(2),DeltaYu
        
        IF (step<1E-10 .AND. Flag==0) THEN
            !WRITE(*,*) 'DeltaY, DeltaYmax = '
            !WRITE(*,*) X0(2)+step-Xp(2),DeltaYu
            !WRITE(14,'1X,E25.18,1X,E25.18') X0(2)+step-Xp(2),DeltaYu
            EXIT
        ENDIF
            
        
        IF (step>1E-10 .AND. Flag==0) THEN
            step = step/2D0
            X0 = Xmid6
            GOTO 201
        ENDIF
        
    ENDDO
    
    
    END SUBROUTINE
    
    SUBROUTINE DeltaY_Bound(YHC,Xp,TF,DeltaY1,Num)
    ! 这个子函数可以找到在y方向改变初始状态量后轨迹在TF时间积分内的相较于平衡点最大y轴偏移量
    USE GlobalDef
    EXTERNAL :: YHC
    REAL(8) :: ERR,HH,T0,TF,Xp(6),EE,XX(6)
    REAL(8) :: step,DeltaY1,DeltaYu,DeltaYd,Xmid
    INTEGER(4) :: Num,Flag
103 format(1X, E25.18,1X,E25.18,1X,E25.18)      
    err = 1D-14
    hh = 1D-5
    t0 = 0D0
    
    !Num = 100
    step = DeltaY1/Num
    OPEN(14,file='DeltaY.DAT', status='REPLACE')
    DO I = -Num,Num
        XX = Xp
        XX(2) = Xp(2) + I*step
        Xmid = XX(1)
        DeltaYu = 0D0
        DeltaYd = 0D0
        T0 = 0D0
        DO WHILE (t0<TF)
            Call RKF78(YHC,hh,t0,XX,EE,err,6)
            IF (Xmid*XX(1)<0 .AND. (XX(2)-rS)>0 .AND. DeltaYu<(XX(2)-rS) )THEN
                DeltaYu = XX(2)-rS
            ENDIF
            IF (Xmid*XX(1)<0 .AND. (XX(2)-rS)<0 .AND. DeltaYd>(XX(2)-rS) )THEN
                DeltaYd = XX(2)-rS
            ENDIF
            Xmid = XX(1)
            
            CALL CheckTotal(T0,XX,Flag)
            IF (Flag==0) THEN
                DeltaYu = 0D0
                DeltaYd = 0D0
                EXIT
            ENDIF 
        ENDDO
        WRITE(*,*) I*step,DeltaYu,DeltaYd
        WRITE(14,103) I*step,DeltaYu,DeltaYd
        !Call RKF78(YHC,TF-t0,t0,XX,EE,1D0,6)
    ENDDO
    CLOSE(14)
    
    END SUBROUTINE