      SUBROUTINE LCE_FLI(YHC,X0,FLI,TTOL,KK)
	IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL :: YHC
	PARAMETER(PI=3.14159265358979323846D0,PI2=PI*2D0)
	PARAMETER(AU=1.49597870691D011,EM=3.84747981D008)
	PARAMETER(AM=1.738D006,AE=6.378137D006)
	DIMENSION X0(6),XX(42),AA(6,6),AT(6,6),BB(6,6)
	DIMENSION ST1(6),ST2(6),DST(6)
      INTEGER(4) :: FLAG
C	COMMON /CN/ U
	OPEN(3,FILE="ORBIT.DAT")
	KK=0
	DD=1D-012
	AA=0D0
	DO I=1,6
	  AA(I,I)=1D0
	ENDDO
	DO I=1,6
	DO J=1,6
	  XX(6*(I-1)+J)=AA(I,J)
	ENDDO
	ENDDO
	XX(37:42)=X0(1:6)
	T0=0D0
	TT=TTOL
	HH=1D-004
      FLAG=1
	DO WHILE(T0.LT.TT)
C	  CALL RKF2(HH,T0,XX,EE,DD,42)
        CALL RKF78(YHC,HH,T0,XX,EE,DD,42)
	  !R1=DSQRT((XX(37)+U)**2+XX(38)**2+XX(39)**2)
	  !R2=DSQRT((XX(37)+U-1D0)**2+XX(38)**2+XX(39)**2)
	  !IF(R1.LT.1D-002.OR.R2.LT.1D-002)THEN
	  !  KK=1
	  !  RETURN
	  !ENDIF
        CALL CheckImpact2(T0,XX(37:42),FLAG)
        IF (FLAG==0) THEN
            KK = 1
            GOTO 201
        ENDIF
        CALL CheckEscape(XX(37:42),FLAG)
        IF (FLAG==0) THEN
            KK = 2
            GOTO 201
        ENDIF
        CALL CheckOrbit(XX(37:42),FLAG)
        IF (FLAG==0) THEN
            KK = 3
            GOTO 201
        ENDIF
	ENDDO
C	CALL RKF2(TT-T0,T0,XX,EE,1D0,42)
      CALL RKF78(YHC,TT-T0,T0,XX,EE,1D0,42)
	!R1=DSQRT((XX(37)+U)**2+XX(38)**2+XX(39)**2)
	!R2=DSQRT((XX(37)+U-1D0)**2+XX(38)**2+XX(39)**2)
	!IF(R1.LT.1D-002.OR.R2.LT.1D-002)THEN
	!  KK=1
	!  RETURN
	!ENDIF
      CALL CheckImpact2(T0,XX(37:42),FLAG)
      IF (FLAG==0) THEN
          KK = 1
          GOTO 201
      ENDIF
      CALL CheckEscape(XX(37:42),FLAG)
      IF (FLAG==0) THEN
          KK = 2
          GOTO 201
      ENDIF
      CALL CheckOrbit(XX(37:42),FLAG)
      IF (FLAG==0) THEN
          KK = 3
          GOTO 201
      ENDIF
	DO I=1,36
	  J=INT((I+5)/6)
	  AA(J,I-6*(J-1))=XX(I)
	ENDDO
	DO I=1,6
	  AT(I,1:6)=AA(1:6,I)
	ENDDO
	CALL MATRIXMATRIX(AT,AA,BB,6,6,6)
	ST1(1:6)=1D0/DSQRT(6D0)
	EPS=1D0
	DO WHILE(EPS.GT.1D-012)
	  CALL MATRIXVECTOR(BB,ST1,ST2,6,6)
	  ST2_MOD=DSQRT(ST2(1)**2+ST2(2)**2+ST2(3)**2
     &               +ST2(4)**2+ST2(5)**2+ST2(6)**2)
	  ST2(1:6)=ST2(1:6)/ST2_MOD
	  DST=ST2-ST1
	  EPS=DSQRT(DST(1)**2+DST(2)**2+DST(3)**2
     &           +DST(4)**2+DST(5)**2+DST(6)**2)
	  ST1=ST2
	ENDDO
	FLI=DSQRT(ST2_MOD)
	RETURN
201	END