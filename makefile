main_FLI: main_FLI.o Global.o LCE_FLI.o YHC.o YHC2_v2.o YHC3_v2.o RKF78.o \
	GetSEP.o YHCN.o CheckImpact.o OTHERS.o PeriodicOrbit.o GetFre.o EIGENS.o

	gfortran -o main_FLI main_FLI.o Global.o LCE_FLI.o YHC.o RKF78.o \
	YHC2_v2.o GetSEP.o YHC3_v2.o YHCN.o CheckImpact.o OTHERS.o PeriodicOrbit.o \
	GetFre.o EIGENS.o
	
#Test_FLI.o: Test_FLI.f90
#	gfortran -c Test_FLI.f90
	
main_FLI.o: main_FLI.f90
	gfortran -c main_FLI.f90
	
Global.o GlobalDef.mod: Global.f90
	gfortran -c Global.f90
	
LCE_FLI.o: LCE_FLI.for
	gfortran -c LCE_FLI.for
	
RKF78.o: GlobalDef.mod RKF78.f90
	gfortran -c RKF78.f90 -std=legacy
	
GetSEP.o: GetSEP.f90
	gfortran -c GetSEP.f90
	
YHC.o: YHC.f90
	gfortran -c YHC.f90
	
YHC2_v2.o: YHC2_v2.f90
	gfortran -c YHC2_v2.f90 -ffree-line-length-none
	
YHC3_v2.o: YHC3_v2.f90
	gfortran -c YHC3_v2.f90 -ffree-line-length-none
	
YHCN.o: YHCN.f90
	gfortran -c YHCN.f90 -ffree-line-length-none
	
CheckImpact.o: CheckImpact.f90
	gfortran -c CheckImpact.f90
	
OTHERS.o: OTHERS.FOR
	gfortran -c OTHERS.FOR
	
EIGENS.o: EIGENS.FOR
	gfortran -c EIGENS.FOR -std=legacy
	
PeriodicOrbit.o: PeriodicOrbit.f90
	gfortran -c PeriodicOrbit.f90 -ffree-line-length-none
	
GetFre.o: GetFre.f90
	gfortran -c GetFre.f90
	
clean:
	rm Test_FLI.o Global.o LCE_FLI.o YHC.o RKF78.o \
	YHC2_v2.o GetSEP.o YHC3_v2.o YHCN.o CheckImpact.o OTHERS.o PeriodicOrbit.o \
	GetFre.o EIGENS.o
	
	