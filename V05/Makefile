FFLAGS =  -O  #-fopenmp 
FF = -c
FC = gfortran

LIBBLAS = -lblas 
LIBLAPACK = -llapack
LIBFFT = # -lfftw3

LIBS = $(LIBFFT) $(LIBLAPACK) $(LIBBLAS)

OBJ= declarations.o declfft.o OpenFiles.o algebra.o io.o spin_parameters.o KrondProd.o \
     H_QD.o QME.o Transport.o SpinTools.o TimeESR.o

TimeESR.out: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o TimeESR05.out

declarations.o: declarations.f90
	$(FC) $(FFLAGS) $(FF)  declarations.f90

declfft.o: declfft.f90
	$(FC) $(FFLAGS) $(FF) declfft.f90 

algebra.o: algebra.f90
	$(FC) $(FFLAGS) $(FF) algebra.f90

io.o: io.f90
	$(FC) $(FFLAGS) $(FF) io.f90

OpenFiles.o: OpenFiles.f90
	$(FC) $(FFLAGS) $(FF) OpenFiles.f90

spin_parameters.o: spin_parameters.f90
	$(FC) $(FFLAGS) $(FF) spin_parameters.f90

KrondProd.o: KrondProd.f90
	$(FC) $(FFLAGS) $(FF) KrondProd.f90

H_QD.o: H_QD.f90
	$(FC) $(FFLAGS) $(FF) H_QD.f90

QME.o: QME.f90
	$(FC) $(FFLAGS) $(FF) QME.f90

Transport.o: Transport.f90
	$(FC) $(FFLAGS) $(FF) Transport.f90

SpinTools.o: SpinTools.f90
	$(FC) $(FFLAGS) $(FF) SpinTools.f90

TimeESR.o: TimeESR.f90
	$(FC) $(FFLAGS) $(FF) TimeESR.f90

clean:
	@echo "CLEANING up!"
	rm -f *.o *.mod

