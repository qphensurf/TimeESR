# -=== ===================== ===-
# - DO NOT EDIT BELOW THIS LINE -
# -=== ===================== ===-

# Makefile Includes
include ./make.sys
include ./objects
LIBS = $(LIBFFT) $(LIBLAPACK) $(LIBBLAS) $(LIBINTEL)

# -------------
# - COMMANDS  -
# -------------
.SUFFIXES :
.SUFFIXES : .o .c .f .f90

all : TimeESR

.f90.o:
	$(FC) $(FFLAGS) -cpp $(DFLAGS) -c $*.f90

.f.o:
	$(FC) $(FFLAGS) -cpp $(DFLAGS) -c $*.f

TimeESR: $(SOURCE)
	rm -f TimeESR.x
	$(FC) $(FFLAGS) -cpp $(DFLAGS) $(SOURCE) $(LIBS) -o TimeESR.x

clean:
	@echo "CLEANING up!"
	rm -f *.mod *.o *.__genmod.f90 *.x *.a; touch *.*
