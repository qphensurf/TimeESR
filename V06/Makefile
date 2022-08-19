# -=== ===================== ===-
# - DO NOT EDIT BELOW THIS LINE -
# -=== ===================== ===-

# Suffix
SUFFIX = .f90

# Makefile Includes
include ./make.sys
include ./objects
LIBS = $(LIBFFT) $(LIBLAPACK) $(LIBBLAS) $(LIBINTEL)

# -------------
# - COMMANDS  -
# -------------
.SUFFIXES : $(SUFFIX)

all : TimeESR

.f90.o:
	$(FC) $(FFLAGS) -cpp $(DFLAGS) -c $*$(SUFFIX)

TimeESR: $(SOURCE)
	rm -f TimeESR06.x
	$(FC) $(FFLAGS) -cpp $(DFLAGS) $(SOURCE) $(LIBS) -o TimeESR06.x

clean:
	@echo "CLEANING up!"
	rm -f *.mod *.o TimeESR06.x; touch *.*
