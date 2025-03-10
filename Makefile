#
# This program needs blas and hdf5 libraries
#
### Fortran compiler. We tested gfortran on Linux x86_64
FC     = gfortran
#FC     = ifx
FC90   = $(FC)
LINKER = $(FC)
# Debug flags (used for 'make debug')
# gfortran:
FDBG = -g -fimplicit-none -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall
# Intel fortran:
# FDBG = -g -traceback -C -fpe0

### Disable/Enable OpenMP support
OMP = 

### Compilation/Linking options
FOPT = -O2 

### Extra libraries
FLIB = -lblas -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_fortran

### Include dir(s)
INC = -I/usr/include/hdf5/serial

# Fortran options
FFLAGS = $(OMP) $(INC)  

.SUFFIXES: .f .F .F90 .f90

BIN = molcasto47

OBJ = 

OBJ90 = definitions.o \
        h5extractor.o \
        orbitalextractor.o \
	functions.o

%.mod: %.f90
	$(FC90) $(FOPT) $(FFLAGS) -c $< -o $@

%.o : %.mod

.F.o:
	$(FC) $(FOPT) $(FFLAGS) -c $< -o $@
.F90.o: 
	$(FC90) $(FOPT) $(FFLAGS) -c $< -o $@
.f90.o: 
	$(FC90) $(FOPT) $(FFLAGS) -c $< -o $@

all: $(BIN)

molcasto47: $(OBJ) $(OBJ90) main.o
	$(LINKER) -o $(BIN).exe main.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

clean:
	rm -f *.o *.mod

realclean:
	rm -f $(BIN).exe *.o *.mod

debug: FFLAGS += $(FDBG)
debug: all

