#
#  Makefile
#

## for Intel oneAPI
FCC = ifx
#FCC = ifort
OPT = -r8 -zero -O2 -mcmodel=large
OMP = -qopenmp
LIB = -qmkl=parallel

# ## for GNU Fortran
# FCC = gfortran
# OPT=-fdefault-real-8 -O3
# OMP = -fopenmp
# LIB = -lblas -llapack

FLG = -c


.SUFFIXES : .f .f90 .o

reppu_emulator: mtfort90.o reppu_par.o srbf.o polemap_l.o reservoir.o reppu_emulator.o
	$(FCC) $^ -o $@ $(OPT) $(LIB)

.f90.o:
	$(FCC) $(FLG) $(OPT) $(OMP) $*.f90

clean:
	rm -f *.o *.mod
