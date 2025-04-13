#
#  Makefile for Intel oneAPI
#

FCC = ifx
FP = mpiifx

OPT = -r8 -zero -O2 -mcmodel=large
FLG = -c
OMP = -qopenmp
LIB = -qmkl=parallel

.SUFFIXES : .f .f90 .o

all: fitreppu reppu_learn_srbf reppu_reconst_srbf

fitreppu: reppu_par.o srbf.o polemapMC.o fitreppu.o
	$(FCC) $(LIB) $(OPT) $^ -o $@

reppu_learn_srbf: mtfort90_2.o reppu_par.o polemapMC.o reservoir.o reppu_learn_srbf.o
	$(FCC) $^ -o $@ $(LIB) $(OMP)

reppu_reconst_srbf: mtfort90_2.o reppu_par.o srbf.o polemapMC.o reservoir.o reppu_reconst_srbf.o
	$(FCC) $^ -o $@ $(LIB) $(OMP)

.f90.o:
	$(FCC) $(FLG) $(OPT) $(OMP) $*.f90

.f.o:
	$(FCC) $(FLG) $(OPT) $*.f

clean:
	rm -f *.o *.mod
