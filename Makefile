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

all: fitreppu reppu_learn_secs reppu_reconst_secs

fitreppu: reppu_par.o secs_parallel.o polemapMC.o fitreppu.o
	$(FCC) $(LIB) $(OPT) $^ -o $@

reppu_learn_secs: mtfort90_2.o reppu_par.o polemapMC.o reservoir.o reppu_learn_secs.o
	$(FCC) $^ -o $@ $(LIB) $(OMP)

reppu_reconst_secs: mtfort90_2.o reppu_par.o secs_parallel.o polemapMC.o reservoir.o reppu_reconst_secs.o
	$(FCC) $^ -o $@ $(LIB) $(OMP)

.f90.o:
	$(FCC) $(FLG) $(OPT) $(OMP) $*.f90

.f.o:
	$(FCC) $(FLG) $(OPT) $*.f

clean:
	rm -f *.o *.mod
