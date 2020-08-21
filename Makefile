# Makefile for RAMSES utils
F90     = gfortran -ffree-form
#F90     = ifort
#DEBUG   = -g -fbounds-check -fbacktrace
FFLAGS  = -O3 -DWITHOUTMPI -DNPRE=8 -DNPACS=3 ${DEBUG}
CC      = cc 3

# Make all targets
.SUFFIXES: .o .f90  .c
.f90.o:
	${F90} ${FFLAGS} -c $<
.c.o:
	$(CC)   -c $<

all: amr2Intcube intcube2velcube

amr2Intcube: ramses_info.o amr2Intcube.o
	$(F90) ${FFLAGS} $(SOURCES1) ramses_info.o amr2Intcube.o -o amr2Intcube

intcube2velcube:
	$(F90) ${FFLAGS} intcube2velcube.f90 -o intcube2velcube

# Make a specific f90 file
#%: %.f90
#	$(F90) $< -o $(BINDIR)/$@

#############################################################################
%.o:%.f90
	$(F90)  $(FFLAGS) -c $^ -o $@
#############################################################################
clean :
	rm *.o *.mod
	rm amr2Intcube
	rm intcube2velcube
#############################################################################
