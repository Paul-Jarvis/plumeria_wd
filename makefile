###############################################################################
#  Makefile for plumeria
#
#  Author is Larry Mastin (lgmastin@usgs.gov)
#
###############################################################################

#  SYSTEM specifies which compiler to use
#    Current available options are:
#      gfortran , f95_dwarf , g95 , open64 , ifort
#    This variable cannot be left blank
#      
SYSTEM = gfortran
#SYSTEM = g95
#SYSTEM = ifort

#  RUN specifies which collection of compilation flags that should be run
#    Current available options are:
#      DEBUG : includes debugging info and issues warnings
#      OPT   : includes optimizations flags for fastest runtime
RUN=OPT

ifeq ($(RUN), DEBUG)
    FFLAGS =  -g3 -pg -Wall -fbounds-check -pedantic -fimplicit-none -Wunderflow -Wuninitialized -ffpe-trap=invalid,zero,overflow -fdefault-real-8
endif
ifeq ($(RUN), OPT)
    FFLAGS =
endif

#Specified compiler.  If you prefer some compiler other than gfortran, modify this line.
ifeq ($(SYSTEM), gfortran)
    FC = gfortran -ffree-form -ffree-line-length-none
endif
ifeq ($(SYSTEM), ifort)
    FC = ifort
endif

OBJECTS =        \
FindT.o          \
FindT_init.o     \
enthfunctions.o  \
zfunctions.o     \
read_input.o     \
Metreader.o      \
Module1.o        \
cashkarpqs.o     \
derivs.o

plume3:  $(OBJECTS) main.f90 zfunctions.o makefile
	$(FC) main.f90 -o plume3 $(OBJECTS) $(FFLAGS)

FindT.o: FindT.f90 enthfunctions.o zfunctions.o Module1.o makefile
	$(FC) -c FindT.f90 $(FFLAGS)

FindT_init.o: FindT_init.f90 enthfunctions.o zfunctions.o Module1.o makefile
	$(FC) -c FindT_init.f90 $(FFLAGS)

enthfunctions.o: enthfunctions.f90 Module1.o makefile
	$(FC) -c enthfunctions.f90 $(FFLAGS)

zfunctions.o: zfunctions.f90 enthfunctions.o Module1.o makefile
	$(FC) -c zfunctions.f90 $(FFLAGS)

read_input.o: read_input.f90 Module1.o Metreader.o makefile
	$(FC) -c read_input.f90 $(FFLAGS)

Metreader.o:  Metreader.f90 Module1.o enthfunctions.o makefile
	$(FC) -c Metreader.f90 $(FFLAGS)

Module1.o:  Module1.f90 makefile
	$(FC) $(FPPFLAGS) -c Module1.f90 $(FFLAGS)

derivs.o:  derivs.f90 Module1.o makefile
	$(FC) -c derivs.f90 $(FFLAGS)

cashkarpqs.o:    cashkarpqs.f90 Module1.o derivs.o makefile
	$(FC) -c cashkarpqs.f90 $(FFLAGS)

clean:
	rm  $(OBJECTS) plume3 *.mod 

