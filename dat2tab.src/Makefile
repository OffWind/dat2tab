#----------------------------------------------------------------------#
# Makefile for DAT2TAB                                                 #
#----------------------------------------------------------------------#
# This Makefile was prepared for use with Intel Fortran compiler.
# Change 'ifort' to your Fortran compile of choice and adjust options

# Intel fortran
CMF   =  ifort -c -O3
CMF_EX   = ifort

dat2tab: dat2tab.o  mjseries_io.o
	$(CMF_EX) -o dat2tab dat2tab.o mjseries_io.o -lm ; \mv dat2tab ..

all:
	 make dat2tab

clean:
	rm -f *.o *.mod dat2tab *~

dat2tab.o :  dat2tab.f90
	$(CMF)	dat2tab.f90 

mjseries_io.o :  mjseries_io.f90
	$(CMF)	mjseries_io.f90
