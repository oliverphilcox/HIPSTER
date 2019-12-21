## MAKEFILE FOR grid_power.cpp. This compiles the grid_power.cpp file into the ./power exececutable.

## Periodic flag from command line input?
PERIODIC_FLAG=$(Periodic)

CC = gcc-9.1
CFLAGS = -O3 -Wall
CXXFLAGS = -DPOWER -Wall -O3  $(PERIODIC_FLAG)
#-DOPENMP
# disable OPENMP to run single threaded
#-DPERIODIC # use this to enable periodic behavior

CXX = g++-9.1 -fopenmp -lgomp -std=c++0x -ffast-math

AUNTIE	= power
AOBJS	= grid_power.o

LD	= g++-9.1
LFLAGS	= -L/usr/local/lib -lgsl -lgslcblas
# -lgomp -L/usr/lib/x86_64-linux-gnu 


main: $(AUNTIE)

$(AUNTIE):	$(AOBJS)
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm power
