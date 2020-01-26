## MAKEFILE FOR grid_power.cpp. This compiles the grid_power.cpp file into the ./power exececutable.

## Periodic flag from command line input?
PERIODIC_FLAG=$(Periodic)

CC = gcc
CFLAGS = -g -Wall
#-O3 -Wall
CXXFLAGS = -DPOWER -Wall -O3 $(PERIODIC_FLAG) -DBISPECTRUM -DPERIODIC
#-DOPENMP
# disable OPENMP to run single threaded
# use -DPERIODIC to run in periodic mode, or set Periodic=-DPERIODIC on the command line

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math

AUNTIE	= power
AOBJS	= grid_power.o

LD	= g++
LFLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lgomp -L/usr/lib/x86_64-linux-gnu

main: $(AUNTIE)

$(AUNTIE):	$(AOBJS)
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm power
