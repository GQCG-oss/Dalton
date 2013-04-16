F95     = gfortran
#FFLAGS  = 
#FFLAGS  = -O3 -xW -ip -w
#LIBS    = -lgoto -lpthread
FFLAGS  = -O3
#LIBS    = -L/home/taylor/src/GotoBLAS -lgoto -lpthread

XOBJS    = blocks.o xyz2dalton.o

AOBJS    = blocks.o aces2dalton.o

DOBJS    = blocks.o distances.o

all:		xyz2dalton aces2dalton distances

xyz2dalton:	$(XOBJS)
		$(F95) $(FFLAGS) -o $@ $(XOBJS) $(LIBS)

aces2dalton:	$(AOBJS)
		$(F95) $(FFLAGS) -o $@ $(AOBJS) $(LIBS)

distances:	$(DOBJS)
		$(F95) $(FFLAGS) -o $@ $(DOBJS) $(LIBS)

.SUFFIXES : .f90 .o

.f90.o:
	$(F95) $(FFLAGS) -c $*.f90

