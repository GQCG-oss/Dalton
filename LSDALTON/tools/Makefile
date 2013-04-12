include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../linears -I../geomopt -I../ddynam -I../deccc -I../lsdaltonsrc

.PHONY: clean

clean:
	rm -f *.o *.mod

# Gao: Fortran compiler creates a file .mod for each module name.
# This confuses the gnumake, which thinks that .mod files are
# Modula2 source code files. We must override this built-in rule
# by adding the following line:
%.o: %.mod
