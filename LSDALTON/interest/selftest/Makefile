include ../../Makefile.config

f90:=${F90}
normalflags:=-c -O2
highflags:=${F90OPTFLAGS}
linkflags:=-03
LIB_PATH:=../lib
LIB:=interest

all: libinterest_selftest.x
	@echo "Interest selftest part is compiled!"

clean:
	@set -e; \
	echo -n "Cleaning Interest selftest part ... "; \
	rm -f *.mod *.o *.a *.x mak.log 2>&1 && echo "Sucessfully done!" \
        || (echo "Error occured!" && exit 1);

ROUTINES = module_interest_selftest.o main_interest_selftest.o
           

libinterest_selftest.x: ${ROUTINES}
	${F77} ${FFLAGS} -o libinterest_selftest.x  $(ROUTINES) -L$(LIB_PATH) -l$(LIB)
#
# ROUTINES
#
main_interest_selftest.o: main_interest_selftest.f90
	${f90} ${highflags} -c main_interest_selftest.f90 
module_interest_selftest.o: module_interest_selftest.f90
	${f90} ${highflags} -c module_interest_selftest.f90 
