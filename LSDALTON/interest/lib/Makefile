include ../../Makefile.config

f90:=${F90}
normalflags:=-c -O2
highflags:=${F90OPTFLAGS}

all: libinterest.a
	@echo "Interest library is compiled!"

clean:
	@set -e; \
	echo -n "Cleaning Interest library ... "; \
	rm -f *.mod *.o *.a mak.log 2>&1 && echo "Sucessfully done!" \
        || (echo "Error occured!" && exit 1);

depend:


ROUTINES = module_interest_osr.o module_interest_hrr.o module_interest.o
           

libinterest.a: ${ROUTINES}
	$(AR) $(ARFLAGS) $@ $?
#
# ROUTINES
#
module_interest_osr.o: ../src/module_interest_osr.f90
	${f90} ${highflags} -c ../src/module_interest_osr.f90 
module_interest_hrr.o: ../src/module_interest_hrr.f90
	${f90} ${highflags} -c ../src/module_interest_hrr.f90 
module_interest.o: ../src/module_interest.f90
	${f90} ${highflags} -c ../src/module_interest.f90 
