include ../Makefile.config

# INTEREST LIBS, MODULES and SUBDIRS must be modified at the same time!
INTEREST_LIBS = \
-Llib -linterest

MODULES = LIB_OBJ ST_MOD
SUBDIRS = lib selftest

#
#     Most common build of Dalton
#
dalton.x: $(MODULES)
depend  :
	for i in $(SUBDIRS); do ( cd $$i  && touch Makefile.depend && $(MAKE) $@ ); done

clean :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done
	$(RM) -f *~

LIB_OBJ	:
	cd lib && $(MAKE) all

ST_MOD :
	cd selftest && $(MAKE) all

ftnchek: pre
	ftnchek $(CHEKFLAGS) */*.i > dalton.ftnchek

pre :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done

check : default
	cd test && ./TEST -y linsca
