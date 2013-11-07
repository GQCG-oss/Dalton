all: libIchorIntegrals.a

include ../Makefile.config

INCLUDES := $(INCLUDES) 

# Keep the file lists sorted for easier modification.

F90_OBJS = \
        IchorPresicion.o \
        IchorCommon.o \
        IchorMem.o \
        IchorGammaTabulation.o \
        IchorEri_CoulombIntegral_McM_seg_seg_SSSS.o \
        IchorEri_CoulombIntegral_OBS_general.o \
        IchorEri.o \
        AGC_VerticalRecurrenceQPA.o \
        AGC_VerticalRecurrenceQPB.o \
        AGC_VerticalRecurrenceQPC.o \
        AGC_VerticalRecurrenceQPD.o \
        AGC_VerticalRecurrenceQPASegQ.o \
        AGC_VerticalRecurrenceQPBSegQ.o \
        AGC_VerticalRecurrenceQPCSegQ.o \
        AGC_VerticalRecurrenceQPDSegQ.o \
        AGC_VerticalRecurrenceQPASegP.o \
        AGC_VerticalRecurrenceQPBSegP.o \
        AGC_VerticalRecurrenceQPCSegP.o \
        AGC_VerticalRecurrenceQPDSegP.o \
        AGC_VerticalRecurrenceQPASeg.o \
        AGC_VerticalRecurrenceQPBSeg.o \
        AGC_VerticalRecurrenceQPCSeg.o \
        AGC_VerticalRecurrenceQPDSeg.o \
        AGC_VerticalRecurrenceQPASeg1Prim.o \
        AGC_VerticalRecurrenceQPBSeg1Prim.o \
        AGC_VerticalRecurrenceQPCSeg1Prim.o \
        AGC_VerticalRecurrenceQPDSeg1Prim.o \
        AGC_TransferRecurrenceQPAtoC.o \
        AGC_TransferRecurrenceQPAtoD.o \
        AGC_TransferRecurrenceQPDtoA.o \
        AGC_TransferRecurrenceQPCtoA.o \
        AGC_TransferRecurrenceQPBtoC.o \
        AGC_TransferRecurrenceQPBtoD.o \
        AGC_TransferRecurrenceQPDtoB.o \
        AGC_TransferRecurrenceQPCtoB.o \
        AGC_HorizontalRecurrencePAtoB.o \
        AGC_HorizontalRecurrencePBtoA.o \
        AGC_HorizontalRecurrenceQCtoD.o \
        AGC_HorizontalRecurrenceQDtoC.o \
        AGC_SphContractOBS1.o \
        AGC_SphContractOBS2.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
libIchorIntegrals.a: $(F_OBJS) $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 
