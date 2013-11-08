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
        IchorEri_CoulombIntegral_OBS_Gen.o \
        IchorEri_CoulombIntegral_OBS_SegQ.o \
        IchorEri_CoulombIntegral_OBS_SegP.o \
        IchorEri_CoulombIntegral_OBS_Seg.o \
        IchorEri_CoulombIntegral_OBS_Seg1Prim.o \
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
        AGC_TransferRecurrenceAtoCGen.o \
        AGC_TransferRecurrenceAtoDGen.o \
        AGC_TransferRecurrenceDtoAGen.o \
        AGC_TransferRecurrenceCtoAGen.o \
        AGC_TransferRecurrenceBtoCGen.o \
        AGC_TransferRecurrenceBtoDGen.o \
        AGC_TransferRecurrenceDtoBGen.o \
        AGC_TransferRecurrenceCtoBGen.o \
        AGC_TransferRecurrenceAtoCSegQ.o \
        AGC_TransferRecurrenceAtoDSegQ.o \
        AGC_TransferRecurrenceDtoASegQ.o \
        AGC_TransferRecurrenceCtoASegQ.o \
        AGC_TransferRecurrenceBtoCSegQ.o \
        AGC_TransferRecurrenceBtoDSegQ.o \
        AGC_TransferRecurrenceDtoBSegQ.o \
        AGC_TransferRecurrenceCtoBSegQ.o \
        AGC_TransferRecurrenceAtoCSegP.o \
        AGC_TransferRecurrenceAtoDSegP.o \
        AGC_TransferRecurrenceDtoASegP.o \
        AGC_TransferRecurrenceCtoASegP.o \
        AGC_TransferRecurrenceBtoCSegP.o \
        AGC_TransferRecurrenceBtoDSegP.o \
        AGC_TransferRecurrenceDtoBSegP.o \
        AGC_TransferRecurrenceCtoBSegP.o \
        AGC_TransferRecurrenceAtoCSeg.o \
        AGC_TransferRecurrenceAtoDSeg.o \
        AGC_TransferRecurrenceDtoASeg.o \
        AGC_TransferRecurrenceCtoASeg.o \
        AGC_TransferRecurrenceBtoCSeg.o \
        AGC_TransferRecurrenceBtoDSeg.o \
        AGC_TransferRecurrenceDtoBSeg.o \
        AGC_TransferRecurrenceCtoBSeg.o \
        AGC_TransferRecurrenceAtoCSeg1Prim.o \
        AGC_TransferRecurrenceAtoDSeg1Prim.o \
        AGC_TransferRecurrenceDtoASeg1Prim.o \
        AGC_TransferRecurrenceCtoASeg1Prim.o \
        AGC_TransferRecurrenceBtoCSeg1Prim.o \
        AGC_TransferRecurrenceBtoDSeg1Prim.o \
        AGC_TransferRecurrenceDtoBSeg1Prim.o \
        AGC_TransferRecurrenceCtoBSeg1Prim.o \
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
