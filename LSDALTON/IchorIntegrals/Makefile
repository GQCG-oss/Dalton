all: libIchorIntegrals.a

include ../Makefile.config

INCLUDES := $(INCLUDES) 

# Keep the file lists sorted according to order of compilation.

F90_OBJS = \
        IchorPresicion.o \
        AGC_CPU_BuildRJ000Gen.o\
        AGC_CPU_BuildRJ000Seg1Prim.o\
        AGC_CPU_VerticalRecurrenceQPAGen.o \
        AGC_CPU_VerticalRecurrenceQPBGen.o \
        AGC_CPU_VerticalRecurrenceQPCGen.o \
        AGC_CPU_VerticalRecurrenceQPDGen.o \
        AGC_CPU_VerticalRecurrenceQPASegQ.o \
        AGC_CPU_VerticalRecurrenceQPBSegQ.o \
        AGC_CPU_VerticalRecurrenceQPCSegQ.o \
        AGC_CPU_VerticalRecurrenceQPDSegQ.o \
        AGC_CPU_VerticalRecurrenceQPASegP.o \
        AGC_CPU_VerticalRecurrenceQPBSegP.o \
        AGC_CPU_VerticalRecurrenceQPCSegP.o \
        AGC_CPU_VerticalRecurrenceQPDSegP.o \
        AGC_CPU_VerticalRecurrenceQPASeg.o \
        AGC_CPU_VerticalRecurrenceQPBSeg.o \
        AGC_CPU_VerticalRecurrenceQPCSeg.o \
        AGC_CPU_VerticalRecurrenceQPDSeg.o \
        AGC_CPU_VerticalRecurrenceQPASeg1Prim.o \
        AGC_CPU_VerticalRecurrenceQPBSeg1Prim.o \
        AGC_CPU_VerticalRecurrenceQPCSeg1Prim.o \
        AGC_CPU_VerticalRecurrenceQPDSeg1Prim.o \
        AGC_CPU_TransferRecurrenceAtoCGen.o \
        AGC_CPU_TransferRecurrenceAtoDGen.o \
        AGC_CPU_TransferRecurrenceDtoAGen.o \
        AGC_CPU_TransferRecurrenceCtoAGen.o \
        AGC_CPU_TransferRecurrenceBtoCGen.o \
        AGC_CPU_TransferRecurrenceBtoDGen.o \
        AGC_CPU_TransferRecurrenceDtoBGen.o \
        AGC_CPU_TransferRecurrenceCtoBGen.o \
        AGC_CPU_TransferRecurrenceAtoCSegQ.o \
        AGC_CPU_TransferRecurrenceAtoDSegQ.o \
        AGC_CPU_TransferRecurrenceDtoASegQ.o \
        AGC_CPU_TransferRecurrenceCtoASegQ.o \
        AGC_CPU_TransferRecurrenceBtoCSegQ.o \
        AGC_CPU_TransferRecurrenceBtoDSegQ.o \
        AGC_CPU_TransferRecurrenceDtoBSegQ.o \
        AGC_CPU_TransferRecurrenceCtoBSegQ.o \
        AGC_CPU_TransferRecurrenceAtoCSegP.o \
        AGC_CPU_TransferRecurrenceAtoDSegP.o \
        AGC_CPU_TransferRecurrenceDtoASegP.o \
        AGC_CPU_TransferRecurrenceCtoASegP.o \
        AGC_CPU_TransferRecurrenceBtoCSegP.o \
        AGC_CPU_TransferRecurrenceBtoDSegP.o \
        AGC_CPU_TransferRecurrenceDtoBSegP.o \
        AGC_CPU_TransferRecurrenceCtoBSegP.o \
        AGC_CPU_TransferRecurrenceAtoCSeg.o \
        AGC_CPU_TransferRecurrenceAtoDSeg.o \
        AGC_CPU_TransferRecurrenceDtoASeg.o \
        AGC_CPU_TransferRecurrenceCtoASeg.o \
        AGC_CPU_TransferRecurrenceBtoCSeg.o \
        AGC_CPU_TransferRecurrenceBtoDSeg.o \
        AGC_CPU_TransferRecurrenceDtoBSeg.o \
        AGC_CPU_TransferRecurrenceCtoBSeg.o \
        AGC_CPU_TransferRecurrenceAtoCSeg1Prim.o \
        AGC_CPU_TransferRecurrenceAtoDSeg1Prim.o \
        AGC_CPU_TransferRecurrenceDtoASeg1Prim.o \
        AGC_CPU_TransferRecurrenceCtoASeg1Prim.o \
        AGC_CPU_TransferRecurrenceBtoCSeg1Prim.o \
        AGC_CPU_TransferRecurrenceBtoDSeg1Prim.o \
        AGC_CPU_TransferRecurrenceDtoBSeg1Prim.o \
        AGC_CPU_TransferRecurrenceCtoBSeg1Prim.o \
        AGC_HorizontalRecurrencePAtoB.o \
        AGC_HorizontalRecurrencePBtoA.o \
        AGC_HorizontalRecurrenceQCtoD.o \
        AGC_HorizontalRecurrenceQDtoC.o \
        AGC_SphContractOBS1.o \
        AGC_SphContractOBS2.o \
        IchorParameters.o \
        IchorCommon.o \
        IchorBatchTools.o \
        IchorMem.o \
	IchorSaveGab.o \
        IchorGammaTabulation.o \
        IchorEri_CoulombIntegral_McM_seg_seg_SSSS.o \
        IchorEri_CoulombIntegral_OBS_Gen.o \
        IchorEri_CoulombIntegral_OBS_SegQ.o \
        IchorEri_CoulombIntegral_OBS_SegP.o \
        IchorEri_CoulombIntegral_OBS_Seg.o \
        IchorEri_CoulombIntegral_OBS_Seg1Prim.o \
        IchorEri_CoulombIntegral_OBS_general.o \
	IchorEri_GabIntegral_OBS_Gen.o \
	IchorEri_GabIntegral_OBS_Seg.o \
	IchorEri_GabIntegral_OBS_general.o \
        IchorEri.o \
        IchorGab.o 

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
