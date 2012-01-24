#
# FILE    : Makefile for dalton version of Gen1Int interface
#           Oct. 2, 2011, by Bin Gao
#

all : depend libgen1int.a

include ../Makefile.config

# directory for object files
OBJDIR := $(PWD)

# Gen1Int library source code
GEN1INT_SRCS := error_stop.F90         \
                xtimer.F90             \
                norm_contr_cgto.F90    \
                norm_contr_sgto.F90    \
                reorder_ints.F90       \
                trace_ints.F90         \
                dump_info.F90          \
                binom_coeff.F90        \
                next_permutation.F90   \
                geom_total.F90         \
                sort_cents.F90         \
                shell_scatter.F90      \
                geom_part_zero.F90     \
                geom_part_one.F90      \
                hgto_to_cgto.F90       \
                hgto_to_lcgto.F90      \
                hgto_to_sgto.F90       \
                london_mom_hgto.F90    \
                const_contr_ints.F90   \
                carmom_deriv.F90       \
                carmom_init.F90        \
                carmom_to_hbra.F90     \
                carmom_hrr_ket.F90     \
                prim_hgto_carmom.F90   \
                contr_cgto_carmom.F90  \
                contr_sgto_carmom.F90  \
                contr_csgto_carmom.F90 \
                aux_boys_vec.F90       \
                nucpot_nderv.F90       \
                nucpot_hket.F90        \
                nucpot_hbra.F90        \
                nucpot_to_carmom.F90   \
                nucpot_to_hgto.F90     \
                prim_hgto_nucpot.F90   \
                contr_cgto_nucpot.F90  \
                contr_sgto_nucpot.F90  \
                gaupot_nderv.F90       \
                prim_hgto_gaupot.F90   \
                contr_cgto_gaupot.F90  \
                contr_sgto_gaupot.F90  \
                gen1int.F90
GEN1INT_OBJS := $(addprefix $(OBJDIR)/,$(GEN1INT_SRCS:.F90=.o))
#GEN1INT_SRCS := $(addprefix $(LIB_GEN1INT_PATH)/src/,$(GEN1INT_SRCS))

$(OBJDIR)/%.o: $(LIB_GEN1INT_PATH)/src/%.F90
	$(F90) $(INCLUDES) $(MODULES) $(CPPFLAGS) $(FFLAGS) -c $<

# Gen1Int interface source code
#
# F90_MODS should contain all files with F90 modules in proper order
#          (i.e. if b.F90 calls a module in k.F90, then the proper order is k.F90 b.F90)
F90_MODS = gen1int_matrix.F90 \
           gen1int_shell.F90  \
           gen1int_interface.F90

F90_MOD_OBJS = $(F90_MODS:.F90=.o)

libgen1int.a : $(GEN1INT_OBJS) $(F90_MOD_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a *.mod

Makefile.depend depend: ../include/*.h *.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend

#mod:
#	$(MAKE) --no-print-directory $(F90_MOD_OBJS)

# Fortran compiler creates a file .mod for each module name.
# This confuses the gnumake, which thinks that .mod files are
# Modula2 source code files. We must override this built-in rule
# by adding the following line:
%.o: %.mod
