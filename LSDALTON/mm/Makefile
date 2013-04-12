all: libmm.a

include ../Makefile.config 

INCLUDES := $(INCLUDES) -I../lsutil -I.

# Keep the file lists sorted for easier modification.

F90_OBJS = \
	mm_box_builder.o	\
	mm_box_processor.o	\
	mm_buffers.o		\
	mm_driver.o		\
	mm_global_paras.o	\
	mm_grid_searcher.o	\
	mm_input.o		\
	mm_interface.o		\
	mm_map_builder.o	\
	mm_memory.o		\
	mm_memory_manager.o	\
	mm_overlap_tree.o	\
	mm_sort.o		\
	mm_stats.o		\
	mm_T_contractors.o	\
	mm_T_pair_builder.o	\
	mm_T_pair_tests.o	\
	mm_tree_T_buffer.o	\
	mm_T_worker.o		\
	mm_Vff_driver.o		\
	mm_Vff_processor.o	\
	mm_W_contractors.o	\
	mm_W_pair_builder.o	\
	mm_W_worker.o           \
        mm_densfit.o

C_OBJS = \
	mm_proc_selector.o


F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

libmm.a: $(F_OBJS) $(F90_OBJS) $(C_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F_OBJS:.o=.F) $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 
