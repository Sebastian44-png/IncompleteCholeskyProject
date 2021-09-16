CC = gcc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include
LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib PCG_ICF_demo 
	- ./PCG_ICF_demo
#	- ./icholesky_demo
#	- ./hpc_demo 4
#	- ./gem_lu_demo 2
#	- ./sed_lu_demo 2
#	- ./sed_icne0 2
#	- ./cg_demo 2
#	- ./linked_list_demo

lib:
	( cd ../Lib ; $(MAKE) )

PCG_ICF_demo: lib PCG_ICF_demo.c Makefile

	$(CC) $(CF) $(I) -o  PCG_ICF_demo PCG_ICF_demo.c $(CS)


clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r PCG_ICF_demo \
    *.a *.dSYM *.obj *.dll
