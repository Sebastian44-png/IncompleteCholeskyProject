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

<<<<<<< HEAD
all: lib ILU_dem
#	- ./hpc_demo 4
=======
all: lib hpc_demo
	- ./sed_lu_demo 2
>>>>>>> d32f0aa98ae08c3cb23931b0a6aad0f35d1a4c5e

lib:
	( cd ../Lib ; $(MAKE) )

<<<<<<< HEAD
ILU_dem: lib ILU_dem.c Makefile
	$(CC) $(CF) $(I) -o ILU_dem ILU_dem.c $(CS)
=======
hpc_demo: lib hpc_demo.c sed_lu_demo.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo hpc_demo.c $(CS)
	$(CC) $(CF) $(I) -o sed_lu_demo sed_lu_demo.c $(CS)
>>>>>>> d32f0aa98ae08c3cb23931b0a6aad0f35d1a4c5e

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r ILU_dem \
    *.a *.dSYM *.obj *.dll
