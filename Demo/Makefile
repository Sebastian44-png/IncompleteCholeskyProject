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

all: lib ILU_dem
#	- ./hpc_demo 4

lib:
	( cd ../Lib ; $(MAKE) )

ILU_dem: lib ILU_dem.c Makefile
	$(CC) $(CF) $(I) -o ILU_dem ILU_dem.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r ILU_dem \
    *.a *.dSYM *.obj *.dll
