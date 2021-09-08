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
all: lib cholesky_demo
	- ./cholesky_demo
#	- ./hpc_demo 4
=======
all: lib ILU_dem
#	- ./hpc_demo 4
	- ./ILU_dem
>>>>>>> 9ec61b2b41b1f142e9e65f11a90e5294f06f660d

lib:
	( cd ../Lib ; $(MAKE) )

<<<<<<< HEAD
cholesky_demo: lib cholesky_demo.c Makefile
	$(CC) $(CF) $(I) -o  cholesky_demo cholesky_demo.c $(CS)
=======
ILU_dem: lib ILU_dem.c Makefile
	$(CC) $(CF) $(I) -o ILU_dem ILU_dem.c $(CS)
>>>>>>> 9ec61b2b41b1f142e9e65f11a90e5294f06f660d

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r cholesky_demo \
    *.a *.dSYM *.obj *.dll
