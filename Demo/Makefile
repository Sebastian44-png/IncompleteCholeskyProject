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

all: lib cholesky_demo
	- ./cholesky_demo
#	- ./hpc_demo 4

lib:
	( cd ../Lib ; $(MAKE) )

cholesky_demo: lib cholesky_demo.c Makefile
	$(CC) $(CF) $(I) -o  cholesky_demo cholesky_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r cholesky_demo \
    *.a *.dSYM *.obj *.dll
