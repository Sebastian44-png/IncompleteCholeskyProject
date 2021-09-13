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

all: lib icholesky_demo
	- ./icholesky_demo
#	- ./hpc_demo 4

lib:
	( cd ../Lib ; $(MAKE) )

icholesky_demo: lib icholesky_demo.c Makefile
	$(CC) $(CF) $(I) -o  icholesky_demo icholesky_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r icholesky_demo \
    *.a *.dSYM *.obj *.dll
