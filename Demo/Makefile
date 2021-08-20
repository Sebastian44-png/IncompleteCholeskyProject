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

all: lib hpc_demo
	- ./hpc_demo 4

lib:
	( cd ../Lib ; $(MAKE) )

hpc_demo: lib hpc_demo.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo hpc_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo \
    *.a *.dSYM *.obj *.dll
