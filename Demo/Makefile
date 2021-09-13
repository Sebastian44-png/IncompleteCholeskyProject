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
	#- ./gem_lu_demo 2
	#- ./sed_lu_demo 2
	#- ./sed_icne0 2
	- ./cg_demo 2
	# ./linked_list_demo

lib:
	( cd ../Lib ; $(MAKE) )

hpc_demo: lib hpc_demo.c Makefile
	$(CC) $(CF) $(I) -o gem_lu_demo gem_lu_demo.c $(CS)
	$(CC) $(CF) $(I) -o hpc_demo hpc_demo.c $(CS)
	$(CC) $(CF) $(I) -o linked_list_demo linked_list_demo.c $(CS)
	$(CC) $(CF) $(I) -o sed_icne0 sed_icne0_demo.c $(CS)
	$(CC) $(CF) $(I) -o cg_demo cg_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo \
    *.a *.dSYM *.obj *.dll
