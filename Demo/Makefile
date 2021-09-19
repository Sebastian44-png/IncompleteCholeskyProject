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

all: lib demo
#	- ./PCG_ICF_demo
#	- ./icholesky_demo
#	- ./hpc_demo 4
#- ./gem_lu_demo 2
#- ./sed_lu_demo j
#- ./sed_icne0 2
#- ./sed_lu_demo 2
#	- ./sed_icne0 3
#	- ./cg_demo 2
# ./linked_list_demo
#   - ./mg_demo1 problem3 2
	- ./pcg_mg_demo1 problem1 1 
#	- ./mg_demo1 problem4 1

lib:
	( cd ../Lib ; $(MAKE) )

demo: lib pcg_mg_demo1.c Makefile
#$(CC) $(CF) $(I) -o  icholesky_demo icholesky_demo.c $(CS)
#$(CC) $(CF) $(I) -o gem_lu_demo gem_lu_demo.c $(CS)
#$(CC) $(CF) $(I) -o hpc_demo hpc_demo.c $(CS)
#$(CC) $(CF) $(I) -o linked_list_demo linked_list_demo.c $(CS)
#$(CC) $(CF) $(I) -o sed_icne0 sed_icne0_demo.c $(CS)
# 	$(CC) $(CF) $(I) -o mg_demo1 mg_demo1.c $(CS)
	$(CC) $(CF) $(I) -o pcg_mg_demo1 pcg_mg_demo1.c $(CS)
#	$(CC) $(CF) $(I) -o  PCG_ICF_demo PCG_ICF_demo.c $(CS)
#	$(CC) $(CF) $(I) -o sed_icne0 sed_icne0_demo.c $(CS)
#	$(CC) $(CF) $(I) -o cg_demo cg_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r pcg_mg_demo \
    *.a *.dSYM *.obj *.dll
