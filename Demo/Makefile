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
#	- ./gem_lu_demo 2
#	- ./linked_list_demo
#	- ./sed_icne0_demo 2
#	- ./ILU_demo 6 0.7
#	- ./ILU_demo
	- ./comparison_demo problem1 1 
	- ./comparison_demo problem1 2 
	- ./comparison_demo problem2 1 
	- ./comparison_demo problem2 2 
	- ./comparison_demo problem3 1 
	- ./comparison_demo problem3 2 
	- ./comparison_demo problem4 1 
	- ./comparison_demo problem4 2

lib:
	( cd ../Lib ; $(MAKE) )

#demo: lib gem_lu_demo.c Makefile
#	$(CC) $(CF) $(I) -o gem_lu_demo gem_lu_demo.c $(CS)
#demo: lib linked_list_demo.c Makefile
#	$(CC) $(CF) $(I) -o linked_list_demo linked_list_demo.c $(CS)
#demo: lib  sed_icne0_demo.c Makefile
#	$(CC) $(CF) $(I) -o sed_icne0_demo sed_icne0_demo.c $(CS)
#demo: lib  ILU_demo.c Makefile
#	$(CC) $(CF) $(I) -o ILU_demo ILU_demo.c $(CS)
demo: lib comparison_demo.c Makefile
	$(CC) $(CF) $(I) -o comparison_demo comparison_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
#	- $(RM) -r gem_lu_demo \
    *.a *.dSYM *.obj *.dll
#	- $(RM) -r linked_list_demo \
    *.a *.dSYM *.obj *.dll
# 	- $(RM) -r sed_icne0_demo \
    *.a *.dSYM *.obj *.dll
#	- $(RM) -r ILU_demo \
    *.a *.dSYM *.obj *.dll
	- $(RM) -r comparison_demo \
    *.a *.dSYM *.obj *.dll
