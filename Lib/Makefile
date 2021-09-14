CC = gcc

LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 bnd_compress.o  \
 bnd_print.o  \
 bnd_util.o  \
 cs_compress.o \
 cs_entry.o \
 cs_spmv.o \
 cs_load.o \
 cs_print.o \
 cs_util.o \
 cs_lapmat_p1_square.o \
 gem_compress.o \
 gem_gauss.o \
 gem_gausssol.o \
 gem_lu.o \
 gem_print.o \
 gem_util.o \
 gem_spmv.o \
 hpc_malloc.o \
 jds_spmv.o \
 jds_print.o \
 jds_util.o  \
 sed_compress.o \
 sed_icholesky.o \
 sed_isoccupied.o \
 sed_spmv.o \
 sed_print.o \
 sed_util.o  \
 sed_gauss_seidel.o \
 sed_ILU.o \
 sed_MILU.o \
 sed_jacobi.o \
 sed_richardson.o \
 sed_lu.o \
 sky_cholesky.o \
 sed_cg.o \
 sed_icne0.o \
 sky_compress.o \
 sky_print.o  \
 sky_util.o \
 sky_load.o \
 sky_spmv.o \
 hpc_cumsum.o \
 hpc_malloc.o \
 s_list.o \
 utils.o 

$(HPC): ../Include/hpc.h Makefile

%.o: ../Source/%.c ../Include/hpc.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(HPC)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

