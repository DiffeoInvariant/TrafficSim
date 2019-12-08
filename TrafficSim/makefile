TARGETS = TSimLib
PETSC_DIR=/home/diffeoinvariant/petsc-3.12.1
PETSC_ARCH = arch-linux2-c-debug

TS_C_DIR = $(shell pwd)
TS_C_INCL = $(TS_C_DIR)/include
TS_C_SRC = $(TS_C_DIR)/src

TS_C_IMPLS = $(TS_C_SRC)/tsdensity.c $(TS_C_SRC)/tshighway.c $(TS_C_SRC)/tsnetwork.c

include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

TS_CC = clang
TS_COPTS = -std=c99 -O3 -Wall -g -shared -fpic



TSimLib: $(TS_C_IMPLS)
	$(TS_CC) $(TS_COPTS) -I $(TS_C_INCL) $(TS_C_IMPLS) -o libtsim.so $(PETSC_CC_INCLUDES) -L $(PETSC_DIR)/$(PETSC_ARCH)/lib/ -lpetsc -lmpi -lm -ldl

all: $(TARGETS)
