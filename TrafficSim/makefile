TARGETS = TSimLib

TS_C_DIR = $(shell pwd)
TS_C_INCL = $(TS_C_DIR)/include
TS_C_SRC = $(TS_C_DIR)/src

TS_C_IMPLS = network.c

TS_CC = clang
TS_COPTS = -std=C99 -O3 -Wall -g -shared -fpic

TSimLib:
	$(TS_CC) $(TS_COPTS) $(TS_C_IMPLS) -o libtsim

all: $(TARGETS)
