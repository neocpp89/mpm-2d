# Project: material models for use with mpm_2d.

# include definitions such as CC, CXX
include ../makefile.common

# override CFLAGS to make shared objects
CFLAGS = -std=gnu99 -g -O3 -Wall -march=native -fPIC -shared -nostartfiles -I../src

# Need to fix this; test each material separately
MATERIAL_TEST_SRC = \
    ../src/tensor.c \
	test_material.c \
	g_local_mu2.c

REL_SOBJ = $(addprefix ../,$(MATERIAL_SRC:.c=.so))
REL_MATERIAL_TEST_BIN = $(addprefix ../,$(MATERIAL_TEST_BIN))

.PHONY: clean

# Default target(s)
all: $(REL_SOBJ) $(REL_MATERIAL_TEST_BIN)

clean:
	-rm $(REL_SOBJ) $(REL_MATERIAL_TEST_BIN);

$(REL_MATERIAL_TEST_BIN): $(MATERIAL_TEST_SRC)
	$(CC) -march=native -I../src -g $^ -o $@ -lm
	./$@

../%.so : %.c
	$(CC) $(CFLAGS) $< -o $@ 

