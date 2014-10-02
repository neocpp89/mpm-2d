# Project: boundary conditions for use with mpm_2d.

# include definitions such as CC, CXX
include ../makefile.common

# override CFLAGS to make shared objects
CFLAGS = -std=gnu99 -g -O3 -Wall -march=native -fPIC -shared -nostartfiles -I../src

REL_SOBJ = $(addprefix ../,$(BC_SRC:.c=.so))

.PHONY: clean

# Default target(s)
all: $(REL_SOBJ)
 
clean:
	-rm $(REL_SOBJ);

../%.so : %.c
	$(CC) $(CFLAGS) $< -o $@ 

