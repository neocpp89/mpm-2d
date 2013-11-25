# Project: mpm_2d.
CC = gcc
BIN = mpm_2d
VIZBIN = viz
CFLAGS = -c -march=native -std=gnu99 -O3 -Wall -Wstrict-prototypes -pedantic -g -funroll-loops -I/usr/lib/openmpi/include/ -rdynamic
#takes form '-ldl -lpthreads' etc.
LIB = -lrt -lm -lhdf5 -pthread -lblas -llapack -lcxsparse -lconfuse -ldl
# modified CFLAGS for libraries
LDFLAGS = $(LIB) -g
# CFLAGS = -o $(BIN) -I<dir> -L<dir> -Wall -Wstrict-prototypes -ansi -pedantic
SRC = \
	main.c \
	bc.c \
	loading.c \
	element.c \
	particle.c \
	process.c \
	interpolate.c \
	visualization.c \
	reader.c \
	writer.c \
	rtsafe.c \
	material.c \
	map.c

MATERIAL_SRC = \
	builtin_material.c \
	dp_ri.c

OBJ = $(SRC:.c=.o)

.PHONY: clean
.PHONY: backup
all: $(BIN) $(VIZBIN)

clean:
	-rm $(OBJ) $(BIN) $(VIZBIN);

backup-nogit:
	tar --exclude-vcs --exclude=jobs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

backup:
	tar --exclude=jobs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

distbackup: clean
	make backup-nogit

doc: Doxyfile
	doxygen

Doxyfile:
	doxygen -g

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(VIZBIN): viz.cpp
	g++ -Wall $< -o $@ -lm -lSDL -lGL -lGLU -lftgl -lpng -lSDL_image -lconfuse

.c.o:
	$(CC) $(CFLAGS) $< -o $@
