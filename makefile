# Project: mpm_2d.
CC = gcc
BIN = mpm_2d
VIZBIN = viz
CFLAGS = -c -std=gnu99 -O3 -Wall -Wstrict-prototypes -pedantic -g -pg -fopenmp -funroll-loops -I/usr/lib/openmpi/include/
#takes form '-ldl -lpthreads' etc.
LIB = -lrt -lm -lhdf5 -pthread -lblas -llapack -lcxsparse
# modified CFLAGS for libraries
LDFLAGS = $(LIB) -pg -g
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
	material.c

OBJ = $(SRC:.c=.o)

.PHONY: clean
.PHONY: backup
all: $(BIN) $(VIZBIN)

clean:
	rm $(OBJ); \
	rm $(BIN); \
	rm $(VIZBIN);

backup:
	tar --exclude=jobs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

distbackup: clean
	backup

doc: Doxyfile
	doxygen

Doxyfile:
	doxygen -g

strip_tabs: $(SRC)
	for file in `echo $(SRC) | sed -e 's/ /\n/g'`; \
	do { sed -e 's/{\t/{       /g' < $$file | \
	sed -e 's/\t/        /g' > $$file.notabs; }; \
	done;

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(VIZBIN): viz.cpp
	g++ -Wall $< -o $@ -lm -lSDL -lGL -lGLU -lftgl -lpng -lSDL_image

.c.o:
	$(CC) $(CFLAGS) $< -o $@
