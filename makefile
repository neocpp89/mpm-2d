# automatically parallelize build
# MAKEFLAGS+="-j 4"

# Project: mpm_2d.
CC = gcc
CXX = g++
BIN = mpm_2d
VIZBIN = viz
CFLAGS = -fno-omit-frame-pointer -c -march=native -std=gnu99 -O3 -Wall -Wstrict-prototypes -pedantic -g -I/usr/lib/openmpi/include/ -rdynamic
#CFLAGS = -fno-omit-frame-pointer -c -std=gnu99 -O3 -Wall -Wstrict-prototypes -pedantic -g -funroll-loops -I/usr/lib/openmpi/include/ -rdynamic
#takes form '-ldl -lpthread' etc.
LIB = -lrt -lm -pthread -lcxsparse -lconfuse -ldl
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
	tensor.c \
	map.c

MATERIAL_SRC = \
	dp_ri.c \
	isolin.c \
	g_local_mu2.c

# Need to fix this; test each material separately
MATERIAL_TEST_SRC = \
	test_material.c \
	g_local_mu2.c

MATERIAL_TEST_BIN = material-test

OBJ = $(SRC:.c=.o)
SOBJ = $(MATERIAL_SRC:.c=.so)

.PHONY: clean
.PHONY: backup

# Default target(s)
all: $(BIN) $(VIZBIN) $(SOBJ)

clean:
	-rm $(SOBJ) $(OBJ) $(BIN) $(VIZBIN) $(MATERIAL_TEST_BIN);

$(MATERIAL_TEST_BIN): $(MATERIAL_TEST_SRC)
	$(CC) -march=native -g $^ -o $@ -lm
	./$@

backup-nogit:
	tar --exclude-vcs --exclude=jobs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

backup:
	tar --exclude=jobs --exclude=figs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

distbackup: clean
	make backup-nogit

doc: Doxyfile
	doxygen

Doxyfile:
	doxygen -g

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(VIZBIN): viz.cpp viz_colormap.cpp viz_reader.cpp
	$(CXX) -std=c++11 -g -O3 -march=native -I/usr/include/freetype2/ -Wall $^ -o $@ -lm -lSDL -lGL -lGLU -lftgl -lpng -lSDL_image -lconfuse

%.o : %.c
	$(CC) $(CFLAGS) $< -o $@

%.so : %.c
	$(CC) -g -O3 -Wall -march=native -fPIC -shared -nostartfiles $< -o $@
