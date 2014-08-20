# Project: visualization program for mpm_2d output.

# include definitions such as CC, CXX. In upper level directory.
include ../makefile.common

# override CXXFLAGS
CXXFLAGS = -std=c++11 -g -O3 -march=native -I/usr/include/freetype2/ -Wall
CXXFLAGS += -I../src

# override needed libraries
LDFLAGS = -lm -lSDL -lGL -lGLU -lftgl -lpng -lSDL_image -lconfuse

SRC = \
    tokenizer.cpp \
    viz.cpp \
    viz_colormap.cpp \
    viz_reader.cpp

OBJ = $(SRC:.cpp=.o)

REL_VIZBIN = $(addprefix ../,$(VIZBIN))

.PHONY: clean

# Default target(s)
all: $(REL_VIZBIN)

clean:
	-rm $(OBJ) $(REL_VIZBIN);

$(REL_VIZBIN): $(OBJ)
	$(CXX) $(OBJ) $(CXXFLAGS) -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

