# Project: actual binary for mpm_2d.

# include definitions such as CC, CXX
include ../makefile.common

OBJ = $(SRC:.c=.o)

REL_BIN = $(addprefix ../, $(BIN))

.PHONY: clean

# Default target(s)
all: $(REL_BIN)

clean:
	-rm $(OBJ) $(REL_BIN);

$(REL_BIN): $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o $@ $(LDFLAGS)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

