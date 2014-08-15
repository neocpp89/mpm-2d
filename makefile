# Project: mpm_2d.
# automatically parallelize build
# MAKEFLAGS+="-j 4"

# include definitions such as CC, CXX
include ./makefile.common

VIZDIR = vizsrc
MATERIALDIR = materialsrc
MPMDIR = src

SOBJ = $(MATERIAL_SRC:.c=.so)

.PHONY: clean
.PHONY: backup
.PHONY: $(VIZBIN)
.PHONY: $(BIN)
.PHONY: $(SOBJ)

# Default target(s)
all: $(BIN) $(VIZBIN) $(SOBJ) $(MATERIAL_TEST_BIN)

clean:
	make -f makefile -C $(VIZDIR) clean
	make -f makefile -C $(MATERIALDIR) clean
	make -f makefile -C $(MPMDIR) clean

$(MATERIAL_TEST_BIN): 
	make -f makefile -C $(MATERIALDIR) ../$@

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

$(BIN):
	make -f makefile -C $(MPMDIR) ../$@

$(VIZBIN):
	make -f makefile -C $(VIZDIR)

materials: $(SOBJ)

$(SOBJ):
	make -f makefile -C $(MATERIALDIR) ../$@
