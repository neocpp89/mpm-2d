# Project: mpm_2d.
# automatically parallelize build
# MAKEFLAGS+="-j 4"

# include definitions such as CC, CXX
include ./makefile.common

VIZDIR = vizsrc
MATERIALDIR = materialsrc
MPMDIR = src
BCDIR = bcsrc

SOBJ = $(MATERIAL_SRC:.c=.so)

.PHONY: clean
.PHONY: backup
.PHONY: $(VIZBIN)
.PHONY: $(BIN)
.PHONY: $(SOBJ)
.PHONY: $(BCSOBJ)

# Default target(s)
all: $(BIN) $(VIZBIN) $(SOBJ) $(BCSOBJ) $(MATERIAL_TEST_BIN)

clean:
	${MAKE} -f makefile -C $(VIZDIR) clean
	${MAKE} -f makefile -C $(MATERIALDIR) clean
	${MAKE} -f makefile -C $(MPMDIR) clean
	${MAKE} -f makefile -C $(BCDIR) clean

$(MATERIAL_TEST_BIN): 
	${MAKE} -f makefile -C $(MATERIALDIR) ../$@

backup-nogit:
	tar --exclude-vcs --exclude=jobs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

backup:
	tar --exclude=jobs --exclude=figs -cvzf ../mpm_2d_`date +%Y%m%d`.tar.gz ../`basename $(CURDIR)`;

distbackup: clean
	${MAKE} backup-nogit

doc: Doxyfile
	doxygen

Doxyfile:
	doxygen -g

$(BIN):
	${MAKE} -f makefile -C $(MPMDIR) ../$@

$(VIZBIN):
	${MAKE} -f makefile -C $(VIZDIR)

materials: $(SOBJ)

boundaries: $(BCSOBJ)

$(SOBJ):
	${MAKE} -f makefile -C $(MATERIALDIR) ../$@

$(BCSOBJ):
	${MAKE} -f makefile -C $(BCDIR) ../$@
