SHELL=/bin/sh
CLASS=A
VERSION=
SFILE=config/suite.def

#default: header
#	@ sys/print_instructions

all: BT SP FT

BT: bt
bt: header
	cd BT; $(MAKE) CLASS=$(CLASS) VERSION=$(VERSION)

SP: sp		       
sp: header	       
	cd SP; $(MAKE) CLASS=$(CLASS)

FT: ft		       
ft: header	       
	cd FT; $(MAKE) CLASS=$(CLASS)

# Awk script courtesy cmg@cray.com, modified by Haoqiang Jin
suite:
	@ awk -f sys/suite.awk SMAKE=$(MAKE) $(SFILE) | $(SHELL)


# It would be nice to make clean in each subdirectory (the targets
# are defined) but on a really clean system this will won't work
# because those makefiles need config/make.def
clean:
	- rm -f core 
	- rm -f *~ */core */*~ */*.o */npbparams.h */*.obj */*.exe
	- rm -f sys/setparams sys/makesuite sys/setparams.h
	- rm -rf */rii_files
	- rm -f bin/sp.* bin/lu.* bin/mg.* bin/ft.* bin/bt.* bin/is.*
	- rm -f bin/ep.* bin/cg.* bin/ua.* bin/dc.*

veryclean: clean
	- rm -f bin/sp.* bin/lu.* bin/mg.* bin/ft.* bin/bt.* bin/is.*
	- rm -f bin/ep.* bin/cg.* bin/ua.* bin/dc.*

header:
	#@ sys/print_header



