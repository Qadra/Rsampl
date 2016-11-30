.PHONY: help copy

# Base files that includes utility functions
SAMPLING_FILES := xutil.c xutil.h

# Utility files for handling generic sampling
SAMPLING_FILES += sample.c sample.h

# Files implementing the RSTree weighted random sampling functions
SAMPLING_FILES += rstree.c rstree.h

# Files implementing binary heap weighted random sampling functions
SAMPLING_FILES += exact.c exact.h

# Files supporting the S+R Method
SAMPLING_FILES += walkers.c walkers.h exact_walkers.c exact_walkers.h

help:
	echo "Help incoming"

# Copy files from subsetsampling into src/ direction for building.
copy:
	cp $(addprefix subsetsampling/src/, $(SAMPLING_FILES)) src/
