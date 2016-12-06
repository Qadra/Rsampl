.PHONY: help copy

# Utility files for handling generic sampling
SAMPLING_FILES += sample.c sample.h

						  # Weighted Random Sampling #

# Files implementing the RSTree weighted random sampling functions
SAMPLING_FILES += rstree.c rstree.h

# Files implementing binary heap weighted random sampling functions
SAMPLING_FILES += exact.c exact.h

# Files supporting the S+R Method
SAMPLING_FILES += walkers.c walkers.h exact_walkers.c exact_walkers.h

							  # Subset Sampling #

SAMPLING_FILES += sss.c sss.h uss.c uss.h uss_base.c uss_base.h

help:
	echo "Help incoming"

# Copy files from subsetsampling into src/ direction for building.
copy:
	cp $(addprefix subsetsampling/src/, $(SAMPLING_FILES)) src/
