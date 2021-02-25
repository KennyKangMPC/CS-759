#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the makefile for compiling the EZ_PARALLEL module. It contains two 
# sets of options: Command-Line System options and Fortran Compiler (FC) 
# options. The Command-Line System options determine the 'make clean' and 
# 'make tidy' commands that will be executed, (e.g., 'rm' for Linux and 'del'
# for Windows. The FC options determine which compiler will be used, and
# it is assumed that msmpi will be used with gfortran. File paths used in 
# this makefile are either standard installation paths.
#
# After running 'make', the archive file libEZ_PARALLEL.a produced can be moved
# to other directories to be used with other programs.
#
# Created by - Jason Turner
# Last Edited - April 2nd, 2020
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MOD = EZ_PARALLEL

# Command-Line System options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran.
FC = gfortran

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Path to msmpi include.
  MPIPTHINC = $(MSMPI_INC)
  # Path to msmpi lib.
  MPIPTHLIB = $(MSMPI_LIB64)
  # Compiler flags.
  FCFLAGS = -g -Wall
  # Include file search path, e.g., C:/lib/mpi/include.
  INCPTHS = $(MPIPTHINC) ./DFFTPACK
  # Library search path, e.g., C:/lib/mpi/lib.
  INCLIBS = $(MPIPTHLIB) ./DFFTPACK
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = msmpi DFFTPACK

else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = 
  # Include file search paths, e.g., /usr/include/.
  INCPTHS = ./DFFTPACK
  # Library search path, e.g., ../../EZ_PARALLEL.
  INCLIBS = ./DFFTPACK
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = DFFTPACK
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# Souce file names.
EZPARSRCS = $(wildcard *.f90)
# EZ_PARALLEL object file names.
EZPAROBJS = $(patsubst %.f90, %.o, $(EZPARSRCS))

# Compilation rules.
# Compile the EZ_PARALLEL source code
$(MOD).o: ez_parallel.f90 ez_parallel_structs.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

ez_parallel_structs.o: ez_parallel_structs.f90
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
ifeq ($(CMD), Windows)
	del *.o 
else ifeq ($(CMD), Linux)
	rm -rf *.o 
endif

# Tidy deletes all .o files, keeping the .a file.
tidy:
ifeq ($(CMD), Windows)
	del *.o 
else ifeq ($(CMD), Linux)
	rm -rf *.o 
endif
