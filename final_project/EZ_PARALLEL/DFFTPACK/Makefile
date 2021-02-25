#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the makefile for compiling the DFFTPACK library. It contains two 
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

LIB = DFFTPACK

# Command-Line System options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran (MPI not required)
FC = gfortran

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Compiler flags.
  FCFLAGS = -O3 -funroll-loops -fexpensive-optimizations -std=f95 -w
  # Include file search path, e.g., C:/lib/mpi/include.
  INCPTHS = 
  # Library search path, e.g., C:/lib/mpi/lib.
  INCLIBS = 
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = 

else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -O3 -funroll-loops -fexpensive-optimizations -std=f95 -w
  # Include file search paths, e.g., /usr/include/.
  INCPTHS = 
  # Library search path, e.g., ../../EZ_PARALLEL.
  INCLIBS =  
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = 
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# Souce file names.
SRCS = $(wildcard *.f)
# Object file names.
OBJS = $(patsubst %.f, %.o, $(SRCS))

# Compilation rules.
# Archive all .o files to library.
lib$(LIB).a: $(OBJS)
	ar -rcs lib$(LIB).a $(OBJS)

# Compile all .o files.
%.o: %.f
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
ifeq ($(CMD), Windows)
	del *.o *.a
else ifeq ($(CMD), Linux)
	rm -rf *.o *.a
endif

# Tidy deletes all .o files, keeping the .a file.
tidy:
ifeq ($(CMD), Windows)
	del *.o 
else ifeq ($(CMD), Linux)
	rm -rf *.o 
endif
