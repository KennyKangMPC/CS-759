#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    This is a makefile that will give instructions on how to link the other
#    files in this folder.
#    Please see http://genius2k.is-programmer.com/posts/40301.html for a 
#    tutorial on makefiles. In short, each stanza is of the form
#
#    TARGET: PREREQUISITES
#    	COMMAND
#
#    The -c commands compile without linking every source file. The -o commands
#    link all of the object files.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filename separator character.
SEP = /

# Path to EZ_PARALLEL module.
EZPPTH = ..$(SEP)..$(SEP)EZ_PARALLEL

# Path to DFFTPACK
DFFTPTH = ..$(SEP)..$(SEP)EZ_PARALLEL$(SEP)DFFTPACK

# Command-Line System Options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran (with msmpi).
FC = gfortran

# Executable name.
TRGT = execute_scheme_ifft_unit_test.exe
# .f90 files.
SRCS = $(wildcard *.f90)
# .o files.
OBJS = $(patsubst %.f90, %.o, $(SRCS)) $(EZPPTH)$(SEP)ez_parallel.o $(EZPPTH)$(SEP)ez_parallel_structs.o
# .mod files.
MODS = $(patsubst %.f90, %.mod, $(SRCS))

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Path to msmpi include.
  MPIPTHINC = $(MSMPI_INC)
  # Path to msmpi lib.
  MPIPTHLIB = $(MSMPI_LIB64)
  # Compiler flags.
  FCFLAGS = -g
  # Include file search path, e.g., C:\lib\mpi\include.
  INCPTHS = $(EZPPTH) $(MPIPTHINC) $(DFFTPTH)
  # Library search path, e.g., C:\lib\mpi\lib.
  INCLIBS =  $(MPIPTHLIB) $(DFFTPTH)
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = msmpi dfftpack

else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -g
  # Include file search paths, e.g., \usr\include\.
  INCPTHS = $(EZPPTH) $(DFFTPTH)
  # Library search path, e.g., ..\..\EZ_PARALLEL.
  INCLIBS = $(DFFTPTH)
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = DFFTPACK
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# Makes the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $(OBJS) -o $(TRGT) $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile all .f90 files. 
%.o: %.f90
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)


.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
ifeq ($(CMD), Windows)
	del *.o *.mod *.exe
else ifeq ($(CMD), Linux)
	rm -rf *.o *.mod *.exe
endif

# Tidy deletes all .o and .mod files, keeping the executable.
tidy:
ifeq ($(CMD), Windows)
	del *.o *.mod
else ifeq ($(CMD), Linux)
	rm -rf *.o *.mod
endif
