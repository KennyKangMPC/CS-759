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

# Path to EZDFFTPACK library.
EZDFFTPTH = .$(SEP)ezdfftpack

# Command-Line System Options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran.
FC = gfortran

# Executable name.
TRGT = two_level_qg_solver_serial.exe
# .f90 files.
SRCS = $(wildcard *.f90)
# .o files.
OBJS = $(patsubst %.f90, %.o, $(SRCS))
# .mod files.
MODS = $(patsubst %.f90, %.mod, $(SRCS))

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Compiler flags.
  FCFLAGS = -O3 -Wall
  # Include file search path, e.g., C:\lib\mpi\include.
  INCPTHS = $(EZDFFTPTH) 
  # Library search path, e.g., C:\lib\mpi\lib.
  INCLIBS = $(EZDFFTPTH)
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = EZDFFTPACK 

else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -O3 -Wall
  # Include file search paths, e.g., /usr/include/.
  INCPTHS = $(EZDFFTPTH)
  # Library search path, e.g., ../../sEZ_PARALLEL.
  INCLIBS = $(EZDFFTPTH)
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = EZDFFTPACK
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# Links objects to the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $(OBJS) -o $(TRGT) $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Make main.o.
main.o: main.f90 time_stepper.o output.o initialize.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile time_stepper.o.
time_stepper.o: time_stepper.f90 output.o initialize.o \
jacobian_ekman_shear_solve.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile jacobian_ekman_shear_solve.o.
jacobian_ekman_shear_solve.o: jacobian_ekman_shear_solve.f90 initialize.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile output.o.
output.o: output.f90 initialize.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile all independent .f90 files. 
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
