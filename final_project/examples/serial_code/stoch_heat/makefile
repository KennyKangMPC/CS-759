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

# Command-Line System Options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran (with msmpi).
FC = gfortran

# Executable name.
TRGT = stoch_heat_solver_serial.exe
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

else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -O3 -Wall
endif

# Makes the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $(OBJS) -o $(TRGT) $(FCFLAGS)

# Make main.o.
main.o: main.f90 time_stepper.o output.o initialize.o
	$(FC) -c $< $(FCFLAGS)

# Make time_stepper.o.
time_stepper.o: time_stepper.f90 output.o initialize.o
	$(FC) -c $< $(FCFLAGS)

# Make output.o.
output.o: output.f90 initialize.o
	$(FC) -c $< $(FCFLAGS)

# Compile all independent .f90 files. 
%.o: %.f90
	$(FC) -c $< $(FCFLAGS)


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
