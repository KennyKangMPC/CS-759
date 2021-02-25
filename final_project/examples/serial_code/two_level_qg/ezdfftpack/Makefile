# Several changes made by HCP so this would build without trouble 
# on a Linux/gfortran system.
# (1) changed step to build library to use ar instead of update
#     (update must mean something different on someone elses Unix.)
# (2) Added make clean step 
# (3) In test step, changed a.out to ./a.out for cautious folk who don't 
#     have "." in their PATH.
# (4) Change FFLAGS from -O to -O2 -funroll-loops
# (5) Specify FC=gcc in case /usr/bin/f77 is not a link to g77 
#     (as it won't be if you have f77reorder installed)
# (6) Added targets shared and installshared to make and install a shared 
#     version of the library. You need /usr/local/lib in /etc/ld.so.conf 
#     for this to work   
# (7) Modified names for dble prec version
LIB = EZDFFTPACK

# Command-Line System Options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran (with msmpi).
FC = mpifort

# .f files.
SRCS = $(wildcard *.f)
# .o files.
OBJS = $(patsubst %.f, %.o, $(SRCS))

# Compiler flags.
FC = gfortran
FCFLAGS = -O3 -funroll-loops -fexpensive-optimizations -w

# Archive all .o files into the .a file.
lib$(LIB).a: $(OBJS)
	ar -rcs lib$(LIB).a $(OBJS)

# Compile all EZDFFTPACK source code.
%.o: %.f
	$(FC) -c $< $(FCFLAGS)

.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
ifeq ($(CMD), Windows)
	del *.o *.mod *.a
else ifeq ($(CMD), Linux)
	rm -rf *.o *.mod *.a
endif

# Tidy deletes all .o and .mod files, keeping the executable.
tidy:
ifeq ($(CMD), Windows)
	del *.o *.mod
else ifeq ($(CMD), Linux)
	rm -rf *.o *.mod
endif