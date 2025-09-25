# THIS MAKEFILE USES THE FOLLOWING ENVIRONMENTAL VARIABLE:
#     APlusPlus
#
# The APlusPlus env variable should normally be set when you login.
# If not you can set it using:
#     setenv APlusPlus /home/henshw/software/AppPpp-0.8.3/A++/install
#
# If you type:
#     ls $APlusPlus/include
# you should see A++.h as well as other include files.

# First target (default)
all: heat1dImp

# Compiler setup
CC   = gcc
CXX  = g++
CCFLAGS = -fPIC -O3 -I$(APlusPlus)/include -I../../include

# Libraries for A++
AppLibraries = -Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lApp -lApp_static
LIBS = $(AppLibraries)

# Rule to compile .C files
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# 1D heat equation, implicit time-stepping with A++ arrays
heat1dImpFiles = heat1dImp.o tridiagonal.o
heat1dImp: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dImpFiles) $(LIBS)

clean:
	rm -f *.o heat1dImp
