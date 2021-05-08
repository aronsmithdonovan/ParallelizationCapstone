# Makefile for fungi simulations
# created by Aron Smith-Donovan

# compilers + flags
CXX=g++
OMP=-fopenmp
DEBUG=-DDEBUG  # show numerical DEBUG prints
COLOR=-DCOLOR  # show colorful grid in DEBUG prints (DEBUG must also be enabled)

# trng library
INCLUDE=/usr/local/include/trng
LIB=trng4

# executables
EXECUTABLES={omp.fungi,seq.fungi}

# make rules
seq.fungi: fungi-seq.cpp
	$(CXX) $(DEBUG) $(COLOR) -o seq.fungi fungi-seq.cpp -I$(INCLUDE) -l$(LIB)

omp.fungi: fungi-omp.cpp
	$(CXX) $(DEBUG) $(COLOR) ${OMP} -o omp.fungi fungi-omp.cpp -I$(INCLUDE) -l$(LIB)

clean:
	rm -f $(EXECUTABLES) *.o

all: $(EXECUTABLES)