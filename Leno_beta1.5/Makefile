#Hand coded Makefile

# Source files. Change this to the name of the correct file.
SRC = main.cpp

#Object file. This is the name of the .o file as a result of compiling $(SRC)
OBJ = main.o

# Target. This is the name of the ELF executable.
target = Leno

# Source directory
SRCDIR = ./src/

# include directory
INCLUDEDIR = ./include/

# Bin directory
BINDIR = ./bin/

#==================================================================================================#

# Compiler details
CC = g++ 
# CFLAGS = -std=c++11 -O2 # add -O2, delete -g
CFLAGS = -g -std=c++11
# The location of include files of various libraries
# Only work with ".lib/...", not "lib/...", why?
# On some computer, may only work without "."
ARMA_INCLUDE = .lib/armadillo-6.100.0/include 
SL_INCLUDE = lib/SuperLU_4.3/SRC
GP_INCLUDE = lib/getpot-c++
MP_INCLUDE = lib/matplotlib-cpp
CU_INCLUDE = lib/cubature-1.0.2

INCLUDES = -I$(ARMA_INCLUDE) -I $(SL_INCLUDE) -I$(GP_INCLUDE) -I$(MP_INCLUDE) -I$(CU_INCLUDE) -I$(SRCDIR) -I$(INCLUDEDIR)

# The location of libraries
SL_LIB = lib/SuperLU_4.3/lib -lsuperlu_4.3
ARMA_LIB = lib/armadillo-6.100.0 -larmadillo
BLAS_LIB = -lblas 
LAPACK_LIB = -llapack
PYTHON_LIB = -lpython2.7

MYLIBS = -L $(ARMA_LIB) -L $(SL_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(PYTHON_LIB)

# The location of external .o file
CU_O = lib/cubature-1.0.2/pcubature.o

#==================================================================================================#

all : $(target)

$(target): $(OBJ)  MathFunct.o Params.o Material.o Device1D.o Device2D.o Poisson1D.o Poisson2D.o Transport.o Tunnelling.o ExtractData.o InOut.o InOut2D.o Run1D.o Run2D.o Capacitance.o GateEfficiency.o $(CU_O)
	$(CC) $(CFLAGS) $(OBJ) MathFunct.o Params.o Material.o Device1D.o Device2D.o Poisson1D.o Poisson2D.o Transport.o Tunnelling.o ExtractData.o InOut.o InOut2D.o Run1D.o Run2D.o Capacitance.o GateEfficiency.o $(CU_O) -o $(target)  $(MYLIBS)
	
$(OBJ) : $(SRCDIR)$(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)$(SRC)

MathFunct.o : $(SRCDIR)MathFunct.cpp   	 $(INCLUDEDIR)MathFunct.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)MathFunct.cpp
Params.o : $(SRCDIR)Params.cpp	       	 $(INCLUDEDIR)Params.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Params.cpp
Material.o : $(SRCDIR)Material.cpp     	 $(INCLUDEDIR)Material.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Material.cpp
Device1D.o : $(SRCDIR)Device1D.cpp     	 $(INCLUDEDIR)Device1D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Device1D.cpp
Device2D.o : $(SRCDIR)Device2D.cpp     	 $(INCLUDEDIR)Device2D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Device2D.cpp
Poisson1D.o : $(SRCDIR)Poisson1D.cpp   	 $(INCLUDEDIR)Poisson1D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Poisson1D.cpp
Poisson2D.o : $(SRCDIR)Poisson2D.cpp   	 $(INCLUDEDIR)Poisson2D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Poisson2D.cpp
Transport.o : $(SRCDIR)Transport.cpp   	 $(INCLUDEDIR)Transport.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Transport.cpp
Tunnelling.o : $(SRCDIR)Tunnelling.cpp 	 $(INCLUDEDIR)Tunnelling.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Tunnelling.cpp
ExtractData.o : $(SRCDIR)ExtractData.cpp $(INCLUDEDIR)ExtractData.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)ExtractData.cpp
InOut.o :     $(SRCDIR)InOut.cpp         $(INCLUDEDIR)InOut.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)InOut.cpp
InOut2D.o :     $(SRCDIR)InOut2D.cpp     $(INCLUDEDIR)InOut2D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)InOut2D.cpp
Run1D.o :     $(SRCDIR)Run1D.cpp         $(INCLUDEDIR)Run1D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Run1D.cpp
Run2D.o :     $(SRCDIR)Run2D.cpp         $(INCLUDEDIR)Run2D.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Run2D.cpp
Capacitance.o : $(SRCDIR)Capacitance.cpp $(INCLUDEDIR)Capacitance.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)Capacitance.cpp
GateEfficiency.o : $(SRCDIR)GateEfficiency.cpp $(INCLUDEDIR)GateEfficiency.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)GateEfficiency.cpp

.phony : clean

clean : 
	rm *.o 


