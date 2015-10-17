g++ -std=c++11 *.cpp lib/cubature-1.0.2/pcubature.o -o poissonSolver -O2 -larmadillo -lgfortran -L/$HOME/Documents/lib_usr -lsuperlu_4.3 -lblas -llapack -lpython2.7

