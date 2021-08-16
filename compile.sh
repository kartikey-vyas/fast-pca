rm myEigenFunctions.o myEigenMain 
g++ -c -o myEigenFunctions.o myEigenFunctions.cpp
g++ -o myEigenMain myEigenMain.cpp myEigenFunctions.o -lm
./myEigenMain
