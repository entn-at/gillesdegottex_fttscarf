#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanFFTW3>("fftw3");

    return 0;
}
