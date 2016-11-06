#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanImplementationPFFFT>("pffft", 5);

    return 0;
}
