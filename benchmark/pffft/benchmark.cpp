#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanPFFFT>("pffft", 5);

    return 0;
}
