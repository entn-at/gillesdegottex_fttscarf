#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanImplementationIPP>("ipp");

    return 0;
}
