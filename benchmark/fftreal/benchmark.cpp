#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanImplementationFFTReal>("fftreal");

    return 0;
}
