#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanFFTReal>("fftreal");

    return 0;
}
