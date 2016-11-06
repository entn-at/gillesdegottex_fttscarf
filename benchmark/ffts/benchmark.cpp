#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanImplementationFFTS>("ffts");

    return 0;
}
