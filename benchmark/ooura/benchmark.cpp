#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanImplementationOoura>("ooura");

    return 0;
}
