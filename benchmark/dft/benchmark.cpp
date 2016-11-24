#include "../benchmark.h"

int main(int argc, char** argcc){

    benchmark<fftscarf::FFTPlanDFT>("dft", 2, 10); // Skip higher orders as it is too slow

    return 0;
}
