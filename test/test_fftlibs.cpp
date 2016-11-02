#include <iostream>
#include <iterator>
#include <deque>
#include <fstream>
#include <sstream>
using namespace std;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestFFTLibs
#include <boost/test/unit_test.hpp>

#include <fftscarf.h>
using namespace fftscarf;

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "../benchmark/stream.h"

template<typename FFTPlanType>
static void test_lib(){

    std::cout << "Testing " << FFTPlanType::libraryName() << " ..." << std::endl;

    int N = 8;
    double accthresh = 100*fftscarf::eps;

    std::cout << "    N=" << N << " accuracy threshold=" << accthresh << std::endl;

    // fftscarf::FFTPlanImplementationFFTW3 fft_ref(true); // TODO
    fftscarf::FFTPlanImplementationOoura fft_ref(true); // TODO
    FFTPlanType fft(true);
    FFTPlanType ifft(false);

    // Prepare the random number generator
    boost::mt19937 rnd_engine(std::time(0));
    boost::normal_distribution<FFFLOAT> rnd_normal_distrib;

    std::vector<std::complex<FFFLOAT> > spec;
    std::vector<std::complex<FFFLOAT> > spec_ref;
    std::vector<FFFLOAT> inframe, outframe;

    // Fill an input frame
    inframe.resize(N);
    for(size_t n=0; n<N; ++n)
        inframe[n] = rnd_normal_distrib(rnd_engine);

    std::cout << "inframe=" << inframe << endl;

    // Run the "reference" implementation
    fft_ref.dft(inframe, spec_ref, N);

    // Run the tested implementation
    fft.dft(inframe, spec, N);


    // Verify: sig->spec == specref
    double spec_err = 0.0;
    for(size_t i=0; i<spec.size(); ++i)
        spec_err += std::abs(spec_ref[i]-spec[i])*std::abs(spec_ref[i]-spec[i]);
    spec_err = sqrt(spec_err/spec_ref.size());

    std::cout << "    spec err=" << spec_err << std::endl;

//     if(spec_err>accthresh){
        std::cout << "spec_err=" << spec_err << " accuracy threshold=" << accthresh << std::endl;
        std::cout << "spec_ref=" << spec_ref << endl;
        std::cout << "spec=" << spec << endl;
//     }
    assert(spec_err<accthresh);


    // Verify: sig->spec->sig' == sig

    // First reverse the spec
    ifft.idft(spec, outframe, N);

    // Then measure relative RMS
    double sqerr = 0.0;
    double sqin = 0.0;
    for(size_t i=0; i<inframe.size(); ++i){
        sqerr += (inframe[i]-outframe[i])*(inframe[i]-outframe[i]);
        sqin += inframe[i]*inframe[i];
    }
    double sig_rerr = sqrt(sqerr/sqin);
    assert(sig_rerr<accthresh);

    std::cout << "    sig err=" << sig_rerr << std::endl;
}

#ifdef FFTSCARF_FFT_OOURA
BOOST_AUTO_TEST_CASE( test_fftlibs_ooura )
{
    test_lib<fftscarf::FFTPlanImplementationOoura>();
}
#endif

#ifdef FFTSCARF_FFT_FFTREAL
BOOST_AUTO_TEST_CASE( test_fftlibs_fftreal )
{
    test_lib<fftscarf::FFTPlanImplementationFFTReal>();
}
#endif

#ifdef FFTSCARF_FFT_FFTW3
BOOST_AUTO_TEST_CASE( test_fftlibs_fftw3 )
{
    test_lib<fftscarf::FFTPlanImplementationFFTW3>();
}
#endif
