#include <ctime>
#include <iostream>
#include <iterator>
#include <deque>
#include <fstream>
#include <sstream>
using namespace std;

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestFFTLibs
#include <boost/test/unit_test.hpp>

#include <fftscarf.h>
using namespace fftscarf;

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "../benchmark/stream.h"

template<typename FFTPlanType>
static void test_lib(){

    std::cout << "Testing " << FFTPlanType::libraryName() << " ..." << std::endl;

    int N = 1024;
    double accthresh = 100*fftscarf::eps<typename FFTPlanType::FloatType>();

    std::cout << "    N=" << N << " accuracy threshold=" << accthresh << std::endl;

//    fftscarf::FFTPlanImplementationOoura fft_ref(true); // TODO Need a reference for both types!
    FFTPlanType fft(true);
    FFTPlanType ifft(false);

    // Prepare the random number generator
    boost::mt19937 rnd_engine(std::time(0));
    boost::normal_distribution<typename FFTPlanType::FloatType> rnd_normal_distrib;
    boost::variate_generator<boost::mt19937&, 
        boost::normal_distribution<typename FFTPlanType::FloatType> > generator(rnd_engine, rnd_normal_distrib);

    std::vector<std::complex<typename FFTPlanType::FloatType> > spec;
    std::vector<std::complex<typename FFTPlanType::FloatType> > spec_ref;
    std::vector<typename FFTPlanType::FloatType> inframe, outframe;

    // Fill an input frame
    inframe.resize(N);
    for(size_t n=0; n<N; ++n)
        inframe[n] = generator();

    // Run the "reference" implementation
//    fft_ref.dft(inframe, spec_ref, N);

    // Run the tested implementation
    fft.dft(inframe, spec, N);


//    // Verify: sig->spec == specref
//    std::cout << "    ref implementation for sepc err: " << fft_ref.libraryName() << std::endl;
//    double spec_err = 0.0;
//    for(size_t i=0; i<spec.size(); ++i)
//        spec_err += std::abs(spec_ref[i]-spec[i])*std::abs(spec_ref[i]-spec[i]);
//    spec_err = sqrt(spec_err/spec_ref.size());
//    std::cout << "    spec err=" << spec_err << std::endl;
//    if(spec_err>accthresh){
//        std::cout << "spec_err=" << spec_err << " accuracy threshold=" << accthresh << std::endl;
//        std::cout << "spec_ref=" << spec_ref << endl;
//        std::cout << "spec=" << spec << endl;
//    }
//    assert(spec_err<accthresh);


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
    std::cout << "    sig err=" << sig_rerr << std::endl;
    assert(sig_rerr<accthresh);
}

#ifdef FFTSCARF_FFT_OOURA
BOOST_AUTO_TEST_CASE( test_fftlibs_ooura )
{
    test_lib<fftscarf::FFTPlanOoura>();
}
#endif

#ifdef FFTSCARF_FFT_FFTREAL
#ifdef FFTSCARF_PRECISION_SINGLE
BOOST_AUTO_TEST_CASE( test_fftlibs_fftreal_single )
{
    test_lib<fftscarf::FFTPlanSingleFFTReal>();
}
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
BOOST_AUTO_TEST_CASE( test_fftlibs_fftreal_double )
{
    test_lib<fftscarf::FFTPlanDoubleFFTReal>();
}
#endif
// #8 The precision is currently not that of long double
//BOOST_AUTO_TEST_CASE( test_fftlibs_fftreal_long_double )
//{
//    test_lib<fftscarf::FFTPlanImplementationFFTReal<long double> >();
//}
#endif

#ifdef FFTSCARF_FFT_PFFFT
#ifdef FFTSCARF_PRECISION_SINGLE
BOOST_AUTO_TEST_CASE( test_fftlibs_pffft )
{
    test_lib<fftscarf::FFTPlanPFFFT>();
}
#endif
#endif

#ifdef FFTSCARF_FFT_FFTS
#ifdef FFTSCARF_PRECISION_SINGLE
BOOST_AUTO_TEST_CASE( test_fftlibs_ffts )
{
    test_lib<fftscarf::FFTPlanFFTS>();
}
#endif
#endif

#ifdef FFTSCARF_FFT_FFTW3
#ifdef FFTSCARF_PRECISION_SINGLE
BOOST_AUTO_TEST_CASE( test_fftlibs_fftw3_single )
{
    test_lib<fftscarf::FFTPlanSingleFFTW3>();
}
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
BOOST_AUTO_TEST_CASE( test_fftlibs_fftw3_double )
{
    test_lib<fftscarf::FFTPlanDoubleFFTW3>();
}
#endif
#endif

#ifdef FFTSCARF_FFT_IPP
#ifdef FFTSCARF_PRECISION_SINGLE
BOOST_AUTO_TEST_CASE( test_fftlibs_ipp_single)
{
    test_lib<fftscarf::FFTPlanSingleIPP>();
}
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
BOOST_AUTO_TEST_CASE( test_fftlibs_ipp_double)
{
    test_lib<fftscarf::FFTPlanDoubleIPP>();
}
#endif
#endif
