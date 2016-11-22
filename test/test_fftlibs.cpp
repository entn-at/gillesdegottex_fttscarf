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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

template<typename FFTPlanType>
static void test_lib(){

    // Default argument
    int N = 4096;

    // Arguments parsing
    po::options_description desc("Options");
    desc.add_options()
        ("size", po::value<uint32_t>(), "The FFT size")
        ("specverif", "Verify the spectrum by comparison with a reference")
        ("specprint", "If verification fails, print the erroneous vectors")
    ;
    //    po::positional_options_description posopt;
    //    posopt.add("size", 1);
    po::variables_map vm;
//    po::store(po::command_line_parser(boost::unit_test::framework::master_test_suite().argc, boost::unit_test::framework::master_test_suite().argv).options(desc).positional(posopt).run(), vm);
    po::store(po::parse_command_line(boost::unit_test::framework::master_test_suite().argc, boost::unit_test::framework::master_test_suite().argv, desc), vm);
    po::notify(vm);
    if (vm.count("size"))
        N = vm["size"].as<uint32_t>();

    std::cout << "Testing " << FFTPlanType::libraryName() << "  N=" << N << " ..." << std::endl;

    long double accthresh = 100*fftscarf::eps<typename FFTPlanType::FloatType>();
    long double specaccthresh = accthresh; // TODO 10*

    #ifdef FFTSCARF_PRECISION_LONGDOUBLE
        FFTPlanLongDoubleOoura fft_ref(true);
    #else
        #ifdef FFTSCARF_PRECISION_DOUBLE
            FFTPlanDoubleOoura fft_ref(true);
        #else
            FFTPlanSingleOoura fft_ref(true);
        #endif
    #endif
    FFTPlanType fft(true);
    FFTPlanType ifft(false);

    // Prepare the random number generator
    boost::mt19937 rnd_engine((uint32_t)std::time(0));
    boost::normal_distribution<typename FFTPlanType::FloatType> rnd_normal_distrib;
    boost::variate_generator<boost::mt19937&, 
        boost::normal_distribution<typename FFTPlanType::FloatType> > generator(rnd_engine, rnd_normal_distrib);

    std::vector<std::complex<typename FFTPlanType::FloatType> > spec;
    std::vector<std::complex<typename FFTPlanType::FloatType> > spec_ref;
    std::vector<typename FFTPlanType::FloatType> inframe, outframe;

    // Fill an input frame
    inframe.resize(N);
    for(int n=0; n<N; ++n)
        inframe[n] = generator();

    // Run the tested implementation
    fft.dft(inframe, spec, N);

    if(vm.count("specverif")){
        // Run the "reference" implementation
        fft_ref.dft(inframe, spec_ref, N);

        // Verify: sig->spec == specref
        long double spec_err = 0.0;
        for(size_t i=0; i<spec.size(); ++i)
            spec_err += std::abs(spec_ref[i]-spec[i])*std::abs(spec_ref[i]-spec[i]);
        spec_err = sqrt(spec_err/spec_ref.size());
        std::cout << "    spec err=" << spec_err << " (threshold=" << specaccthresh << "; ref implementation: " << fft_ref.libraryName() << ")" << std::endl;
        if(spec_err>specaccthresh && vm.count("specprint")){
            std::cout << "spec_ref=" << spec_ref << endl;
            std::cout << "spec=" << spec << endl;
        }
        BOOST_CHECK(spec_err<specaccthresh);
    }

    // Verify: sig->spec->sig' == sig

    // First reverse the spec
    ifft.idft(spec, outframe, N);

    // Then measure relative RMS
    long double sqerr = 0.0;
    long double sqin = 0.0;
    for(size_t i=0; i<inframe.size(); ++i){
        sqerr += (inframe[i]-outframe[i])*(inframe[i]-outframe[i]);
        sqin += inframe[i]*inframe[i];
    }
    long double sig_rerr = sqrt(sqerr/sqin);
    std::cout << "    sig err=" << sig_rerr << " (threshold=" << accthresh << ")" << std::endl;
    BOOST_CHECK(sig_rerr<accthresh);
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
    #ifdef FFTSCARF_PRECISION_LONGDOUBLE
    BOOST_AUTO_TEST_CASE( test_fftlibs_fftreal_longdouble )
    {
        test_lib<fftscarf::FFTPlanLongDoubleFFTReal>();
    }
    #endif
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
    #ifdef FFTSCARF_PRECISION_LONGDOUBLE
    BOOST_AUTO_TEST_CASE( test_fftlibs_fftw3_longdouble )
    {
        test_lib<fftscarf::FFTPlanLongDoubleFFTW3>();
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

#ifdef FFTSCARF_FFT_DFT
    #ifdef FFTSCARF_PRECISION_SINGLE
    BOOST_AUTO_TEST_CASE( test_fftlibs_dft_single)
    {
        test_lib<fftscarf::FFTPlanSingleDFT>();
    }
    #endif
    #ifdef FFTSCARF_PRECISION_DOUBLE
    BOOST_AUTO_TEST_CASE( test_fftlibs_dft_double)
    {
        test_lib<fftscarf::FFTPlanDoubleDFT>();
    }
    #endif
    #ifdef FFTSCARF_PRECISION_LONGDOUBLE
    BOOST_AUTO_TEST_CASE( test_fftlibs_dft_longdouble)
    {
        test_lib<fftscarf::FFTPlanLongDoubleDFT>();
    }
    #endif
#endif
