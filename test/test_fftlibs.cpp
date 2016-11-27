#include <cmath>
#include <ctime>
#include <iostream>
#include <iterator>
#include <deque>
#include <fstream>
#include <sstream>
using namespace std;

#include <boost/random.hpp>

#define BOOST_TEST_MAIN
// #define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestFFTLibs
#include <boost/test/unit_test.hpp>

#include <fftscarf.h>
using namespace fftscarf;

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

    boost::mt19937 rnd_engine((uint32_t)std::time(0));
    std::vector<std::complex<typename FFTPlanType::FloatType> > spec;
    std::vector<std::complex<typename FFTPlanType::FloatType> > spec_ref;
    std::vector<typename FFTPlanType::FloatType> inframe, outframe;

    long double accthresh = 100*fftscarf::eps<typename FFTPlanType::FloatType>();

    #ifdef FFTSCARF_PRECISION_LONGDOUBLE
        FFTPlanLongDoubleFFTReal fft_ref(true);
    #else
        #ifdef FFTSCARF_PRECISION_DOUBLE
            FFTPlanDoubleFFTReal fft_ref(true);
        #else
            FFTPlanSingleFFTReal fft_ref(true);
        #endif
    #endif
    FFTPlanType fft(true);
    FFTPlanType ifft(false);

    // Test transforms of sinusoids --------------------------------------------
    boost::random::uniform_int_distribution<> binrnd(0,N/2); // For random frequency
    boost::random::uniform_real_distribution<typename FFTPlanType::FloatType> phirnd(0.0, 2*fftscarf::pi); // For random phase
    for(size_t b=0; b<10; ++b){
        int binref = binrnd(rnd_engine); // Frequency
        int ampref = N/2;            // Amplitude // TODO Randomize!
        if(binref==0 || binref==N/2) ampref *= 2;
        long double phiref = phirnd(rnd_engine); // Phase

        // Fill an input frame
        // Test with a simple sinusoid centered on an exact bin
        inframe.resize(N);
        for(int n=0; n<N; ++n)
            inframe[n] = cosl(binref*2*fftscarf::pi*n/((long double)N) + phiref);
        
        // Run the tested implementation
        fft.dft(inframe, spec, N);

        // Check the amplitude
        long double ampmeas = std::abs(spec[binref]);
        if(abs(ampref-ampmeas)>10*N*accthresh)
            std::cout << "    spec err=" << abs(ampref-ampmeas) << " (threshold=" << 10*N*accthresh << ")" << std::endl;
        BOOST_CHECK(abs(ampref-ampmeas)<10*N*accthresh);

        // Check the phase
        long double phimeas = std::arg(spec[binref]);
        if(abs(wrap(phiref-phimeas))>10*N*accthresh)
            std::cout << "    spec err=" << abs(wrap(phiref-phimeas)) << " (threshold=" << 10*N*accthresh << ")" << std::endl;
        BOOST_CHECK(abs(wrap(phiref-phimeas))<10*N*accthresh);

        // Check the zeros
        long double spec_err = 0.0;
        for(size_t k=0; k<N/2; ++k)
            if(k!=binref)
                spec_err += abs(spec[k])*abs(spec[k]);
         spec_err = sqrt(spec_err/spec.size());
         if(spec_err>10*N*accthresh)
             std::cout << "    spec err=" << spec_err << " (threshold=" << 10*N*accthresh << ")" << std::endl;
//         BOOST_CHECK(spec_err<10*N*accthresh);
    }

    // Test transforms of Gaussian noise ---------------------------------------
    boost::normal_distribution<typename FFTPlanType::FloatType> rnd_normal_distrib;
    boost::variate_generator<boost::mt19937&, 
        boost::normal_distribution<typename FFTPlanType::FloatType> > generator(rnd_engine, rnd_normal_distrib);

    for(size_t b=0; b<100; ++b){
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
            if(spec_err>accthresh){
                std::cout << "    spec err=" << spec_err << " (threshold=" << accthresh << "; ref implementation: " << fft_ref.libraryName() << ")" << std::endl;
                if(vm.count("specprint")){
                    std::cout << "spec_ref=" << spec_ref << endl;
                    std::cout << "spec=" << spec << endl;
                }
            }
            BOOST_CHECK(spec_err<accthresh);
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
        if(sig_rerr>accthresh)
            std::cout << "    sig err=" << sig_rerr << " (threshold=" << accthresh << ")" << std::endl;
        BOOST_CHECK(sig_rerr<accthresh);
    }
}

//#ifdef FFTSCARF_FFT_DJBFFT
//BOOST_AUTO_TEST_CASE( test_fftlibs_djbfft )
//{
//    test_lib<fftscarf::FFTPlanDJBFFT>();
//}
//#endif

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
