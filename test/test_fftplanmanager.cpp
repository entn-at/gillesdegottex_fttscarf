#include <ctime>
#include <iostream>
#include <limits>
using namespace std;

#define BOOST_TEST_MAIN
// #define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestPlanManager
#include <boost/test/unit_test.hpp>

#include <boost/random.hpp>

#include <fftscarf.h>
using namespace fftscarf;

BOOST_AUTO_TEST_CASE( test_fftplanmanager )
{
    int p_min = 5;
    int p_max = 16; // TODO Going above 13 crashes with FFTS on Travis (!?!?!)

    #ifdef FFTSCARF_PRECISION_SINGLE
        FFTPlanManager<FFTPlanSingle> pm_single;
        pm_single.write_status(std::cout);
        pm_single.prepare_osf(5, 11);
        pm_single.prepare(8192);
        pm_single.prepare(4096);
        pm_single.write_status(std::cout);
    #endif
    #ifdef FFTSCARF_PRECISION_DOUBLE
        FFTPlanManager<FFTPlanDouble> pm_double;
        pm_double.write_status(std::cout);
        pm_double.prepare_osf(p_min, p_max);
        pm_double.write_status(std::cout);
    #endif

    FFTPlanManager<FFTPlan> pm_default_fwd(true);
    FFTPlanManager<FFTPlan> pm_default_bck(false);
    pm_default_fwd.write_status(std::cout);
    pm_default_bck.write_status(std::cout);

    boost::mt19937 rnd_engine((uint32_t)std::time(0));
    boost::normal_distribution<typename FFTPlan::FloatType> rnd_normal_distrib;
    boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<typename FFTPlan::FloatType> > generator(rnd_engine, rnd_normal_distrib);

    std::cout << "testing: " << std::flush;

    std::vector<typename FFTPlan::FloatType> inframe;
    for(int p=p_min; p<=p_max; ++p){
        int dftlen = std::pow(2, p);
        std::cout << dftlen << " " << std::flush;

        // Fill an input frame
        inframe.resize(dftlen);
        for(int n=0; n<dftlen; ++n)
            inframe[n] = generator();

        // Run the tested implementation
        std::vector<std::complex<typename FFTPlan::FloatType> > spec;
        pm_default_fwd.fft(inframe, spec, dftlen);
        pm_default_bck.ifft(spec, inframe);
    }
    std::cout << std::endl;
    pm_default_fwd.write_status(std::cout);
    pm_default_bck.write_status(std::cout);


//    BOOST_CHECK(!fftscarf::isPow2(p)); // The reminder should be non-power of 2^a
}
