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
    std::cout << __LINE__ << std::endl;
    #ifdef FFTSCARF_PRECISION_SINGLE
        FFTPlanManager<FFTPlanSingle> pm_single;
        std::cout << __LINE__ << std::endl;
        pm_single.prepare_osf(5, 12);
        std::cout << __LINE__ << std::endl;
        pm_single.prepare(16384);
        std::cout << __LINE__ << std::endl;
        pm_single.prepare(8192);
        std::cout << __LINE__ << std::endl;
        pm_single.write_status(std::cout);
        std::cout << __LINE__ << std::endl;
    #endif
    #ifdef FFTSCARF_PRECISION_DOUBLE
        std::cout << __LINE__ << std::endl;
        FFTPlanManager<FFTPlanDouble> pm_double;
        std::cout << __LINE__ << std::endl;
        pm_double.prepare_osf(5, 8);
        std::cout << __LINE__ << std::endl;
        pm_double.write_status(std::cout);
        std::cout << __LINE__ << std::endl;
    #endif

    int p_min = 5;
    int p_max = 20;
    std::cout << __LINE__ << std::endl;
    FFTPlanManager<FFTPlan> pm_default_fwd(true);
    std::cout << __LINE__ << std::endl;
    FFTPlanManager<FFTPlan> pm_default_bck(false);
    std::cout << __LINE__ << std::endl;
    pm_default_fwd.write_status(std::cout);
    pm_default_bck.write_status(std::cout);

    std::cout << __LINE__ << std::endl;
    boost::mt19937 rnd_engine((uint32_t)std::time(0));
    boost::normal_distribution<typename FFTPlan::FloatType> rnd_normal_distrib;
    boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<typename FFTPlan::FloatType> > generator(rnd_engine, rnd_normal_distrib);

    std::cout << __LINE__ << std::endl;
    std::cout << "testing: " << std::flush;

    std::vector<typename FFTPlan::FloatType> inframe;
    std::cout << __LINE__ << std::endl;
    for(int p=p_min; p<=p_max; ++p){
        int dftlen = std::pow(2, p);
        std::cout << dftlen << " " << std::flush;

        // Fill an input frame
        std::cout << __LINE__ << std::endl;
        inframe.resize(dftlen);
        std::cout << __LINE__ << std::endl;
        for(int n=0; n<dftlen; ++n)
            inframe[n] = generator();

        std::cout << __LINE__ << std::endl;
        // Run the tested implementation
        std::cout << __LINE__ << std::endl;
        std::vector<std::complex<typename FFTPlan::FloatType> > spec;
        std::cout << __LINE__ << std::endl;
        pm_default_fwd.fft(inframe, spec, dftlen);
        std::cout << __LINE__ << std::endl;
        pm_default_bck.ifft(spec, inframe);
        std::cout << __LINE__ << std::endl;
    }
    std::cout << __LINE__ << std::endl;
    std::cout << std::endl;
    pm_default_fwd.write_status(std::cout);
    pm_default_bck.write_status(std::cout);


//    BOOST_CHECK(!fftscarf::isPow2(p)); // The reminder should be non-power of 2^a
}
