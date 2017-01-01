#include <ctime>
#include <iostream>
#include <limits>
#include <iomanip>
using namespace std;

#define BOOST_TEST_MAIN
// #define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestMisc
#include <boost/test/unit_test.hpp>

#include <boost/chrono/system_clocks.hpp>
#include <boost/random.hpp>

#include <fftscarf.h>

int g_nb_tests = 100000;

BOOST_AUTO_TEST_CASE( test_ispow2 )
{
    // Test power of 2 up to 2^30
    int p=2;
    for(int l=1; l<30; ++l, p*=2)
        BOOST_CHECK(fftscarf::isPow2(p));

    // Test non-power of 2
    boost::random::mt19937 rng((const uint32_t)std::time(0));
    boost::random::uniform_int_distribution<int> intrnd(1,std::numeric_limits<int>::max());
    for(int n=0; n<g_nb_tests; ++n){
        p = intrnd(rng);
        while(p%2) p /= 2; // Remove potential powers of 2
        if(p>1)
            BOOST_CHECK(!fftscarf::isPow2(p)); // The reminder should be non-power of 2^a
    }
}

BOOST_AUTO_TEST_CASE( test_ispow235 )
{
    // Test power of 2^a * 3^b * 5^c
    boost::random::mt19937 rng((const uint32_t)std::time(0));
    boost::random::uniform_int_distribution<int> powsrnd(1,30);
    int p;
    for(int n=0; n<g_nb_tests; ++n){
        p = std::pow(2, powsrnd(rng)) * std::pow(3, powsrnd(rng)) * std::pow(5, powsrnd(rng));
        BOOST_CHECK(fftscarf::isPow235(p));
    }

    // Test non-power of 2^a * 3^b * 5^c
    boost::random::uniform_int_distribution<int> intrnd(1,std::numeric_limits<int>::max());
    for(int n=0; n<g_nb_tests; ++n){
        p = intrnd(rng);
        while(p%2==0) p /= 2; // Remove potential powers of 2
        while(p%3==0) p /= 3; // Remove potential powers of 3
        while(p%5==0) p /= 5; // Remove potential powers of 5
        if(p>1)
            BOOST_CHECK(!fftscarf::isPow235(p)); // The reminder should be non-power of 2^a * 3^b * 5^c
    }
}

// Test wrapq -----------------------------------------------------------------

template<typename FloatType>
inline FloatType refwrap(FloatType value){
    return std::arg(std::complex<FloatType>(std::cos(value),std::sin(value)));
}

template<typename ValueType>
void check_wrap(ValueType value){
    long double err = refwrap(value) - fftscarf::wrapq(value);
    if(!(std::abs(err)<100*fftscarf::eps<ValueType>())){
        std::cout << std::setprecision(std::numeric_limits<ValueType>::digits10+2);
        std::cout << "value=" << value << " refwrap=" << refwrap(value) << " wrap=" << fftscarf::wrapq(value) << " err=" << err << " eps=" << fftscarf::eps<ValueType>() << std::endl;
    }
    BOOST_CHECK(std::abs(err)<100*fftscarf::eps<ValueType>());
}

template<typename ValueType>
void check_multi_wrap(){

    std::cout << "test " << sizeof(ValueType) << std::endl;

    boost::mt19937 rnd_engine((uint32_t)std::time(0));
    boost::random::uniform_real_distribution<ValueType> phirnd(0.0, 2*fftscarf::pi); // For random phase

    check_wrap<ValueType>(0.0);
    check_wrap<ValueType>(fftscarf::pi/2);
    check_wrap<ValueType>(fftscarf::pi);
    check_wrap<ValueType>(-fftscarf::pi/2);
    check_wrap<ValueType>(-fftscarf::pi);
    check_wrap<ValueType>(2*fftscarf::pi);
    for(int N=-16; N<=16; ++N)
        check_wrap<ValueType>(phirnd(rnd_engine) + N*2*fftscarf::pi);


    // Check speed
    int Nmax = 500000;
    boost::chrono::system_clock::time_point tstart;
    boost::chrono::system_clock::time_point tend;

    tstart = boost::chrono::system_clock::now();
    for(int N=-Nmax; N<=Nmax; ++N)
        fftscarf::wrap(phirnd(rnd_engine) + N*2*fftscarf::pi);
    tend = boost::chrono::system_clock::now();
    boost::chrono::nanoseconds wrap_nanosec = tend-tstart;

    tstart = boost::chrono::system_clock::now();
    for(int N=-Nmax; N<=Nmax; ++N)
        fftscarf::wrapq(phirnd(rnd_engine) + N*2*fftscarf::pi);
    tend = boost::chrono::system_clock::now();
    boost::chrono::nanoseconds wrapq_nanosec = tend-tstart;

    std::cout << (double(wrap_nanosec.count())/wrapq_nanosec.count()) << std::endl;
}

BOOST_AUTO_TEST_CASE( test_wrap )
{
    check_multi_wrap<float>();
    check_multi_wrap<double>();
    check_multi_wrap<long double>();
}
