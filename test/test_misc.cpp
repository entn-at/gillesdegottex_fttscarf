#include <iostream>
#include <limits>
using namespace std;

#define BOOST_TEST_MODULE TestMisc
#include <boost/test/unit_test.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <fftscarf.h>

int g_nb_tests = 1e5;

BOOST_AUTO_TEST_CASE( test_ispow2 )
{
    // Test power of 2
    int p=2;
    for(size_t l=1; l<30; ++l, p*=2)
        BOOST_CHECK(fftscarf::isPow2(p));

    // Test non-power of 2
    boost::random::mt19937 rng(std::time(0));
    boost::random::uniform_int_distribution<int> intrnd(1,std::numeric_limits<int>::max());
    for(size_t n=0; n<g_nb_tests; ++n){
        p = intrnd(rng);
        while(p%2) p /= 2;
        if(p>1)
            BOOST_CHECK(!fftscarf::isPow2(p));
    }
}

BOOST_AUTO_TEST_CASE( test_ispow235 )
{
    // Test power of 2^a * 3^b * 5^c
    boost::random::mt19937 rng(std::time(0));
    boost::random::uniform_int_distribution<int> powsrnd(1,30);
    int p;
    for(size_t n=0; n<g_nb_tests; ++n){
        p = std::pow(2, powsrnd(rng)) * std::pow(3, powsrnd(rng)) * std::pow(5, powsrnd(rng));
        BOOST_CHECK(fftscarf::isPow235(p));
    }

    // Test non-power of 2^a * 3^b * 5^c
    boost::random::uniform_int_distribution<int> intrnd(1,std::numeric_limits<int>::max());
    for(size_t n=0; n<g_nb_tests; ++n){
        p = intrnd(rng);
        while(p%2==0) p /= 2;
        while(p%3==0) p /= 3;
        while(p%5==0) p /= 5;
        if(p>1)
            BOOST_CHECK(!fftscarf::isPow235(p));
    }
}
