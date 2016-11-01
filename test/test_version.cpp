#include <iostream>
using namespace std;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestRangedArray
#include <boost/test/unit_test.hpp>

#include <fftscarf.h>

BOOST_AUTO_TEST_CASE( test_version )
{
    std::cout << "Testing FFScarf compilation..." << std::endl;

    // Print some info related to versions and compilation
    fftscarf::write_compile_info(cout);
}
