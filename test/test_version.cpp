#include <iostream>
using namespace std;

#define BOOST_TEST_MODULE TestVersion
#include <boost/test/unit_test.hpp>

#include <fftscarf.h>

BOOST_AUTO_TEST_CASE( test_version )
{
    std::cout << "Testing FFScarf version and compilation information..." << std::endl;

    // Print some info related to versions and compilation
    fftscarf::write_compile_info(cout);
}
