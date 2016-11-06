#include "fftscarf_pffft.h"

#include <cmath>
#include <string>
using namespace std;

#include <boost/lexical_cast.hpp>

namespace fftscarf {

std::string FFTPlanImplementationPFFFT::version(){
    return string("2014-08-10"); // This is the current built-in version
}
std::string FFTPlanImplementationPFFFT::libraryName(){
    return string("Pretty Fast FFT (PFFFT) ")+version()+string(" (SIMD size ")+boost::lexical_cast<std::string>(pffft_simd_size())+string(")"); // This is the current built-in version
}

FFTPlanImplementationPFFFT::FFTPlanImplementationPFFFT(bool forward)
    : FFTPlanImplementation(forward)
{
    m_setup = NULL;
    m_input = NULL;
    m_output = NULL;
    m_work = NULL;
}
FFTPlanImplementationPFFFT::FFTPlanImplementationPFFFT(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_setup = NULL;
    m_input = NULL;
    m_output = NULL;
    m_work = NULL;

    resize(n);
}

void FFTPlanImplementationPFFFT::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    m_size = n;

    if(m_setup || m_input || m_output){
        pffft_destroy_setup(m_setup);
        m_setup = NULL;
        pffft_aligned_free(m_input);
        m_input = NULL;
        pffft_aligned_free(m_output);
        m_output = NULL;
    }
    if(m_work){
        pffft_aligned_free(m_work);
        m_work = NULL;
    }

    m_setup = pffft_new_setup(n, PFFFT_REAL);

    m_input = (float*)pffft_aligned_malloc(n*sizeof(float));
    m_output = (float*)pffft_aligned_malloc(n*sizeof(float));

    // Following the doc
    if(n>=16384)
        m_work = (float*)pffft_aligned_malloc(n*sizeof(float));
}

FFTPlanImplementationPFFFT::~FFTPlanImplementationPFFFT() {
    if(m_setup || m_input || m_output){
        pffft_destroy_setup(m_setup);
        m_setup = NULL;
        pffft_aligned_free(m_input);
        m_input = NULL;
        pffft_aligned_free(m_output);
        m_output = NULL;
    }
    if(m_work){
        pffft_aligned_free(m_work);
        m_work = NULL;
    }
}

}
