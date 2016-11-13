#include "fftscarf_pffft.h"

#include <cmath>
#include <string>
using namespace std;

namespace fftscarf {

std::string FFTPlanPFFFT::version(){
    return string("2014-08-10"); // This is the current built-in version
}
std::string FFTPlanPFFFT::libraryName(){
    stringstream result;
    result << "Pretty Fast FFT (PFFFT) " << version() << " (SIMD size " << pffft_simd_size() << ")" << " (precision " << 8*sizeof(FloatType) << "b)"; // This is the current built-in version
    return result.str();
}

FFTPlanPFFFT::FFTPlanPFFFT(bool forward)
    : FFTPlanImplementation(forward)
{
    m_setup = NULL;
    m_input = NULL;
    m_output = NULL;
    m_work = NULL;
}
FFTPlanPFFFT::FFTPlanPFFFT(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_setup = NULL;
    m_input = NULL;
    m_output = NULL;
    m_work = NULL;

    resize(n);
}

void FFTPlanPFFFT::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    FFTSCARF_PLAN_ACCESS_LOCK
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

    m_input = (FloatType*)pffft_aligned_malloc(n*sizeof(FloatType));
    m_output = (FloatType*)pffft_aligned_malloc(n*sizeof(FloatType));

    // Following the PFFFT's documentation
    if(n>=16384)
        m_work = (FloatType*)pffft_aligned_malloc(n*sizeof(FloatType));

    FFTSCARF_PLAN_ACCESS_UNLOCK
}

FFTPlanPFFFT::~FFTPlanPFFFT() {
    FFTSCARF_PLAN_ACCESS_LOCK
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
    FFTSCARF_PLAN_ACCESS_UNLOCK
}

}
