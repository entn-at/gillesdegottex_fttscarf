#include <fftscarf_ffts.h>

// memalign
//#include <malloc.h>
// posix_memalign
//#include <stdlib.h>
//#include <mm.h>
//#define HAVE_SSE

#include <cmath>
#include <string>
using namespace std;

#include <boost/lexical_cast.hpp>

#include <fftscarf.h>

namespace fftscarf {

std::string FFTPlanImplementationFFTS::version(){
    return string("0.9.0"); // This is the current built-in version
}
std::string FFTPlanImplementationFFTS::libraryName(){
    return string("FFTS ")+version()+string(" (precision ")+boost::lexical_cast<std::string>(8*sizeof(FFFLOAT))+string("b)");
}

FFTPlanImplementationFFTS::FFTPlanImplementationFFTS(bool forward)
    : FFTPlanImplementation(forward)
{
    m_p = NULL;
    m_signal = NULL;
    m_spec = NULL;

    m_p = NULL;
}
FFTPlanImplementationFFTS::FFTPlanImplementationFFTS(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_p = NULL;
    m_signal = NULL;
    m_spec = NULL;

    resize(n);
}

void FFTPlanImplementationFFTS::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    m_size = n;

    if(m_signal || m_spec){
    #ifdef HAVE_SSE
        _mm_free(m_signal);
        _mm_free(m_spec);
//        free(m_signal);
//        free(m_spec);
    #else
        free(m_signal);
        free(m_spec);
    #endif
        m_signal = NULL;
        m_spec = NULL;
    }

    // See http://www.delorie.com/gnu/docs/glibc/libc_31.html
    #ifdef HAVE_SSE
        m_signal = _mm_malloc(n * sizeof(FFFLOAT), 32);
        m_spec = _mm_malloc(2*(m_size/2+1) * sizeof(FFFLOAT), 32);
//        posix_memalign((void**)&m_signal, n * sizeof(FFFLOAT), 32); // CRASHES
//        posix_memalign((void**)&m_spec, 2*(m_size/2+1) * sizeof(FFFLOAT), 32); // CRASHES
//        m_signal = (FFFLOAT FFTS_ALIGN(32) *) posix_memalign(n * sizeof(FFFLOAT), 32);
//        m_spec = (FFFLOAT FFTS_ALIGN(32) *) posix_memalign(2*(m_size/2+1) * sizeof(FFFLOAT), 32);
    #else
        m_signal = (FFFLOAT FFTS_ALIGN(32) *) valloc(n * sizeof(FFFLOAT));
        m_spec = (FFFLOAT FFTS_ALIGN(32) *) valloc((2*n+1) * sizeof(FFFLOAT));
//        m_signal = (FFFLOAT*) valloc(m_size * sizeof(FFFLOAT));
//        m_spec = (FFFLOAT*) valloc(2*(m_size/2+1) * sizeof(FFFLOAT));
    #endif

    if(m_p){
        ffts_free(m_p);
        m_p = NULL;
    }

    if(m_forward)
        m_p = ffts_init_1d_real(n, FFTS_FORWARD); // TODO sign for the backward I suppose
    else
        m_p = ffts_init_1d_real(n, FFTS_BACKWARD); // TODO sign for the backward I suppose
}

void FFTPlanImplementationFFTS::dft() {
    if (m_forward) // DFT
        throw std::string("FFTPlanImplementationFFTS::dft(): Not implemented.");
    else
        throw std::string("FFTPlanImplementationFFTS::dft(): DFT using iDFT plan is not implemented.");
}
void FFTPlanImplementationFFTS::idft() {
    throw string("FFTPlanImplementationFFTS::idft(): Not implemented.");
}

FFTPlanImplementationFFTS::~FFTPlanImplementationFFTS() {

    if(m_signal || m_spec){
    #ifdef HAVE_SSE
        _mm_free(m_signal);
        _mm_free(m_spec);
//        free(m_signal);
//        free(m_spec);
    #else
        free(m_signal);
        free(m_spec);
    #endif
        m_signal = NULL;
        m_spec = NULL;
    }

    if(m_p){
        ffts_free(m_p);
        m_p = NULL;
    }
}

}
