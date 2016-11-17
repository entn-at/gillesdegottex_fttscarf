#include <fftscarf_ffts.h>

// memalign
#include <malloc.h>
// posix_memalign
#include <stdlib.h>
//#include <mm.h>
//#define HAVE_SSE

#include <cmath>
#include <string>
using namespace std;

#include <fftscarf.h>

namespace fftscarf {

std::string FFTPlanFFTS::version(){
    return string("0.9.0"); // This is the current built-in version
}
std::string FFTPlanFFTS::libraryName(){
    stringstream result;
    result << "FFTS " << version() << " (precision " << 8*sizeof(FloatType) << "b)"; // This is the current built-in version
    return result.str();
}

FFTPlanFFTS::FFTPlanFFTS(bool forward)
    : FFTPlanImplementation(forward)
{
    m_p = NULL;
    m_signal = NULL;
    m_spec = NULL;

    m_p = NULL;
}
FFTPlanFFTS::FFTPlanFFTS(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_p = NULL;
    m_signal = NULL;
    m_spec = NULL;

    resize(n);
}

void FFTPlanFFTS::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    FFTSCARF_PLAN_ACCESS_LOCK
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

    #if (defined(_WIN32) || defined(WIN32))
        m_signal = (FloatType FFTS_ALIGN(32) *) _aligned_malloc(n * sizeof(FloatType), 32);
        m_spec = (FloatType FFTS_ALIGN(32) *) _aligned_malloc((2*n+1) * sizeof(FloatType), 32);
    #else
        // See http://www.delorie.com/gnu/docs/glibc/libc_31.html
        #ifdef HAVE_SSE
            m_signal = _mm_malloc(n * sizeof(FloatType), 32);
            m_spec = _mm_malloc(2*(m_size/2+1) * sizeof(FloatType), 32);
    //        posix_memalign((void**)&m_signal, n * sizeof(FloatType), 32); // CRASHES
    //        posix_memalign((void**)&m_spec, 2*(m_size/2+1) * sizeof(FloatType), 32); // CRASHES
    //        m_signal = (FloatType FFTS_ALIGN(32) *) posix_memalign(n * sizeof(FloatType), 32);
    //        m_spec = (FloatType FFTS_ALIGN(32) *) posix_memalign(2*(m_size/2+1) * sizeof(FloatType), 32);
        #else
            m_signal = (FloatType FFTS_ALIGN(32) *) valloc(n * sizeof(FloatType));
            m_spec = (FloatType FFTS_ALIGN(32) *) valloc((2*n+1) * sizeof(FloatType)); // TODO (2*n+1) ?????
    //        m_signal = (FloatType*) valloc(m_size * sizeof(FloatType));
    //        m_spec = (FloatType*) valloc(2*(m_size/2+1) * sizeof(FloatType));
        #endif
    #endif

    if(m_p){
        ffts_free(m_p);
        m_p = NULL;
    }

    if(m_forward)
        m_p = ffts_init_1d_real(n, FFTS_FORWARD); // TODO sign for the backward I suppose
    else
        m_p = ffts_init_1d_real(n, FFTS_BACKWARD); // TODO sign for the backward I suppose

    FFTSCARF_PLAN_ACCESS_UNLOCK
}

FFTPlanFFTS::~FFTPlanFFTS() {

    FFTSCARF_PLAN_ACCESS_LOCK

    if(m_signal || m_spec){
        #if (defined(_WIN32) || defined(WIN32))
            _aligned_free(m_signal);
            _aligned_free(m_spec);
        #else
            #ifdef HAVE_SSE
                _mm_free(m_signal);
                _mm_free(m_spec);
            //        free(m_signal);
            //        free(m_spec);
            #else
                free(m_signal);
                free(m_spec);
            #endif
        #endif
        m_signal = NULL;
        m_spec = NULL;
    }

    if(m_p){
        ffts_free(m_p);
        m_p = NULL;
    }

    FFTSCARF_PLAN_ACCESS_UNLOCK
}

}
