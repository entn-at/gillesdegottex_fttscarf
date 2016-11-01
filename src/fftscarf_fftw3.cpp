#include <fftscarf_fftw3.h>

#include <cmath>
using namespace std;

namespace fftscarf {

std::string FFTPlanImplementationFFTW3::version(){
    return string("3");
}
std::string FFTPlanImplementationFFTW3::libraryName(){
    return string("FFTW ")+version(); // Used interface's version
}

void FFTPlanImplementationFFTW3::setTimeLimitForPlanPreparation(float t){
    #ifdef FFTW3RESIZINGMAXTIMESPENT
        fftwg_set_timelimit(t); // From FFTW 3.1 only, no means exists to check version at compile time ...
    #else
        // Do nothing
    #endif
}

// To protect the access to the FFT and external variables
//boost::signals2::mutex FFTPlan::m_fftw3_planner_access;

FFTPlanImplementationFFTW3::FFTPlanImplementationFFTW3(bool forward)
{
    m_size = -1;
    m_forward = forward;

    m_fftw3_plan = NULL;
    m_fftw3_sig = NULL;
    m_fftw3_spec = NULL;
//    #ifdef FFTW3RESIZINGMAXTIMESPENT
//        fftwg_set_timelimit(1.0); // From FFTW 3.1. No means exist to check version at compile time ...
//    #endif
}
FFTPlanImplementationFFTW3::FFTPlanImplementationFFTW3(int n, bool forward)
{
    m_size = -1;
    m_forward = forward;

    m_fftw3_plan = NULL;
    m_fftw3_sig = NULL;
    m_fftw3_spec = NULL;

    resize(n);
}
void FFTPlanImplementationFFTW3::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    m_size = n;

//         m_fftw3_planner_access.lock();
    if(m_fftw3_sig)
        fftwg_free((void*)m_fftw3_sig);
    m_fftw3_sig = NULL;
    if(m_fftw3_spec)
        fftwg_free((void*)m_fftw3_spec);
    m_fftw3_spec = NULL;

    m_fftw3_sig = (FFFLOAT*) fftwg_malloc(sizeof(FFFLOAT) * m_size);
    m_fftw3_spec = (fftwg_complex*) fftwg_malloc(sizeof(fftwg_complex) * m_size);
    //  | FFTW_PRESERVE_INPUT
//         unsigned int flags = FFTW_ESTIMATE;
    unsigned int flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT;
    // The following is likely to generate non-deterministic runs !
    // See: http://www.fftw.org/faq/section3.html#nondeterministic
    // unsigned int flags = FFTW_MEASURE;
    if(m_forward){
        m_fftw3_plan = fftwg_plan_dft_r2c_1d(m_size, m_fftw3_sig, m_fftw3_spec, flags);
//            m_fftw3_plan = fftw_plan_dft_1d(m_size, m_fftw3_sig, m_fftw3_spec, FFTW_FORWARD, FFTW_MEASURE);
    }
    else{
        m_fftw3_plan = fftwg_plan_dft_c2r_1d(m_size, m_fftw3_spec, m_fftw3_sig, flags);
//            m_fftw3_plan = fftw_plan_dft_1d(m_size, m_fftw3_sig, m_fftw3_spec, FFTW_BACKWARD, FFTW_MEASURE);
    }
//         m_fftw3_planner_access.unlock();
}

FFTPlanImplementationFFTW3::~FFTPlanImplementationFFTW3()
{
//         m_fftw3_planner_access.lock();
    if(m_fftw3_plan) fftwg_destroy_plan(m_fftw3_plan);
    if(m_fftw3_sig)  delete[] m_fftw3_sig;
    if(m_fftw3_spec) delete[] m_fftw3_spec;
//         m_fftw3_planner_access.unlock();
}

}
