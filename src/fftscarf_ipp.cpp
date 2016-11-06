#include "fftscarf_ipp.h"

#include <cmath>
#include <string>
using namespace std;

#include <boost/lexical_cast.hpp>

static bool s_ipp_initialized = false;

namespace fftscarf {

std::string FFTPlanImplementationIPP::version(){
    const IppLibraryVersion *lib = ippGetLibVersion();
    return string(lib->Name)+" "+string(lib->Version); // This is the current used version
}
std::string FFTPlanImplementationIPP::libraryName(){
    return string("Intel Integrated Performance Primitives (Intel IPP) ")+version();
}

void FFTPlanImplementationIPP::setTimeLimitForPlanPreparation(float t){
    throw string("FFTPlanImplementationPFFFT::setTimeLimitForPlanPreparation: Does nothing for FFTReal implementation");
}

void FFTPlanImplementationIPP::init(){
    if(!s_ipp_initialized){
        ippInit();
        s_ipp_initialized = true;
    }
    m_pFFTSpec = NULL;
    m_pFFTSpecBuf=NULL;
    m_pFFTWorkBuf=NULL;
    m_pSrc=NULL;
    m_pDst=NULL;
}

FFTPlanImplementationIPP::FFTPlanImplementationIPP(bool forward)
    : FFTPlanImplementation(forward)
{
    init();
}
FFTPlanImplementationIPP::FFTPlanImplementationIPP(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    init();

    resize(n);
}

void FFTPlanImplementationIPP::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    m_size = n;

    const int order=(int)(log((double)m_size)/log(2.0));

    if(m_pSrc)  ippFree(m_pSrc);
    m_pSrc = ippsMalloc_32f(m_size);
    if(m_pDst)  ippFree(m_pDst);
    m_pDst = ippsMalloc_32f(m_size);

    // Query to get buffer sizes
    int sizeFFTSpecBuf;
    int sizeFFTInitBuf;
    int sizeFFTWorkBuf;
    ippsFFTGetSize_R_32f(order, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeFFTSpecBuf, &sizeFFTInitBuf, &sizeFFTWorkBuf);

    // Alloc FFT buffers
    Ipp8u *pFFTInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
    if(m_pFFTSpecBuf)   ippFree(m_pFFTSpecBuf);
    m_pFFTSpecBuf = ippsMalloc_8u(sizeFFTSpecBuf);
    if(m_pFFTWorkBuf)   ippFree(m_pFFTWorkBuf);
    m_pFFTWorkBuf = ippsMalloc_8u(sizeFFTWorkBuf);

    // Initialize FFT
    ippsFFTInit_R_32f(&m_pFFTSpec, order, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, m_pFFTSpecBuf, pFFTInitBuf);

    if(pFFTInitBuf) ippFree(pFFTInitBuf);
}

void FFTPlanImplementationIPP::dft() {
    if (m_forward) // DFT
        throw std::string("FFTPlanImplementationIPP::dft(): Not implemented.");
    else
        throw std::string("FFTPlanImplementationIPP::dft(): DFT using iDFT plan is not implemented.");
}
void FFTPlanImplementationIPP::idft() {
    throw string("FFTPlanImplementationIPP::idft(): Not implemented.");
}

FFTPlanImplementationIPP::~FFTPlanImplementationIPP() {
    if(m_pSrc)  ippFree(m_pSrc);
    if(m_pDst)  ippFree(m_pDst);
    if(m_pFFTSpecBuf)   ippFree(m_pFFTSpecBuf);
    if(m_pFFTWorkBuf)   ippFree(m_pFFTWorkBuf);
}

}
