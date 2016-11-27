#ifndef __FFTSCARF_IPP_H__
#define __FFTSCARF_IPP_H__

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <ipp.h>

namespace fftscarf {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern bool FFTPlanIPP__s_ipp_initialized;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

template<typename _FloatType, typename IppsFFTSpec_R_NNf, typename IppNNf, IppNNf* (*ippsMalloc_NNf)(int len), IppStatus (*ippsFFTGetSize_R_NNf)(int, int, IppHintAlgorithm, int*, int*, int*), IppStatus (*ippsFFTInit_R_NNf)(IppsFFTSpec_R_NNf**, int, int, IppHintAlgorithm, Ipp8u*, Ipp8u*),  IppStatus (*ippsFFTFwd_RToPerm_NNf)(const IppNNf*, IppNNf*, const IppsFFTSpec_R_NNf*, Ipp8u*), IppStatus (*ippsFFTInv_PermToR_NNf)(const IppNNf*, IppNNf*, const IppsFFTSpec_R_NNf*, Ipp8u*) >
class FFTPlanIPPTemplate : public FFTPlanImplementation
{
public:
    typedef _FloatType FloatType;

private:
    IppsFFTSpec_R_NNf *m_pFFTSpec;
    Ipp8u *m_pFFTSpecBuf;
    Ipp8u *m_pFFTWorkBuf;
    IppNNf *m_pSrc;
    IppNNf *m_pDst;
    FFTSCARF_PLAN_ACCESS_DECLARE

    void init(){
        if(!FFTPlanIPP__s_ipp_initialized){
            ippInit();
            FFTPlanIPP__s_ipp_initialized = true;
        }
        m_pFFTSpec = NULL;
        m_pFFTSpecBuf=NULL;
        m_pFFTWorkBuf=NULL;
        m_pSrc=NULL;
        m_pDst=NULL;
    }

public:
    static std::string version(){
        const IppLibraryVersion *lib = ippGetLibVersion();
        return std::string(lib->Name)+" "+std::string(lib->Version); // This is the current used version
    }
    static std::string libraryName(){
        std::stringstream result;
        result << std::string("Intel IPP ") << version() << " (precision " << 8*sizeof(FloatType) << "b)";
        return result.str();
    }

    FFTPlanIPPTemplate(bool forward=true)
        : FFTPlanImplementation(forward)
    {
        init();
    }
    FFTPlanIPPTemplate(int n, bool forward=true)
        : FFTPlanImplementation(n, forward)
    {
        init();

        resize(n);
    }

    void resize(int n)
    {
        if(n==m_size) return;

        assert(n>0);
        assert(isPow2(n));

        FFTSCARF_PLAN_ACCESS_LOCK
        m_size = n;

        const int order=(int)(log((double)m_size)/log(2.0));

        if(m_pSrc)  ippFree(m_pSrc);
        if(m_pDst)  ippFree(m_pDst);
        m_pSrc = ippsMalloc_NNf(m_size);
        m_pDst = ippsMalloc_NNf(m_size);

        // Query to get buffer sizes
        int sizeFFTSpecBuf;
        int sizeFFTInitBuf;
        int sizeFFTWorkBuf;
        ippsFFTGetSize_R_NNf(order, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeFFTSpecBuf, &sizeFFTInitBuf, &sizeFFTWorkBuf);

        // Alloc FFT buffers
        Ipp8u *pFFTInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
        if(m_pFFTSpecBuf)   ippFree(m_pFFTSpecBuf);
        m_pFFTSpecBuf = ippsMalloc_8u(sizeFFTSpecBuf);
        if(m_pFFTWorkBuf)   ippFree(m_pFFTWorkBuf);
        m_pFFTWorkBuf = ippsMalloc_8u(sizeFFTWorkBuf);

        // Initialize FFT
        ippsFFTInit_R_NNf(&m_pFFTSpec, order, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, m_pFFTSpecBuf, pFFTInitBuf);

        if(pFFTInitBuf) ippFree(pFFTInitBuf);

        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    template<typename TypeInContainer, typename TypeOutContainer>
    void fft(const TypeInContainer& in, TypeOutContainer& out, int dftlen=-1){
        if (!m_forward)
            throw std::string("A backward IDFT FFTPlan cannot compute the forward DFT");

        if(dftlen>0 && m_size!=dftlen)
            resize(dftlen);
        dftlen = m_size;

        // Fill pSrc
        size_t u = 0;
        for(; u<in.size(); ++u)
            m_pSrc[u] = in[u];
        for(; u<m_size; ++u)
            m_pSrc[u] = 0.0;

        // Do the FFT
        FFTSCARF_PLAN_ACCESS_LOCK
        ippsFFTFwd_RToPerm_NNf(m_pSrc, m_pDst, m_pFFTSpec, m_pFFTWorkBuf);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        // Unpack pDst
        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        if(int(out.size())!=neededoutsize)
            out.resize(neededoutsize);

        out[0] = m_pDst[0]; // DC
        out[m_size/2] = m_pDst[1]; // Nyquist
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_pDst[2*f], m_pDst[2*f+1]);
    }

    template<typename TypeInContainer, typename TypeOutContainer>
    void ifft(const TypeInContainer& in, TypeOutContainer& out, int winlen=-1){
        if(m_forward)
            throw std::string("A forward DFT FFTPlan cannot compute the backward IDFT");

        int dftlen = (in.size()-1)*2;
        if(m_size!=dftlen)
            resize(dftlen);

        if(winlen==-1)
            winlen = m_size;

        m_pSrc[0] = in[0].real(); // DC
        m_pSrc[1] = in[m_size/2].real(); // Nyquist
        for(int f=1; f<m_size/2; f++){
            m_pSrc[2*f] = in[f].real();
            m_pSrc[2*f+1] = in[f].imag();
        }

        // Do the inverse FFT
        FFTSCARF_PLAN_ACCESS_LOCK
        ippsFFTInv_PermToR_NNf(m_pSrc, m_pDst, m_pFFTSpec, m_pFFTWorkBuf);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        // Unpack pDst
        if(int(out.size())!=winlen)
            out.resize(winlen);

        for(size_t u=0; u<winlen; ++u)
            out[u] = m_pDst[u];
    }

    ~FFTPlanIPPTemplate() {
        FFTSCARF_PLAN_ACCESS_LOCK
        if(m_pSrc)  ippFree(m_pSrc);
        if(m_pDst)  ippFree(m_pDst);
        if(m_pFFTSpecBuf)   ippFree(m_pFFTSpecBuf);
        if(m_pFFTWorkBuf)   ippFree(m_pFFTWorkBuf);
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }
};

#ifdef FFTSCARF_PRECISION_SINGLE
    typedef FFTPlanIPPTemplate<float, IppsFFTSpec_R_32f, Ipp32f, ippsMalloc_32f, ippsFFTGetSize_R_32f, ippsFFTInit_R_32f, ippsFFTFwd_RToPerm_32f, ippsFFTInv_PermToR_32f > FFTPlanSingleIPP;
    #ifndef FFTSCARF_FFTPLANSINGLE
        #define FFTSCARF_FFTPLANSINGLE
        typedef FFTPlanSingleIPP FFTPlanSingle;
    #endif
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
    typedef FFTPlanIPPTemplate<double, IppsFFTSpec_R_64f, Ipp64f, ippsMalloc_64f, ippsFFTGetSize_R_64f, ippsFFTInit_R_64f, ippsFFTFwd_RToPerm_64f, ippsFFTInv_PermToR_64f > FFTPlanDoubleIPP;
    #ifndef FFTSCARF_FFTPLANDOUBLE
        #define FFTSCARF_FFTPLANDOUBLE
        typedef FFTPlanDoubleIPP FFTPlanDouble;
    #endif
#endif
        
#if FFTSCARF_PRECISION_DEFAULT == 32
    typedef FFTPlanSingleIPP FFTPlanIPP;
#elif FFTSCARF_PRECISION_DEFAULT == 64
    typedef FFTPlanDoubleIPP FFTPlanIPP;
#endif

#ifndef FFTSCARF_FFTPLAN
    #if FFTSCARF_PRECISION_DEFAULT == 32
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanSingleIPP FFTPlan;
    #elif FFTSCARF_PRECISION_DEFAULT == 64
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanDoubleIPP FFTPlan;
    #endif
#endif
}

#endif // __FFTSCARF_IPP_H__
