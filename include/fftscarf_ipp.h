#ifndef __FFTSCARF_IPP_H__
#define __FFTSCARF_IPP_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>

#include <iostream>

#include <fftscarf.h>
#ifndef FFTScarf
#define FFTScarf FFTPlanImplementationIPP
#endif


#include <ipp.h>

namespace fftscarf {

class FFTPlanImplementationIPP : public FFTPlanImplementation
{
    IppsFFTSpec_R_32f *m_pFFTSpec;
    Ipp8u *m_pFFTSpecBuf;
    Ipp8u *m_pFFTWorkBuf;

    Ipp32f *m_pSrc;
    Ipp32f *m_pDst;

    void init();

public:
    static std::string version();
    static std::string libraryName();
    static void setTimeLimitForPlanPreparation(float t); // t[s]

    FFTPlanImplementationIPP(bool forward=true);
    FFTPlanImplementationIPP(int n, bool forward=true);
    ~FFTPlanImplementationIPP();

    virtual void resize(int n);

    template<typename TypeInContainer>
    void dft(const TypeInContainer& in, std::vector<std::complex<FFFLOAT> >& out, int dftlen=-1){
        if (!m_forward)
            throw std::string("A backward IDFT FFTPlan cannot compute the forward DFT");

        if(dftlen>0 && m_size!=dftlen)
            resize(dftlen);
        dftlen = m_size;

        // Fill pSrc
        size_t u = 0;
        for(; u<in.size(); ++u)
            m_pSrc[u] = in[u];
        for(; u<dftlen; ++u)
            m_pSrc[u] = 0.0;

        // Do the FFT
        ippsFFTFwd_RToPerm_32f(m_pSrc, m_pDst, m_pFFTSpec, m_pFFTWorkBuf);

        // Unpack pDst
        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        if(int(out.size())!=neededoutsize)
            out.resize(neededoutsize);

        out[0] = m_pDst[0]; // DC
        out[m_size/2] = m_pDst[1]; // Nyquist
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_pDst[2*f], m_pDst[2*f+1]);
    }

    template<typename TypeOutContainer>
    void idft(const std::vector<std::complex<FFFLOAT> >& in, TypeOutContainer& out, int winlen=-1){
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
        ippsFFTInv_PermToR_32f(m_pSrc, m_pDst, m_pFFTSpec, m_pFFTWorkBuf);

        // Unpack pDst
        if(int(out.size())!=winlen)
            out.resize(winlen);

        for(size_t u=0; u<winlen; ++u)
            out[u] = m_pDst[u];
    }

    virtual void dft();
    virtual void idft();

//    // This works only for forward DFT !
//    inline void setInput(size_t n, FFFLOAT value){m_signal[n] = value;}
//    inline std::complex<FFFLOAT> getOutput(size_t n){
//        if(n==0)
//            return fftscarf::make_complex(m_fftreal_spec[0], FFFLOAT(0.0));
//        if(int(n)==m_size/2)
//            return fftscarf::make_complex(m_fftreal_spec[m_size/2], FFFLOAT(0.0));
//        return make_complex(m_fftreal_spec[n], -m_fftreal_spec[m_size/2+n]);
//    }
//    inline std::complex<FFFLOAT> getDCOutput(){return fftscarf::make_complex(m_fftreal_spec[0], FFFLOAT(0.0));} // Avoid index checking
//    inline std::complex<FFFLOAT> getMidOutput(size_t n){return fftscarf::make_complex(m_fftreal_spec[n], -m_fftreal_spec[m_size/2+n]);} // Avoid index checking
//    inline std::complex<FFFLOAT> getNyquistOutput(){return fftscarf::make_complex(m_fftreal_spec[m_size/2],FFFLOAT(0.0));}// Avoid index checking
};

}

#endif // __FFTSCARF_IPP_H__
