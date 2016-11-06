#ifndef __FFTSCARF_FFT_PFFFT_H__
#define __FFTSCARF_FFT_PFFFT_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>

#include <iostream>

#include <fftscarf.h>
#ifndef FFTScarf
#define FFTScarf FFTPlanImplementationPFFFT
#endif


#include "../fftlibs/pffft/pffft.h"

namespace fftscarf {

class FFTPlanImplementationPFFFT : public FFTPlanImplementation
{
    PFFFT_Setup *m_setup;

    float *m_input;
    float *m_output;
    float *m_work;

public:
    static std::string version();
    static std::string libraryName();

    FFTPlanImplementationPFFFT(bool forward=true);
    FFTPlanImplementationPFFFT(int n, bool forward=true);
    ~FFTPlanImplementationPFFFT();

    virtual void resize(int n);

    template<typename TypeInContainer>
    void dft(const TypeInContainer& in, std::vector<std::complex<FFFLOAT> >& out, int dftlen=-1){
        if (!m_forward)
            throw std::string("A backward IDFT FFTPlan cannot compute the forward DFT");

        if(dftlen>0 && m_size!=dftlen)
            resize(dftlen);

        dftlen = m_size;

        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        if(int(out.size())!=neededoutsize)
            out.resize(neededoutsize);

        size_t u = 0;
        for(; u<in.size(); ++u)
            m_input[u] = in[u];
        for(; u<dftlen; ++u)
            m_input[u] = 0.0;

        pffft_transform_ordered(m_setup, m_input, m_output, m_work, PFFFT_FORWARD);

        out[0] = m_output[0]; // DC
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_output[2*f], m_output[2*f+1]);
        out[m_size/2] = m_output[1]; // Nyquist
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

        if(int(out.size())!=winlen)
            out.resize(winlen);

        m_input[0] = in[0].real(); // DC
        m_input[1] = in[m_size/2].real(); // Nyquist
        for(int f=1; f<m_size/2; f++){
            m_input[2*f] = in[f].real();
            m_input[2*f+1] = in[f].imag();
        }

        pffft_transform_ordered(m_setup, m_input, m_output, m_work, PFFFT_BACKWARD);

        FFFLOAT oneoverdftlen = 1.0/m_size;
        for(size_t u=0; u<winlen; ++u)
            out[u] = m_output[u]*oneoverdftlen;
    }
};

}

#endif // __FFTSCARF_FFT_PFFFT_H__
