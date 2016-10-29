#ifndef __FFTSCARF_FFTREAL_H__
#define __FFTSCARF_FFTREAL_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>

#include <fftscarf.h>
#ifndef FFTScarf
#define FFTScarf FFTPlanImplementationFFTReal
#endif

#include "../fftlibs/FFTReal/FFTReal.h"

namespace fftscarf {

class FFTPlanImplementationFFTReal : public FFTPlanImplementation
{
    ffft::FFTReal<FFFLOAT> *m_fftreal_fft;
    std::vector<FFFLOAT> m_signal;
    FFFLOAT* m_fftreal_spec;

public:
    static std::string version();
    static std::string libraryName();

    FFTPlanImplementationFFTReal(bool forward=true);
    FFTPlanImplementationFFTReal(int n, bool forward=true);
    ~FFTPlanImplementationFFTReal();

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

        if(in.size()==m_size)
            m_fftreal_fft->do_fft(m_fftreal_spec, &(in[0]));
        else{
            size_t u = 0;
            if(m_signal.size()!=m_size)
                m_signal.resize(m_size);
            for(; u<in.size(); ++u)
                m_signal[u] = in[u];
            for(; u<dftlen; ++u)
                m_signal[u] = 0.0;
            m_fftreal_fft->do_fft(m_fftreal_spec, &(m_signal[0]));
        }

        out[0] = m_fftreal_spec[0]; // DC
        // TODO manage odd size
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_fftreal_spec[f], -m_fftreal_spec[m_size/2+f]);
        out[m_size/2] = m_fftreal_spec[m_size/2]; // Nyquist
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
        
        // TODO SOMETHING CRASHING IN IDFT

        if(int(out.size())!=winlen)
            out.resize(winlen);


        // TODO manage odd size
        FFFLOAT oneoverdftlen = 1.0/m_size;
        m_fftreal_spec[0] = in[0].real()*oneoverdftlen; // DC
        for(int f=1; f<m_size/2; f++){
            m_fftreal_spec[f] = in[f].real()*oneoverdftlen;
            m_fftreal_spec[m_size/2+f] = -in[f].imag()*oneoverdftlen;
        }
        m_fftreal_spec[m_size/2] = in[m_size/2].real()*oneoverdftlen; // Nyquist

        if(winlen==m_size)
            m_fftreal_fft->do_ifft(m_fftreal_spec, &(out[0])); // IDFT
        else{
            if(m_signal.size()!=m_size)
                m_signal.resize(m_size);

            m_fftreal_fft->do_ifft(m_fftreal_spec, &(m_signal[0])); // IDFT

            for(size_t i=0; i<winlen; i++)
                out[i] = m_signal[i];
        }
    }
};

}

#endif // __FFTSCARF_FFTREAL_H__
