#ifndef __FFTSCARF_FFTS_H__
#define __FFTSCARF_FFTS_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>

#include <iostream>

#include <fftscarf.h>
#ifndef FFTScarf
#define FFTScarf FFTPlanImplementationFFTS
#endif


#include "../fftlibs/ffts/include/ffts.h"
#include "../fftlibs/ffts/src/ffts_attributes.h"

namespace fftscarf {

class FFTPlanImplementationFFTS : public FFTPlanImplementation
{
    FFFLOAT FFTS_ALIGN(32) *m_signal;
    FFFLOAT FFTS_ALIGN(32) *m_spec;

    ffts_plan_t *m_p;

public:
    static std::string version();
    static std::string libraryName();
    static void setTimeLimitForPlanPreparation(float t); // t[s]

    FFTPlanImplementationFFTS(bool forward=true);
    FFTPlanImplementationFFTS(int n, bool forward=true);
    ~FFTPlanImplementationFFTS();

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
            m_signal[u] = in[u];
        for(; u<dftlen; ++u)
            m_signal[u] = 0.0;

        ffts_execute(m_p, m_signal, m_spec);

        for(int f=0; f<m_size/2+1; f++)
            out[f] = make_complex(m_spec[2*f], m_spec[2*f+1]);
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

        for(int f=0; f<m_size/2+1; f++){
            m_spec[2*f] = in[f].real();
            m_spec[2*f+1] = in[f].imag();
        }

        ffts_execute(m_p, m_spec, m_signal);

        FFFLOAT oneoverdftlen = 1.0/m_size;
        for(size_t u=0; u<winlen; ++u)
            out[u] = m_signal[u]*oneoverdftlen;
    }

    virtual void dft();
    virtual void idft();
};

}

#endif // __FFTSCARF_FFTS_H__
