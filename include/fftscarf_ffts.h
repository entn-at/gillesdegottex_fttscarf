#ifndef __FFTSCARF_FFTS_H__
#define __FFTSCARF_FFTS_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <fftscarf.h>

#include "../fftlibs/ffts/include/ffts.h"
#include "../fftlibs/ffts/src/ffts_attributes.h"

namespace fftscarf {

class FFTPlanFFTS : public FFTPlanImplementation
{
public:
    typedef float FloatType;

private:
    FloatType FFTS_ALIGN(32) *m_signal;
    FloatType FFTS_ALIGN(32) *m_spec;
    ffts_plan_t *m_p;
    FFTSCARF_PLAN_ACCESS_DECLARE

public:
    static std::string version();
    static std::string libraryName();

    FFTPlanFFTS(bool forward=true);
    FFTPlanFFTS(int n, bool forward=true);
    ~FFTPlanFFTS();

    virtual void resize(int n);

    template<typename TypeInContainer>
    void dft(const TypeInContainer& in, std::vector<std::complex<FloatType> >& out, int dftlen=-1){
        if (!m_forward)
            throw std::string("A backward IDFT FFTPlan cannot compute the forward DFT");

        if(dftlen>0 && m_size!=dftlen)
            resize(dftlen);

        dftlen = m_size;

        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        if(int(out.size())!=neededoutsize)
            out.resize(neededoutsize);

        int u = 0;
        for(; u<int(in.size()); ++u)
            m_signal[u] = in[u];
        for(; u<dftlen; ++u)
            m_signal[u] = 0.0;

        FFTSCARF_PLAN_ACCESS_LOCK
        ffts_execute(m_p, m_signal, m_spec);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        for(int f=0; f<m_size/2+1; f++)
            out[f] = make_complex(m_spec[2*f], m_spec[2*f+1]);
    }

    template<typename TypeOutContainer>
    void idft(const std::vector<std::complex<FloatType> >& in, TypeOutContainer& out, int winlen=-1){
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

        FFTSCARF_PLAN_ACCESS_LOCK
        ffts_execute(m_p, m_spec, m_signal);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        FloatType oneoverdftlen = FloatType(1.0)/m_size;
        for(int u=0; u<winlen; ++u)
            out[u] = m_signal[u]*oneoverdftlen;
    }
};

typedef FFTPlanFFTS FFTPlanSingleFFTS;
#ifndef FFTSCARF_FFTPLANSINGLE
    #define FFTSCARF_FFTPLANSINGLE
    typedef FFTPlanSingleFFTS FFTPlanSingle;
#endif
#ifndef FFTSCARF_FFTPLAN
    #define FFTSCARF_FFTPLAN
    typedef FFTPlanFFTS FFTPlan;
#endif
}

#endif // __FFTSCARF_FFTS_H__
