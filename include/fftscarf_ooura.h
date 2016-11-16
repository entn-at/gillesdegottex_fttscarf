#ifndef __FFTSCARF_OOURA_H__
#define __FFTSCARF_OOURA_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>

#include <fftscarf.h>
#ifdef FFTSCARF_PLAN_PROTECTACCESS
#include <boost/thread/mutex.hpp>
#endif

extern "C" {
#include "../fftlibs/ooura/fftsg.h"
}

namespace fftscarf {

class FFTPlanOoura : public FFTPlanImplementation
{
public:
    typedef OOFLOAT FloatType;

private:
    FloatType* m_ooura_a;
    int* m_ooura_ip;
    FloatType* m_ooura_w;
    FFTSCARF_PLAN_ACCESS_DECLARE

public:
    static std::string version();
    static std::string libraryName();

    FFTPlanOoura(bool forward=true);
    FFTPlanOoura(int n, bool forward=true);
    ~FFTPlanOoura();

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
        for(; u<in.size(); ++u)
            m_ooura_a[u] = in[u];
        for(; u<dftlen; ++u)
            m_ooura_a[u] = 0.0;

        FFTSCARF_PLAN_ACCESS_LOCK
        rdft(m_size, 1, m_ooura_a, m_ooura_ip, m_ooura_w);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        out[0] = m_ooura_a[0]; // DC
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_ooura_a[2*f], -m_ooura_a[2*f+1]);
        out[m_size/2] = m_ooura_a[1]; // Nyquist
    }

    template<typename TypeOutContainer>
    void idft(const std::vector<std::complex<FloatType> >& in, TypeOutContainer& out, int winlen=-1){
        if(m_forward)
            throw std::string("A forward DFT FFTPlan cannot compute the backward IDFT");

        int dftlen = int((in.size()-1)*2);
        if(m_size!=dftlen)
            resize(dftlen);

        if(winlen==-1)
            winlen = m_size;
        
        if(int(out.size())!=winlen)
            out.resize(winlen);

        m_ooura_a[0] = in[0].real(); // DC
        for(int f=1; f<m_size/2; f++){
            m_ooura_a[2*f] = in[f].real();
            m_ooura_a[2*f+1] = -in[f].imag();
        }
        m_ooura_a[1] = in[m_size/2].real(); // Nyquist

        FFTSCARF_PLAN_ACCESS_LOCK
        rdft(m_size, -1, m_ooura_a, m_ooura_ip, m_ooura_w);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        FloatType oneoverdftlen = 2.0/m_size;
        for(int u=0; u<winlen; ++u)
            out[u] = m_ooura_a[u]*oneoverdftlen;
    }
};

#ifdef FFTSCARF_PRECISION_DEFAULTSINGLE
    typedef FFTPlanOoura FFTPlanSingleOoura;
    #ifndef FFTSCARF_FFTPLANSINGLE
        #define FFTSCARF_FFTPLANSINGLE
        typedef FFTPlanSingleOoura FFTPlanSingle;
    #endif
#else
    typedef FFTPlanOoura FFTPlanDoubleOoura;
    #ifndef FFTSCARF_FFTPLANDOUBLE
        #define FFTSCARF_FFTPLANDOUBLE
        typedef FFTPlanDoubleOoura FFTPlanDouble;
    #endif
#endif
#ifndef FFTSCARF_FFTPLAN
    #define FFTSCARF_FFTPLAN
    typedef FFTPlanOoura FFTPlan;
#endif
}

#endif // __FFTSCARF_OOURA_H__
