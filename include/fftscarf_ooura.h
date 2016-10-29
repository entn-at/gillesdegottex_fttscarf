#ifndef __FFTSCARF_OOURA_H__
#define __FFTSCARF_OOURA_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>

#include <fftscarf.h>
#ifndef FFTScarf
#define FFTScarf FFTPlanImplementationOoura
#endif

extern "C" {
#include "../fftlibs/ooura/fftsg.h"
}


namespace fftscarf {

class FFTPlanImplementationOoura : public FFTPlanImplementation
{
    FFFLOAT* m_ooura_a;
    int* m_ooura_ip;
    FFFLOAT* m_ooura_w;

public:
    static std::string version();
    static std::string libraryName();

    FFTPlanImplementationOoura(bool forward=true);
    FFTPlanImplementationOoura(int n, bool forward=true);
    ~FFTPlanImplementationOoura();

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
            m_ooura_a[u] = in[u];
        for(; u<dftlen; ++u)
            m_ooura_a[u] = 0.0;

        rdft(m_size, 1, m_ooura_a, m_ooura_ip, m_ooura_w);

        out[0] = m_ooura_a[0]; // DC
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_ooura_a[2*f], -m_ooura_a[2*f+1]);
        out[m_size/2] = m_ooura_a[1]; // Nyquist
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

        m_ooura_a[0] = in[0].real(); // DC
        for(int f=1; f<m_size/2; f++){
            m_ooura_a[2*f] = in[f].real();
            m_ooura_a[2*f+1] = -in[f].imag();
        }
        m_ooura_a[1] = in[m_size/2].real(); // Nyquist

        rdft(m_size, -1, m_ooura_a, m_ooura_ip, m_ooura_w);

        FFFLOAT oneoverdftlen = 2.0/m_size;
        for(size_t u=0; u<winlen; ++u)
            out[u] = m_ooura_a[u]*oneoverdftlen;
    }
};

}

#endif // __FFTSCARF_OOURA_H__
