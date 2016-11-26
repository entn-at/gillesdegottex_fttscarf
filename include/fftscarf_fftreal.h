#ifndef __FFTSCARF_FFTREAL_H__
#define __FFTSCARF_FFTREAL_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <sstream>

#include <FFTReal/FFTReal.h>

namespace fftscarf {

template<typename _FloatType>
class FFTPlanFFTRealTemplate : public FFTPlanImplementation
{
public:
    typedef _FloatType FloatType;

private:
    ffft::FFTReal<FloatType> *m_fftreal_fft;
    std::vector<FloatType> m_signal;
    FloatType* m_fftreal_spec;
    FFTSCARF_PLAN_ACCESS_DECLARE

public:
    static std::string version() {
        return std::string("2.11"); // This is the current built-in version
    }
    static std::string libraryName() {
        std::stringstream result;
        result << "FFTReal " << version() << " (precision " << 8*sizeof(FloatType) << "b)";
        return result.str();
    }

    FFTPlanFFTRealTemplate(bool forward=true)
        : FFTPlanImplementation(forward)
    {
        m_fftreal_spec = NULL;
        m_fftreal_fft = NULL;
    }
    FFTPlanFFTRealTemplate(int n, bool forward=true)
        : FFTPlanImplementation(n, forward)
    {
        m_fftreal_spec = NULL;
        m_fftreal_fft = NULL;

        resize(n);
    }
    virtual void resize(int n)
    {
        assert(n>0);

        if(n==m_size) return;

        FFTSCARF_PLAN_ACCESS_LOCK
        m_size = n;

        delete[] m_fftreal_spec;
        m_fftreal_spec = NULL;

        m_fftreal_fft = new ffft::FFTReal<FloatType>(m_size);
        m_fftreal_spec = new FloatType[m_size];

        if(m_forward){
            m_signal.resize(m_size);
        }
        else{
            m_signal.resize((m_size%2==1)?(m_size-1)/2+1:m_size/2+1);
        }
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    template<typename TypeInContainer, typename TypeOutContainer>
    void dft(const TypeInContainer& in, TypeOutContainer& out, int dftlen=-1){
        if (!m_forward)
            throw std::string("A backward IDFT FFTPlan cannot compute the forward DFT");

        if(dftlen>0 && m_size!=dftlen)
            resize(dftlen);

        dftlen = m_size;

        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        if(int(out.size())!=neededoutsize)
            out.resize(neededoutsize);

        FFTSCARF_PLAN_ACCESS_LOCK
//         if(in.size()==m_size)
//             m_fftreal_fft->do_fft(m_fftreal_spec, &(in[0]));
//         else{
            int u = 0;
            if(m_signal.size()!=m_size)
                m_signal.resize(m_size);
            for(; u<int(in.size()); ++u)
                m_signal[u] = in[u];
            for(; u<dftlen; ++u)
                m_signal[u] = 0.0;
            m_fftreal_fft->do_fft(m_fftreal_spec, &(m_signal[0]));
//         }
        FFTSCARF_PLAN_ACCESS_UNLOCK

        out[0] = m_fftreal_spec[0]; // DC
        // TODO manage odd size
        for(int f=1; f<m_size/2; f++)
            out[f] = make_complex(m_fftreal_spec[f], -m_fftreal_spec[m_size/2+f]);
        out[m_size/2] = m_fftreal_spec[m_size/2]; // Nyquist
    }

    template<typename TypeInContainer, typename TypeOutContainer>
    void idft(const TypeInContainer& in, TypeOutContainer& out, int winlen=-1){
        if(m_forward)
            throw std::string("A forward DFT FFTPlan cannot compute the backward IDFT");

        int dftlen = int((in.size()-1)*2);
        if(m_size!=dftlen)
            resize(dftlen);

        if(winlen==-1)
            winlen = m_size;

        if(int(out.size())!=winlen)
            out.resize(winlen);


        // TODO manage odd size
        FloatType oneoverdftlen = FloatType(1.0)/m_size;
        m_fftreal_spec[0] = in[0].real()*oneoverdftlen; // DC
        for(int f=1; f<m_size/2; f++){
            m_fftreal_spec[f] = in[f].real()*oneoverdftlen;
            m_fftreal_spec[m_size/2+f] = -in[f].imag()*oneoverdftlen;
        }
        m_fftreal_spec[m_size/2] = in[m_size/2].real()*oneoverdftlen; // Nyquist

        FFTSCARF_PLAN_ACCESS_LOCK
//         if(winlen==m_size)
//             m_fftreal_fft->do_ifft(m_fftreal_spec, &(out[0])); // IDFT
//         else{
            if(m_signal.size()!=m_size)
                m_signal.resize(m_size);

            m_fftreal_fft->do_ifft(m_fftreal_spec, &(m_signal[0])); // IDFT

            for(int i=0; i<winlen; i++)
                out[i] = m_signal[i];
//         }
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    ~FFTPlanFFTRealTemplate()
    {
        if(m_fftreal_fft)	delete m_fftreal_fft;
        if(m_fftreal_spec)	delete[] m_fftreal_spec;
    }
};

#ifdef FFTSCARF_PRECISION_SINGLE
    typedef FFTPlanFFTRealTemplate<float> FFTPlanSingleFFTReal;
    #ifndef FFTSCARF_FFTPLANSINGLE
        #define FFTSCARF_FFTPLANSINGLE
        typedef FFTPlanSingleFFTReal FFTPlanSingle;
    #endif
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
    typedef FFTPlanFFTRealTemplate<double> FFTPlanDoubleFFTReal;
    #ifndef FFTSCARF_FFTPLANDOUBLE
        #define FFTSCARF_FFTPLANDOUBLE
        typedef FFTPlanDoubleFFTReal FFTPlanDouble;
    #endif
#endif
#ifdef FFTSCARF_PRECISION_LONGDOUBLE
    typedef FFTPlanFFTRealTemplate<long double> FFTPlanLongDoubleFFTReal;
    #ifndef FFTSCARF_FFTPLANLONGDOUBLE
        #define FFTSCARF_FFTPLANLONGDOUBLE
        typedef FFTPlanLongDoubleFFTReal FFTPlanLongDouble;
    #endif
#endif

#if FFTSCARF_PRECISION_DEFAULT == 32
    typedef FFTPlanSingleFFTReal FFTPlanFFTReal;
#elif FFTSCARF_PRECISION_DEFAULT == 64
    typedef FFTPlanDoubleFFTReal FFTPlanFFTReal;
#elif FFTSCARF_PRECISION_DEFAULT == 128
    typedef FFTPlanLongDoubleFFTReal FFTPlanFFTReal;
#endif

#ifndef FFTSCARF_FFTPLAN
    #if FFTSCARF_PRECISION_DEFAULT == 32
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanSingleFFTReal FFTPlan;
    #elif FFTSCARF_PRECISION_DEFAULT == 64
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanDoubleFFTReal FFTPlan;
    #elif FFTSCARF_PRECISION_DEFAULT == 128
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanLongDoubleFFTReal FFTPlan;
    #endif
#endif
}

#endif // __FFTSCARF_FFTREAL_H__
