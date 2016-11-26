#ifndef __FFTSCARF_DFT_H__
#define __FFTSCARF_DFT_H__

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <iostream>

namespace fftscarf {

//# define M_PIl          3.141592653589793238462643383279502884L /* pi */

template<typename _FloatType>
class FFTPlanDFTTemplate : public FFTPlanImplementation
{
public:
    typedef _FloatType FloatType;

private:
    std::vector<std::vector<std::complex<FloatType> > > m_FF;
    FFTSCARF_PLAN_ACCESS_DECLARE

public:
    static std::string version(){
        return std::string("1.0");
    }
    static std::string libraryName(){
        std::stringstream result;
        result << std::string("Standard DFT ") << version() << " (precision " << 8*sizeof(FloatType) << "b)";
        return result.str();
    }

    FFTPlanDFTTemplate(bool forward=true) {
        m_size = -1;
        m_forward = forward;
    }
    FFTPlanDFTTemplate(int n, bool forward=true) {
        m_size = -1;
        m_forward = forward;
        resize(n);
    }

    void resize(int n) {
        assert(n>0);
        if(n==m_size) return;

        FFTSCARF_PLAN_ACCESS_LOCK
        m_size = n;

        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        m_FF.resize(neededoutsize);
        for(size_t k=0; k<m_FF.size(); ++k){
            m_FF[k].resize(n);
            for(size_t n=0; n<m_FF[k].size(); ++n){
                long double phi = -2.0L*M_PIl*k*n/m_size;
                m_FF[k][n].real() = std::cos(phi);
                m_FF[k][n].imag() = std::sin(phi);
            }
        }
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    template<typename TypeInContainer, typename TypeOutContainer>
    void dft(const TypeInContainer& in, TypeOutContainer& out, int dftlen=-1) {
        if (!m_forward)
            throw std::string("A backward IDFT FFTPlan cannot compute the forward DFT");

        if(dftlen>0 && m_size!=dftlen)
            resize(dftlen);
        dftlen = m_size;

        FFTSCARF_PLAN_ACCESS_LOCK
        int neededoutsize = (m_size%2==1)?(m_size-1)/2+1:m_size/2+1;
        if(int(out.size())!=neededoutsize)
            out.resize(neededoutsize);

        for(size_t k=0; k<out.size(); ++k){
            out[k] = 0.0;
            for(size_t n=0; n<in.size(); ++n)
                out[k] += in[n]*m_FF[k][n];
        }
        out[0].imag() = 0.0; // Ensure spectrum of real data
        out[out.size()-1].imag() = 0.0; // Ensure spectrum of real data
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    template<typename TypeInContainer, typename TypeOutContainer>
    void idft(const TypeInContainer& in, TypeOutContainer& out, int winlen=-1) {
        if(m_forward)
            throw std::string("A forward DFT FFTPlan cannot compute the backward IDFT");

        int dftlen = int((in.size()-1)*2);
        if(m_size!=dftlen)
            resize(dftlen);

        if(winlen==-1)
            winlen = m_size;
        
        FFTSCARF_PLAN_ACCESS_LOCK
        if(int(out.size())!=winlen)
            out.resize(winlen);

        FloatType oneoverdftlen = FloatType(1.0)/m_size;
        for(size_t n=0; n<winlen; ++n){
            out[n] = in[0].real()*m_FF[0][n].real();
            for(size_t k=1; k<in.size()-1; ++k){
                std::complex<FloatType> f = m_FF[k][n];
                FloatType c = in[k].real()*f.real()+in[k].imag()*f.imag();
                out[n] += c+c; // Add positive and negative components
            }
            out[n] += in[(in.size()-1)].real()*m_FF[(in.size()-1)][n].real();
            out[n] *= oneoverdftlen;
        }
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    ~FFTPlanDFTTemplate() {
    }
};

#ifdef FFTSCARF_PRECISION_SINGLE
    typedef FFTPlanDFTTemplate<float> FFTPlanSingleDFT;
    #ifndef FFTSCARF_FFTPLANSINGLE
        #define FFTSCARF_FFTPLANSINGLE
        typedef FFTPlanSingleDFT FFTPlanSingle;
    #endif
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
    typedef FFTPlanDFTTemplate<double> FFTPlanDoubleDFT;
    #ifndef FFTSCARF_FFTPLANDOUBLE
        #define FFTSCARF_FFTPLANDOUBLE
        typedef FFTPlanDoubleDFT FFTPlanDouble;
    #endif
#endif
#ifdef FFTSCARF_PRECISION_LONGDOUBLE
    typedef FFTPlanDFTTemplate<long double> FFTPlanLongDoubleDFT;
    #ifndef FFTSCARF_FFTPLANLONGDOUBLE
        #define FFTSCARF_FFTPLANLONGDOUBLE
        typedef FFTPlanLongDoubleDFT FFTPlanLongDouble;
    #endif
#endif

#if FFTSCARF_PRECISION_DEFAULT == 32
    typedef FFTPlanSingleDFT FFTPlanDFT;
    #ifndef FFTSCARF_FFTPLAN
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanSingleDFT FFTPlan;
    #endif
#elif FFTSCARF_PRECISION_DEFAULT == 64
    typedef FFTPlanDoubleDFT FFTPlanDFT;
    #ifndef FFTSCARF_FFTPLAN
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanDoubleDFT FFTPlan;
    #endif
#elif FFTSCARF_PRECISION_DEFAULT == 128
    typedef FFTPlanLongDoubleDFT FFTPlanDFT;
    #ifndef FFTSCARF_FFTPLAN
        #define FFTSCARF_FFTPLAN
        typedef FFTPlanLongDoubleDFT FFTPlan;
    #endif
#endif
}

#endif // __FFTSCARF_DFT_H__
