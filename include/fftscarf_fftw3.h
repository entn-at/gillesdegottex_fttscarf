#ifndef __FFTSCARF_FFTW3_H__
#define __FFTSCARF_FFTW3_H__

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <fftscarf.h>
#ifdef FFTSCARF_PLAN_PROTECTACCESS
#include <boost/thread/mutex.hpp>
#endif

#include <fftw3.h>

namespace fftscarf {

template<typename _FloatType, typename fftwg_plan, typename fftwg_complex, fftwg_plan (*fftwg_plan_dft_r2c_1d)(int, _FloatType*, fftwg_complex*, unsigned), fftwg_plan (*fftwg_plan_dft_c2r_1d)(int, fftwg_complex*, _FloatType*, unsigned), void (*fftwg_execute)(const fftwg_plan), void (*fftwg_destroy_plan)(fftwg_plan), void* (*fftwg_malloc)(size_t), void (*fftwg_free)(void*)>
class FFTPlanFFTW3Template : public FFTPlanImplementation
{
public:
    typedef _FloatType FloatType;

private:
    fftwg_plan m_fftw3_plan;
    FloatType *m_fftw3_sig;
    fftwg_complex *m_fftw3_spec;
    FFTSCARF_PLAN_ACCESS_DECLARE

public:
    static std::string version(){
        return std::string("3");
    }
    static std::string libraryName(){
        std::stringstream result;
        result << std::string("FFTW ") << version() << " (precision " << 8*sizeof(FloatType) << "b)";
        return result.str();
    }

    #ifdef FFTW3RESIZINGMAXTIMESPENT
    static void setTimeLimitForPlanPreparation(float t){
        fftwg_set_timelimit(t); // From FFTW 3.1 only, no means to check version at compile time...
    }
    #endif

    FFTPlanFFTW3Template(bool forward=true) {
        m_size = -1;
        m_forward = forward;

        m_fftw3_plan = NULL;
        m_fftw3_sig = NULL;
        m_fftw3_spec = NULL;
    //    #ifdef FFTW3RESIZINGMAXTIMESPENT
    //        fftwg_set_timelimit(1.0); // From FFTW 3.1. No means exist to check version at compile time ...
    //    #endif
    }
    FFTPlanFFTW3Template(int n, bool forward=true) {
        m_size = -1;
        m_forward = forward;

        m_fftw3_plan = NULL;
        m_fftw3_sig = NULL;
        m_fftw3_spec = NULL;

        resize(n);
    }

    void resize(int n) {
        assert(n>0);

        if(n==m_size) return;

        FFTSCARF_PLAN_ACCESS_LOCK
        m_size = n;

        if(m_fftw3_sig)
            fftwg_free((void*)m_fftw3_sig);
        m_fftw3_sig = NULL;
        if(m_fftw3_spec)
            fftwg_free((void*)m_fftw3_spec);
        m_fftw3_spec = NULL;

        m_fftw3_sig = (FloatType*) fftwg_malloc(sizeof(FloatType) * m_size);
        m_fftw3_spec = (fftwg_complex*) fftwg_malloc(sizeof(fftwg_complex) * m_size);
        //  | FFTW_PRESERVE_INPUT
    //         unsigned int flags = FFTW_ESTIMATE;
        unsigned int flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT;
        // The following is likely to generate non-deterministic runs !
        // See: http://www.fftw.org/faq/section3.html#nondeterministic
        // unsigned int flags = FFTW_MEASURE;
        if(m_forward){
            m_fftw3_plan = fftwg_plan_dft_r2c_1d(m_size, m_fftw3_sig, m_fftw3_spec, flags);
    //            m_fftw3_plan = fftw_plan_dft_1d(m_size, m_fftw3_sig, m_fftw3_spec, FFTW_FORWARD, FFTW_MEASURE);
        }
        else{
            m_fftw3_plan = fftwg_plan_dft_c2r_1d(m_size, m_fftw3_spec, m_fftw3_sig, flags);
    //            m_fftw3_plan = fftw_plan_dft_1d(m_size, m_fftw3_sig, m_fftw3_spec, FFTW_BACKWARD, FFTW_MEASURE);
        }
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }

    template<typename TypeInContainer>
    void dft(const TypeInContainer& in, std::vector<std::complex<FloatType> >& out, int dftlen=-1) {
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
            m_fftw3_sig[u] = in[u];
        for(; u<dftlen; ++u)
            m_fftw3_sig[u] = 0.0;

        FFTSCARF_PLAN_ACCESS_LOCK
        fftwg_execute(m_fftw3_plan);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        for(size_t i=0; i<out.size(); ++i)
            out[i] = make_complex(m_fftw3_spec[i]);
    }

    template<typename TypeOutContainer>
    void idft(const std::vector<std::complex<FloatType> >& in, TypeOutContainer& out, int winlen=-1) {
        if(m_forward)
            throw std::string("A forward DFT FFTPlan cannot compute the backward IDFT");

        int dftlen = int((in.size()-1)*2);
        if(m_size!=dftlen)
            resize(dftlen);

        if(winlen==-1)
            winlen = m_size;
        
        if(int(out.size())!=winlen)
            out.resize(winlen);

        for(int i=0; i<m_size; i++){
            m_fftw3_spec[i][0] = in[i].real();
            m_fftw3_spec[i][1] = in[i].imag();
        }
        FFTSCARF_PLAN_ACCESS_LOCK
        fftwg_execute(m_fftw3_plan);
        FFTSCARF_PLAN_ACCESS_UNLOCK

        FloatType oneoverdftlen = FloatType(1.0)/m_size;
        for(int i=0; i<winlen; i++)
            out[i] = m_fftw3_sig[i]*oneoverdftlen;
    }

    ~FFTPlanFFTW3Template() {
        FFTSCARF_PLAN_ACCESS_LOCK
        if(m_fftw3_plan) fftwg_destroy_plan(m_fftw3_plan);
        if(m_fftw3_sig)  delete[] m_fftw3_sig;
        if(m_fftw3_spec) delete[] m_fftw3_spec;
        FFTSCARF_PLAN_ACCESS_UNLOCK
    }
};

#ifdef FFTSCARF_PRECISION_SINGLE
    typedef FFTPlanFFTW3Template<float, fftwf_plan, fftwf_complex, fftwf_plan_dft_r2c_1d, fftwf_plan_dft_c2r_1d, fftwf_execute, fftwf_destroy_plan, fftwf_malloc, fftwf_free> FFTPlanSingleFFTW3;
    #ifndef FFTSCARF_FFTPLANSINGLE
        #define FFTSCARF_FFTPLANSINGLE
        typedef FFTPlanSingleFFTW3 FFTPlanSingle;
    #endif
#endif
#ifdef FFTSCARF_PRECISION_DOUBLE
    typedef FFTPlanFFTW3Template<double, fftw_plan, fftw_complex, fftw_plan_dft_r2c_1d, fftw_plan_dft_c2r_1d, fftw_execute, fftw_destroy_plan, fftw_malloc, fftw_free> FFTPlanDoubleFFTW3;
    #ifndef FFTSCARF_FFTPLANDOUBLE
        #define FFTSCARF_FFTPLANDOUBLE
        typedef FFTPlanDoubleFFTW3 FFTPlanDouble;
    #endif
#endif

#ifdef FFTSCARF_PRECISION_DEFAULTSINGLE
    typedef FFTPlanSingleFFTW3 FFTPlanFFTW3;
#else
    typedef FFTPlanDoubleFFTW3 FFTPlanFFTW3;
#endif

#ifndef FFTSCARF_FFTPLAN
    #define FFTSCARF_FFTPLAN
    #ifdef FFTSCARF_PRECISION_DEFAULTSINGLE
        typedef FFTPlanSingleFFTW3 FFTPlan;
    #else
        typedef FFTPlanDoubleFFTW3 FFTPlan;
    #endif
#endif
}

#endif // __FFTSCARF_FFTW3_H__
