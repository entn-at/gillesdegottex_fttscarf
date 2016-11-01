#ifndef __FFTSCARF_FFTW3_H__
#define __FFTSCARF_FFTW3_H__

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cmath>

#include <iostream>

#include <fftscarf.h>
#ifndef FFTScarf
#define FFTScarf FFTPlanImplementationFFTW3
#endif

// #include <boost/signals2/mutex.hpp>

#include <fftw3.h>
#ifdef FFTSCARF_PRECISION_FLOAT32
    #define fftwg_plan fftwf_plan
    #define fftwg_complex fftwf_complex
    #define fftwg_set_timelimit fftwf_set_timelimit
    #define fftwg_plan_dft_r2c_1d fftwf_plan_dft_r2c_1d
    #define fftwg_plan_dft_c2r_1d fftwf_plan_dft_c2r_1d
    #define fftwg_execute fftwf_execute
    #define fftwg_destroy_plan fftwf_destroy_plan
    #define fftwg_malloc fftwf_malloc
    #define fftwg_free fftwf_free
#else
    #define fftwg_plan fftw_plan
    #define fftwg_complex fftw_complex
    #define fftwg_set_timelimit fftw_set_timelimit
    #define fftwg_plan_dft_r2c_1d fftw_plan_dft_r2c_1d
    #define fftwg_plan_dft_c2r_1d fftw_plan_dft_c2r_1d
    #define fftwg_execute fftw_execute
    #define fftwg_destroy_plan fftw_destroy_plan
    #define fftwg_malloc fftw_malloc
    #define fftwg_free fftw_free
#endif

namespace fftscarf {

class FFTPlanImplementationFFTW3 : public FFTPlanImplementation
{
    fftwg_plan m_fftw3_plan;
    FFFLOAT *m_fftw3_sig;
    fftwg_complex *m_fftw3_spec;
//    static boost::signals2::mutex m_fftw3_planner_access; // To protect the access to the FFT and external variables

public:
    static std::string version();
    static std::string libraryName();
    static void setTimeLimitForPlanPreparation(float t); // t[s]

    FFTPlanImplementationFFTW3(bool forward=true);
    FFTPlanImplementationFFTW3(int n, bool forward=true);
    ~FFTPlanImplementationFFTW3();

    void resize(int n);

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
            m_fftw3_sig[u] = in[u];
        for(; u<dftlen; ++u)
            m_fftw3_sig[u] = 0.0;

//        m_fftw3_planner_access.lock();
        fftwg_execute(m_fftw3_plan);
//        m_fftw3_planner_access.unlock();

        for(size_t i=0; i<out.size(); ++i)
            out[i] = make_complex(m_fftw3_spec[i]);
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

        for(int i=0; i<m_size; i++){
            m_fftw3_spec[i][0] = in[i].real();
            m_fftw3_spec[i][1] = in[i].imag();
        }
//        m_fftw3_planner_access.lock();
        fftwg_execute(m_fftw3_plan);
//        m_fftw3_planner_access.unlock();

        FFFLOAT oneoverdftlen = 1.0/m_size;
        for(size_t i=0; i<winlen; i++)
            out[i] = m_fftw3_sig[i]*oneoverdftlen;
    }
};
    
}

#endif // __FFTSCARF_FFTW3_H__
