#include <fftscarf.h>

#include <string>
#include <sstream>

const struct AvailableFFTLibraries : public std::list<std::string> {
    AvailableFFTLibraries(){
    #ifdef FFTSCARF_FFT_IPP
        #ifdef FFTSCARF_PRECISION_SINGLE
            push_back(fftscarf::FFTPlanSingleIPP::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_DOUBLE
            push_back(fftscarf::FFTPlanDoubleIPP::libraryName());
        #endif
    #endif
    #ifdef FFTSCARF_FFT_FFTS
        push_back(fftscarf::FFTPlanFFTS::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_PFFFT
        push_back(fftscarf::FFTPlanPFFFT::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_FFTW3
        #ifdef FFTSCARF_PRECISION_SINGLE
            push_back(fftscarf::FFTPlanSingleFFTW3::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_DOUBLE
            push_back(fftscarf::FFTPlanDoubleFFTW3::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_LONGDOUBLE
            push_back(fftscarf::FFTPlanLongDoubleFFTW3::libraryName());
        #endif
    #endif
    #ifdef FFTSCARF_FFT_OOURA
        push_back(fftscarf::FFTPlanOoura::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_FFTREAL
        #ifdef FFTSCARF_PRECISION_SINGLE
            push_back(fftscarf::FFTPlanSingleFFTReal::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_DOUBLE
            push_back(fftscarf::FFTPlanDoubleFFTReal::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_LONGDOUBLE
            push_back(fftscarf::FFTPlanLongDoubleFFTReal::libraryName());
        #endif
    #endif
//    #ifdef FFTSCARF_FFT_DJBFFT
//        push_back(fftscarf::FFTPlanDJBFFT::libraryName());
//    #endif
    #ifdef FFTSCARF_FFT_DFT
        #ifdef FFTSCARF_PRECISION_SINGLE
            push_back(fftscarf::FFTPlanSingleDFT::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_DOUBLE
            push_back(fftscarf::FFTPlanDoubleDFT::libraryName());
        #endif
        #ifdef FFTSCARF_PRECISION_LONGDOUBLE
            push_back(fftscarf::FFTPlanLongDoubleDFT::libraryName());
        #endif
    #endif
    }
} s_available_fft_libraries;

namespace fftscarf {

    std::string version(){
        std::stringstream ver;
        ver << FFTSCARF_VERSION_MAJOR << "." << FFTSCARF_VERSION_MINOR << "." << FFTSCARF_VERSION_REVISION;
        return ver.str();
    }

    void write_compile_info(std::ostream& out){
        out << "    FFScarf version: " << fftscarf::version()  << std::endl;
        #ifdef FFTSCARF_LICENSE_GPLENFORCED
        out << "    License: ATTENTION: GPL is enforced on this software" << std::endl;
        #endif
        #ifdef FFTSCARF_PLAN_PROTECTACCESS
        out << "    Plan access is protected" << std::endl;
        #else
        out << "    Plan access is NOT protected" << std::endl;
        #endif
        std::list<std::string> fftlibs = availableLibraries();
        if(fftlibs.size()==0)
            out << "    No FFT implementation available";
        else if(fftlibs.size()==1)
            out << "    FFT implementation: " << *(fftlibs.begin()) << std::endl;
        else{
            out << "    FFT implementations: " << std::endl;
            for(std::list<std::string>::iterator it=fftlibs.begin(); it!=fftlibs.end(); ++it)
                out << "        " << *it << std::endl;
        }
        #if FFTSCARF_PRECISION_DEFAULT == 32
        out << "    Default precision: single (float 32b)" << std::endl;
        #elif FFTSCARF_PRECISION_DEFAULT == 64
        out << "    Default precision: double (float 64b)" << std::endl;
        #elif FFTSCARF_PRECISION_DEFAULT == 128
        out << "    Default precision: long double (float 128b)" << std::endl;
        #endif
        #ifdef FFTSCARF_FFTPLANSINGLE
        out << "    Default FFTPlanSingle: " << FFTPlanSingle::libraryName() << std::endl;
        #endif
        #ifdef FFTSCARF_FFTPLANDOUBLE
        out << "    Default FFTPlanDouble: " << FFTPlanDouble::libraryName() << std::endl;
        #endif
        #ifdef FFTSCARF_FFTPLANLONGDOUBLE
        out << "    Default FFTPlanLongDouble: " << FFTPlanLongDouble::libraryName() << std::endl;
        #endif
        #ifdef FFTSCARF_FFTPLAN
        out << "    Default FFTPlan: " << FFTPlan::libraryName() << std::endl;
        #endif
    }

    std::list<std::string> availableLibraries() {
        return s_available_fft_libraries;
    }

    FFTPlanImplementation::FFTPlanImplementation(bool forward) {
        m_forward = forward;
        m_size = -1;
    }

    FFTPlanImplementation::FFTPlanImplementation(int n, bool forward) {
        m_forward = forward;
        // resize(n); // Don't call it here the vtable might not be ready. Let the child class deal with it.
    }

    #ifdef FFTSCARF_FFT_IPP
    bool FFTPlanIPP__s_ipp_initialized = false;
    #endif
}
