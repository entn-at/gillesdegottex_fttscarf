#include <fftscarf.h>

#include <string>
#include <sstream>

const struct AvailableFFTLibraries : public std::list<std::string> {
    AvailableFFTLibraries(){
    #ifdef FFTSCARF_FFT_OOURA
        push_back(fftscarf::FFTPlanImplementationOoura::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_FFTREAL
        push_back(fftscarf::FFTPlanImplementationFFTReal::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_PFFFT
        push_back(fftscarf::FFTPlanImplementationPFFFT::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_DJBFFT
        push_back(fftscarf::FFTPlanImplementationDJBFFT::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_FFTS
        push_back(fftscarf::FFTPlanImplementationFFTS::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_FFTW3
        push_back(fftscarf::FFTPlanImplementationFFTW3::libraryName());
    #endif
    #ifdef FFTSCARF_FFT_IPP
        push_back(fftscarf::FFTPlanImplementationIPP::libraryName());
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
        out << "    Precision: " << 8*sizeof(FFFLOAT) << "b eps=" << fftscarf::eps;
        #ifdef FFTSCARF_PRECISION_FLOAT32
        out << " (float32)" << std::endl;
        #else
        out << " (double64)" << std::endl;
        #endif
        std::list<std::string> fftlibs = FFTPlanImplementation::availableLibraries();
        if(fftlibs.size()==0)
            out << "    No FFT implementation available";
        else if(fftlibs.size()==1)
            out << "    FFT implementation: " << *(fftlibs.begin()) << std::endl;
        else{
            out << "    FFT implementations: " << std::endl;
            for(std::list<std::string>::iterator it=fftlibs.begin(); it!=fftlibs.end(); ++it)
                out << "        " << *it << std::endl;
        }
    }

    std::list<std::string> FFTPlanImplementation::availableLibraries() {
        return s_available_fft_libraries;
    }

    FFTPlanImplementation::FFTPlanImplementation(bool forward) {
        m_forward = forward;
        m_size = -1;
    }

    FFTPlanImplementation::FFTPlanImplementation(int n, bool forward) {
        m_forward = forward;
        // resize(n); // Don't call it here the vtable might not be ready
    }
}
