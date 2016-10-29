#include "fftscarf_fftreal.h"

#include <cmath>
#include <string>
using namespace std;

namespace fftscarf {

std::string FFTPlanImplementationFFTReal::version(){
    return string("2.11"); // This is the current built-in version
}
std::string FFTPlanImplementationFFTReal::libraryName(){
    return string("FFTReal ")+version();
}

FFTPlanImplementationFFTReal::FFTPlanImplementationFFTReal(bool forward)
    : FFTPlanImplementation(forward)
{
    m_fftreal_spec = NULL;
    m_fftreal_fft = NULL;
}
FFTPlanImplementationFFTReal::FFTPlanImplementationFFTReal(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_fftreal_spec = NULL;
    m_fftreal_fft = NULL;

    resize(n);
}
void FFTPlanImplementationFFTReal::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    m_size = n;

    delete[] m_fftreal_spec;
    m_fftreal_spec = NULL;

    m_fftreal_fft = new ffft::FFTReal<FFFLOAT>(m_size);
    m_fftreal_spec = new FFFLOAT[m_size];

    if(m_forward){
        m_signal.resize(m_size);
    }
    else{
        m_signal.resize((m_size%2==1)?(m_size-1)/2+1:m_size/2+1);
    }
}

FFTPlanImplementationFFTReal::~FFTPlanImplementationFFTReal()
{
    if(m_fftreal_fft)	delete m_fftreal_fft;
    if(m_fftreal_spec)	delete[] m_fftreal_spec;
}

}
