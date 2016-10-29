#include "fftscarf_ooura.h"

#include <cmath>
#include <string>
#include <sstream>
using namespace std;

namespace fftscarf {

std::string FFTPlanImplementationOoura::version(){
    return string("2006.12"); // This is the current built-in version
}
std::string FFTPlanImplementationOoura::libraryName(){
    stringstream result;
    result << "Ooura " << version() << " (precision " << 8*sizeof(OOFLOAT) << "b)"; // This is the current built-in version
    return result.str();
}

FFTPlanImplementationOoura::FFTPlanImplementationOoura(bool forward)
    : FFTPlanImplementation(forward)
{
    m_ooura_a = NULL;
    m_ooura_ip = NULL;
    m_ooura_w = NULL;
}
FFTPlanImplementationOoura::FFTPlanImplementationOoura(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_ooura_a = NULL;
    m_ooura_ip = NULL;
    m_ooura_w = NULL;

    resize(n);
}

void FFTPlanImplementationOoura::resize(int n)
{
    assert(n>0);

    if(n==m_size) return;

    m_size = n;

    delete[] m_ooura_a;
    m_ooura_a = NULL;
    delete[] m_ooura_ip;
    m_ooura_ip = NULL;
    delete[] m_ooura_w;
    m_ooura_w = NULL;

    m_ooura_a = new OOFLOAT[m_size];
    m_ooura_ip = new int[2+(1<<(int)(std::log(m_size/2+0.5)/std::log(2.0))/2)];
    m_ooura_w = new OOFLOAT[m_size/2];
    m_ooura_ip[0] = 0; // first time only
    rdft(m_size, 1, m_ooura_a, m_ooura_ip, m_ooura_w); // init cos/sin table
}

FFTPlanImplementationOoura::~FFTPlanImplementationOoura() {
    delete[] m_ooura_a;
    m_ooura_a = NULL;
    delete[] m_ooura_ip;
    m_ooura_ip = NULL;
    delete[] m_ooura_w;
    m_ooura_w = NULL;
}

}
