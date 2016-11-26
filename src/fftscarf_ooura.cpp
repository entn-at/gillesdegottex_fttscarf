#include <fftscarf.h>

#include <assert.h>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;

namespace fftscarf {

std::string FFTPlanOoura::version(){
    return string("2006.12"); // This is the current built-in version
}
std::string FFTPlanOoura::libraryName(){
    stringstream result;
    result << "Ooura " << version() << " (precision " << 8*sizeof(FloatType) << "b)"; // This is the current built-in version
    return result.str();
}

FFTPlanOoura::FFTPlanOoura(bool forward)
    : FFTPlanImplementation(forward)
{
    m_ooura_a = NULL;
    m_ooura_ip = NULL;
    m_ooura_w = NULL;
}
FFTPlanOoura::FFTPlanOoura(int n, bool forward)
    : FFTPlanImplementation(n, forward)
{
    m_ooura_a = NULL;
    m_ooura_ip = NULL;
    m_ooura_w = NULL;

    resize(n);
}

void FFTPlanOoura::resize(int n)
{
    if(n==m_size) return;

    assert(n>0);
    assert(isPow2(n));

    FFTSCARF_PLAN_ACCESS_LOCK
    m_size = n;

    delete[] m_ooura_a;
    m_ooura_a = NULL;
    delete[] m_ooura_ip;
    m_ooura_ip = NULL;
    delete[] m_ooura_w;
    m_ooura_w = NULL;

    m_ooura_a = new FloatType[m_size];
    m_ooura_ip = new int[2+(1<<(int)(std::log(m_size/2+0.5)/std::log(2.0))/2)];
    m_ooura_w = new FloatType[m_size/2];
    m_ooura_ip[0] = 0; // first time only
    rdft(m_size, 1, m_ooura_a, m_ooura_ip, m_ooura_w); // init cos/sin table

    FFTSCARF_PLAN_ACCESS_UNLOCK
}

FFTPlanOoura::~FFTPlanOoura() {
    FFTSCARF_PLAN_ACCESS_LOCK
    delete[] m_ooura_a;
    m_ooura_a = NULL;
    delete[] m_ooura_ip;
    m_ooura_ip = NULL;
    delete[] m_ooura_w;
    m_ooura_w = NULL;
    FFTSCARF_PLAN_ACCESS_UNLOCK
}

}
