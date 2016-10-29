#ifndef __FFTSCARF_STREAM_H__
#define __FFTSCARF_STREAM_H__

#include <iostream>
#include <vector>
#include <complex>
#include <string>

namespace fftscarf {

// Print data in the Python format

template<typename TypeValue>
inline std::ostream& operator<<(std::ostream& stream, const std::vector<std::complex<TypeValue> >& vec) {
    if(vec.empty()){
        stream << "[]";
        return stream;
    }

    stream << "[";
    size_t i=0;
    for(; i<vec.size()-1; ++i)
        stream << "(" << vec[i].real() << "," << vec[i].imag() << "), ";
    stream << "(" << vec[i].real() << "," << vec[i].imag() << ")]";
    return stream;
}

template<typename TypeValue>
inline std::ostream& operator<<(std::ostream& stream, const std::vector<TypeValue>& vec) {
    if(vec.empty()){
        stream << "[]";
        return stream;
    }

    stream << "[";
    size_t i = 0;
    for(; i<vec.size()-1; ++i)
        stream << vec[i] << ", ";
    stream << vec[i] << "]";
    return stream;
}

}

#endif // __FFTSCARF_STREAM_H__
