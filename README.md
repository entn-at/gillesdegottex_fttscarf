[![Build Status](https://travis-ci.org/gillesdegottex/fftscarf.svg?branch=master)](https://travis-ci.org/gillesdegottex/fftscarf)
[![Build status](https://ci.appveyor.com/api/projects/status/ehsedi0p9xl5r663/branch/master?svg=true)](https://ci.appveyor.com/project/gillesdegottex/fftscarf/branch/master)

# FFTScarf
### A wrapper for FFT implementations dedicated to audio processing


## Purpose
* Depending on usage or license issue, we might want to switch to another FFT
implementation when coding softwares using signal processing. Though, we surely
don't want to change our code for this. FFTScarf provides a wrapper, as light
and simple as possible, so that we can switch from one FFT implementation to
another using a simple compilation flag.
* The benchmark/ directory provides also a comparison pipeline for FFT
implementations.
* Depending on the documentation of the FFT implementation, the used vectors
representation and its usage might not be obvious. FFTScarf can also be 
considered as a collection of FFT recipes.
* Multi dimensional FFTs are currently omitted since audio processing is the
main target.


## Legal
Each FFT implementation has obviously its own license. If you intend to use
FFTScarf, you first have to have a look at the license of the implementation you
intend to use.

The code of the wrapper itself, the tests and the benchmark pipeline is in the
public domain (see UNLICENSE.md file).


## Compilation

Go into FFTScarf root directory, then:
```
$ mkdir build; cd build
$ cmake ..
```

Then you can use `fftscarf.h` and `libfftscarf.a` in your own project.

When using external libraries (FFTW3 or IPP), and only in this case, your project has to link with the proper lib files.

If using FFTReal, its code must be made available during compilation of your project.
The other FFT libraries are merged into fftscarf.h and libfftscarf.a (except for the external libraries).


## Benchmark

To run the benchmark, first compile it:
```
$ mkdir build-benchmark; cd build-benchmark
$ cmake ../benchmark; make
```

and run it:
```
$ make benchmark_run
$ make plot
```


## Author
Feel free to throw rotten tomatoes to:
Gilles Degottex <gilles.degottex@gmail.com>

However, if you do so, please first raise a polite and well written issue on GitHub.com
