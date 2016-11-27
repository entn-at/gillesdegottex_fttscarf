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

## Available implementations

Name         | Size          | single | double | long dbl | Availability | License | Notes
------------ | ------------- | ------ | ------ | -------- | ------------ | ------- | -----
[IPP][1]     | 2^a %         | Yes    | Yes    | No       | external$    | "Community License"  | .
[FFTS][2]    | 2^a %         | Yes    | No     | No       | built-in     | 3-clause BSD  | Limited efficiency on 32b architectures
[PFFFT][3]   | (2^a)*(3^b)*(5^c) [min=32] | Yes | No | No | built-in     | 3-clause BSD  | Initially unordered
[FFTW3][4]   | Any           | Yes    | Yes    | Yes      | external$    | GPL (version >=2)  | MIT License can be obtained for a charge
[Ooura][5]   | 2^a           | Yes&   | Yes&   | Yes&     | built-in     | None~   | .
[FFTReal][6] | 2^a           | Yes    | Yes    | Yes      | shipped with | WTFPL~  | Sources need to be made available
[DFT][7]     | Any           | Yes    | Yes    | Yes      | built-in     | The Unlicense~  | This is the O(N^2) DFT impl. (for comparison)

% Currently limited by FFTScarf.

& Selected during FFTScarf compilation.

^ An unordered FFT is always re-ordered by the FFTScarf wrapper.

$ When using external libraries, the end software needs to be linked with the appropriate library files (i.e. linking with fftscarf.a is not enough).

~ i.e. Public Domain.

Note that on Windows, long double (so-called 128b) might actually be a 64b or might not work at all.

[1]: https://software.intel.com/en-us/articles/how-to-use-intel-ipp-s-1d-fourier-transform-functions
[2]: https://github.com/linkotec/ffts
[3]: https://bitbucket.org/jpommier/pffft
[4]: http://www.fftw.org/
[5]: http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
[6]: http://ldesoras.free.fr/prod.html
[7]: https://en.wikipedia.org/wiki/Fourier_transform


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

When using external libraries (FFTW3 or IPP), and only in this case, your project has to access the corresponding header files (e.g. fftw3.h) and link with the lib files (e.g. -lfftw3).

When using FFTReal, its code must be made available during compilation of your project.
The other FFT libraries are fully merged into fftscarf.h and libfftscarf.a (except for the external libraries).


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
