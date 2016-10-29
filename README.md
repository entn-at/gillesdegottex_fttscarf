# FFTScarf
### A wrapper for FFT implementations dedicated to audio processing

## Purpose
* Depending on usage or license issue, we might want to switch to another FFT
implementation. Though we surely don't want to change our code for this.
FFTScarf provides a wrapper, as light and simple as possible, so that we can
switch from one FFT implementation to another using a simple compilation flag.
* The benchmark/ directory provides also a comparison pipeline for FFT impl.
* Depending on the documentation of the FFT impl. its usage and the vectors
representation used might not be obvious. So, FFTScarf can also me considered
as a collection of FFT recipes.
* Multi dimensional FFTs are currently omitted since FFTScarf is mainly intended
for audio processing.

## Legal
Each FFT implementation has obviously its own license. If you intend to use
FFTScarf, you first have to have a look at the license of the implementation you
intend to use.

The code of the wrapper itself, the tests and the benchmark pipeline is in the
public domain (see UNLICENSE.md file).

## Author
Feel free to throw rotten tomatoes to:
Gilles Degottex <gilles.degottex@gmail.com>

However, if you do so, please first raise a polite and well written issue on GitHub.com
