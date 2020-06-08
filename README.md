# Filters
A library for calculating the transfer function of continuous-time filters.
Currently the following filter types are supported:
* Lowpass Butterworth
* Highpass Butterworth
* Lowpass Chebyshev
* Lowpass Inverse Chebyshev

To compile, simply include the filters.h file in your source file and call the compiler as normal.
Using g++ this would look something like:
 `g++ output.cpp filters.cpp -o output`
