# Filters
A library for calculating the transfer function of continuous-time filters.
Currently the following filter types are supported:
* Lowpass Butterworth
* Highpass Butterworth
* Lowpass Chebyshev
* Lowpass Inverse Chebyshev

This allows the tranfer function of common filter classes to be calculated given specifications without needing to consult tables.

# Features
* Supports arbitrary filter order
* Lowpass and highpass versions of common filter classes

# Limitations
* Currently has no error handling for invalid specifications (non-integer or negative order, etc.)
* No intention to support circuit design

# Planned Features
* Bandpass filters
* Arbitrary transmission zeros for lowpass and highpass
* Highpass Chebyshev and Inverse Chebyshev filters

# Usage

An example for a fifth order lowpass Butterworth filter with a cutoff frequency of 1 rad/s and a passband gain of 0 dB would be as follows:

```
#include "filters.h"

int
main ()
{
	ButterworthLP myFilter (5, 1, 0);

	myFilter.calcCoefficients ();
	myFilter.filterPrintf ();

	return 0;
}
```

To compile, simply include the filters.h file in your source file and call the compiler as normal.
Using g++ this would look something like:
 `g++ output.cpp filters.cpp -o output`
