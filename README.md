# Filters
A library for calculating the transfer function of continuous-time filters.
Currently the following filter types are supported:
* Lowpass Butterworth
* Highpass Butterworth
* Lowpass Chebyshev
* Highpass Chebyshev
* Lowpass Inverse Chebyshev

This allows the tranfer function of common filter classes to be calculated much faster than could be done by hand and without the need to consult tables.

The results can be easily verified by plotting the frequency response of the transfer function.

# Features
* Supports arbitrary filter order
* Lowpass and highpass versions of common filter classes
* Simple programming interface consisting of two functions
* Only requires standard math library and vectors

# Limitations
* Currently has no error handling for invalid specifications (non-integer or negative order, etc.)
* No intention to support circuit design
* No way to specify passband gain, so calculated gain might not match desired specification

# Planned Features
* Bandpass filters
* Arbitrary transmission zeros for lowpass and highpass
* Highpass Inverse Chebyshev filters
* Error handling
* Specify desired passband gain

# Usage
The classes for the available filter types are as follows:
* ButterworthLP
	* Lowpass Butterworth
* ButterworthHP
	* Highpass Butterworth
* ChebyshevLP
	* Lowpass Chebyshev
* ChebyshevHP
	* Highpass Chebyshev
* InverseChebyshevLP
	* Lowpass Inverse Chebyshev

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

This will give the following output:
```
Numerator: 1

Denominator:
0 1 1
1 1.61803 1
1 0.618034 1
```

Another example with a fourth order highpass Chebyshev filter with a cutoff frequency of 15 rad/s and maximum passband ripple of 20 dB would be as follows:

```
#include "filters.h"

int
main ()
{
	ChebyshevHP myFilter (6, 15, 20);

	myFilter.calcCoefficients ();
	myFilter.filterPrintf ();

	return 0;
}
```

This will give the following output:

```
Numerator:
1 0 0
1 0 0
Gain = 5062.5

Denominator:
1 0.337171 263.41 
1 4.72753 1529.82
```

To compile, simply include the filters.h file in your source file and call the compiler as normal.
Using g++ this would look something like:
 `g++ output.cpp filters.cpp -o output`

