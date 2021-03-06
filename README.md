# Filters
A library for calculating the transfer function of continuous-time filters.
Currently the following filter types are supported:
* Lowpass Butterworth
* Highpass Butterworth
* Lowpass Chebyshev
* Highpass Chebyshev
* Lowpass Inverse Chebyshev
* Highpass Inverse Chebyshev

This allows the tranfer function of common filter classes to be calculated much faster than could be done by hand and without the need to consult tables.

The results can be easily verified by plotting the frequency response of the transfer function or checking tables.

This library came about after taking a university class on analog signal processing and not wanting to have to look up transfer function coefficients in tables.
It started as a simple program to calculate the pole locations of Butterworth filters which expanded to calculating the transfer function. From there it has expanded
to include Chebyshev and Inverse Chebyshev filters.

The reference for this project is *Design of Analog Filters, 2nd Edition* by Rolf Schaumann, Haiqiao Xiao, and Mac E. Van Valkenburg. The class variable names are made to match the syntax
in this book.

# Features
* Supports arbitrary filter order
* Supports specifying passband gain
* Lowpass and highpass versions of common filter classes
* Filter stages are sorted by increasing Q for optimal cascade design
* Simple programming interface consisting of two functions
	* Calculate the transfer function
	* Print the transfer function
* Only requires standard math library and vectors
* Compatible with C++03

# Limitations
* Currently has no error handling for invalid specifications (non-integer or negative order, etc.)
* No intention to support circuit design
* Can't specify seperate passband and stopband frequencies for Chebyshev filters
* Filters cannot be cascaded

# Planned Features
* Bandpass filters
* Arbitrary transmission zeros for lowpass and highpass
* Highpass Inverse Chebyshev filters
* Error handling
* Ability to cascade filters
* Ability to build as a shared library

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
* InverseChebyshevHP
	* Highpass Inverse Chebyshev

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

