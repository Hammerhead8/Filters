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
* Error handling
* Ability to cascade filters

# Compilation, Installation, and Removal
The library can be built from source by downloading the latest release under the Releases section.
Use the following commands after unpacking the archive in the desired directory:
1. `./configure` to configure the package for your system
2. `make` to compile the package
3. `make install` to install the package

To remove the library, use `make uninstall` in the source directory.

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

Here the arguments for the constructor function are the order of the filter, cutoff frequency, and passband gain. This will give the following output:
```
Numerator: 1

Denominator:
0 1 1
1 1.61803 1
1 0.618034 1
```

Another example with a sixth order highpass Chebyshev filter with a cutoff frequency of 15 rad/s and maximum passband ripple of 20 dB would be as follows:

```
#include "filters.h"

int
main ()
{
	ChebyshevHP myFilter (6, 15, 0, 20);

	myFilter.calcCoefficients ();
	myFilter.filterPrintf ();

	return 0;
}
```

Note that the arguments for the constructor function are the order of the filter, cutoff frequency, passband gain, and passband ripple. This will give the following output:

```
Chebyshev highpass:

Gain = 0.1

Numerator:
1 0 0 
1 0 0 
1 0 0 
Denominator:
1 0.139131 241.082 
1 7.20421 3344.88 
1 0.709114 449.748
```

The if the program is called `example.cc`, it can be compiled using g++ as follows:  `g++ example.cc -o example -lfilters`.
