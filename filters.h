/* filters.h
 * Class definitions for the filters library.
 * This library provides a resource for calculating the transfer function
 * of continuous-time filters. The supported filters are
 * Lowpass Butterworth, Highpass Butterworth, Lowpass Chebyshev,
 * Highpass Chebyshev and Lowpass Inverse Chebyshev.
 *
 * TODO:  Implement highpass Inverse Chebyshev filters */
#ifndef FILTERS_H
#define FILTERS_H

/* Class for Butterworth lowpass filter.
 * n is an integer holding the order of the filter
 * cutoffFreq is the cutoff frequency of the filter in rad/s
 * max is the gain of the filter in dB */
class
ButterworthLP
{
	public:
		ButterworthLP (int n, double cutoffFreq, double max); /* Constructor */
		~ButterworthLP (); /* Destructor */
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order; /* The order of the filter */
		int quads; /* Number of linear and quadratic factors */
		double w0; /* Cutoff frequency of the filter */
		double gain; /* Gain of the filter in dB */
		double *poleAngles; /* Vector for the pole angles in radians */
		double *Q; /* Vector for the Q values */
		double **coefficients; /* The factored transfer function */
};

/* Class for Butterworth highpass filter
 * n is an integer holding the order of the filter
 * cutoffFreq is the cutoff frequency of the filter in rad/s
 * max is the gain of the filter in dB */
class
ButterworthHP
{
	public:
		ButterworthHP (int n, double cutoffFreq, double max); /* Constructor */
		~ButterworthHP (); /* Destructor */
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order;
		int quads;
		double w0;
		double gain;
		double *poleAngles;
		double *Q;
		double **coefficients;
		double **numerator;
};

/* Class for Chebyshev lowpass filter
 * n is the order of the filter
 * cutoffFreq is the cutoff frequency of the filter in rad/s
 * max is the maximum passband gain in dB */
class
ChebyshevLP
{
	public:
		ChebyshevLP (int n, double cutoffFreq, double max);
		~ChebyshevLP ();
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order; /* The order of the filter */
		int quads; /* The number of linear and quadratic factors */
		double epsilon;
		double aMax; /* Maximum atenuation in dB */
		double w0;
		double numerator; /* Numerator of the transfer function */
		double *sigma;
		double *omega;
		double *poleFreq; /* The pole frequencies */
		double *Q; /* The Q values */
		double **coefficients; /* The coefficients of the transfer function */
};

/* Class for Chebyshev highpass filter
 * n is the order of the filter
 * cutoffFreq is the cutoff frequency of the filter
 * max is the maximum passband attenuation of the filter */
class
ChebyshevHP
{
	public:
		ChebyshevHP (int n, double cutoffFreq, double max);
		~ChebyshevHP ();
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order; /* Order of the filter */
		int quads; /* Number of linear and quadratic factors */
		double epsilon; /* Damping factor */
		double aMax; /* Maximum attenuation in dB */
		double w0; /* Cutoff frequency */
		double gain; /* Gain of the filter */
		double **numerator; /* Numerator of the transfer function */
		double *sigma; /* Real part of the poles */
		double *omega; /* Imaginary part of the poles */
		double *poleFreq; /* The pole frequencies */
		double *Q; /* The Q values for the poles */
		double **coefficients; /* The coefficients of the transfer function */
};

/* Class for Inverse Chebyshev filter
 * n is the order of the filter
 * cutoffFreq is the cutoff frequency of the filter
 * max is the maximum stopband gain */
class
InverseChebyshevLP
{
	public:
		InverseChebyshevLP (int n, double cutoffFreq, double max);
		~InverseChebyshevLP ();
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order;
		int quads;
		double epsilon;
		double aMin; /* minimum stopband gain */
		double w0; /* Cutoff frequency */
		double *sigma; /* Used to calculate pole frequencies */
		double *omega; /* Used to calculate pole frequencies */
		double *poleFreq; /* Pole frequencies */
		double *Q; /* Q values of the poles */
		double *M; /* Maximum values of each stage's output */
		double *k; /* Gain */
		double *zeroFreq; /* Zero frequencies */
		double **numerator; /* Numerator coefficients */
		double **coefficients; /* Denominator coefficients */
};

#endif
