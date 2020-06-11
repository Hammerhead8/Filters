/* filters.h
 * Class definitions for the filters library.
 *
 * Copyright (C) 2020 Joshua Cates
 * hammerhead810@gmail.com
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>
 *
 * TODO:  Implement highpass Inverse Chebyshev filters.
	  Implement bandpass butterworth filters.
	  Implement positive passband gain for Chebyshev filters. */
#ifndef FILTERS_H
#define FILTERS_H

#include <vector>

/* Class for Butterworth lowpass filter.
 * n is an integer holding the order of the filter
 * cutoffFreq is the cutoff frequency of the filter in rad/s
 * max is the gain of the filter in dB */
class
ButterworthLP
{
	public:
		ButterworthLP (int n, double cutoffFreq, double max); /* Constructor */
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order; /* The order of the filter */
		int quads; /* Number of linear and quadratic factors */
		double w0; /* Cutoff frequency of the filter */
		double gain; /* Gain of the filter in dB */
		std::vector<double> poleAngles; /* Pole angles in radians */
		std::vector<double> Q; /* Q values for the poles */
		std::vector <std::vector<double> > coefficients; /* The factored denominator of the transfer function */
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
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order;
		int quads;
		double w0;
		double gain;
		std::vector<double> poleAngles;
		std::vector<double> Q;
		std::vector< std::vector<double> > coefficients;
		std::vector< std::vector<double> > numerator;
};

/* Class for Chebyshev lowpass filter
 * n is the order of the filter
 * cutoffFreq is the cutoff frequency of the filter in rad/s
 * passGain is the passband gain in dB
 * max is the maximum passband ripple in dB */
class
ChebyshevLP
{
	public:
		ChebyshevLP (int n, double cutoffFreq, double passGain, double max);
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order; /* The order of the filter */
		int quads; /* The number of linear and quadratic factors */
		double epsilon; /* Damping factor */
		double aMax; /* Maximum atenuation in dB */
		double w0; /* Cutoff frequency */
//		double passbandGain; /* Passband gain of the function before ripple starts */
		double passbandGain; /* Gain in the passband */
		double numerator; /* Numerator of the transfer function */
		std::vector<double> sigma; /* Real part of the poles */
		std::vector<double> omega; /* Imaginary part of the poles */
		std::vector<double> poleFreq; /* The pole frequencies */
		std::vector<double> Q; /* Q values of the poles */
		std::vector< std::vector<double> > coefficients; /* Factors of the denominator */
};

/* Class for Chebyshev highpass filter
 * n is the order of the filter
 * cutoffFreq is the cutoff frequency of the filter
 * passGain is the desired passband gain of the filter
 * max is the maximum passband ripple in dB */
class
ChebyshevHP
{
	public:
		ChebyshevHP (int n, double cutoffFreq, double passGain, double max);
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order; /* Order of the filter */
		int quads; /* Number of linear and quadratic factors */
		double epsilon; /* Damping factor */
		double aMax; /* Maximum attenuation in dB */
		double w0; /* Cutoff frequency */
		double gain; /* Gain of the filter */
		double passbandGain; /* Factor by which to divide numerator to get desired gain */
		std::vector<double> sigma; /* Real part of the poles */
		std::vector<double> omega; /* Imaginary part of the poles */
		std::vector<double> poleFreq; /* The pole frequencies */
		std::vector<double> Q; /* The Q values for the poles */
		std::vector < std::vector<double> > coefficients; /* The coefficients of the transfer function */
		std::vector < std::vector<double> > numerator; /* Numerator of the transfer function */
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
		void filterPrintf ();
		void calcCoefficients ();

	private:
		int order;
		int quads;
		double epsilon;
		double aMin; /* minimum stopband gain */
		double w0; /* Cutoff frequency */
		std::vector<double> sigma; /* Real part of the poles */
		std::vector<double> omega; /* Imaginary part of the poles */
		std::vector<double> poleFreq; /* The pole frequencies */
		std::vector<double> Q; /* Q values of the poles */
		std::vector<double> M; /* Maximum values of each stage's output */
		std::vector<double> k; /* Gain */
		std::vector<double> zeroFreq; /* The frequencies of the zeros */
		std::vector< std::vector<double> > numerator; /* Numerator in factored form */
		std::vector< std::vector<double> > coefficients; /* Denominator in factored form */
};

#endif
