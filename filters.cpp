/* filters.cpp
 * Function definitions for the filters library.
 *
 * Copyright (C) 2020 Joshua Cates
 * hammerhead810@gmail.com
 *
 * This program is free software; you can redistribute it and/or modify1
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
 *	  Implement bandpass Butterworth filters.
	  Implement positive passband gain for Chebyshev filters. */
#include <iostream>
#include <cmath>
#include <vector>
#include "filters.h"

/***********************
 * Butterworth Lowpass *
******** ***************/

/* Create and delete a filter */

/* Class constructor */
ButterworthLP::ButterworthLP (int n, double cutoff_freq, double max)
{
	int i;

	/* Create the vectors for the pole angles, Q values, and coefficients */
	/* If the order of the filter is even */
	if (n % 2 == 0) {
		this->quads = n / 2;

		this->pole_angles.resize (quads);
		this->Q.resize (quads);

		this->coefficients.resize (quads);

		for (i = 0; i < quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	/* If the order of the filter is odd */
	else {
		this->quads = (n + 1) / 2;

		this->pole_angles.resize (quads);
		this->Q.resize (quads);

		this->coefficients.resize (quads);

		for (i = 0; i < quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	/* Initialize the constants */
	this->order = n;
	this->w0 = cutoff_freq;
	this->gain = pow (10, max / 20);
}

/* Print the numerator and denominator coefficients of the transfer function */
void
ButterworthLP::filterPrintf ()
{
	int i, j;

	std::cout << "Butterworth Lowpass:" << std::endl;

	/* Print the gain */
	std::cout << "Gain = " << this->gain << std::endl;

	std::cout <<std::endl;

	/* Print the denominators */
	std::cout << "Denominator:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->coefficients[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}

/* Calculate the coefficients of the transfer function */
void
ButterworthLP::calcCoefficients ()
{
	int n;
	int i, j;
	std::vector<double> temp (3);

	/* If the order of the filter is even */
	if (this->order % 2 == 0) {
		n = this->quads;

		/* Calcuate the pole angles and Q values */
//		this->pole_angles[0] = M_PI / (2 * this->order);

		/* Calculate the Q values of the factors */
		this->pole_angles[0] = M_PI / (2 * this->order);
		this->Q[0] = 1 / (2 * cos (this->pole_angles[0]));

		if (n > 1) {
			for (i = 1; i < n; ++i) {
				this->pole_angles[i] = this->pole_angles[i - 1] +
								M_PI / this->order;

				this->Q[i] = 1 / (2 * cos (this->pole_angles[i]));
			}
		}

		/* Now calculate the coefficients.
		 * The factors are of the form (s^2 + w0 / Q_i) * s + w0^2 */
		for (i = 0; i < n; ++i) {
			this->coefficients[i][0] = 1;
			this->coefficients[i][1] = this->w0 / this->Q[i];
			this->coefficients[i][2] = pow (this->w0, 2);
		}
	}

	/* Otherwise the order is odd.
	 * This makes the calculations slightly trickier. */
	else {
		n = this->quads;

		/* Calculate the pole angles, Q values, and the linear factor */
		this->pole_angles[0] = 0;
		this->Q[0] = .5;

		this->coefficients[0][0] = 0;
		this->coefficients[0][1] = 1;
		this->coefficients[0][2] = this->w0;

		/* Now calculate the quadratic factors */
		if (n > 1) {
			for (i = 1; i < n; ++i) {
				this->pole_angles[i] = i * M_PI / this->order;
				this->Q[i] = 1 / (2 * cos (this->pole_angles[i]));

				this->coefficients[i][0] = 1;
				this->coefficients[i][1] = this->w0 /
									this->Q[i];
				this->coefficients[i][2] = pow (this->w0, 2);
			}
		}

		/* Now calculate the coefficients.
		 * We will deal with the linear factor first.
		 * This factor is of the form (s + w0) */
//		this->coefficients[0][0] = 0;
//		this->coefficients[0][1] = 1;
//		this->coefficients[0][2] = this->w0;

		/* Now deal with the quadratic factors.
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
//		for (i = 1; i < n; ++i) {
//			this->coefficients[i][0] = 1;
//			this->coefficients[i][1] = this->w0 / this->Q[i];
//			this->coefficients[i][2] = pow (this->w0, 2);
//		}
	}

	/* Finally sort the stages in order of increasing Q */
	for (i = 0; i < this->quads - 1; ++i) {
		/* If a stage with a smaller Q is found after the current stage */
		if (this->Q[i + 1] < this->Q[i]) {
			temp[0] = this->coefficients[i][0];
			temp[1] = this->coefficients[i][1];
			temp[2] = this->coefficients[i][2];

			this->coefficients[i][0] = this->coefficients[i+1][0];
			this->coefficients[i][1] = this->coefficients[i+1][1];
			this->coefficients[i][2] = this->coefficients[i+1][2];

			this->coefficients[i+1][0] = temp[0];
			this->coefficients[i+1][1] = temp[1];
			this->coefficients[i+1][2] = temp[2];
		}
	}
}

/************************
 * Butterworth highpass *
 ************************
 *
 * The coefficients for a highpass Butterworth filter are found the same way
 * as for a lowpass Butterworth filter. The only difference is the numerator.
 * Rather than just having a constant, there will be an s term whose exponent
 * is equal to the order of the filter */

/* Class constructor */
ButterworthHP::ButterworthHP (int n, double cutoff_freq, double max)
{
	int i;

	/* Create the vectors for the pole angles, Q values, and coefficients */
	if (n % 2 == 0) {
		this->quads = n / 2;

		this->pole_angles.resize (quads);
		this->Q.resize (quads);

		this->coefficients.resize (quads);
		this->numerator.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}

	}
	else {
		this->quads = (n + 1) / 2;

		this->pole_angles.resize (quads);
		this->Q.resize (quads);

		this->coefficients.resize (quads);
		this->numerator.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	/* Initialize the constants */
	this->order = n;
	this->w0 = cutoff_freq;
	this->gain = pow (10, max / 20);
}

/* Print the numerator and denominator coefficients of the transfer function */
void
ButterworthHP::filterPrintf ()
{
	int i, j;

	std::cout << "Butterworth Highpass:" << std::endl;

	/* Print the numerators */
	std::cout << "Gain = " << this->gain << std::endl;

	std::cout << "Numerator:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->numerator[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout <<std::endl;

	/* Print the denominators */
	std::cout << "Denominator:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->coefficients[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}

/* Calculate the coefficients of the transfer function */
void
ButterworthHP::calcCoefficients ()
{
	int n;
	int i, j;
	std::vector<double> temp (3);

	/* If the order of the filter is even */
	if (this->order % 2 == 0) {
		n = this->quads;

		/* Calcuate the pole angles and Q values */
		this->pole_angles[0] = M_PI / (2 * this->order);
		this->Q[0] = 1 / (2 * cos (this->pole_angles[0]));

		if (n > 1) {
			for (i = 1; i < n; ++i) {
				this->pole_angles[i] = this->pole_angles[i - 1] +
								M_PI / this->order;
				this->Q[i] = 1 / (2 * cos (this->pole_angles[i]));
			}
		}

		/* Now calculate the coefficients.
		 * The factors are of the form (s^2 + w0 / Q_i) * s + w0^2 */
		for (i = 0; i < n; ++i) {
			this->coefficients[i][0] = 1;
			this->coefficients[i][1] = this->w0 / this->Q[i];
			this->coefficients[i][2] = pow (this->w0, 2);

			this->numerator[i][0] = 1;
			this->numerator[i][1] = 0;
			this->numerator[i][2] = 0;
		}
	}

	/* Otherwise the order is odd.
	 * This makes the calculations slightly trickier. */
	else {
		n = this->quads;

		/* Calculate the pole angles, Q values, and numerator for the linear factor
		 * This factor is of the form (s + w0) */
		this->pole_angles[0] = 0;
		this->Q[0] = .5;

		this->coefficients[0][0] = 0;
		this->coefficients[0][1] = 1;
		this->coefficients[0][2] = this->w0;

		this->numerator[0][0] = 0;
		this->numerator[0][1] = 1;
		this->numerator[0][2] = 0;

		/* If the order is greater than one then calculate the quadratic factors
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
		if (n > 1) {
			for (i = 1; i < n; ++i) {
				this->pole_angles[i] = i * M_PI / this->order;
				this->Q[i] = 1 / (2 * cos (this->pole_angles[i]));

				this->coefficients[i][0] = 1;
				this->coefficients[i][1] = this->w0 /
									this->Q[i];
				this->coefficients[i][2] = pow (this->w0, 2);

				this->numerator[i][0] = 1;
				this->numerator[i][1] = 0;
				this->numerator[i][2] = 0;
			}
		}

		/* Now calculate the coefficients.
		 * We will deal with the linear factor first.
		 * This factor is of the form (s + w0) */
//		this->coefficients[0][0] = 0;
//		this->coefficients[0][1] = 1;
//		this->coefficients[0][2] = this->w0;

//		this->numerator[0][0] = 0;
//		this->numerator[0][1] = 1;
//		this->numerator[0][2] = 0;

		/* Now deal with the quadratic factors.
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
//		for (i = 1; i < n; ++i) {
//			this->coefficients[i][0] = 1;
//			this->coefficients[i][1] = this->w0 / this->Q[i];
//			this->coefficients[i][2] = pow (this->w0, 2);

//			this->numerator[i][0] = 1;
//			this->numerator[i][1] = 0;
//			this->numerator[i][2] = 0;
//		}
	}

	/* Finally sort the stages in order of increasing Q */
	for (i = 0; i < this->quads - 1; ++i) {
		/* If a stage with a smaller Q is found after the current stage */
		if (this->Q[i + 1] < this->Q[i]) {
			temp[0] = this->coefficients[i][0];
			temp[1] = this->coefficients[i][1];
			temp[2] = this->coefficients[i][2];

			this->coefficients[i][0] = this->coefficients[i+1][0];
			this->coefficients[i][1] = this->coefficients[i+1][1];
			this->coefficients[i][2] = this->coefficients[i+1][2];

			this->coefficients[i+1][0] = temp[0];
			this->coefficients[i+1][1] = temp[1];
			this->coefficients[i+1][2] = temp[2];
		}
	}
}

/*********************
 * Chebyshev Lowpass *
 *********************/

/* Note for ChebyshevLP filters:
 * The coefficients are calculated for unity cutoff frequency.
 * To find the transfer function of a filter for cutoff frequencies other than
 * unity, simply multiply the calculated pole frequencies by the desired cutoff
 * frequency.
 *
 * For example:  For a fifth order filter with a_max = .5 dB
 * T(s) = .1789 / [(s + .36232)(s^2 + .5862s + .47677)(s^2 + .22393s + 1.0358)]
 * For a cutoff frequency of 1000 rad/s, the transfer function would be
 * T(s) = .1789 / [(s + 362.32)(s^2 + 586.2s + 476.77)(s^2 + 223.93s + 1035.8)] */

/* Class constructor */
ChebyshevLP::ChebyshevLP (int n, double cutoff_freq, double pass_gain, double max)
{
	int i;

	if (n % 2 == 0) {
		this->quads = n / 2;

		/* Create the array for the coefficients */
		this->coefficients.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	else {
		this->quads = (n + 1) / 2;

		/* Create the array for the coefficients */
		this->coefficients.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	this->sigma.resize (quads);
	this->omega.resize (quads);
	this->pole_freq.resize (quads);
	this->Q.resize (quads);

	this->order = n;
	this->w0 = cutoff_freq;
	this->a_max = max;

	this->epsilon = sqrt (pow (10, a_max / 10) - 1);
	this->passband_gain = pow (10, pass_gain / 20) / pow (w0, order);
//	passband_gain = 1 / pow (w0, order);
}

/* Print the coefficients */
void
ChebyshevLP::filterPrintf ()
{
	int i, j;

	std::cout << "Chebyshev lowpass:\n\n";

	/* Print the numerator */
	std::cout << "Numerator: " << this->numerator << "\n\n";

	/* Now print the denominator */
	std::cout << "Denominator:\n";

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->coefficients[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}


/* Calculate the coefficients for a ChebyshevLP filter */
void
ChebyshevLP::calcCoefficients ()
{
	int i, j;
	int n;
	double a, sinhA, coshA;
	double temp[3];

	/* Pole locations:  s_k = o_k + jw_k
	 * where o_k = - sinh (a) * sin ((2 * k - 1) / (2 * n) * pi)
	 * and w_k = cosh (a) * cos ((2 * k - 1) / (2 * n) * pi)
	 * k = 1, 2, ..., n
	 *
	 * The coefficients of the linear factors are of the form
	 * s - s_1
	 * The quadratic factors are of the form
	 * s^2 + (w0_i / Q_i)s + w0^2
	 *
	 * w0_i = sqrt (o_k^2 + w_k^2)
	 * Q_i = w0_k / (2 * o_k) */

	n = this->order;

	/* The first step is to calculate the value of a */
	a = asinh (1 / this->epsilon) / this->order;

	/* Next calculate the values of sinh(a) and cosh(a) */
	sinhA = this->w0 * sinh (a);
	coshA = this->w0 * cosh (a);

	/* Caculate the locations of the poles as well as the pole frequencies and Q values */
	for (i = 0; i < this->quads; ++i) {
		this->sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		this->omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		this->pole_freq[i] = sqrt (pow (this->sigma[i], 2) +
					       pow (this->omega[i], 2));

		this->Q[i] = this->pole_freq[i] / (2 * this->sigma[i]);
	}

	/* Now we can calculate the coefficients. */
	if (n % 2 == 0) { /* If the order is even then there is no linear factor */
		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i][0] = 1;
			this->coefficients[i][1] = this->pole_freq[i] / this->Q[i];
			this->coefficients[i][2] = pow (this->pole_freq[i], 2);
		}
	}
	else { /* If the order is odd then there will be a linear factor */
		/* We first need to find which pole corresponds to the linear factor */
		for (i = 0; i < this->quads; ++i) {
			if (this->omega[i] < .001) {
				this->coefficients[i][0] = 0;
				this->coefficients[i][1] = 1;
				this->coefficients[i][2] = this->pole_freq[i];
			}
			else {
				this->coefficients[i][0] = 1;
				this->coefficients[i][1] = this->pole_freq[i] / this->Q[i];
				this->coefficients[i][2] = pow (this->pole_freq[i], 2);
			}
		}

		/* If the linear factor is not the first term then we need to make it first */
		for (i = 0; i < this->quads; ++i) {
			if (this->coefficients[i][0] == 0) {
				temp[0] = this->coefficients[0][0];
				temp[1] = this->coefficients[0][1];
				temp[2] = this->coefficients[0][2];

				this->coefficients[0][0] = this->coefficients[i][0];
				this->coefficients[0][1] = this->coefficients[i][1];
				this->coefficients[0][2] = this->coefficients[i][2];

				this->coefficients[i][0] = temp[0];
				this->coefficients[i][1] = temp[1];
				this->coefficients[i][2] = temp[2];

				break;
			}
		}
	}

	/* Now sort the stages in order of increasing Q */
/*	for (i = 0; i < this->quads - 1; ++i) {
		for (j = 1; j < this->quads; ++j) { */

	for (i = 1; i < this->quads - 1; ++i) {
		if (this->Q[i+1] < this->Q[i]) { /* If a stage with a smaller Q is found after the current stage */
				temp[0] = this->coefficients[i][0];
				temp[1] = this->coefficients[i][1];
				temp[2] = this->coefficients[i][2];

				this->coefficients[i][0] = this->coefficients[i+1][0];
				this->coefficients[i][1] = this->coefficients[i+1][1];
				this->coefficients[i][2] = this->coefficients[i+1][2];

				this->coefficients[i+1][0] = temp[0];
				this->coefficients[i+1][1] = temp[1];
				this->coefficients[i+1][2] = temp[2];
		}
	}

	/* The coefficients should now be fully calculated.
	 * The only thing left to find is the numerator.
	 * This is found using the following:
	 * numerator = w_0 / (2^(n - 1) * epsilon) */
	this->numerator = 1 / (pow (2, n - 1) * this->epsilon);
	this->numerator /= this->passband_gain;

	/* Now everything should be fully calculated */
}

/**********************
 * Chebyshev Highpass *
 **********************/

/* Note for ChebyshevHP filters:
 * The coefficients are calculated for unity cutoff frequency.
 * To find the transfer function of a filter for cutoff frequencies other than
 * unity, simply multiply the calculated pole frequencies by the desired cutoff
 * frequency.
 *
 * For example:  For a fifth order filter with a_max = .5 dB
 * T(s) = .1789 / [(s + .36232)(s^2 + .5862s + .47677)(s^2 + .22393s + 1.0358)]
 * For a cutoff frequency of 1000 rad/s, the transfer function would be
 * T(s) = .1789 / [(s + 362.32)(s^2 + 586.2s + 476.77)(s^2 + 223.93s + 1035.8)] */

/* Class constructor */
ChebyshevHP::ChebyshevHP (int n, double cutoff_freq, double pass_gain, double max)
{
	int i;

	/* If the order of the filter is even */
	if (n % 2 == 0) {
		this->quads = n  / 2;

		this->order = n;
		this->w0 = cutoff_freq;
		this->a_max = max;

		this->coefficients.resize (quads);
		this->numerator.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	/* Otherwise the order of the filter is even */
	else {
		this->quads = (n + 1) / 2;

		this->order = n;
		this->w0 = cutoff_freq;
		this->a_max = max;

		this->coefficients.resize (quads);
		this->numerator.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	this->pole_freq.resize (quads);
	this->sigma.resize (quads);
	this->omega.resize (quads);
	this->Q.resize (quads);

	/* Initialize the constants */
	this->a_max = max;
	this->epsilon = sqrt (pow (10, a_max / 10) - 1);
	this->gain = 1;

	/* Calculate the passband gain of the filter.
	 * This is done by converting the passband attenuation
	 * from dB to linear and dividing by w0^order.
	 * The cutoff frequency raised to the order is used
	 * to scale the gain with respect to the cutoff frequency.
	 * This is done because otherwise the gain will be
	 * greater than the desired value for cutoff frequencies
	 * greater than 1 rad/s and less than the desired value
	 * for cutoff frequencies less than 1 rad/s. */
	this->passband_gain = pow (10, pass_gain / 20) / pow (w0, order);
}

/* Print the coefficients */
void
ChebyshevHP::filterPrintf ()
{
	int i, j;

	std::cout << "Chebyshev highpass:\n\n";

	/* Print the numerator */
	std::cout << "Gain = " << this->gain << "\n\n";

	std::cout << "Numerator:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->numerator[i][j] << " ";
		}

		std::cout << std::endl;
	}

	/* Print the denominator */
	std::cout << "Denominator:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->coefficients[i][j] << " ";
		}

		std::cout << std::endl;
	}
}

/* Calculate the transfer function */
void
ChebyshevHP::calcCoefficients ()
{
	int i, j;
	int n;
	double a, sinhA, coshA;
	double quadTerm, gainMult = 1;
	double tempNum;
	std::vector<double> temp (3);

	/* Pole locations:  s_k = o_k + jw_k
	 * where o_k = - sinh (a) * sin ((2 * k - 1) / (2 * n) * pi)
	 * and w_k = cosh (a) * cos ((2 * k - 1) / (2 * n) * pi)
	 * k = 1, 2, ..., n
	 *
	 * The coefficients of the linear factors are of the form
	 * s - s_1
	 * The quadratic factors are of the form
	 * s^2 + (w0_i / Q_i)s + w0^2
	 *
	 * w0_i = sqrt (o_k^2 + w_k^2)
	 * Q_i = w0_k / (2 * o_k) */

	n = this->order;

	/* The first step is to calculate the value of a */
	a = asinh (1 / this->epsilon) / n;

	/* Next calculate the values of sinh(a) and cosh(a)
	 * We also divide each by the cutoff frequency.
	 * This is because performing the lowpass to highpass transformation
	 * inverts the poles of the filter. */
	sinhA = sinh (a) / this->w0;
	coshA = cosh (a) / this->w0;

	/* Caculate the locations of the poles as well as the pole frequencies and Q values */
	for (i = 0; i < this->quads; ++i) {
		this->sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		this->omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		this->pole_freq[i] = sqrt (pow (this->sigma[i], 2) +
					       pow (this->omega[i], 2));

		this->Q[i] = this->pole_freq[i] / (2 * this->sigma[i]);
	}

	/* Now we can calculate the coefficients. */
	if (n % 2 == 0) { /* If the order is even then there is no linear factor */
		for (i = 0; i < this->quads; ++i) {
			this->numerator[i][0] = 1;
			this->numerator[i][1] = 0;
			this->numerator[i][2] = 0;

			this->coefficients[i][0] = 1;
			this->coefficients[i][1] = this->pole_freq[i] / this->Q[i];
			this->coefficients[i][2] = pow (this->pole_freq[i], 2);
		}
	}
	else { /* If the order is odd then there will be a linear factor */
		/* We first need to find which pole corresponds to the linear factor */
		for (i = 0; i < this->quads; ++i) {
			if (this->omega[i] < .001) {
				this->numerator[i][0] = 0;
				this->numerator[i][1] = 1;
				this->numerator[i][2] = 0;

				this->coefficients[i][0] = 0;
				this->coefficients[i][1] = 1;
				this->coefficients[i][2] = this->pole_freq[i];
			}
			else {
				this->numerator[i][0] = 1;
				this->numerator[i][1] = 0;
				this->numerator[i][2] = 0;

				this->coefficients[i][0] = 1;
				this->coefficients[i][1] = this->pole_freq[i] / this->Q[i];
				this->coefficients[i][2] = pow (this->pole_freq[i], 2);
			}
		}

		/* If the linear factor is not the first term then we need to make it first */
		for (i = 0; i < this->quads; ++i) {
			if (this->coefficients[i][0] == 0) {
				temp[0] = this->numerator[0][0];
				temp[1] = this->numerator[0][1];
				temp[2] = this->numerator[0][2];

				this->numerator[0][0] = this->numerator[i][0];
				this->numerator[0][1] = this->numerator[i][1];
				this->numerator[0][2] = this->numerator[i][2];

				this->numerator[i][0] = temp[0];
				this->numerator[i][1] = temp[1];
				this->numerator[i][2] = temp[2];

				temp[0] = this->coefficients[0][0];
				temp[1] = this->coefficients[0][1];
				temp[2] = this->coefficients[0][2];

				this->coefficients[0][0] = this->coefficients[i][0];
				this->coefficients[0][1] = this->coefficients[i][1];
				this->coefficients[0][2] = this->coefficients[i][2];

				this->coefficients[i][0] = temp[0];
				this->coefficients[i][1] = temp[1];
				this->coefficients[i][2] = temp[2];

				break;
			}
		}
	}

	/* Now sort the stages in order of increasing Q */
/*	for (i = 0; i < this->quads - 1; ++i) {
		for (j = 1; j < this->quads; ++j) { */
	
	for (i = 1; i < this->quads - 1; ++i) {
		if (this->Q[i+1] < this->Q[i]) { /* If a stage with a smaller Q is found after the current stage */
				temp[0] = this->coefficients[i][0];
				temp[1] = this->coefficients[i][1];
				temp[2] = this->coefficients[i][2];

				this->coefficients[i][0] = this->coefficients[i+1][0];
				this->coefficients[i][1] = this->coefficients[i+1][1];
				this->coefficients[i][2] = this->coefficients[i+1][2];

				this->coefficients[i+1][0] = temp[0];
				this->coefficients[i+1][1] = temp[1];
				this->coefficients[i+1][2] = temp[2];
		}
	}

	/* The last step is to perform the lowpass-to-highpass transformation.
	 * This is done by taking the lowpass transfer function T_L(S) and
	 * substituting s = 1/S to get T_H(s). That is:  T_H(s) = T_L(1/s).
	 *
	 * After making the substitution the quadratic factors are of the form
	 * 	a * s^-2 + b * s^-1 + c
	 * and the linear factors are of the form
	 * 	a * s^-1 + b
	 *
	 * We then swap the first and last terms in each factor since we multiply
	 * the quadratic factors by s^2 and the linear factors by s to remove the
	 * negative exponent. Then the quadratic factors will be of the form
	 * 	c * s^2 + b * s + a
	 * and the linear factors will be of the form
	 * 	b * s + a
	 *
	 * The last step will be to divide each factor by its leading coefficient
	 * to normalize the leading term. This will leave the quadratic factors in the form of
	 * 	s^2 + (b / c) * s + (a / c)
	 * and the linear factors in the form of
	 * 	s + (a / b)
	 *
	 * We also need to divide the gain by the product of the leading terms in each factor */

	/* If the order is even then there is no linear term to deal with */
	if (n % 2 == 0) {
		for (i = 0; i < this->quads; ++i) {
			tempNum = this->coefficients[i][0];
			this->coefficients[i][0] = this->coefficients[i][2];
			this->coefficients[i][2] = tempNum;

			quadTerm = this->coefficients[i][0];

			this->coefficients[i][0] /= quadTerm;
			this->coefficients[i][1] /= quadTerm;
			this->coefficients[i][2] /= quadTerm;

			gainMult *= quadTerm;
		}
	}

	/* Otherwise the order is odd and we need to deal with the linear factor */
	else {
		tempNum = this->coefficients[0][1];
		this->coefficients[0][1] = this->coefficients[0][2];
		this->coefficients[0][2] = tempNum;

		quadTerm = this->coefficients[0][1];

		this->coefficients[0][1] /= quadTerm;
		this->coefficients[0][2] /= quadTerm;

		gainMult *= quadTerm;

		for (i = 1; i < this->quads; ++i) {
			tempNum = this->coefficients[i][0];
			this->coefficients[i][0] = this->coefficients[i][2];
			this->coefficients[i][2] = tempNum;

			quadTerm = this->coefficients[i][0];

			this->coefficients[i][0] /= quadTerm;
			this->coefficients[i][1] /= quadTerm;
			this->coefficients[i][2] /= quadTerm;

			gainMult *= quadTerm;
		}
	}

	/* The only thing left to find is the numerator.
	 * We will use the passband gain calculated in the
	 * class constructor and the value of gainMult we
	 * calculated above when we performed the lowpass-highpass
	 * transform. The expression for the numerator is:
	 * gain = passband_gain / (2^(n - 1) * epsilon * gainMult) */
	this->gain = this->passband_gain;
	this->gain /= (pow (2, n - 1) * this->epsilon);
	this->gain /= gainMult;

	/* Now everything should be fully calculated */
}

/*****************************
 * Inverse Chebyshev Lowpass *
 *****************************
 *
 * The poles for an Inverse ChebyshevLP filter are calculated the same way as
 * for a ChebyshevLP filter, but then their reciprocal is taken.
 * Also, for a ChebyshevLP filter, the maximum passband attenuation is used for
 * calculating epsilon. For an Inverse ChebyshevLP filter, however, the minimum
 * stopband attenuation is used. Another thing is that since an Inverse Chebyshev
 * filter has ripple in the stopband, there will be zeros in the numerator of
 * the transfer function. There are always as many zeros as there are
 * quadratic factors in the denominator.
 */

/* Class constructor */
InverseChebyshevLP::InverseChebyshevLP (int n, double cutoff_freq, double pass_gain, double min)
{
	int i;

	if (n % 2 == 0) {
		this->quads = n / 2;

		this->coefficients.resize (quads);
		this->numerator.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}
	else {
		this->quads = (n + 1) / 2;

		this->coefficients.resize (quads);
		this->numerator.resize (quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	this->sigma.resize (quads);
	this->omega.resize (quads);
	this->pole_freq.resize (quads);
	this->Q.resize (quads);
	this->M.resize (quads);
	this->zero_freq.resize (quads);

	/* Initialize the constants */
	this->order = n;
	this->a_min = min;
	this->epsilon = 1 / sqrt (pow (10, a_min / 10) - 1);
	this->w0 = cutoff_freq;

	/* Convert the desired passband gain from dB to linear */
	this->pass_gain = pow (10, pass_gain / 20);
	this->passband_gain = pass_gain;
	this->K = 1;
}

/* Print the coefficients */
void
InverseChebyshevLP::filterPrintf ()
{
	int i, j;

	std::cout << "Inverse ChebyshevLP:" << std::endl;

	std::cout << "Gain = " << this->passband_gain << std::endl;

	std::cout << "Numerators:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->numerator[i][j] << " ";
		}
		std::cout << "\n";
	}

	std::cout << "\nDenominators:" << std::endl;

	for (i = 0; i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->coefficients[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;
}

/* Calculate the coefficients */
void
InverseChebyshevLP::calcCoefficients ()
{
	int i, j;
	double a, sinhA, coshA;
	int n;
	int term;
	double temp;
	std::vector<double> temp_array (3);

	a = asinh (1 / this->epsilon) / this->order;
	sinhA = sinh (a);
	coshA = cosh (a);
	n = this->order;

	/* Calculate the values of sigma_i, omega_i, pole frequencies, Q value
	 * and zero frequencies.= */
	for (i = 0; i < this->quads; ++i) {
		this->sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		this->omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		this->pole_freq[i] = sqrt (pow (this->sigma[i], 2) +
					       	      pow (this->omega[i], 2));

		this->Q[i] = this->pole_freq[i] / (2 * this->sigma[i]);

		/* sigma + j*omega is the pole of a chebyshev filter
		 * so we need to take the reciprocal to get the pole
		 * for an inverse chebyshev filter.
		 * p = 1 / s = 1 / (sigma + j*omega)
		 * p = (sigma - j * omega) / (sigma^2 + omega^2) */
		this->sigma[i] /= pow (this->pole_freq[i], 2);
		this->omega[i] /= pow (this->pole_freq[i], 2);

		/* Also invert the pole frequency and scale it by the cutoff frequency */
		this->pole_freq[i] = this->w0 / this->pole_freq[i];
	}

	/* Now sort the Q values in ascending order */
	/* If the order is even there are no linear factors */
	if (n % 2 == 0) {
		/* Sort the Q values in ascending order
		 * Also sort the pole frequencies with their corresponding Q */
		for (i = 0; i <  this->quads - 1; ++i) {
			if (this->Q[i + 1] < this->Q[i]) {
				temp = this->Q[i];
				this->Q[i] = this->Q[i + 1];
				this->Q[i + 1] = temp;

				temp = this->pole_freq[i];
				this->pole_freq[i] = this->Q[i + 1];
				this->Q[i + 1] = temp;
			}
		}

		/* Now calculate the coefficients */
		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i][0] = 1;
			this->coefficients[i][1] = this->pole_freq[i] /
								this->Q[i];
			this->coefficients[i][2] = this->pole_freq[i] * this->pole_freq[i];
		}
	}

	/* Otherwise the order is odd and there will be a linear factor */
	else {
		/* Make sure that the linear factor (Q = .5) is first */
		for (i = 0; i < this->quads; ++i) {
			if (this->Q[i] == .5) {
				if (i == 0) {
					break;
				}
				else {
					temp = this->Q[i];
					this->Q[i] = this->Q[0];
					this->Q[0] = temp;

					temp = this->pole_freq[i];
					this->pole_freq[i] = this->pole_freq[0];
					this->pole_freq[0] = temp;
					break;
				}
			}
		}

		/* Now calculate the coefficients */
		this->coefficients[0][0] = 0;
		this->coefficients[0][1] = 1;
		this->coefficients[0][2] = this->pole_freq[0];

		/* Now calculate the quadratic factors */
		for (i = 1; i < this->quads; ++i) {
			this->coefficients[i][0] = 1;
			this->coefficients[i][1] = this->pole_freq[i] /
								this->Q[i];
			this->coefficients[i][2] = this->pole_freq[i] * this->pole_freq[i];
		}
	}

	/* Now calculate the zero frequencies.
	 * If the order of the filter is even, then there will be as many
	 * zeros as there are poles.
	 * If the order is odd then there will be one fewer zeros than poles
	 * since the linear factor in the denominator does not have
	 * a corresponding zero. */
	if (n % 2 == 0) {
		for (i = 0; i < this->quads; ++i) {
			this->zero_freq[i] = this->w0 / cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		}
	}

	else {
		this->zero_freq[0] = 0;

		for (i = 1; i < this->quads; ++i) {
			this->zero_freq[i] = this->w0 / cos (((2 * i - 1) / (2 * (double)n)) * M_PI);
		}
	}

	/* Now calculate the numerators of the transfer function.
	 * Since each quadratic factor has a zero, if the order of the
	 * filter is even then there will be as many zeros as there
	 * are factors. If the order is odd, then the linear factor
	 * will not have a zero associated with it, so there will
	 * be one less zero than there are factors */
	if (n % 2 == 0) {
		for (i = 0; i < this->quads; i++) {
			this->numerator[i][0] = 1;
			this->numerator[i][1] = 0;
			this->numerator[i][2] = this->zero_freq[i] * this->zero_freq[i];
		}
	}

	else {
		this->numerator[0][0] = 0;
		this->numerator[0][1] = 0;
		this->numerator[0][2] = 1;

		for (i = 1; i < this->quads; i++) {
			this->numerator[i][0] = 1;
			this->numerator[i][1] = 0;
			this->numerator[i][2] = this->zero_freq[i] * this->zero_freq[i];

		}
	}

	/* The last step is to sort the factors by increasing Q */
	for (i = 1; i < this->quads - 1; ++i) {
		if (this->Q[i+1] < this->Q[i]) { // If a stage with a smaller Q is found after the current stage
				temp_array[0] = this->coefficients[i][0];
				temp_array[1] = this->coefficients[i][1];
				temp_array[2] = this->coefficients[i][2];

				this->coefficients[i][0] = this->coefficients[i+1][0];
				this->coefficients[i][1] = this->coefficients[i+1][1];
				this->coefficients[i][2] = this->coefficients[i+1][2];

				this->coefficients[i+1][0] = temp_array[0];
				this->coefficients[i+1][1] = temp_array[1];
				this->coefficients[i+1][2] = temp_array[2];

				temp_array[0] = this->numerator[i][0];
				temp_array[1] = this->numerator[i][1];
				temp_array[2] = this->numerator[i][2];

				this->numerator[i][0] = this->numerator[i+1][0];
				this->numerator[i][1] = this->numerator[i+1][1];
				this->numerator[i][2] = this->numerator[i+1][2];

				this->numerator[i+1][0] = temp_array[0];
				this->numerator[i+1][1] = temp_array[1];
				this->numerator[i+1][2] = temp_array[2];
		}
	}

	/* The last step is to calculate the gain.
	 * To do this we will start with K = 1
	 * We will then multiply all of the constant terms in the numerator factors
	 * and divide by the product of all of the constant terms in the denominator factors
	 * Finally we will divide this result by the scaling factor for the cutoff frequency */
	for (i = 0; i < this->quads; ++i) {
		this->K *= this->numerator[i][2];
		this->K /= this->coefficients[i][2];
	}

	this->passband_gain /= this->K;

	/* Now the transfer function should be fully calculated */
}
