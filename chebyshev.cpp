/* chebyshev.cpp
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
 * TODO:  Implement bandpass Butterworth filters.
	  Implement positive passband gain for Chebyshev filters. */
#include <iostream>
#include <cmath>
#include <vector>
#include "filters.h"

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
		this->coefficients.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	else {
		this->quads = (n + 1) / 2;

		/* Create the array for the coefficients */
		this->coefficients.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	this->sigma.resize (this->quads);
	this->omega.resize (this->quads);
	this->pole_freq.resize (this->quads);
	this->Q.resize (this->quads);

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

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

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

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	this->pole_freq.resize (this->quads);
	this->sigma.resize (this->quads);
	this->omega.resize (this->quads);
	this->Q.resize (this->quads);

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
	double quad_term, gain_mult = 1;
	double temp_num;
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
			temp_num = this->coefficients[i][0];
			this->coefficients[i][0] = this->coefficients[i][2];
			this->coefficients[i][2] = temp_num;

			quad_term = this->coefficients[i][0];

			this->coefficients[i][0] /= quad_term;
			this->coefficients[i][1] /= quad_term;
			this->coefficients[i][2] /= quad_term;

			gain_mult *= quad_term;
		}
	}

	/* Otherwise the order is odd and we need to deal with the linear factor */
	else {
		temp_num = this->coefficients[0][1];
		this->coefficients[0][1] = this->coefficients[0][2];
		this->coefficients[0][2] = temp_num;

		quad_term = this->coefficients[0][1];

		this->coefficients[0][1] /= quad_term;
		this->coefficients[0][2] /= quad_term;

		gain_mult *= quad_term;

		for (i = 1; i < this->quads; ++i) {
			temp_num = this->coefficients[i][0];
			this->coefficients[i][0] = this->coefficients[i][2];
			this->coefficients[i][2] = temp_num;

			quad_term = this->coefficients[i][0];

			this->coefficients[i][0] /= quad_term;
			this->coefficients[i][1] /= quad_term;
			this->coefficients[i][2] /= quad_term;

			gain_mult *= quad_term;
		}
	}

	/* The only thing left to find is the numerator.
	 * We will use the passband gain calculated in the
	 * class constructor and the value of gain_mult we
	 * calculated above when we performed the lowpass-highpass
	 * transform. The expression for the numerator is:
	 * gain = passband_gain / (2^(n - 1) * epsilon * gain_mult) */
	this->gain = this->passband_gain;
	this->gain /= (pow (2, n - 1) * this->epsilon);
	this->gain /= gain_mult;

	/* Now everything should be fully calculated */
}
