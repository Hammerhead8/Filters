/* butterworth.cpp
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

		this->pole_angles.resize (this->quads);
		this->Q.resize (this->quads);

		this->coefficients.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	/* If the order of the filter is odd */
	else {
		this->quads = (n + 1) / 2;

		this->pole_angles.resize (this->quads);
		this->Q.resize (this->quads);

		this->coefficients.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
		}
	}

	/* Initialize the constants */
	this->order = n;
	this->w0 = cutoff_freq;
	this->gain = pow (10, max / 20) * pow (this->w0, n);
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

/*********************************
 * Butterworth transfer function *
 *********************************/

/* Class constructor */
ButterworthLP_TF::ButterworthLP_TF (int n, double numerator, double *denominator)
{
	int i;

	/* If the order is even */
	if (n % 2 == 0) {
		this->order = n;
		this->quads = n / 2;
	}

	/* Otherwise the order is odd */
	else {
		this->order = n;
		this->quads = (n + 1) / 2;
	}

	/* Set the numerator of the filter */
	this->num = numerator;

	/* Allocate the denoninator */
	this->coefficients.resize (this->quads);

	for (i = 0; i < this->quads; ++i) {
		this->coefficients[i].resize (3);

		this->coefficients[i][0] = denominator[3 * i];
		this->coefficients[i][1] = denominator[3 * i + 1];
		this->coefficients[i][2] = denominator[3 * i + 2];
	}
}

/* Calculate the cutoff frequency of the filter using the constant term of the first factor in the denominator */
double
ButterworthLP_TF::calcCutoffFreq ()
{
	/* For a Butterworth filter, the constant term of each factor is either the cutoff frequency for a linear term or
	 * the square of the cutoff frequency for a quadratic factor */
	if (this->order % 2 == 0) {
		this->w0 = sqrt (this->coefficients[0][2]);
	}

	else {
		this->w0 = this->coefficients[0][2];
	}

	return this->w0;
}

/* Calculate the passband gain of the filter.
 * G = 20 * log_10 (N / w^n)
 *
 * where G is the gain of the filter in dB
 * N is the numerator of the transfer function
 * w is the cutoff frequency of the filter
 * n is the order of the filter */
double
ButterworthLP_TF::calcGain ()
{
	this->gain = 20 * log (this->num / pow (this->w0, this->order)) / log (10);
	return this->gain;
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

		this->pole_angles.resize (this->quads);
		this->Q.resize (this->quads);

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}

	}
	else {
		this->quads = (n + 1) / 2;

		this->pole_angles.resize (this->quads);
		this->Q.resize (this->quads);

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

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
