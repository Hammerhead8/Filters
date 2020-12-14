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
ButterworthLP::ButterworthLP (int n, double cutoffFreq, double max)
{
	int i;

	/* Create the vectors for the pole angles, Q values, and coefficients */
	/* If the order of the filter is even */
	if (n % 2 == 0) {
		quads = n / 2;

		poleAngles.resize (quads);
		Q.resize (quads);

		coefficients.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
		}
	}

	/* If the order of the filter is odd */
	else {
		quads = (n + 1) / 2;

		poleAngles.resize (quads);
		Q.resize (quads);

		coefficients.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
		}
	}

	/* Initialize the constants */
	order = n;
	w0 = cutoffFreq;
	gain = pow (10, max / 20);
}

/* Print the numerator and denominator coefficients of the transfer function */
void
ButterworthLP::filterPrintf ()
{
	int i, j;

	std::cout << "Butterworth Lowpass:" << std::endl;

	/* Print the gain */
	std::cout << "Gain = " << ButterworthLP::gain << std::endl;

	std::cout <<std::endl;

	/* Print the denominators */
	std::cout << "Denominator:" << std::endl;

	for (i = 0; i < ButterworthLP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << ButterworthLP::coefficients[i][j] << " ";
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
	if (ButterworthLP::order % 2 == 0) {
		n = ButterworthLP::quads;

		/* Calcuate the pole angles and Q values */
//		ButterworthLP::poleAngles[0] = M_PI / (2 * ButterworthLP::order);

		/* Calculate the Q values of the factors */
		ButterworthLP::poleAngles[0] = M_PI / (2 * ButterworthLP::order);
		ButterworthLP::Q[0] = 1 / (2 * cos (ButterworthLP::poleAngles[0]));

		if (n > 1) {
			for (i = 1; i < n; ++i) {
				ButterworthLP::poleAngles[i] = ButterworthLP::poleAngles[i - 1] +
								M_PI / ButterworthLP::order;

				ButterworthLP::Q[i] = 1 / (2 * cos (ButterworthLP::poleAngles[i]));
			}
		}

		/* Now calculate the coefficients.
		 * The factors are of the form (s^2 + w0 / Q_i) * s + w0^2 */
		for (i = 0; i < n; ++i) {
			ButterworthLP::coefficients[i][0] = 1;
			ButterworthLP::coefficients[i][1] = ButterworthLP::w0 / ButterworthLP::Q[i];
			ButterworthLP::coefficients[i][2] = pow (ButterworthLP::w0, 2);
		}
	}

	/* Otherwise the order is odd.
	 * This makes the calculations slightly trickier. */
	else {
		n = ButterworthLP::quads;

		/* Calculate the pole angles, Q values, and the linear factor */
		ButterworthLP::poleAngles[0] = 0;
		ButterworthLP::Q[0] = .5;

		ButterworthLP::coefficients[0][0] = 0;
		ButterworthLP::coefficients[0][1] = 1;
		ButterworthLP::coefficients[0][2] = ButterworthLP::w0;

		/* Now calculate the quadratic factors */
		if (n > 1) {
			for (i = 1; i < n; ++i) {
				ButterworthLP::poleAngles[i] = i * M_PI / ButterworthLP::order;
				ButterworthLP::Q[i] = 1 / (2 * cos (ButterworthLP::poleAngles[i]));

				ButterworthLP::coefficients[i][0] = 1;
				ButterworthLP::coefficients[i][1] = ButterworthLP::w0 /
									ButterworthLP::Q[i];
				ButterworthLP::coefficients[i][2] = pow (ButterworthLP::w0, 2);
			}
		}

		/* Now calculate the coefficients.
		 * We will deal with the linear factor first.
		 * This factor is of the form (s + w0) */
//		ButterworthLP::coefficients[0][0] = 0;
//		ButterworthLP::coefficients[0][1] = 1;
//		ButterworthLP::coefficients[0][2] = ButterworthLP::w0;

		/* Now deal with the quadratic factors.
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
//		for (i = 1; i < n; ++i) {
//			ButterworthLP::coefficients[i][0] = 1;
//			ButterworthLP::coefficients[i][1] = ButterworthLP::w0 / ButterworthLP::Q[i];
//			ButterworthLP::coefficients[i][2] = pow (ButterworthLP::w0, 2);
//		}
	}

	/* Finally sort the stages in order of increasing Q */
	for (i = 0; i < ButterworthLP::quads - 1; ++i) {
		/* If a stage with a smaller Q is found after the current stage */
		if (ButterworthLP::Q[i + 1] < ButterworthLP::Q[i]) {
			temp[0] = ButterworthLP::coefficients[i][0];
			temp[1] = ButterworthLP::coefficients[i][1];
			temp[2] = ButterworthLP::coefficients[i][2];

			ButterworthLP::coefficients[i][0] = ButterworthLP::coefficients[i+1][0];
			ButterworthLP::coefficients[i][1] = ButterworthLP::coefficients[i+1][1];
			ButterworthLP::coefficients[i][2] = ButterworthLP::coefficients[i+1][2];

			ButterworthLP::coefficients[i+1][0] = temp[0];
			ButterworthLP::coefficients[i+1][1] = temp[1];
			ButterworthLP::coefficients[i+1][2] = temp[2];
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
ButterworthHP::ButterworthHP (int n, double cutoffFreq, double max)
{
	int i;

	/* Create the vectors for the pole angles, Q values, and coefficients */
	if (n % 2 == 0) {
		quads = n / 2;

		poleAngles.resize (quads);
		Q.resize (quads);

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}

	}
	else {
		quads = (n + 1) / 2;

		poleAngles.resize (quads);
		Q.resize (quads);

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}

	/* Initialize the constants */
	order = n;
	w0 = cutoffFreq;
	gain = pow (10, max / 20);
}

/* Print the numerator and denominator coefficients of the transfer function */
void
ButterworthHP::filterPrintf ()
{
	int i, j;

	std::cout << "Butterworth Highpass:" << std::endl;

	/* Print the numerators */
	std::cout << "Numerator:" << std::endl;

	for (i = 0; i < ButterworthHP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << ButterworthHP::numerator[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << "Gain = " << ButterworthHP::gain << std::endl;

	std::cout <<std::endl;

	/* Print the denominators */
	std::cout << "Denominator:" << std::endl;

	for (i = 0; i < ButterworthHP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << ButterworthHP::coefficients[i][j] << " ";
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
	if (ButterworthHP::order % 2 == 0) {
		n = ButterworthHP::quads;

		/* Calcuate the pole angles and Q values */
		ButterworthHP::poleAngles[0] = M_PI / (2 * ButterworthHP::order);
		ButterworthHP::Q[0] = 1 / (2 * cos (ButterworthHP::poleAngles[0]));

		if (n > 1) {
			for (i = 1; i < n; ++i) {
				ButterworthHP::poleAngles[i] = ButterworthHP::poleAngles[i - 1] +
								M_PI / ButterworthHP::order;
				ButterworthHP::Q[i] = 1 / (2 * cos (ButterworthHP::poleAngles[i]));
			}
		}

		/* Now calculate the coefficients.
		 * The factors are of the form (s^2 + w0 / Q_i) * s + w0^2 */
		for (i = 0; i < n; ++i) {
			ButterworthHP::coefficients[i][0] = 1;
			ButterworthHP::coefficients[i][1] = ButterworthHP::w0 / ButterworthHP::Q[i];
			ButterworthHP::coefficients[i][2] = pow (ButterworthHP::w0, 2);

			ButterworthHP::numerator[i][0] = 1;
			ButterworthHP::numerator[i][1] = 0;
			ButterworthHP::numerator[i][2] = 0;
		}
	}

	/* Otherwise the order is odd.
	 * This makes the calculations slightly trickier. */
	else {
		n = ButterworthHP::quads;

		/* Calculate the pole angles, Q values, and numerator for the linear factor
		 * This factor is of the form (s + w0) */
		ButterworthHP::poleAngles[0] = 0;
		ButterworthHP::Q[0] = .5;

		ButterworthHP::coefficients[0][0] = 0;
		ButterworthHP::coefficients[0][1] = 1;
		ButterworthHP::coefficients[0][2] = ButterworthHP::w0;

		ButterworthHP::numerator[0][0] = 0;
		ButterworthHP::numerator[0][1] = 1;
		ButterworthHP::numerator[0][2] = 0;

		/* If the order is greater than one then calculate the quadratic factors
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
		if (n > 1) {
			for (i = 1; i < n; ++i) {
				ButterworthHP::poleAngles[i] = i * M_PI / ButterworthHP::order;
				ButterworthHP::Q[i] = 1 / (2 * cos (ButterworthHP::poleAngles[i]));

				ButterworthHP::coefficients[i][0] = 1;
				ButterworthHP::coefficients[i][1] = ButterworthHP::w0 /
									ButterworthHP::Q[i];
				ButterworthHP::coefficients[i][2] = pow (ButterworthHP::w0, 2);

				ButterworthHP::numerator[i][0] = 1;
				ButterworthHP::numerator[i][1] = 0;
				ButterworthHP::numerator[i][2] = 0;
			}
		}

		/* Now calculate the coefficients.
		 * We will deal with the linear factor first.
		 * This factor is of the form (s + w0) */
//		ButterworthHP::coefficients[0][0] = 0;
//		ButterworthHP::coefficients[0][1] = 1;
//		ButterworthHP::coefficients[0][2] = ButterworthHP::w0;

//		ButterworthHP::numerator[0][0] = 0;
//		ButterworthHP::numerator[0][1] = 1;
//		ButterworthHP::numerator[0][2] = 0;

		/* Now deal with the quadratic factors.
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
//		for (i = 1; i < n; ++i) {
//			ButterworthHP::coefficients[i][0] = 1;
//			ButterworthHP::coefficients[i][1] = ButterworthHP::w0 / ButterworthHP::Q[i];
//			ButterworthHP::coefficients[i][2] = pow (ButterworthHP::w0, 2);

//			ButterworthHP::numerator[i][0] = 1;
//			ButterworthHP::numerator[i][1] = 0;
//			ButterworthHP::numerator[i][2] = 0;
//		}
	}

	/* Finally sort the stages in order of increasing Q */
	for (i = 0; i < ButterworthHP::quads - 1; ++i) {
		/* If a stage with a smaller Q is found after the current stage */
		if (ButterworthHP::Q[i + 1] < ButterworthHP::Q[i]) {
			temp[0] = ButterworthHP::coefficients[i][0];
			temp[1] = ButterworthHP::coefficients[i][1];
			temp[2] = ButterworthHP::coefficients[i][2];

			ButterworthHP::coefficients[i][0] = ButterworthHP::coefficients[i+1][0];
			ButterworthHP::coefficients[i][1] = ButterworthHP::coefficients[i+1][1];
			ButterworthHP::coefficients[i][2] = ButterworthHP::coefficients[i+1][2];

			ButterworthHP::coefficients[i+1][0] = temp[0];
			ButterworthHP::coefficients[i+1][1] = temp[1];
			ButterworthHP::coefficients[i+1][2] = temp[2];
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
 * For example:  For a fifth order filter with aMax = .5 dB
 * T(s) = .1789 / [(s + .36232)(s^2 + .5862s + .47677)(s^2 + .22393s + 1.0358)]
 * For a cutoff frequency of 1000 rad/s, the transfer function would be
 * T(s) = .1789 / [(s + 362.32)(s^2 + 586.2s + 476.77)(s^2 + 223.93s + 1035.8)] */

/* Class constructor */
ChebyshevLP::ChebyshevLP (int n, double cutoffFreq, double passGain, double max)
{
	int i;

	if (n % 2 == 0) {
		quads = n / 2;

		/* Create the array for the coefficients */
		coefficients.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
		}
	}

	else {
		quads = (n + 1) / 2;

		/* Create the array for the coefficients */
		coefficients.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
		}
	}

	sigma.resize (quads);
	omega.resize (quads);
	poleFreq.resize (quads);
	Q.resize (quads);

	order = n;
	w0 = cutoffFreq;
	aMax = max;

	epsilon = sqrt (pow (10, aMax / 10) - 1);
	passbandGain = pow (10, passGain / 20) / pow (w0, order);
//	passbandGain = 1 / pow (w0, order);
}

/* Print the coefficients */
void
ChebyshevLP::filterPrintf ()
{
	int i, j;

	std::cout << "Chebyshev lowpass:\n\n";

	/* Print the numerator */
	std::cout << "Numerator: " << ChebyshevLP::numerator << "\n\n";

	/* Now print the denominator */
	std::cout << "Denominator:\n";

	for (i = 0; i < ChebyshevLP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << ChebyshevLP::coefficients[i][j] << " ";
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

	n = ChebyshevLP::order;

	/* The first step is to calculate the value of a */
	a = asinh (1 / ChebyshevLP::epsilon) / ChebyshevLP::order;

	/* Next calculate the values of sinh(a) and cosh(a) */
	sinhA = ChebyshevLP::w0 * sinh (a);
	coshA = ChebyshevLP::w0 * cosh (a);

	/* Caculate the locations of the poles as well as the pole frequencies and Q values */
	for (i = 0; i < ChebyshevLP::quads; ++i) {
		ChebyshevLP::sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		ChebyshevLP::omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		ChebyshevLP::poleFreq[i] = sqrt (pow (ChebyshevLP::sigma[i], 2) +
					       pow (ChebyshevLP::omega[i], 2));

		ChebyshevLP::Q[i] = ChebyshevLP::poleFreq[i] / (2 * ChebyshevLP::sigma[i]);
	}

	/* Now we can calculate the coefficients. */
	if (n % 2 == 0) { /* If the order is even then there is no linear factor */
		for (i = 0; i < ChebyshevLP::quads; ++i) {
			ChebyshevLP::coefficients[i][0] = 1;
			ChebyshevLP::coefficients[i][1] = ChebyshevLP::poleFreq[i] / ChebyshevLP::Q[i];
			ChebyshevLP::coefficients[i][2] = pow (ChebyshevLP::poleFreq[i], 2);
		}
	}
	else { /* If the order is odd then there will be a linear factor */
		/* We first need to find which pole corresponds to the linear factor */
		for (i = 0; i < ChebyshevLP::quads; ++i) {
			if (ChebyshevLP::omega[i] < .001) {
				ChebyshevLP::coefficients[i][0] = 0;
				ChebyshevLP::coefficients[i][1] = 1;
				ChebyshevLP::coefficients[i][2] = ChebyshevLP::poleFreq[i];
			}
			else {
				ChebyshevLP::coefficients[i][0] = 1;
				ChebyshevLP::coefficients[i][1] = ChebyshevLP::poleFreq[i] / ChebyshevLP::Q[i];
				ChebyshevLP::coefficients[i][2] = pow (ChebyshevLP::poleFreq[i], 2);
			}
		}

		/* If the linear factor is not the first term then we need to make it first */
		for (i = 0; i < ChebyshevLP::quads; ++i) {
			if (ChebyshevLP::coefficients[i][0] == 0) {
				temp[0] = ChebyshevLP::coefficients[0][0];
				temp[1] = ChebyshevLP::coefficients[0][1];
				temp[2] = ChebyshevLP::coefficients[0][2];

				ChebyshevLP::coefficients[0][0] = ChebyshevLP::coefficients[i][0];
				ChebyshevLP::coefficients[0][1] = ChebyshevLP::coefficients[i][1];
				ChebyshevLP::coefficients[0][2] = ChebyshevLP::coefficients[i][2];

				ChebyshevLP::coefficients[i][0] = temp[0];
				ChebyshevLP::coefficients[i][1] = temp[1];
				ChebyshevLP::coefficients[i][2] = temp[2];

				break;
			}
		}
	}

	/* Now sort the stages in order of increasing Q */
/*	for (i = 0; i < ChebyshevLP::quads - 1; ++i) {
		for (j = 1; j < ChebyshevLP::quads; ++j) { */

	for (i = 1; i < ChebyshevLP::quads - 1; ++i) {
		if (ChebyshevLP::Q[i+1] < ChebyshevLP::Q[i]) { /* If a stage with a smaller Q is found after the current stage */
				temp[0] = ChebyshevLP::coefficients[i][0];
				temp[1] = ChebyshevLP::coefficients[i][1];
				temp[2] = ChebyshevLP::coefficients[i][2];

				ChebyshevLP::coefficients[i][0] = ChebyshevLP::coefficients[i+1][0];
				ChebyshevLP::coefficients[i][1] = ChebyshevLP::coefficients[i+1][1];
				ChebyshevLP::coefficients[i][2] = ChebyshevLP::coefficients[i+1][2];

				ChebyshevLP::coefficients[i+1][0] = temp[0];
				ChebyshevLP::coefficients[i+1][1] = temp[1];
				ChebyshevLP::coefficients[i+1][2] = temp[2];
		}
	}

	/* The coefficients should now be fully calculated.
	 * The only thing left to find is the numerator.
	 * This is found using the following:
	 * numerator = w_0 / (2^(n - 1) * epsilon) */
	ChebyshevLP::numerator = 1 / (pow (2, n - 1) * ChebyshevLP::epsilon);
	ChebyshevLP::numerator /= ChebyshevLP::passbandGain;

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
 * For example:  For a fifth order filter with aMax = .5 dB
 * T(s) = .1789 / [(s + .36232)(s^2 + .5862s + .47677)(s^2 + .22393s + 1.0358)]
 * For a cutoff frequency of 1000 rad/s, the transfer function would be
 * T(s) = .1789 / [(s + 362.32)(s^2 + 586.2s + 476.77)(s^2 + 223.93s + 1035.8)] */

/* Class constructor */
ChebyshevHP::ChebyshevHP (int n, double cutoffFreq, double passGain, double max)
{
	int i;

	/* If the order of the filter is even */
	if (n % 2 == 0) {
		quads = n  / 2;

		order = n;
		w0 = cutoffFreq;
		aMax = max;

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}

	/* Otherwise the order of the filter is even */
	else {
		quads = (n + 1) / 2;

		order = n;
		w0 = cutoffFreq;
		aMax = max;

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}

	poleFreq.resize (quads);
	sigma.resize (quads);
	omega.resize (quads);
	Q.resize (quads);

	/* Initialize the constants */
	aMax = max;
	epsilon = sqrt (pow (10, aMax / 10) - 1);
	gain = 1;

	/* Calculate the passband gain of the filter.
	 * This is done by converting the passband attenuation
	 * from dB to linear and dividing by w0^order.
	 * The cutoff frequency raised to the order is used
	 * to scale the gain with respect to the cutoff frequency.
	 * This is done because otherwise the gain will be
	 * greater than the desired value for cutoff frequencies
	 * greater than 1 rad/s and less than the desired value
	 * for cutoff frequencies less than 1 rad/s. */
	passbandGain = pow (10, passGain / 20) / pow (w0, order);
}

/* Print the coefficients */
void
ChebyshevHP::filterPrintf ()
{
	int i, j;

	std::cout << "Chebyshev highpass:\n\n";

	/* Print the numerator */
	std::cout << "Numerator:" << std::endl;

	for (i = 0; i < ChebyshevHP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << ChebyshevHP::numerator[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << "Gain = " << ChebyshevHP::gain << "\n\n";

	/* Print the denominator */
	std::cout << "Denominator:" << std::endl;

	for (i = 0; i < ChebyshevHP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << ChebyshevHP::coefficients[i][j] << " ";
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

	n = ChebyshevHP::order;

	/* The first step is to calculate the value of a */
	a = asinh (1 / ChebyshevHP::epsilon) / n;

	/* Next calculate the values of sinh(a) and cosh(a)
	 * We also divide each by the cutoff frequency.
	 * This is because performing the lowpass to highpass transformation
	 * inverts the poles of the filter. */
	sinhA = sinh (a) / ChebyshevHP::w0;
	coshA = cosh (a) / ChebyshevHP::w0;

	/* Caculate the locations of the poles as well as the pole frequencies and Q values */
	for (i = 0; i < ChebyshevHP::quads; ++i) {
		ChebyshevHP::sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		ChebyshevHP::omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		ChebyshevHP::poleFreq[i] = sqrt (pow (ChebyshevHP::sigma[i], 2) +
					       pow (ChebyshevHP::omega[i], 2));

		ChebyshevHP::Q[i] = ChebyshevHP::poleFreq[i] / (2 * ChebyshevHP::sigma[i]);
	}

	/* Now we can calculate the coefficients. */
	if (n % 2 == 0) { /* If the order is even then there is no linear factor */
		for (i = 0; i < ChebyshevHP::quads; ++i) {
			ChebyshevHP::numerator[i][0] = 1;
			ChebyshevHP::numerator[i][1] = 0;
			ChebyshevHP::numerator[i][2] = 0;

			ChebyshevHP::coefficients[i][0] = 1;
			ChebyshevHP::coefficients[i][1] = ChebyshevHP::poleFreq[i] / ChebyshevHP::Q[i];
			ChebyshevHP::coefficients[i][2] = pow (ChebyshevHP::poleFreq[i], 2);
		}
	}
	else { /* If the order is odd then there will be a linear factor */
		/* We first need to find which pole corresponds to the linear factor */
		for (i = 0; i < ChebyshevHP::quads; ++i) {
			if (ChebyshevHP::omega[i] < .001) {
				ChebyshevHP::numerator[i][0] = 0;
				ChebyshevHP::numerator[i][1] = 1;
				ChebyshevHP::numerator[i][2] = 0;

				ChebyshevHP::coefficients[i][0] = 0;
				ChebyshevHP::coefficients[i][1] = 1;
				ChebyshevHP::coefficients[i][2] = ChebyshevHP::poleFreq[i];
			}
			else {
				ChebyshevHP::numerator[i][0] = 1;
				ChebyshevHP::numerator[i][1] = 0;
				ChebyshevHP::numerator[i][2] = 0;

				ChebyshevHP::coefficients[i][0] = 1;
				ChebyshevHP::coefficients[i][1] = ChebyshevHP::poleFreq[i] / ChebyshevHP::Q[i];
				ChebyshevHP::coefficients[i][2] = pow (ChebyshevHP::poleFreq[i], 2);
			}
		}

		/* If the linear factor is not the first term then we need to make it first */
		for (i = 0; i < ChebyshevHP::quads; ++i) {
			if (ChebyshevHP::coefficients[i][0] == 0) {
				temp[0] = ChebyshevHP::numerator[0][0];
				temp[1] = ChebyshevHP::numerator[0][1];
				temp[2] = ChebyshevHP::numerator[0][2];

				ChebyshevHP::numerator[0][0] = ChebyshevHP::numerator[i][0];
				ChebyshevHP::numerator[0][1] = ChebyshevHP::numerator[i][1];
				ChebyshevHP::numerator[0][2] = ChebyshevHP::numerator[i][2];

				ChebyshevHP::numerator[i][0] = temp[0];
				ChebyshevHP::numerator[i][1] = temp[1];
				ChebyshevHP::numerator[i][2] = temp[2];

				temp[0] = ChebyshevHP::coefficients[0][0];
				temp[1] = ChebyshevHP::coefficients[0][1];
				temp[2] = ChebyshevHP::coefficients[0][2];

				ChebyshevHP::coefficients[0][0] = ChebyshevHP::coefficients[i][0];
				ChebyshevHP::coefficients[0][1] = ChebyshevHP::coefficients[i][1];
				ChebyshevHP::coefficients[0][2] = ChebyshevHP::coefficients[i][2];

				ChebyshevHP::coefficients[i][0] = temp[0];
				ChebyshevHP::coefficients[i][1] = temp[1];
				ChebyshevHP::coefficients[i][2] = temp[2];

				break;
			}
		}
	}

	/* Now sort the stages in order of increasing Q */
/*	for (i = 0; i < ChebyshevHP::quads - 1; ++i) {
		for (j = 1; j < ChebyshevHP::quads; ++j) { */
	
	for (i = 1; i < ChebyshevHP::quads - 1; ++i) {
		if (ChebyshevHP::Q[i+1] < ChebyshevHP::Q[i]) { /* If a stage with a smaller Q is found after the current stage */
				temp[0] = ChebyshevHP::coefficients[i][0];
				temp[1] = ChebyshevHP::coefficients[i][1];
				temp[2] = ChebyshevHP::coefficients[i][2];

				ChebyshevHP::coefficients[i][0] = ChebyshevHP::coefficients[i+1][0];
				ChebyshevHP::coefficients[i][1] = ChebyshevHP::coefficients[i+1][1];
				ChebyshevHP::coefficients[i][2] = ChebyshevHP::coefficients[i+1][2];

				ChebyshevHP::coefficients[i+1][0] = temp[0];
				ChebyshevHP::coefficients[i+1][1] = temp[1];
				ChebyshevHP::coefficients[i+1][2] = temp[2];
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
		for (i = 0; i < ChebyshevHP::quads; ++i) {
			tempNum = ChebyshevHP::coefficients[i][0];
			ChebyshevHP::coefficients[i][0] = ChebyshevHP::coefficients[i][2];
			ChebyshevHP::coefficients[i][2] = tempNum;

			quadTerm = ChebyshevHP::coefficients[i][0];

			ChebyshevHP::coefficients[i][0] /= quadTerm;
			ChebyshevHP::coefficients[i][1] /= quadTerm;
			ChebyshevHP::coefficients[i][2] /= quadTerm;

			gainMult *= quadTerm;
		}
	}

	/* Otherwise the order is odd and we need to deal with the linear factor */
	else {
		tempNum = ChebyshevHP::coefficients[0][1];
		ChebyshevHP::coefficients[0][1] = ChebyshevHP::coefficients[0][2];
		ChebyshevHP::coefficients[0][2] = tempNum;

		quadTerm = ChebyshevHP::coefficients[0][1];

		ChebyshevHP::coefficients[0][1] /= quadTerm;
		ChebyshevHP::coefficients[0][2] /= quadTerm;

		gainMult *= quadTerm;

		for (i = 1; i < ChebyshevHP::quads; ++i) {
			tempNum = ChebyshevHP::coefficients[i][0];
			ChebyshevHP::coefficients[i][0] = ChebyshevHP::coefficients[i][2];
			ChebyshevHP::coefficients[i][2] = tempNum;

			quadTerm = ChebyshevHP::coefficients[i][0];

			ChebyshevHP::coefficients[i][0] /= quadTerm;
			ChebyshevHP::coefficients[i][1] /= quadTerm;
			ChebyshevHP::coefficients[i][2] /= quadTerm;

			gainMult *= quadTerm;
		}
	}

	/* The only thing left to find is the numerator.
	 * We will use the passband gain calculated in the
	 * class constructor and the value of gainMult we
	 * calculated above when we performed the lowpass-highpass
	 * transform. The expression for the numerator is:
	 * gain = passbandGain / (2^(n - 1) * epsilon * gainMult) */
	ChebyshevHP::gain = ChebyshevHP::passbandGain;
	ChebyshevHP::gain /= (pow (2, n - 1) * ChebyshevHP::epsilon);
	ChebyshevHP::gain /= gainMult;

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
InverseChebyshevLP::InverseChebyshevLP (int n, double cutoffFreq, double passGain, double min)
{
	int i;

	if (n % 2 == 0) {
		quads = n / 2;

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}
	else {
		quads = (n + 1) / 2;

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}

	sigma.resize (quads);
	omega.resize (quads);
	poleFreq.resize (quads);
	Q.resize (quads);
	M.resize (quads);
	zeroFreq.resize (quads);

	/* Initialize the constants */
	order = n;
	aMin = min;
	epsilon = 1 / sqrt (pow (10, aMin / 10) - 1);
	w0 = cutoffFreq;

	/* Convert the desired passband gain from dB to linear */
	passGain = pow (10, passGain / 20);
	passbandGain = passGain;
	K = 1;
}

/* Print the coefficients */
void
InverseChebyshevLP::filterPrintf ()
{
	int i, j;

	std::cout << "Inverse ChebyshevLP:" << std::endl;

	std::cout << "Numerators:" << std::endl;

	for (i = 0; i < InverseChebyshevLP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << InverseChebyshevLP::numerator[i][j] << " ";
		}
		std::cout << "\n";
	}

	std::cout << "Gain = " << InverseChebyshevLP::passbandGain << std::endl;

	std::cout << "\nDenominators:" << std::endl;

	for (i = 0; i < InverseChebyshevLP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << InverseChebyshevLP::coefficients[i][j] << " ";
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
	std::vector<double> tempArray (3);

	a = asinh (1 / InverseChebyshevLP::epsilon) / InverseChebyshevLP::order;
	sinhA = sinh (a);
	coshA = cosh (a);
	n = InverseChebyshevLP::order;

	/* Calculate the values of sigma_i, omega_i, pole frequencies, Q value
	 * and zero frequencies.= */
	for (i = 0; i < InverseChebyshevLP::quads; ++i) {
		InverseChebyshevLP::sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		InverseChebyshevLP::omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		InverseChebyshevLP::poleFreq[i] = sqrt (pow (InverseChebyshevLP::sigma[i], 2) +
					       	      pow (InverseChebyshevLP::omega[i], 2));

		InverseChebyshevLP::Q[i] = InverseChebyshevLP::poleFreq[i] / (2 * InverseChebyshevLP::sigma[i]);

		/* sigma + j*omega is the pole of a chebyshev filter
		 * so we need to take the reciprocal to get the pole
		 * for an inverse chebyshev filter.
		 * p = 1 / s = 1 / (sigma + j*omega)
		 * p = (sigma - j * omega) / (sigma^2 + omega^2) */
		InverseChebyshevLP::sigma[i] /= pow (InverseChebyshevLP::poleFreq[i], 2);
		InverseChebyshevLP::omega[i] /= pow (InverseChebyshevLP::poleFreq[i], 2);

		/* Also invert the pole frequency and scale it by the cutoff frequency */
		InverseChebyshevLP::poleFreq[i] = InverseChebyshevLP::w0 / InverseChebyshevLP::poleFreq[i];
	}

	/* Now sort the Q values in ascending order */
	/* If the order is even there are no linear factors */
	if (n % 2 == 0) {
		/* Sort the Q values in ascending order
		 * Also sort the pole frequencies with their corresponding Q */
		for (i = 0; i <  InverseChebyshevLP::quads - 1; ++i) {
			if (InverseChebyshevLP::Q[i + 1] < InverseChebyshevLP::Q[i]) {
				temp = InverseChebyshevLP::Q[i];
				InverseChebyshevLP::Q[i] = InverseChebyshevLP::Q[i + 1];
				InverseChebyshevLP::Q[i + 1] = temp;

				temp = InverseChebyshevLP::poleFreq[i];
				InverseChebyshevLP::poleFreq[i] = InverseChebyshevLP::Q[i + 1];
				InverseChebyshevLP::Q[i + 1] = temp;
			}
		}

		/* Now calculate the coefficients */
		for (i = 0; i < InverseChebyshevLP::quads; ++i) {
			InverseChebyshevLP::coefficients[i][0] = 1;
			InverseChebyshevLP::coefficients[i][1] = InverseChebyshevLP::poleFreq[i] /
								InverseChebyshevLP::Q[i];
			InverseChebyshevLP::coefficients[i][2] = InverseChebyshevLP::poleFreq[i] * InverseChebyshevLP::poleFreq[i];
		}
	}

	/* Otherwise the order is odd and there will be a linear factor */
	else {
		/* Make sure that the linear factor (Q = .5) is first */
		for (i = 0; i < InverseChebyshevLP::quads; ++i) {
			if (InverseChebyshevLP::Q[i] == .5) {
				if (i == 0) {
					break;
				}
				else {
					temp = InverseChebyshevLP::Q[i];
					InverseChebyshevLP::Q[i] = InverseChebyshevLP::Q[0];
					InverseChebyshevLP::Q[0] = temp;

					temp = InverseChebyshevLP::poleFreq[i];
					InverseChebyshevLP::poleFreq[i] = InverseChebyshevLP::poleFreq[0];
					InverseChebyshevLP::poleFreq[0] = temp;
					break;
				}
			}
		}

		/* Now calculate the coefficients */
		InverseChebyshevLP::coefficients[0][0] = 0;
		InverseChebyshevLP::coefficients[0][1] = 1;
		InverseChebyshevLP::coefficients[0][2] = InverseChebyshevLP::poleFreq[0];

		/* Now calculate the quadratic factors */
		for (i = 1; i < InverseChebyshevLP::quads; ++i) {
			InverseChebyshevLP::coefficients[i][0] = 1;
			InverseChebyshevLP::coefficients[i][1] = InverseChebyshevLP::poleFreq[i] /
								InverseChebyshevLP::Q[i];
			InverseChebyshevLP::coefficients[i][2] = InverseChebyshevLP::poleFreq[i] * InverseChebyshevLP::poleFreq[i];
		}
	}

	/* Now calculate the zero frequencies.
	 * If the order of the filter is even, then there will be as many
	 * zeros as there are poles.
	 * If the order is odd then there will be one fewer zeros than poles
	 * since the linear factor in the denominator does not have
	 * a corresponding zero. */
	if (n % 2 == 0) {
		for (i = 0; i < InverseChebyshevLP::quads; ++i) {
			InverseChebyshevLP::zeroFreq[i] = InverseChebyshevLP::w0 / cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		}
	}

	else {
		InverseChebyshevLP::zeroFreq[0] = 0;

		for (i = 1; i < InverseChebyshevLP::quads; ++i) {
			InverseChebyshevLP::zeroFreq[i] = InverseChebyshevLP::w0 / cos (((2 * i - 1) / (2 * (double)n)) * M_PI);
		}
	}

	/* Now calculate the numerators of the transfer function.
	 * Since each quadratic factor has a zero, if the order of the
	 * filter is even then there will be as many zeros as there
	 * are factors. If the order is odd, then the linear factor
	 * will not have a zero associated with it, so there will
	 * be one less zero than there are factors */
	if (n % 2 == 0) {
		for (i = 0; i < InverseChebyshevLP::quads; i++) {
			InverseChebyshevLP::numerator[i][0] = 1;
			InverseChebyshevLP::numerator[i][1] = 0;
			InverseChebyshevLP::numerator[i][2] = InverseChebyshevLP::zeroFreq[i] * InverseChebyshevLP::zeroFreq[i];
		}
	}

	else {
		InverseChebyshevLP::numerator[0][0] = 0;
		InverseChebyshevLP::numerator[0][1] = 0;
		InverseChebyshevLP::numerator[0][2] = 1;

		for (i = 1; i < InverseChebyshevLP::quads; i++) {
			InverseChebyshevLP::numerator[i][0] = 1;
			InverseChebyshevLP::numerator[i][1] = 0;
			InverseChebyshevLP::numerator[i][2] = InverseChebyshevLP::zeroFreq[i] * InverseChebyshevLP::zeroFreq[i];

		}
	}

	/* The last step is to sort the factors by increasing Q */
	for (i = 1; i < InverseChebyshevLP::quads - 1; ++i) {
		if (InverseChebyshevLP::Q[i+1] < InverseChebyshevLP::Q[i]) { // If a stage with a smaller Q is found after the current stage
				tempArray[0] = InverseChebyshevLP::coefficients[i][0];
				tempArray[1] = InverseChebyshevLP::coefficients[i][1];
				tempArray[2] = InverseChebyshevLP::coefficients[i][2];

				InverseChebyshevLP::coefficients[i][0] = InverseChebyshevLP::coefficients[i+1][0];
				InverseChebyshevLP::coefficients[i][1] = InverseChebyshevLP::coefficients[i+1][1];
				InverseChebyshevLP::coefficients[i][2] = InverseChebyshevLP::coefficients[i+1][2];

				InverseChebyshevLP::coefficients[i+1][0] = tempArray[0];
				InverseChebyshevLP::coefficients[i+1][1] = tempArray[1];
				InverseChebyshevLP::coefficients[i+1][2] = tempArray[2];

				tempArray[0] = InverseChebyshevLP::numerator[i][0];
				tempArray[1] = InverseChebyshevLP::numerator[i][1];
				tempArray[2] = InverseChebyshevLP::numerator[i][2];

				InverseChebyshevLP::numerator[i][0] = InverseChebyshevLP::numerator[i+1][0];
				InverseChebyshevLP::numerator[i][1] = InverseChebyshevLP::numerator[i+1][1];
				InverseChebyshevLP::numerator[i][2] = InverseChebyshevLP::numerator[i+1][2];

				InverseChebyshevLP::numerator[i+1][0] = tempArray[0];
				InverseChebyshevLP::numerator[i+1][1] = tempArray[1];
				InverseChebyshevLP::numerator[i+1][2] = tempArray[2];
		}
	}

	/* The last step is to calculate the gain.
	 * To do this we will start with K = 1
	 * We will then multiply all of the constant terms in the numerator factors
	 * and divide by the product of all of the constant terms in the denominator factors
	 * Finally we will divide this result by the scaling factor for the cutoff frequency */
	for (i = 0; i < InverseChebyshevLP::quads; ++i) {
		InverseChebyshevLP::K *= InverseChebyshevLP::numerator[i][2];
		InverseChebyshevLP::K /= InverseChebyshevLP::coefficients[i][2];
	}

	InverseChebyshevLP::passbandGain /= InverseChebyshevLP::K;

	/* Now the transfer function should be fully calculated */
}

/* Class constructor */
InverseChebyshevHP::InverseChebyshevHP (int n, double cutoffFreq, double passGain, double min)
{
	int i;

	/* If the order of the filter is even */
	if (n % 2 == 0) {
		quads = n / 2;

		coefficients.resize (quads);
		numerator.resize (quads);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}

	/* Otherwise the order of the filter is odd */
	else {
		quads = (n + 1) / 2;

		coefficients.resize (quads);
		numerator.resize (quads + 1);

		for (i = 0; i < quads; ++i) {
			coefficients[i].resize (3);
			numerator[i].resize (3);
		}
	}

	sigma.resize (quads);
	omega.resize (quads);
	poleFreq.resize (quads);
	Q.resize (quads);
	M.resize (quads);
	zeroFreq.resize (quads);

	/* Initialize the constants */
	order = n;
	aMin = min;
	epsilon = 1 / sqrt (pow (10, aMin / 10) - 1);
	w0 = cutoffFreq;

	/* Convert the desired passband gain from dB to linear */
	passGain = pow (10, passGain / 20);
	passbandGain = passGain;
	K = 1;
}

/* Print the coefficients */
void
InverseChebyshevHP::filterPrintf ()
{
	int i, j;

	std::cout << "Inverse Chebyshev highpass:" << std::endl;

	std::cout << "Numerators:" << std::endl;

	for (i = 0; i < InverseChebyshevHP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << InverseChebyshevHP::numerator[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << "Gain = " << InverseChebyshevHP::passbandGain << std::endl << std::endl;

	std::cout << "Denominators:" << std::endl;

	for (i = 0;i < InverseChebyshevHP::quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << InverseChebyshevHP::coefficients[i][j] << " ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}

/* Calculate the transfer function */
void
InverseChebyshevHP::calcCoefficients ()
{
	int i, j;
	double a, sinhA, coshA;
	int n;
	int term;
	double tempNum;
	double quadTerm;
	double gainMult = 1;
	std::vector<double> tempArray (3);

	a = asinh (1 / InverseChebyshevHP::epsilon) / InverseChebyshevHP::order;
	sinhA = sinh (a);
	coshA = cosh (a);
	n = InverseChebyshevHP::order;

	/* Calculate the values of sigma_i, omega_i, pole frequencies, Q value
	 * and zero frequencies.= */
	for (i = 0; i < InverseChebyshevHP::quads; ++i) {
		InverseChebyshevHP::sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		InverseChebyshevHP::omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		InverseChebyshevHP::poleFreq[i] = sqrt (pow (InverseChebyshevHP::sigma[i], 2) +
					       	      pow (InverseChebyshevHP::omega[i], 2));

		InverseChebyshevHP::Q[i] = InverseChebyshevHP::poleFreq[i] / (2 * InverseChebyshevHP::sigma[i]);

		/* sigma + j*omega is the pole of a chebyshev filter
		 * so we need to take the reciprocal to get the pole
		 * for an inverse chebyshev filter.
		 * p = 1 / s = 1 / (sigma + j*omega)
		 * p = (sigma - j * omega) / (sigma^2 + omega^2) */
		InverseChebyshevHP::sigma[i] /= pow (InverseChebyshevHP::poleFreq[i], 2);
		InverseChebyshevHP::omega[i] /= pow (InverseChebyshevHP::poleFreq[i], 2);

		/* Also invert the pole frequency and scale it by the cutoff frequency */
		InverseChebyshevHP::poleFreq[i] = InverseChebyshevHP::w0 / InverseChebyshevHP::poleFreq[i];
	}

	/* Now sort the Q values in ascending order */
	/* If the order is even there are no linear factors */
	if (n % 2 == 0) {
		/* Sort the Q values in ascending order
		 * Also sort the pole frequencies with their corresponding Q */
		for (i = 0; i <  InverseChebyshevHP::quads - 1; ++i) {
			if (InverseChebyshevHP::Q[i + 1] < InverseChebyshevHP::Q[i]) {
				tempNum = InverseChebyshevHP::Q[i];
				InverseChebyshevHP::Q[i] = InverseChebyshevHP::Q[i + 1];
				InverseChebyshevHP::Q[i + 1] = tempNum;

				tempNum = InverseChebyshevHP::poleFreq[i];
				InverseChebyshevHP::poleFreq[i] = InverseChebyshevHP::Q[i + 1];
				InverseChebyshevHP::Q[i + 1] = tempNum;
			}
		}

		/* Now calculate the coefficients */
		for (i = 0; i < InverseChebyshevHP::quads; ++i) {
			InverseChebyshevHP::coefficients[i][0] = 1;
			InverseChebyshevHP::coefficients[i][1] = InverseChebyshevHP::poleFreq[i] /
								InverseChebyshevHP::Q[i];
			InverseChebyshevHP::coefficients[i][2] = InverseChebyshevHP::poleFreq[i] * InverseChebyshevHP::poleFreq[i];
		}
	}

	/* Otherwise the order is odd and there will be a linear factor */
	else {
		/* Make sure that the linear factor (Q = .5) is first */
		for (i = 0; i < InverseChebyshevHP::quads; ++i) {
			if (InverseChebyshevHP::Q[i] == .5) {
				if (i == 0) {
					break;
				}
				else {
					tempNum = InverseChebyshevHP::Q[i];
					InverseChebyshevHP::Q[i] = InverseChebyshevHP::Q[0];
					InverseChebyshevHP::Q[0] = tempNum;

					tempNum = InverseChebyshevHP::poleFreq[i];
					InverseChebyshevHP::poleFreq[i] = InverseChebyshevHP::poleFreq[0];
					InverseChebyshevHP::poleFreq[0] = tempNum;
					break;
				}
			}
		}

		/* Now calculate the coefficients */
		InverseChebyshevHP::coefficients[0][0] = 0;
		InverseChebyshevHP::coefficients[0][1] = 1;
		InverseChebyshevHP::coefficients[0][2] = InverseChebyshevHP::poleFreq[0];

		/* Now calculate the quadratic factors */
		for (i = 1; i < InverseChebyshevHP::quads; ++i) {
			InverseChebyshevHP::coefficients[i][0] = 1;
			InverseChebyshevHP::coefficients[i][1] = InverseChebyshevHP::poleFreq[i] /
								InverseChebyshevHP::Q[i];
			InverseChebyshevHP::coefficients[i][2] = InverseChebyshevHP::poleFreq[i] * InverseChebyshevHP::poleFreq[i];
		}
	}

	/* Now calculate the zero frequencies.
	 * If the order of the filter is even, then there will be as many
	 * zeros as there are poles.
	 * If the order is odd then there will be one fewer zeros than poles
	 * since the linear factor in the denominator does not have
	 * a corresponding zero. */
	if (n % 2 == 0) {
		for (i = 0; i < InverseChebyshevHP::quads; ++i) {
			InverseChebyshevHP::zeroFreq[i] = InverseChebyshevHP::w0 / cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		}
	}

	else {
		InverseChebyshevHP::zeroFreq[0] = 0;

		for (i = 1; i < InverseChebyshevHP::quads; ++i) {
			InverseChebyshevHP::zeroFreq[i] = InverseChebyshevHP::w0 / cos (((2 * i - 1) / (2 * (double)n)) * M_PI);
		}
	}

	/* Now calculate the numerators of the transfer function.
	 * Since each quadratic factor has a zero, if the order of the
	 * filter is even then there will be as many zeros as there
	 * are factors. If the order is odd, then the linear factor
	 * will not have a zero associated with it, so there will
	 * be one less zero than there are factors */
	if (n % 2 == 0) {
		for (i = 0; i < InverseChebyshevHP::quads; i++) {
			InverseChebyshevHP::numerator[i][0] = 1;
			InverseChebyshevHP::numerator[i][1] = 0;
			InverseChebyshevHP::numerator[i][2] = InverseChebyshevHP::zeroFreq[i] * InverseChebyshevHP::zeroFreq[i];
		}
	}

	else {
		InverseChebyshevHP::numerator[0][0] = 0;
		InverseChebyshevHP::numerator[0][1] = 0;
		InverseChebyshevHP::numerator[0][2] = 1;

		for (i = 1; i < InverseChebyshevHP::quads; i++) {
			InverseChebyshevHP::numerator[i][0] = 1;
			InverseChebyshevHP::numerator[i][1] = 0;
			InverseChebyshevHP::numerator[i][2] = InverseChebyshevHP::zeroFreq[i] * InverseChebyshevHP::zeroFreq[i];

		}
	}

	/* The last step is to sort the factors by increasing Q */
	for (i = 1; i < InverseChebyshevHP::quads - 1; ++i) {
		if (InverseChebyshevHP::Q[i+1] < InverseChebyshevHP::Q[i]) { // If a stage with a smaller Q is found after the current stage
				tempArray[0] = InverseChebyshevHP::coefficients[i][0];
				tempArray[1] = InverseChebyshevHP::coefficients[i][1];
				tempArray[2] = InverseChebyshevHP::coefficients[i][2];

				InverseChebyshevHP::coefficients[i][0] = InverseChebyshevHP::coefficients[i+1][0];
				InverseChebyshevHP::coefficients[i][1] = InverseChebyshevHP::coefficients[i+1][1];
				InverseChebyshevHP::coefficients[i][2] = InverseChebyshevHP::coefficients[i+1][2];

				InverseChebyshevHP::coefficients[i+1][0] = tempArray[0];
				InverseChebyshevHP::coefficients[i+1][1] = tempArray[1];
				InverseChebyshevHP::coefficients[i+1][2] = tempArray[2];

				tempArray[0] = InverseChebyshevHP::numerator[i][0];
				tempArray[1] = InverseChebyshevHP::numerator[i][1];
				tempArray[2] = InverseChebyshevHP::numerator[i][2];

				InverseChebyshevHP::numerator[i][0] = InverseChebyshevHP::numerator[i+1][0];
				InverseChebyshevHP::numerator[i][1] = InverseChebyshevHP::numerator[i+1][1];
				InverseChebyshevHP::numerator[i][2] = InverseChebyshevHP::numerator[i+1][2];

				InverseChebyshevHP::numerator[i+1][0] = tempArray[0];
				InverseChebyshevHP::numerator[i+1][1] = tempArray[1];
				InverseChebyshevHP::numerator[i+1][2] = tempArray[2];
		}
	}

	/* The next step is to perform the lowpass-to-highpass transformation.
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

	/* If the order of the filter is even then there is no linear term to deal with */
	if (n % 2 == 0) {
		for (i = 0; i < InverseChebyshevHP::quads; ++i) {
			/* Swap the first and third and coefficients of each factor */
			tempNum = InverseChebyshevHP::coefficients[i][0];
			InverseChebyshevHP::coefficients[i][0] = InverseChebyshevHP::coefficients[i][2];
			InverseChebyshevHP::coefficients[i][2] = tempNum;

			/* Put the factor in the form s^2 + b*s + c */
			quadTerm = InverseChebyshevHP::coefficients[i][0];
			InverseChebyshevHP::coefficients[i][0] = 1;
			InverseChebyshevHP::coefficients[i][1] /= InverseChebyshevHP::coefficients[i][0];
			InverseChebyshevHP::coefficients[i][2] /= InverseChebyshevHP::coefficients[i][0];

			gainMult *= quadTerm;
		}
	}

	/* Otherwise the order is odd and we need to deal with a linear factor */
	else {
		/* First deal with the linear term */
		tempNum = InverseChebyshevHP::coefficients[0][1];
		InverseChebyshevHP::coefficients[0][1] = InverseChebyshevHP::coefficients[0][2];
		InverseChebyshevHP::coefficients[0][2] = tempNum;

		quadTerm = InverseChebyshevHP::coefficients[0][1];

		/* Put the factor in the form s + b */
		InverseChebyshevHP::coefficients[0][1] = 1;
		InverseChebyshevHP::coefficients[0][2] /= InverseChebyshevHP::coefficients[0][1];

		gainMult *= quadTerm;

		/* Now deal with the quadratic terms */
		for (i = 1; i < InverseChebyshevHP::quads; ++i) {
			/* Swap the first and third and coefficients of each factor */
			tempNum = InverseChebyshevHP::coefficients[i][0];
			InverseChebyshevHP::coefficients[i][0] = InverseChebyshevHP::coefficients[i][2];
			InverseChebyshevHP::coefficients[i][2] = tempNum;

			/* Put the factor in the form s^2 + b*s + c */
			quadTerm = InverseChebyshevHP::coefficients[i][0];
			InverseChebyshevHP::coefficients[i][0] = 1;
			InverseChebyshevHP::coefficients[i][1] /= InverseChebyshevHP::coefficients[i][0];
			InverseChebyshevHP::coefficients[i][2] /= InverseChebyshevHP::coefficients[i][0];

			gainMult *= quadTerm;
		}
	}			

	/* The last step is to calculate the gain.
	 * To do this we will start with K = 1
	 * We will then multiply all of the constant terms in the numerator factors
	 * and divide by the product of all of the constant terms in the denominator factors
	 * Finally we will divide this result by the scaling factor for the cutoff frequency */
	for (i = 0; i < InverseChebyshevHP::quads; ++i) {
		InverseChebyshevHP::K *= InverseChebyshevHP::numerator[i][2];
		InverseChebyshevHP::K /= InverseChebyshevHP::coefficients[i][2];
	}

	InverseChebyshevHP::passbandGain /= InverseChebyshevHP::K;

	/* Now the transfer function should be fully calculated */
}
