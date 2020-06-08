/* filters.cpp
 * Function definitions for the filters library.
 *
 * TODO:  Implement highpass Chebyshev and Inverse Chebyshev filters.
 *	  Also implement bandpass Butterworth filters. */
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

/* Print the coefficients */
void
ButterworthLP::filterPrintf ()
{
	int i, j;

	std::cout << "Butterworth Lowpass:" << std::endl;

	std::cout << std::endl;

	std::cout << "Numerator: " << ButterworthLP::gain << std::endl;

	std::cout << std::endl;

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

/*			for (i = 0; i < n; ++i) {
				ButterworthLP::Q[i] = 1 / (2 * cos (ButterworthLP::poleAngles[i]));
			}*/
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

		/* Calculate the pole angles, Q values, and numerator */
		ButterworthLP::poleAngles[0] = 0;
		ButterworthLP::Q[0] = .5;

		if (n > 1) {
			for (i = 1; i < n; ++i) {
				ButterworthLP::poleAngles[i] = i * M_PI / ButterworthLP::order;
				ButterworthLP::Q[i] = 1 / (2 * cos (ButterworthLP::poleAngles[i]));
			}
		}

		/* Now calculate the coefficients.
		 * We will deal with the linear factor first.
		 * This factor is of the form (s + w0) */
		ButterworthLP::coefficients[0][0] = 0;
		ButterworthLP::coefficients[0][1] = 1;
		ButterworthLP::coefficients[0][2] = ButterworthLP::w0;

		/* Now deal with the quadratic factors.
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
		for (i = 1; i < n; ++i) {
			ButterworthLP::coefficients[i][0] = 1;
			ButterworthLP::coefficients[i][1] = ButterworthLP::w0 / ButterworthLP::Q[i];
			ButterworthLP::coefficients[i][2] = pow (ButterworthLP::w0, 2);
		}
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
 * is equal to the degree of the filter */

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

		/* Calculate the pole angles, Q values, and numerator */
		ButterworthHP::poleAngles[0] = 0;
		ButterworthHP::Q[0] = .5;

		if (n > 1) {
			for (i = 1; i < n; ++i) {
				ButterworthHP::poleAngles[i] = i * M_PI / ButterworthHP::order;
				ButterworthHP::Q[i] = 1 / (2 * cos (ButterworthHP::poleAngles[i]));
			}
		}

		/* Now calculate the coefficients.
		 * We will deal with the linear factor first.
		 * This factor is of the form (s + w0) */
		ButterworthHP::coefficients[0][0] = 0;
		ButterworthHP::coefficients[0][1] = 1;
		ButterworthHP::coefficients[0][2] = ButterworthHP::w0;

		ButterworthHP::numerator[0][0] = 0;
		ButterworthHP::numerator[0][1] = 1;
		ButterworthHP::numerator[0][2] = 0;

		/* Now deal with the quadratic factors.
		 * These are of the form s^2 + (w0 / Q_i)s + w0^2 */
		for (i = 1; i < n; ++i) {
			ButterworthHP::coefficients[i][0] = 1;
			ButterworthHP::coefficients[i][1] = ButterworthHP::w0 / ButterworthHP::Q[i];
			ButterworthHP::coefficients[i][2] = pow (ButterworthHP::w0, 2);

			ButterworthHP::numerator[i][0] = 1;
			ButterworthHP::numerator[i][1] = 0;
			ButterworthHP::numerator[i][2] = 0;
		}
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
ChebyshevLP::ChebyshevLP (int n, double cutoffFreq, double max)
{
	int i;

	if (n % 2 == 0) {
		quads = n / 2;

		/* Create the array for the coefficients */
		coefficients = new double *[quads];

		for (i = 0; i < quads; ++i) {
			coefficients[i] = new double [3];
		}

		/* Create the vector for the Q values */
		Q = new double [quads];
	}

	else {
		quads = (n / 2) + 1;

		/* Create the array for the coefficients */
		coefficients = new double *[quads];

		for (i = 0; i < quads; ++i) {
			coefficients[i] = new double [3];
		}

		/* Create the vector for the Q values */
		Q = new double [quads];
	}

	poleFreq = new double [quads];
	sigma = new double [quads];
	omega = new double [quads];

	epsilon = sqrt (pow (10, max / 10) - 1);

	order = n;
	w0 = cutoffFreq;
	aMax = max;
}

/* Class destructor */
ChebyshevLP::~ChebyshevLP ()
{
	int i;

	delete [] ChebyshevLP::sigma;
	delete [] ChebyshevLP::omega;
	delete [] ChebyshevLP::poleFreq; /* Error happens right here */
	delete [] ChebyshevLP::Q;

	for (i = 0; i < ChebyshevLP::quads; ++i) {
		delete [] ChebyshevLP::coefficients[i];
	}

	delete [] ChebyshevLP::coefficients;
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
 * stopband attenuation is used. Another thing is that since an Inverse ChebyshevLP
 * filter has ripple in the stopband, there will be zeros in the numerator of
 * the transfer function. There are always as many zeros as there are
 * quadratic factors in the transfer function.
 */

/* Class constructor */
InverseChebyshevLP::InverseChebyshevLP (int n, double cutoffFreq, double min)
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
	k.resize (quads);
	zeroFreq.resize (quads);

	/* Initialize the constants */
	epsilon = 1 / sqrt (pow (10, min / 10) - 1);
	order = n;
	aMin = min;
	w0 = cutoffFreq;
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

	/* Calculate the values of sigma_i, omega_i, pole frequencies, Q value,
	 * zero frequencies, and gains for the stages */
	for (i = 0; i < InverseChebyshevLP::quads; ++i) {
		InverseChebyshevLP::sigma[i] = sinhA * sin (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		InverseChebyshevLP::omega[i] = coshA * cos (((2 * (i + 1) - 1) / (2 * (double)n)) * M_PI);
		InverseChebyshevLP::poleFreq[i] = sqrt (pow (InverseChebyshevLP::sigma[i], 2) +
					       	      pow (InverseChebyshevLP::omega[i], 2));

		InverseChebyshevLP::Q[i] = InverseChebyshevLP::poleFreq[i] / (2 * InverseChebyshevLP::sigma[i]);

		InverseChebyshevLP::poleFreq[i] = InverseChebyshevLP::w0 / InverseChebyshevLP::poleFreq[i];
	}

	/* Now sort the Q values in ascending order.
	 * Also sort the pole frequencies with their corresponding Q */
	for (i = 0; i < InverseChebyshevLP::quads - 1; ++i) {
		/* If a stage with a smaller Q is found after the current stage */
		if (InverseChebyshevLP::Q[i + 1] < InverseChebyshevLP::Q[i]) {
				temp = InverseChebyshevLP::Q[i];
				InverseChebyshevLP::Q[i] = InverseChebyshevLP::Q[i+1];
				InverseChebyshevLP::Q[i+1] = temp;

				temp = InverseChebyshevLP::poleFreq[i];
				InverseChebyshevLP::poleFreq[i] = InverseChebyshevLP::poleFreq[i+1];
				InverseChebyshevLP::poleFreq[i+1] = temp;
		}
	}

	/* Now calculate the coefficients */
	/* If the order is even there are no linear factors */
	if (n % 2 == 0) {
		for (i = 0; i < InverseChebyshevLP::quads; ++i) {
			InverseChebyshevLP::coefficients[i][0] = 1;
			InverseChebyshevLP::coefficients[i][1] = InverseChebyshevLP::poleFreq[i] /
								InverseChebyshevLP::Q[i];
			InverseChebyshevLP::coefficients[i][2] = InverseChebyshevLP::poleFreq[i] * InverseChebyshevLP::poleFreq[i];
		}
	}

	/* Otherwise the order is odd and there will be a linear factor */
	else {
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
	 * a corresponding pole. */
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

	/* Now the transfer function should be fully calculated */
}
