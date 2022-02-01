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
 * TODO:  Implement bandpass Butterworth filters.
	  Implement positive passband gain for Chebyshev filters. */
#include <iostream>
#include <cmath>
#include <vector>
#include "filters.h"

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

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}
	else {
		this->quads = (n + 1) / 2;

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	this->sigma.resize (this->quads);
	this->omega.resize (this->quads);
	this->pole_freq.resize (this->quads);
	this->Q.resize (this->quads);
	this->M.resize (this->quads);
	this->zero_freq.resize (this->quads);

	/* Initialize the constants */
	this->order = n;
	this->a_min = min;
	this->epsilon = 1 / sqrt (pow (10, a_min / 10) - 1);
	this->w0 = cutoff_freq;

	/* Convert the desired passband gain from dB to linear */
	this->passband_gain = pow (10, pass_gain / 20);
	this->K = 1;
}

/* Print the coefficients */
void
InverseChebyshevLP::filterPrintf ()
{
	int i, j;

	std::cout << "Inverse Chebyshev lowpass:" << std::endl;

	/* Print the numerator */
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

/* Class constructor */
InverseChebyshevHP::InverseChebyshevHP (int n, double cutoff_freq, double pass_gain, double min)
{
	int i;

	/* If the order of the filter is even */
	if (n % 2 == 0) {
		this->quads = n / 2;

		this->coefficients.resize (this->quads);
		this->numerator.resize (this->quads);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}
	}

	/* Otherwise the order of the filter is odd */
	else {
		this->quads = (n + 1) / 2;

		/* Numerator has an extra row because the lowpass to highpass transformation introduces a zero to the transfer function */
		this->coefficients.resize (quads);
		this->numerator.resize (quads + 1);

		for (i = 0; i < this->quads; ++i) {
			this->coefficients[i].resize (3);
			this->numerator[i].resize (3);
		}

		this->numerator[this->quads].resize (3);
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
	this->passband_gain = pow (10, pass_gain / 20);
//	passband_gain = pass_gain;
	this->K = 1;
}

/* Print the coefficients */
void
InverseChebyshevHP::filterPrintf ()
{
	int i, j;

	std::cout << "Inverse Chebyshev highpass:" << std::endl;

	std::cout << "Numerators:" << std::endl;

	/* If the order is even then there is no introduced zero to worry about printing */
	if (this->order % 2 == 0) {
		for (i = 0; i < this->quads; ++i) {
			for (j = 0; j < 3; ++j) {
				std::cout << this->numerator[i][j] << " ";
			}
	
			std::cout << std::endl;
		}
	}

	/* Otherwise the order is odd and we need to print the introduced zero */
	else {
		for (i = 0; i < this->quads + 1; ++i) {
			for (j = 0; j < 3; ++j) {
				std::cout << this->numerator[i][j] << " ";
			}
	
			std::cout << std::endl;
		}
	}

	std::cout << "Gain = " << this->passband_gain << std::endl << std::endl;

	std::cout << "Denominators:" << std::endl;

	for (i = 0;i < this->quads; ++i) {
		for (j = 0; j < 3; ++j) {
			std::cout << this->coefficients[i][j] << " ";
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
	double temp_num;
	double quad_term;
	double gain_mult = 1;
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
				temp_num = this->Q[i];
				this->Q[i] = this->Q[i + 1];
				this->Q[i + 1] = temp_num;

				temp_num = this->pole_freq[i];
				this->pole_freq[i] = this->Q[i + 1];
				this->Q[i + 1] = temp_num;
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
					temp_num = this->Q[i];
					this->Q[i] = this->Q[0];
					this->Q[0] = temp_num;

					temp_num = this->pole_freq[i];
					this->pole_freq[i] = this->pole_freq[0];
					this->pole_freq[0] = temp_num;
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

	/* The next step is to calculate the gain.
	 * To do this we will start with K = 1
	 * We will then multiply all of the constant terms in the numerator factors
	 * and divide by the product of all of the constant terms in the denominator factors
	 * Finally we will divide this result by the scaling factor for the cutoff frequency */
	for (i = 0; i < this->quads; ++i) {
		this->K *= this->numerator[i][2];
		this->K /= this->coefficients[i][2];
	}

	this->passband_gain /= this->K;

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

	/* If the order of the filter is even then there is no linear term to deal with */
	if (n % 2 == 0) {
		this->K = 1;

		for (i = 0; i < this->quads; ++i) {
			/* Swap the first and third and coefficients of each factor in the numerator and denominator */
			temp_num = this->coefficients[i][0];
			this->coefficients[i][0] = this->coefficients[i][2];
			this->coefficients[i][2] = temp_num;

			temp_num = this->numerator[i][0];
			this->numerator[i][0] = this->numerator[i][2];
			this->numerator[i][2] = temp_num;

			/* Put the factors in the form s^2 + b*s + c starting with the denominator */
			quad_term = this->coefficients[i][0];
			this->coefficients[i][0] /= quad_term;
			this->coefficients[i][1] /= quad_term;
			this->coefficients[i][2] /= quad_term;

			this->K /= quad_term;

			/* Now do the same for the numerator */
			quad_term = this->numerator[i][0];
			this->numerator[i][0] /= quad_term;
			this->numerator[i][1] = 0;
			this->numerator[i][2] /= quad_term;

			this->K *= quad_term;
		}
	}

	/* Otherwise the order is odd and we need to deal with a linear factor */
	else {
		this->K = 1;

		/* First deal with the linear term */
		temp_num = this->coefficients[0][1];
		this->coefficients[0][1] = this->coefficients[0][2];
		this->coefficients[0][2] = temp_num;

		quad_term = this->coefficients[0][1];

		/* Put the factor in the form s + b */
		this->coefficients[0][1] = 1;
		this->coefficients[0][2] /= quad_term;

		this->K /= quad_term;

		/* Now deal with the quadratic terms */
		for (i = 1; i < this->quads; ++i) {
			/* Swap the first and third and coefficients of each factor in the numerator and denominator */
			temp_num = this->coefficients[i][0];
			this->coefficients[i][0] = this->coefficients[i][2];
			this->coefficients[i][2] = temp_num;

			temp_num = this->numerator[i][0];
			this->numerator[i][0] = this->numerator[i][2];
			this->numerator[i][2] = temp_num;

			/* Put the factors in the form s^2 + b*s + c starting with the denominator */
			quad_term = this->coefficients[i][0];
			this->coefficients[i][0] /= quad_term;
			this->coefficients[i][1] /= quad_term;
			this->coefficients[i][2] /= quad_term;

			this->K /= quad_term;

			/* Now do the same for the numerator */
			quad_term = this->numerator[i][0];
			this->numerator[i][0] = 1;
			this->numerator[i][1] = 0;
			this->numerator[i][2] /= quad_term;

			this->K *= quad_term;
		}

		this->numerator[i][0] = 0;
		this->numerator[i][1] = 1;
		this->numerator[i][2] = 0;

		this->numerator[i].swap (this->numerator[0]);
	}

	this->passband_gain *= this->K;			

	/* Now the transfer function should be fully calculated */
}
