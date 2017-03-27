/*
 * cordic_algorithm.cpp
 *
 *  Created on: 2017-01-01
 *      Author: oliviercotte
 */

#include "cordic_algorithm.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <numeric>

///
vector<double> compute_hyperbolic_expansion_roc_vector(int M, int pipeline_size);

///
vector<double> compute_hyperbolic_expansion_roc_vector(int M, int pipeline_size) {
	int j, i;
	double sum1, sum2;
	double roc_hyperbolic_original, roc_hyperbolic_extended;
	vector<double> hyperbolic_expansion_roc_vector, atanh_lut(pipeline_size, 0.0f), atanh_expansion_lut;

	roc_hyperbolic_original = atanh(pow(2.0, -pipeline_size));
	for(i = 1; i < pipeline_size; ++i) {
		atanh_lut[i] = atanh(pow(2.0, double(-i)));
	}
	sum1 = std::accumulate(atanh_lut.begin(), atanh_lut.end(), 0.0f);
	roc_hyperbolic_original += sum1;

	for(j = 0; j <= M; ++j) {
		sum2 = roc_hyperbolic_extended = 0.0f;
		atanh_expansion_lut.resize(j+1, 0.0f);
		for(i = -j; i <= 0; ++i) {
			atanh_expansion_lut[-i] = atanh(1.0 - pow(2.0, double(i-2)));
		}

		// range of convergence
		sum2 = std::accumulate(atanh_expansion_lut.begin(), atanh_expansion_lut.end(), 0.0f);
		roc_hyperbolic_extended += roc_hyperbolic_original + sum2;
		hyperbolic_expansion_roc_vector.push_back(roc_hyperbolic_extended);

		atanh_expansion_lut.clear();
	}

	return hyperbolic_expansion_roc_vector;
}

///
cordic_algorithm::cordic_algorithm(int pipeline_size, int hyperbolic_range_expansion)
: m_pipeline_size(pipeline_size),
  m_hyperbolic_range_expansion(hyperbolic_range_expansion),
  m_iteration_num(pipeline_size-2),
  mk_scaling_hyperbolic_original(1.0f),
  mk_scaling_hyperbolic_extended(1.0f),
  mk_scaling_circular(1.0f),
  m_atanh_lut(m_hyperbolic_range_expansion+pipeline_size),
  m_atan_lut(pipeline_size),
  m_power_of_2_lut(pipeline_size),
  mk_phase_roc_hyperbolic_original(atanh(pow(2.0, -m_pipeline_size))),
  mk_phase_roc_hyperbolic_extended(0.0f),
  mk_phase_roc_circular(atan(pow(2.0, -m_pipeline_size))),
  mk_phase_roc_linear(pow(2.0, -m_pipeline_size)),
  mk_arg_roc_hyperbolic_original(0.0f),
  mk_arg_roc_hyperbolic_extended(0.0f)
{
	init_simple_iteration();
	// TODO
	//init_double_iteration();
	//init_triple_iteration();
}

cordic_algorithm::~cordic_algorithm()
{

}

void cordic_algorithm::init_simple_iteration()
{
	int i, j, M = 10;
	double sum1, sum2, K;

	for(i = 0; i < m_pipeline_size; ++i) {
		K = sqrt((1 + pow(2.0, -2.0*i)));
		m_power_of_2_lut[i] = pow(2.0, -1.0*i);
		m_atan_lut[i] = atan(m_power_of_2_lut[i]);
		mk_scaling_circular *= K;
	}

	for(j = -m_hyperbolic_range_expansion; j < m_pipeline_size; ++j) {
		i = j + m_hyperbolic_range_expansion;
		if (j <= 0) {
			K = sqrt(1.0 - pow(1.0 - pow(2.0, j-2.0), 2.0));
			m_atanh_lut[i] = atanh(1.0 - pow(2.0, j-2.0));
			mk_scaling_hyperbolic_extended *= K;
		}
		else {
			K = sqrt((1 - pow(2.0, -2.0*j)));
			m_atanh_lut[i] = atanh(pow(2.0, -j));
			mk_scaling_hyperbolic_original *= K;
		}
	}

	mk_scaling_hyperbolic = mk_scaling_hyperbolic_extended * mk_scaling_hyperbolic_original;

	// range of convergence hyperbolic
	sum1 = std::accumulate(m_atanh_lut.begin(), m_atanh_lut.begin()+m_hyperbolic_range_expansion+1, 0.0f);
	sum2 = std::accumulate(m_atanh_lut.begin()+m_hyperbolic_range_expansion+1, m_atanh_lut.end(), 0.0f);

	mk_phase_roc_hyperbolic_original += sum2;
	mk_phase_roc_hyperbolic_extended += mk_phase_roc_hyperbolic_original + sum1;

	// range of convergence circular & linear
	sum1 = std::accumulate(m_atan_lut.begin(), m_atan_lut.end(), 0.0f);
	sum2 = std::accumulate(m_power_of_2_lut.begin(), m_power_of_2_lut.end(), 0.0f);

	mk_phase_roc_circular += sum1;
	mk_phase_roc_linear += sum2;

	// argument range
	mk_arg_roc_hyperbolic_original = tanh(mk_phase_roc_hyperbolic_original);
	mk_arg_roc_hyperbolic_extended = tanh(mk_phase_roc_hyperbolic_extended);
	mk_arg_roc_circular = tan(mk_phase_roc_circular);
	//mk_arg_roc_circular = std::abs(mk_arg_roc_circular);

	mk_arg_roc_tanh = atanh(mk_phase_roc_linear);
	mk_arg_roc_tan = atan(mk_phase_roc_linear);

	m_vec_range_expansion_roc = compute_hyperbolic_expansion_roc_vector(M, m_pipeline_size);
}

void cordic_algorithm::initial_circular_rotation(
		const double &x_in, const double &y_in, const double &z_in,
		double &x_out, double &y_out, double &z_out) const
{
	// Assume z_in = [0, 2*PI]
	// z_out = [-PI/2, PI/2] (nominal cordic range of convergence in circular mode)
	// (check mk_roc_circular)
	if((z_in > M_PI/2) && (z_in <= M_PI)) {
		x_out = -x_in;
		y_out = -y_in;
		z_out = z_in - M_PI;
	}
	else if((z_in >= M_PI) && (z_in < 3*M_PI/2)) {
		x_out = -x_in;
		y_out = -y_in;
		z_out =  3*M_PI/2 - z_in;
	}
	else if (z_in >= 3*M_PI/2) {
		z_out =  z_in - 2*M_PI;
	}
	else {
		z_out = z_in;
	}
}

void cordic_algorithm::initial_circular_vectoring(
		const double &x_in, const double &y_in, const double &z_in,
		double &x_out, double &y_out, double &z_out) const
{
	// tan is defined for [-PI/2, PI/2] -> x_out >= 0
	if(x_in < 0.0) {
		int d = y_in < 0.0 ? +1 : -1;
		x_out = -d * y_in;
		y_out = d * x_in;
		z_out = z_in + d * M_PI/2;
	}
}

void cordic_algorithm::cfg_hyperbolic_extended(
		const mode  &m,
		const double &x_in, const double &y_in, const double &z_in,
		double &x_out, double &y_out, double &z_out) const
{
	int i, j, d;
	double xtmp, xtmp1, ytmp, ztmp;
	double x_next = x_in, y_next = y_in, z_next = z_in;
	vector <double> x(m_hyperbolic_range_expansion + m_pipeline_size, 0.0);
	vector <double> y(x), z(x);

	x[0] = x_next;
	y[0] = y_next;
	z[0] = z_next;

	for (j = -m_hyperbolic_range_expansion; j <= m_iteration_num; ++j) {
		i = j + m_hyperbolic_range_expansion;

		if(m == ROTATION) {
			d = z[i] < 0.0 ? -1.0 : +1.0;
		}
		else {
			d = ((x[i] * y[i]) >= 0.0) ? -1.0 : +1.0;
		}

		if(j <= 0) {
			x[i+1] = x[i] + d * (1.0 - pow(2.0, j-2)) * y[i];
			y[i+1] = y[i] + d * (1.0 - pow(2.0, j-2)) * x[i];
			z[i+1] = z[i] - d * m_atanh_lut[i];
		}
		else {
			xtmp = x[i] + d * y[i] * pow(2.0, -j);
			ytmp = y[i] + d * x[i] * pow(2.0, -j);
			ztmp = z[i] - d * m_atanh_lut[i];

			// In hyperbolic mode, iterations 4, 13, 40, 121, ..., j, 3j+1,... must be repeated
			if ((j == 4) || (j == 13) ||(j == 40) || (j == 121)) {
				xtmp1 = xtmp;
				if(m == ROTATION) {
					d = ztmp < 0.0 ? -1.0 : +1.0;
				}
				else {
					d = ((xtmp * ytmp) >= 0.0) ? -1.0 : +1.0;
				}

				xtmp = xtmp + d * ytmp * pow(2.0, -j);
				ytmp = ytmp + d * xtmp1 * pow(2.0, -j);
				ztmp = ztmp - d * m_atanh_lut[i];
			}

			x[i+1] = xtmp;
			y[i+1] = ytmp;
			z[i+1] = ztmp;
		}
	}

	x_out = x[m_hyperbolic_range_expansion + m_iteration_num + 1];
	y_out = y[m_hyperbolic_range_expansion + m_iteration_num + 1];
	z_out = z[m_hyperbolic_range_expansion + m_iteration_num + 1];
}

void cordic_algorithm::cfg_circular(
		const mode &m,
		const double &x_in, const double &y_in, const double &z_in,
		double &x_out, double &y_out, double &z_out) const
{
	int i, d;
	double x_next = x_in, y_next = y_in, z_next = z_in;
	vector <double> x(m_pipeline_size, 0.0);
	vector <double> y(x), z(x);

	// assuming a 2pi modulus on the phase input
	// the unified cordic iteration algorithm is defined for [-pi/2,pi/2]
	if(m == ROTATION) {
		initial_circular_rotation(x_in, y_in, z_in, x_next, y_next, z_next);
	}
	else {
		initial_circular_vectoring(x_in, y_in, z_in, x_next, y_next, z_next);
	}

	x[0] = x_next;
	y[0] = y_next;
	z[0] = z_next;

	for(i = 0; i <= m_iteration_num; ++i) {
		if(m == ROTATION) {
			d = z[i] < 0.0 ? -1.0 : +1.0;
		}
		else {
			d = y[i] < 0.0 ? +1.0 : -1.0;
		}

		x[i+1] = x[i] - y[i] * d * pow(2.0, -i);
		y[i+1] = y[i] + x[i] * d * pow(2.0, -i);
		z[i+1] = z[i] - d * m_atan_lut[i];
	}

	x_out = x[m_iteration_num+1];
	y_out = y[m_iteration_num+1];
	z_out = z[m_iteration_num+1];
}

void cordic_algorithm::cfg_linear(
		const mode  &m,
		const double &x_in, const double &y_in, const double &z_in,
		double &x_out, double &y_out, double &z_out) const
{
	int i, d;
	double x_next = x_in, y_next = y_in, z_next = z_in;
	vector <double> x(m_pipeline_size, 0.0);
	vector <double> y(x), z(x);

	x[0] = x_next;
	y[0] = y_next;
	z[0] = z_next;

	for(i = 0; i <= m_iteration_num; ++i) {
		if(m == ROTATION) {
			d = z[i] < 0.0 ? -1.0 : +1.0;
		}
		else {
			d = y[i] < 0.0 ? +1.0 : -1.0;
		}

		x[i+1] = x[i];
		y[i+1] = y[i] + d * x[i] * pow(2.0, -i);
		z[i+1] = z[i] - d * m_power_of_2_lut[i];
	}

	x_out = x[m_iteration_num+1];
	y_out = y[m_iteration_num+1];
	z_out = z[m_iteration_num+1];
}

void cordic_algorithm::unified_cordic(
		const coordinate_system &cs, const mode  &m,
		const double &x_in, const double &y_in, const double &z_in,
		double &x_out, double &y_out, double &z_out) const
{
	switch(cs) {
	case HYPERBOLIC:
		cfg_hyperbolic_extended(m, x_in, y_in, z_in, x_out, y_out, z_out);
		break;
	case LINEAR:
		cfg_linear(m, x_in, y_in, z_in, x_out, y_out, z_out);
		break;
	case CIRCULAR:
		cfg_circular(m, x_in, y_in, z_in, x_out, y_out, z_out);
		break;
	default:
		break;
	}
}

double cordic_algorithm::get_tanh(const double &arg, double &tanh_arg_out) const
{
	double tanh_approx_err;
	_get_tanh(arg, tanh_arg_out, tanh_approx_err);
	return tanh_approx_err;
}

double cordic_algorithm::get_tan(const double arg, double &tan_arg_out) const
{
	double tan_approx_err;
	_get_tan(arg, tan_arg_out, tan_approx_err);
	return tan_approx_err;
}

double cordic_algorithm::get_log(const double &alpha, double &log_alpha_out) const
{
	double log_approx_err;
	_get_log(alpha, log_alpha_out, log_approx_err);
	return log_approx_err;
}

double cordic_algorithm::get_sqrt(const double &alpha, double &sqrt_alpha_out) const
{
	double sqrt_approx_err;
	_get_sqrt(alpha, sqrt_alpha_out, sqrt_approx_err);
	return sqrt_approx_err;
}

double cordic_algorithm::get_cosh_sinh(const double &phase_in, double &cosh_out, double &sinh_out) const
{
	double phase_approx_err;
	_get_cosh_sinh(phase_in, cosh_out, sinh_out, phase_approx_err);
	return phase_approx_err;
}

double cordic_algorithm::get_atanh(const double &ratio, double &normh_out, double &atanh_out) const
{
	double ratio_approx_err;
	_get_atanh(ratio, normh_out, ratio_approx_err, atanh_out);
	return ratio_approx_err;
}

double cordic_algorithm::get_atanh(const double &xh_in, const double &yh_in, double &normh_out, double &atanh_out) const
{
	double ratio_approx_err;
	_get_atanh(xh_in, yh_in, normh_out, ratio_approx_err, atanh_out);
	return ratio_approx_err;
}

double cordic_algorithm::get_exp(const double &phase_in, double &exp_out) const
{
	double phase_approx_err;
	_get_exp(phase_in, exp_out, phase_approx_err);
	return phase_approx_err;
}

double cordic_algorithm::get_cos_sin(const double &phase_in, double &cosine_out, double &sine_out) const
{
	double phase_approx_err;
	_get_cos_sin(phase_in, cosine_out, sine_out, phase_approx_err);
	return phase_approx_err;
}

double cordic_algorithm::get_atan(const double &ratio, double &norm_out, double &atan_out) const
{
	double ratio_approx_err;
	_get_atan(ratio, norm_out, ratio_approx_err, atan_out);
	return ratio_approx_err;
}

double cordic_algorithm::get_atan(const double &x_in, const double &y_in, double &norm_out, double &atan_out) const
{
	double ratio_approx_err;
	_get_atan(x_in, y_in, norm_out, ratio_approx_err, atan_out);
	return ratio_approx_err;
}

double cordic_algorithm::get_mult(const double &x, const double &z, double &mult_out) const
{
	// z = x * z;
	double mult_approx_err;
	_mult(x, z, mult_out, mult_approx_err);
	return mult_approx_err;
}

double cordic_algorithm::get_div(const double &x, const double &y, double &div_out) const
{
	// z = y / x;
	double div_approx_err;
	_div(x, y, div_approx_err, div_out);
	return div_approx_err;
}

void cordic_algorithm::_mult(
		const double &x_in, const double &z_in,
		double &y_out, double &z_out) const
{
	// y_out = x_in * z_in;
	double dummy;
	unified_cordic(
			coordinate_system::LINEAR, mode::ROTATION,
			x_in, 0.0f, z_in,
			dummy, y_out, z_out);
}

void cordic_algorithm::_div(
		const double &x_in, const double &y_in,
		double &y_out, double &z_out) const
{
	// z_out = y_in / x_in;
	double dummy;
	unified_cordic(
			coordinate_system::LINEAR, mode::VECTORING,
			x_in, y_in, 0.0f,
			dummy, y_out, z_out);
}

void cordic_algorithm::_get_tanh(
		const double &xh_in,
		const double &yh_in,
		double &tanh_out,
		double &tanh_approx_err_out) const
{
	// tan_out = y_in / x_in
	_div(xh_in, yh_in, tanh_approx_err_out, tanh_out);
}

void cordic_algorithm::_get_tanh(
		const double &ratio_in,
		double &tanh_out,
		double &tanh_approx_err_out) const
{
	double dummy, x_out, y_out;
	_get_cosh_sinh(ratio_in, x_out, y_out, dummy);
	_get_tanh(x_out, y_out, tanh_out, tanh_approx_err_out);
}

void cordic_algorithm::_get_tan(
		const double &x_in,
		const double &y_in,
		double &tan_out,
		double &tan_approx_err_out) const
{
	// tan_out = y_in / x_in
	_div(x_in, y_in, tan_approx_err_out, tan_out);
}

void cordic_algorithm::_get_tan(
		const double &ratio_in,
		double &tan_out,
		double &tan_approx_err_out) const
{
	// tan_out = y_in / x_in
	double dummy, x_out, y_out;
	_get_cos_sin(ratio_in, x_out, y_out, dummy);
	_get_tan(x_out, y_out, tan_out, tan_approx_err_out);
}

void cordic_algorithm::_get_sech(
		const double &phase_in,
		double &sech_out,
		double &approx_err) const
{
	double cosh_out, sinh_out, dummy;
	_get_cosh_sinh(phase_in, cosh_out, sinh_out, dummy);
	_div(cosh_out, 1.0, approx_err, sech_out);
}

void cordic_algorithm::_get_sech_2(
		const double &phase_in,
		double &sech_2_out,
		double &approx_err) const
{
	double cosh_out, sinh_out, sech_out, dummy;
	_get_cosh_sinh(phase_in, cosh_out, sinh_out, dummy);
	_div(cosh_out, 1.0, approx_err, sech_out);
	_mult(cosh_out, sinh_out, approx_err, sech_2_out);
}

void cordic_algorithm::_get_log(
		const double &alpha_in,
		double &log_out,
		double &approx_err) const
{
	double x_in = alpha_in + 1.0, y_in = alpha_in - 1.0, dummy;
	_get_atanh(x_in, y_in, dummy, approx_err, log_out);
	log_out *= 2.0;
}

void cordic_algorithm::_get_sqrt(
		const double &alpha_in,
		double &sqrt_diff_out,
		double &approx_err) const
{
	double x_in = alpha_in + 0.25, y_in = alpha_in - 0.25, dummy, sqrt_diff_out_unscaled;
	_get_atanh(x_in, y_in, sqrt_diff_out_unscaled, approx_err, dummy);
	_div(mk_scaling_hyperbolic, sqrt_diff_out_unscaled, approx_err, sqrt_diff_out);
}

void cordic_algorithm::_get_cosh_sinh(const double &phase_in, double &cosh_out, double &sinh_out, double &phase_out) const
{
	double x_in = 1.0 / mk_scaling_hyperbolic;
	double y_in = 0.0;
	double z_in = phase_in;

	unified_cordic(
			coordinate_system::HYPERBOLIC, mode::ROTATION,
			x_in, y_in, z_in,
			cosh_out, sinh_out, phase_out);
}

void cordic_algorithm::_get_exp(const double &phase_in, double &exp_out, double &phase_out) const
{
	double cosh_out, sinh_out;
	_get_cosh_sinh(phase_in, cosh_out, sinh_out, phase_out);
	exp_out = cosh_out + sinh_out;
}

void cordic_algorithm::_get_atanh(
		const double &xh_in,
		const double &yh_in,
		double &normh_out,
		double &ratio_out,
		double &phase_out) const
{
	double z_in = 0.0;

	// computed using the vectoring mode CORDIC rotator if the angle
	// accumulator is initialized with zero.
	// the atanh result is taken from the angle accumulator
	unified_cordic(
			coordinate_system::HYPERBOLIC, mode::VECTORING,
			xh_in, yh_in, z_in,
			normh_out, ratio_out, phase_out);
}

void cordic_algorithm::_get_atanh(
		const double &ratio_in,
		double &norm_out,
		double &ratio_out,
		double &phase_out) const
{
	double x_in = 1.0, y_in = ratio_in, z_in = 0.0;

	// computed using the vectoring mode CORDIC rotator if the angle
	// accumulator is initialized with zero.
	// the atanh result is taken from the angle accumulator
	unified_cordic(
			coordinate_system::HYPERBOLIC, mode::VECTORING,
			x_in, y_in, z_in,
			norm_out, ratio_out, phase_out);
}

void cordic_algorithm::_get_cos_sin(
		const double &phase_in,
		double &cosine_out, double &sine_out, double &phase_out) const
{
	const double x_in = 1.0f / mk_scaling_circular;
	const double y_in = 0.0f;
	const double z_in = phase_in;

	unified_cordic(
			coordinate_system::CIRCULAR, mode::ROTATION,
			x_in, y_in, z_in,
			cosine_out, sine_out, phase_out);
}

void cordic_algorithm::_get_atan(
		const double &x_in,
		const double &y_in,
		double &norm_out,
		double &ratio_out,
		double &phase_out) const
{
	double z_in = 0.0;

	// computed using the vectoring mode CORDIC rotator if the angle
	// accumulator is initialized with zero.
	// the atan result is taken from the angle accumulator
	unified_cordic(
			coordinate_system::CIRCULAR, mode::VECTORING,
			x_in, y_in, z_in,
			norm_out, ratio_out, phase_out);
}

void cordic_algorithm::_get_atan(
		const double &ratio_in,
		double &norm_out,
		double &ratio_out,
		double &phase_out) const
{
	double x_in = 1.0, y_in = ratio_in, z_in = 0.0;

	// computed using the vectoring mode CORDIC rotator if the angle
	// accumulator is initialized with zero.
	// the atan result is taken from the angle accumulator
	unified_cordic(
			coordinate_system::CIRCULAR, mode::VECTORING,
			x_in, y_in, z_in,
			norm_out, ratio_out, phase_out);
}
