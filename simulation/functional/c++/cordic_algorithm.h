/*
 * cordic_algorithm.h
 *
 *  Created on: 2017-01-01
 *      Author: oliviercotte
 */

#ifndef _CORDIC_ALGORITHM_H_
#define _CORDIC_ALGORITHM_H_

#include <vector>
using std::vector;

class cordic_algorithm {
public:
	enum coordinate_system {
		HYPERBOLIC = -1,
		LINEAR = 0,
		CIRCULAR = 1
	};

	enum mode {
		ROTATION = 0,
		VECTORING
	};

	enum function_type {
		COSH,
		SINH,
		EXP,
		ATANH,
		COS,
		SIN,
		ATAN,
		SECH,
		SECH_2
	};

	cordic_algorithm(int pipeline_size = 18, int hyperbolic_range_expansion = 2);
	virtual ~cordic_algorithm();

	// Linear configuration
	double get_mult(const double &x, const double &z, double &mult_out) const;
	double get_div(const double &x, const double &y, double &div_out) const;

	// Hyperbolic extended (Linear) configuration
	double get_tanh(const double &arg, double &tanh_arg_out) const;
	double get_tan(const double arg, double &tan_arg_out) const;
	double get_sech(const double &x, double &sech_out) const;
	double get_sech_2(const double &x, double &sech_2_out) const;
	double get_log(const double &alpha, double &log_alpha_out) const;
	double get_sqrt(const double &alpha, double &sqrt_alpha_out) const;

	// Hyperbolic configuration
	double get_cosh_sinh(const double &phase_in, double &cosh_out, double &sinh_out) const;
	double get_exp(const double &phase_in, double &exp_out) const;
	double get_atanh(const double &ratio, double &normh_out, double &atanh_out) const;
	double get_atanh(const double &xh_in, const double &yh_in, double &normh_out, double &atanh_out) const;

	// Circular configuration
	double get_cos_sin(const double &phase_in, double &cosine_out, double &sine_out) const;
	double get_atan(const double &ratio, double &norm_out, double &atan_out) const;
	double get_atan(const double &x_in, const double &y_in, double &norm_out, double &atan_out) const;

	// roc information
	inline double get_roc_hyperbolic_original() const {
		return mk_phase_roc_hyperbolic_original;
	}
	inline double get_roc_hyperbolic_extended() const {
		return mk_phase_roc_hyperbolic_extended;
	}
	inline double get_roc_circular() const {
		return mk_phase_roc_circular;
	}
	inline double get_roc_linear() const {
		return mk_phase_roc_linear;
	}

	// arg information
	inline double get_arg_roc_hyperbolic_extended() const {
		return mk_arg_roc_hyperbolic_extended;
	}
	inline double get_arg_roc_hyperbolic_original() const {
		return mk_arg_roc_hyperbolic_original;
	}
	inline double get_arg_roc_circular() const {
		return mk_arg_roc_circular;
	}
	inline double get_arg_roc_tanh() const {
		return mk_arg_roc_tanh;
	}
	inline double get_arg_roc_tan() const {
		return mk_arg_roc_tan;
	}

protected:
	void init_simple_iteration();

	void initial_circular_rotation(
			const double &x_in, const double &y_in, const double &z_in,
			double &x_out, double &y_out, double &z_out) const;

	void initial_circular_vectoring(
			const double &x_in, const double &y_in, const double &z_in,
			double &x_out, double &y_out, double &z_out) const;

	void cfg_hyperbolic_extended(
			const mode &m,
			const double &x_in, const double &y_in, const double &z_in,
			double &x_out, double &y_out, double &z_out) const;

	void cfg_circular(
			const mode &m,
			const double &x_in, const double &y_in, const double &z_in,
			double &x_out, double &y_out, double &z_out) const;

	void cfg_linear(
			const mode &m,
			const double &x_in, const double &y_in, const double &z_in,
			double &x_out, double &y_out, double &z_out) const;

	void unified_cordic(
			const coordinate_system &cs, const mode &m,
			const double &x_in, const double &y_in, const double &z_in,
			double &x_out, double &y_out, double &z_out) const;

	// Linear configuration
	void _mult(
			const double &x_in, const double &z_in,
			double &y_out, double &z_out) const;
	void _div(
			const double &x_in, const double &y_in,
			double &y_out, double &z_out) const;

	// Hyperbolic extended (Linear) configuration
	void _get_tanh(
			const double &xh_in,
			const double &yh_in,
			double &tanh_out,
			double &tanh_approx_err_out) const;

	void _get_tanh(
			const double &ratio_in,
			double &tanh_out,
			double &tanh_approx_err_out) const;

	void _get_tan(
			const double &x_in,
			const double &y_in,
			double &tan_out,
			double &div_approx_err_out) const;

	void _get_tan(
			const double& ratio_in,
			double& tan_out,
			double& div_approx_err_out) const;

	void _get_sech(
			const double &phase_in,
			double &sech_out,
			double &approx_err) const;

	void _get_sech_2(
			const double &phase_in,
			double &sech_2_out,
			double &approx_err) const;

	void _get_log(
			const double &alpha_in,
			double &ln_out,
			double &approx_err) const;

	void _get_sqrt(
			const double &alpha_in,
			double &sqrt_diff_out,
			double &approx_err) const;

	// Hyperbolic configuration
	void _get_cosh_sinh(
			const double &phase_in,
			double &cosh_out,
			double &sinh_out,
			double &phase_out) const;

	void _get_exp(
			const double &phase_in,
			double &exp_out,
			double &phase_out) const;

	void _get_atanh(
			const double &xh_in,
			const double &yh_in,
			double &normh_out,
			double &ratio_out,
			double &phase_out) const;

	void _get_atanh(
			const double &ratio_in,
			double &ratio_out,
			double &normh_out,
			double &phase_out) const;

	// Circular configuration
	void _get_cos_sin(
			const double &phase_in,
			double &cosine_out,
			double &sine_out,
			double &phase_out) const;

	void _get_atan(
			const double &x_in,
			const double &y_in,
			double &norm_out,
			double &ratio_out,
			double &phase_out) const;

	void _get_atan(
			const double &ratio_in,
			double &norm_out,
			double &ratio_out,
			double &phase_out) const;

private:
	// iteration configuration
	int m_pipeline_size;
	int m_hyperbolic_range_expansion;
	int m_iteration_num;

	// constant scaling factor (the gain for the linear configuration is always 1)
	double mk_scaling_hyperbolic_original;
	double mk_scaling_hyperbolic_extended;
	double mk_scaling_hyperbolic;
	double mk_scaling_circular;

	// look-up table
	vector<double> m_atanh_lut;
	vector<double> m_atan_lut;
	vector<double> m_power_of_2_lut;

	// phase range of convergence (input phase limitation)
	double mk_phase_roc_hyperbolic_original;
	double mk_phase_roc_hyperbolic_extended;
	double mk_phase_roc_circular;
	double mk_phase_roc_linear;

	// argument range of convergence (tan, tanh ratio)
	double mk_arg_roc_hyperbolic_original;
	double mk_arg_roc_hyperbolic_extended;
	double mk_arg_roc_circular;

	double mk_arg_roc_tanh;
	double mk_arg_roc_tan;

	// Additional phase range of convergence from 0..M iteration for hyperbolic extended case
	vector<double> m_vec_range_expansion_roc;
};

#endif // _CORDIC_ALGORITHM_H_
