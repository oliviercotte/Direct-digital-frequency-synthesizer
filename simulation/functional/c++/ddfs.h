/*
 * ddfs.h
 *
 *  Created on: 2017-01-01
 *      Author: oliviercotte
 */

#ifndef _DDFS_H_
#define _DDFS_H_

#include "cordic_algorithm.h"

#include <valarray>
#include <vector>

using std::valarray;
using std::vector;

class ddfs {
public:
	ddfs(double fout = 1.0, double fref = 500.0);
	virtual ~ddfs();

	void reset(double fout = 1.0, double fref = 500.0) {
		set_fout(fout, fref);
	}

	// Phase generators
	// Inverse trigonometric
	inline vector<double> get_phase_asin_acos() const {
		return m_phase_asin_acos;
	}
	// Hyperbolic extended (Linear)
	inline vector<double> get_phase_sech() const {
		return m_phase_sinh_cosh; // ]-inf,inf[
	}
	inline vector<double> get_phase_sech_2() const {
		return m_phase_sinh_cosh; // ]-inf,inf[
	}
	inline vector<double> get_phase_ln() const {
		return m_phase_ln; // [0, m_phase_sinh_cosh]
	}
	inline vector<double> get_phase_sqrt() const {
		return m_phase_sqrt; // [0, 2]
	}
	// Hyperbolic
	inline vector<double> get_phase_sinh_cosh() const {
		return m_phase_sinh_cosh; // (M=2[-5.0994395017623901, 5.0994395017623901]
	}
	inline vector<double> get_phase_exp() const {
		return m_phase_exp; // [-m_phase_sinh_cosh, m_phase_sinh_cosh]
	}
	inline vector<double> get_phase_atanh() const {
		return m_phase_atanh; // [-m_phase_sinh_cosh, m_phase_sinh_cosh]
	}
	inline vector<double> get_phase_tanh() const {
		return m_phase_sinh_cosh; // ]-inf,inf[
	}
	// Circular
	inline vector<double> get_phase_sin_cos() const {
		return m_phase_sin_cos; // [0, 2*M_PI]
	}
	inline vector<double> get_phase_atan() const {
		return m_phase_atan; // [-M_PI_2, M_PI_2]
	}
	inline vector<double> get_phase_tan() const {
		return m_phase_tan; // [-2, 2]
	}

	// Function generators
	// Hyperbolic extended (Linear)
	inline valarray<double> get_sech() const {
		return valarray<double>(m_tan.data(), m_tan.size());
	}
	inline valarray<double> get_sech_2() const {
		return valarray<double>(m_tan.data(), m_tan.size());
	}
	inline valarray<double> get_log() const {
		return valarray<double>(m_log.data(), m_log.size());
	}
	inline valarray<double> get_sqrt() const {
		return valarray<double>(m_tan.data(), m_tan.size());
	}
	// Hyperbolic
	inline valarray<double> get_cosh() const {
		return valarray<double>(m_cosh.data(), m_cosh.size());
	}
	inline valarray<double> get_sinh() const {
		return valarray<double>(m_sinh.data(), m_sinh.size());
	}
	inline valarray<double> get_exp1() const {
		return valarray<double>(m_exp1.data(), m_exp1.size());
	}
	inline valarray<double> get_exp2() const {
		return valarray<double>(m_exp2.data(), m_exp2.size());
	}
	inline valarray<double> get_atanh1() const {
		return valarray<double>(m_atanh1.data(), m_atanh1.size());
	}
	inline valarray<double> get_atanh2() const {
		return valarray<double>(m_atanh2.data(), m_atanh2.size());
	}
	inline valarray<double> get_tanh() const {
		return valarray<double>(m_tanh.data(), m_tanh.size());
	}
	// Circular
	inline valarray<double> get_cos() const {
		return valarray<double>(m_cos.data(), m_cos.size());
	}
	inline valarray<double> get_sin() const {
		return valarray<double>(m_sin.data(), m_sin.size());
	}
	inline valarray<double> get_atan1() const {
		return valarray<double>(m_atan1.data(), m_atan1.size());
	}
	inline valarray<double> get_atan2() const {
		return valarray<double>(m_atan2.data(), m_atan2.size());
	}
	inline valarray<double> get_tan() const {
		return valarray<double>(m_tan.data(), m_tan.size());
	}

protected:
	void set_fout(double fout, double fref);

	void init_phase_gen();
	void init_function_gen();
	void init_waveform_vec();

	void gen_phase_vec(const double &min, const double &max, vector<double> &phase_vec);

	void make_phase_vec();
	void make_waveform_vec();

	void gen_sinh_cosh();
	void gen_exp();
	void gen_atanh();
	void gen_sin_cos();
	void gen_atan();
	void gen_tan();
	void gen_log();
	void gen_sqrt();

private:
	double m_fout, m_fref;
	cordic_algorithm trig_fct_gen;

	double m_phase_increment;
	double mk_min_asin_a_cos, mk_max_asin_acos;
	double mk_min_sinh_cosh, mk_max_sinh_cosh;
	double mk_min_exp, mk_max_exp;
	double mk_min_atanh, mk_max_atanh;
	double mk_min_sin_cos, mk_max_sin_cos;
	double mk_min_atan, mk_max_atan;

	// Extension functions
	//double mk_min_tanh, mk_max_tanh;
	double mk_min_tan, mk_max_tan;
	double mk_min_ln, mk_max_ln;
	double mk_min_sqrt, mk_max_sqrt;

	// Phase generators
	vector<double> m_phase_asinh;
	vector<double> m_phase_acosh;
	vector<double> m_phase_asin_acos;
	vector<double> m_phase_sinh_cosh;
	vector<double> m_phase_exp;
	vector<double> m_phase_atanh;
	vector<double> m_phase_sin_cos;
	vector<double> m_phase_atan;
	vector<double> m_phase_tan;
	vector<double> m_phase_ln;
	vector<double> m_phase_sqrt;
	//vector<double> m_phase_tanh;

	// Hyperbolic extended (Linear)
	vector<double> m_tanh;
	vector<double> m_tan;
	vector<double> m_log;
	vector<double> m_sqrt;
	vector<double> m_sech;
	vector<double> m_sech_2;

	vector<double> m_div_tanh_approx_err;
	vector<double> m_div_tan_approx_err;
	vector<double> m_log_approx_err;
	vector<double> m_sqrt_approx_err;
	vector<double> m_div_sech_approx_err;
	vector<double> m_mult_sech_2_approx_err;

	// Hyperbolic
	vector<double> m_cosh, m_sinh;
	vector<double> m_exp1, m_exp2;
	vector<double> m_atanh1, m_normh1;
	vector<double> m_atanh2, m_normh2;

	vector<double> m_phase_sinh_cosh_approx_err;
	vector<double> m_phase_atanh_approx_err1;
	vector<double> m_phase_atanh_approx_err2;

	// Circular
	vector<double> m_cos, m_sin;
	vector<double> m_atan1, m_norm1;
	vector<double> m_atan2, m_norm2;

	vector<double> m_phase_sin_cos_approx_err;
	vector<double> m_phase_atan_approx_err1;
	vector<double> m_phase_atan_approx_err2;
};

#endif /* _DDFS_H_ */
