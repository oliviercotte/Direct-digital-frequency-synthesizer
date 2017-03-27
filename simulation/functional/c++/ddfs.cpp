/*
 * ddfs.cpp
 *
 *  Created on: 2017-01-01
 *      Author: oliviercotte
 */
#include "ddfs.h"

ddfs::ddfs(double fout, double fref)
: m_phase_increment(fout/fref),
  mk_min_asin_a_cos(-1.0), mk_max_asin_acos(1.0),
  mk_min_sinh_cosh(-trig_fct_gen.get_roc_hyperbolic_extended()), mk_max_sinh_cosh(-mk_min_sinh_cosh),
  mk_min_exp(0.0), mk_max_exp(trig_fct_gen.get_roc_hyperbolic_extended()),
  mk_min_atanh(-trig_fct_gen.get_arg_roc_hyperbolic_extended()), mk_max_atanh(-mk_min_atanh),
  mk_min_sin_cos(0.0), mk_max_sin_cos(2.0 * M_PI),
  mk_min_atan(trig_fct_gen.get_arg_roc_circular()), mk_max_atan(-mk_min_atan),
  mk_min_tan(-trig_fct_gen.get_arg_roc_tan()), mk_max_tan(-mk_min_tan),
  mk_min_ln(1.0), mk_max_ln(mk_max_sinh_cosh),
  mk_min_sqrt(0.0), mk_max_sqrt(4.0)
{
	reset(fout, fref);
}

ddfs::~ddfs()
{

}

void ddfs::set_fout(double fout, double fref)
{
	m_fout = fout;
	m_fref = fref;

	m_phase_increment = m_fout / m_fref; // adjust frequency

	make_phase_vec();
	make_waveform_vec();
}

void ddfs::init_phase_gen()
{
	m_phase_asinh.clear();
	m_phase_acosh.clear();
	m_phase_asin_acos.clear();
	m_phase_sinh_cosh.clear();
	m_phase_exp.clear();
	m_phase_atanh.clear();
	m_phase_sin_cos.clear();
	m_phase_atan.clear();
	m_phase_tan.clear();
	m_phase_ln.clear();
	m_phase_sqrt.clear();
}

void ddfs::init_function_gen()
{
	// Hyperbolic extended (Linear) Configuration
	m_tanh.clear();
	m_tan.clear();

	m_log.clear();
	m_sqrt.clear();

	m_sech.clear();
	m_sech_2.clear();

	m_div_tanh_approx_err.clear();
	m_div_tan_approx_err.clear();

	m_div_sech_approx_err.clear();
	m_mult_sech_2_approx_err.clear();

	m_log_approx_err.clear();
	m_sqrt_approx_err.clear();

	// Hyperbolic
	m_cosh.clear();
	m_sinh.clear();
	m_exp1.clear();
	m_exp2.clear();
	m_normh1.clear();
	m_normh2.clear();
	m_atanh1.clear();
	m_atanh2.clear();

	m_phase_sinh_cosh_approx_err.clear();
	m_phase_atanh_approx_err1.clear();
	m_phase_atanh_approx_err2.clear();

	// Circular
	m_cos.clear();
	m_sin.clear();
	m_norm1.clear();
	m_norm2.clear();
	m_atan1.clear();
	m_atan2.clear();

	m_phase_sin_cos_approx_err.clear();
	m_phase_atan_approx_err1.clear();
	m_phase_atan_approx_err2.clear();
}

void ddfs::gen_phase_vec(const double &min, const double &max, vector<double> &phase_vec)
{
	double phase, phase_step = max * m_phase_increment;
	for(phase = min; phase <= max; phase += phase_step) {
		phase_vec.push_back(phase);
	}
}

void ddfs::gen_tan()
{
	int i = 0;
	for(double& ratio : m_phase_tan) {
		m_div_tan_approx_err[i] = trig_fct_gen.get_tan(ratio, m_tan[i]);
		++i;
	}
}

void ddfs::gen_log()
{
	int i = 0;
	for(double& alpha : m_phase_ln) {
		m_log_approx_err[i] = trig_fct_gen.get_log(alpha, m_log[i]);
		++i;
	}
}

void ddfs::gen_sqrt()
{
	int i = 0;
	for(double& alpha : m_phase_sqrt) {
		m_sqrt_approx_err[i] = trig_fct_gen.get_sqrt(alpha, m_sqrt[i]);
		++i;
	}
}

void ddfs::gen_sinh_cosh()
{
	int i = 0;
	for(double& phase : m_phase_sinh_cosh) {
		m_phase_sinh_cosh_approx_err[i] = trig_fct_gen.get_cosh_sinh(phase, m_cosh[i], m_sinh[i]);
		m_phase_atanh_approx_err1[i] = trig_fct_gen.get_atanh(m_cosh[i], m_sinh[i], m_normh1[i], m_atanh1[i]);
		m_div_tanh_approx_err[i] = trig_fct_gen.get_tanh(phase, m_tanh[i]);
		m_exp1[i] = m_cosh[i] + m_sinh[i];
		m_div_sech_approx_err[i] = trig_fct_gen.get_div(m_cosh[i], 1.0, m_sech[i]);
		m_mult_sech_2_approx_err[i] = trig_fct_gen.get_mult(m_sech[i], m_sech[i], m_sech_2[i]);
		++i;
	}
}

void ddfs::gen_exp()
{
	int i = 0;
	for(double& phase : m_phase_exp) {
		trig_fct_gen.get_exp(phase, m_exp2[i]);
		++i;
	}
}

void ddfs::gen_atanh()
{
	int i = 0;
	for(double& ratio : m_phase_atanh) {
		m_phase_atanh_approx_err2[i] = trig_fct_gen.get_atanh(ratio, m_normh2[i], m_atanh2[i]);
		++i;
	}
}

void ddfs::gen_sin_cos()
{
	int i = 0;
	for(double& phase : m_phase_sin_cos) {
		m_phase_sin_cos_approx_err[i] = trig_fct_gen.get_cos_sin(phase, m_cos[i], m_sin[i]);
		m_phase_atan_approx_err1[i] = trig_fct_gen.get_atan(m_cos[i], m_sin[i], m_norm1[i], m_atan1[i]);
		++i;
	}
}

void ddfs::gen_atan()
{
	int i = 0;
	for(double& ratio : m_phase_atan) {
		m_phase_atan_approx_err2[i] = trig_fct_gen.get_atan(ratio, m_norm2[i], m_atan2[i]);
		++i;
	}
}

void ddfs::make_phase_vec()
{
	// asin, acos
	gen_phase_vec(mk_min_asin_a_cos, mk_max_asin_acos, m_phase_asin_acos);

	// sinh, cosh
	gen_phase_vec(mk_min_sinh_cosh, mk_max_sinh_cosh, m_phase_sinh_cosh);

	// exp
	gen_phase_vec(mk_min_exp, mk_max_exp, m_phase_exp);

	// atanh
	gen_phase_vec(mk_min_atanh, mk_max_atanh, m_phase_atanh);

	// sin, cos
	gen_phase_vec(mk_min_sin_cos, mk_max_sin_cos, m_phase_sin_cos);

	// atan
	gen_phase_vec(mk_min_atan, mk_max_atan, m_phase_atan);

	// tanh
	//phase_generator(mk_min_tanh, mk_max_tanh, m_phase_tanh);

	// tan
	gen_phase_vec(mk_min_tan, mk_max_tan, m_phase_tan);

	// ln
	gen_phase_vec(mk_min_ln, mk_max_ln, m_phase_ln);

	// sqrt
	gen_phase_vec(mk_min_sqrt, mk_max_sqrt, m_phase_sqrt);

	init_waveform_vec();
}

void ddfs::init_waveform_vec()
{
	// sinh, cosh []
	m_phase_sinh_cosh_approx_err.resize(m_phase_sinh_cosh.size());
	m_cosh.resize(m_phase_sinh_cosh.size());
	m_sinh.resize(m_phase_sinh_cosh.size());

	// exp []
	m_exp1.resize(m_phase_sinh_cosh.size());
	m_exp2.resize(m_phase_exp.size());

	// atanh with input xh, yh
	m_normh1.resize(m_phase_sinh_cosh.size());
	m_phase_atanh_approx_err1.resize(m_phase_sinh_cosh.size());
	m_atanh1.resize(m_phase_sinh_cosh.size());

	// atanh with input as a ratio (yh/xh)
	m_normh2.resize(m_phase_atanh.size());
	m_phase_atanh_approx_err2.resize(m_phase_atanh.size());
	m_atanh2.resize(m_phase_atanh.size());

	// tanh []
	m_tanh.resize(m_phase_sinh_cosh.size());
	m_div_tanh_approx_err.resize(m_phase_sinh_cosh.size());

	// tan []
	m_tan.resize(m_phase_tan.size());
	m_div_tan_approx_err.resize(m_phase_tan.size());

	// sech []
	m_div_sech_approx_err.resize(m_phase_sinh_cosh.size());
	m_sech.resize(m_phase_sinh_cosh.size());

	// sech_2 []
	m_mult_sech_2_approx_err.resize(m_phase_sinh_cosh.size());
	m_sech_2.resize(m_phase_sinh_cosh.size());

	// ln []
	m_log_approx_err.resize(m_phase_ln.size());
	m_log.resize(m_phase_ln.size());

	// sqrt []
	m_sqrt_approx_err.resize(m_phase_sqrt.size());
	m_sqrt.resize(m_phase_sqrt.size());

	// sin, cos
	m_phase_sin_cos_approx_err.resize(m_phase_sin_cos.size());
	m_cos.resize(m_phase_sin_cos.size());
	m_sin.resize(m_phase_sin_cos.size());

	// atan with input x, y
	m_norm1.resize(m_phase_sin_cos.size());
	m_phase_atan_approx_err1.resize(m_phase_sin_cos.size());
	m_atan1.resize(m_phase_sin_cos.size());

	// atan with as a ratio (y/x)
	m_norm2.resize(m_phase_atan.size());
	m_phase_atan_approx_err2.resize(m_phase_atan.size());
	m_atan2.resize(m_phase_atan.size());
}

void ddfs::make_waveform_vec()
{
	// Hyperbolic extended (+Linear)
	gen_sinh_cosh(); // gen_sech, gen_sech_2, gen_exp1, gen_atanh1
	gen_exp(); // gen_exp2
	gen_atanh(); // gen_atanh2,
	gen_log();
	gen_sqrt();

	// Circular
	gen_sin_cos(); // gen_atan1
	gen_atan(); // gen_atan2,
	gen_tan();
}
