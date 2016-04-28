// STD headers.
#include <vector>
#include <stdexcept>

// GSL headers.
#include <gsl/gsl_multimin.h>

// Local headers.
#include "regression.h"
#include "BetaBinomialReg.h"

using std::vector;

static double neg_loglik_proxy (const gsl_vector *parameters, void *object) //minimize -loglikelihood= max likelihood
{
	Regression *regression = (Regression *)(object);
	return (-1)*regression->loglik(parameters);
}

static void neg_gradient_proxy (const gsl_vector *parameters, void *object, gsl_vector *d_loglik_val) //first derivative for -loglikelihood function 
{
	Regression *regression = (Regression *)(object);
	regression->gradient(parameters, d_loglik_val);
	gsl_vector_scale(d_loglik_val, -1.0);
}

static void neg_loglik_and_gradient_proxy (const gsl_vector *parameters, void *object, double *loglik_val, gsl_vector *d_loglik_val) //-log(likelihood) and -d(log(likelihood))/d(x)
{
	*loglik_val = neg_loglik_proxy(parameters, object);
	neg_gradient_proxy(parameters, object, d_loglik_val);
}

bool BetaBinomialReg(Regression &r, vector<double> initial_parameters) 
{
	if (initial_parameters.empty()) 
	{
		for(size_t ind = 0; ind < r.num_parameters_ - 1; ++ind) initial_parameters.push_back(0.0);
		initial_parameters.push_back(-2.5);
	}
  
	if (initial_parameters.size() != r.num_parameters_)
		throw std::runtime_error("Wrong number of initial parameters.");

	int status = 0;
	int iter = 0;

	gsl_multimin_function_fdf loglik_bundle;

	loglik_bundle.f = &neg_loglik_proxy;
	loglik_bundle.df = &neg_gradient_proxy;
	loglik_bundle.fdf = &neg_loglik_and_gradient_proxy;
	loglik_bundle.n = r.num_parameters_;
	loglik_bundle.params = (void *)&r;

	gsl_vector *parameters = gsl_vector_alloc(r.num_parameters_);

	for (size_t parameter = 0; parameter < initial_parameters.size(); ++parameter) 
	{
		gsl_vector_set(parameters, parameter, initial_parameters[parameter]);
	}

	const gsl_multimin_fdfminimizer_type *T;
  
	T = gsl_multimin_fdfminimizer_conjugate_fr;
  
	gsl_multimin_fdfminimizer *s; 
	s = gsl_multimin_fdfminimizer_alloc (T, r.num_parameters_);
  
	gsl_multimin_fdfminimizer_set (s, &loglik_bundle, parameters, 0.001, 1e-4);

	do 
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);

		if (status) break;

		status = gsl_multimin_test_gradient (s->gradient, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 700);

	r.maximum_likelihood_ = r.loglik(s->x);
	r.sign_treated = r.fitted_treatment_sign(s->x);
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(parameters);

	return status == GSL_SUCCESS;
}
