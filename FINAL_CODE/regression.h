
#ifndef REGRESSION_H
#define REGRESSION_H

// Std headers.
#include <vector>

// GSL headers.
#include <gsl/gsl_matrix_double.h>
using namespace std;

class Regression 
{
	public:
		Regression() : num_parameters_(0), maximum_likelihood_(0) {}
		void set_response(const std::vector<int> &tot1, const std::vector<int> &met1,const std::vector<int> &tot2,const std::vector<int> &met2,const std::vector< vector<double> > &cov1,const std::vector< vector<double> > &cov2);

		std::vector<int> tot1() const {return tot1_; }
		std::vector<int> tot2() const {return tot2_; }
		std::vector<int> met1() const {return met1_; }
		std::vector<int> met2() const {return met2_; }
		std::vector< vector<double> > cov1() const {return cov1_; }
		std::vector< vector<double> > cov2() const {return cov2_; }

		void FullModel(); //Run null model first, after run FullModel, the treatment effect has been incorporated. 
		double prob(vector <double> Cov, const gsl_vector *v) const;

		double loglik(const gsl_vector *parameters) const; //log likelihood 
		int fitted_treatment_sign(const gsl_vector *parameters) const; //sign of the treatment parameters
		void gradient(const gsl_vector *parameters, gsl_vector *output) const; //first derivative

		double maximum_likelihood() { return maximum_likelihood_; }
		bool low_coverage_check();   //coverage=0;
		bool extreme_coverage_check(); //Completely methylated or completely unmethylated.
		int sign_treated_print();
	friend	bool BetaBinomialReg(Regression &r, std::vector<double> initial_parameters);

	private:
		std::vector<int> tot1_;
		std::vector<int> met1_;
		std::vector<int> tot2_;
		std::vector<int> met2_;

		vector < vector<double> > cov1_; //Covariates, each row represents an individual and each column represents a covariate
		vector < vector<double> > cov2_;

		const size_t num_parameters_;
		double maximum_likelihood_;
		int sign_treated;
};

#endif //REGRESSION_H
