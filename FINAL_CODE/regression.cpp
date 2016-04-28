#include <vector>
#include <algorithm>
#include "regression.h"
#include <gsl/gsl_multimin.h>
#include <iostream>
using std::vector;

 
void Regression::set_response(const std::vector<int> &tot1, const std::vector<int> &met1,const std::vector<int> &tot2,const std::vector<int> &met2,const std::vector< vector<double> > &cov1,const std::vector< vector<double> > &cov2)
{
	tot1_ = tot1;
	met1_ = met1;
	tot2_ = tot2;
	met2_ = met2;
	cov1_ = cov1;
	cov2_ = cov2;
	*((int*)(&num_parameters_)) =cov1_[0].size()+1;
  
}

double Regression::prob(vector <double> Cov, const gsl_vector *parameters) const
{
// The last parameter is the dispersion parameter //
	double dot_prod = 0;
	for(size_t factor = 0; factor < Cov.size(); ++factor)
		dot_prod += Cov[factor]*gsl_vector_get(parameters, factor);
	double p = exp(dot_prod)/(1 + exp(dot_prod));
	return p;
}

int Regression::fitted_treatment_sign(const gsl_vector *parameters) const
{
	int sign;
	const double para = gsl_vector_get(parameters, 0);
	if(para>=0) sign=1;
	else sign=-1;
	return sign;
}
int Regression::sign_treated_print()
{
	return sign_treated;

}
double Regression::loglik(const gsl_vector *parameters) const 
{
	double log_lik = 0;

	//dispersion parameter phi is the last element of parameter vector 
	const double dispersion_param = gsl_vector_get(parameters, num_parameters_ - 1);
	const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));
	//cout<<dispersion_param<<"    "<<phi<<endl;


	for(size_t i=0; i<tot1_.size(); i++)
	{
		int tot1_tmp=tot1_[i];
		int met1_tmp=met1_[i];
		vector < double > cov1_tmp=cov1_[i];
		
		double p_s=prob(cov1_tmp, parameters);

		for(int k=0; k< met1_tmp; k++)
		{
			log_lik=log_lik+ log((1-phi)*p_s+phi*k);
		}

		for(int k=0; k<(tot1_tmp-met1_tmp); k++)
		{		
			log_lik =log_lik + log((1 - phi)*(1 - p_s) + phi*k);
		}
		
		for(int k = 0; k < tot1_tmp; k++) 
		{
			log_lik =log_lik - log(1 + phi*(k - 1));
		}
	}

	for(size_t i=0; i<tot2_.size(); i++)
	{

		int tot2_tmp=tot2_[i];
		int met2_tmp=met2_[i];

		vector < double > cov2_tmp=cov2_[i];

		double p_s=prob(cov2_tmp, parameters);

		for(int k=0; k< met2_tmp; k++)
		{
			log_lik=log_lik+ log((1-phi)*p_s+phi*k);
		}

		for(int k=0; k<(tot2_tmp-met2_tmp); k++)
		{		
			log_lik =log_lik + log((1 - phi)*(1 - p_s) + phi*k);
		}
		
		for(int k = 0; k < tot2_tmp; k++) 
		{
			log_lik =log_lik - log(1 + phi*(k - 1));
		}
	}
	//std::cout<<log_lik<<	std::endl;
	return log_lik;
}

void Regression::gradient(const gsl_vector *parameters, gsl_vector *output) const 
{

	const double dispersion_param = gsl_vector_get(parameters, num_parameters_ - 1);
  
	const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

	for(size_t f = 0; f < num_parameters_; f++) 
	{

		double deriv = 0;

		for(size_t i=0; i<tot1_.size(); i++)
		{
			int tot1_tmp=tot1_[i];
			int met1_tmp=met1_[i];
			vector < double > cov1_tmp=cov1_[i];
			double p_s=prob(cov1_tmp, parameters);

			double term = 0;

			//a parameter linked to p p=exp(xb)/(1+exp(xb));
			if(f < (num_parameters_-1)) 
			{
				double factor = (1 - phi)*p_s*(1 - p_s)*cov1_[i][f];
				if (factor == 0) continue;

				for(int k = 0; k < met1_tmp; ++k)
					term += 1/((1 - phi)*p_s + phi*k);

				for(int k = 0; k < tot1_tmp -met1_tmp; ++k)
					term -= 1/((1 - phi)*(1 - p_s) + phi*k);

				deriv += term*factor;
			} 

			else 
			{ // the parameter linked to phi
				for(int k = 0; k < met1_tmp; ++k)
					term += (k - p_s)/((1 - phi)*p_s + phi*k);
      
				for(int k = 0; k < tot1_tmp - met1_tmp; ++k)
					term += (k - (1 - p_s))/((1 - phi)*(1 - p_s) + phi*k);

				for(int k = 0; k < tot1_tmp; ++k) 
				{
					term -= (k - 1)/(1 + phi*(k - 1));
				}
        
				deriv += term * phi * (1 - phi);
			}

		}


		for(size_t i=0; i<tot2_.size(); i++)
		{
			int tot2_tmp=tot2_[i];
			int met2_tmp=met2_[i];
			vector < double > cov2_tmp=cov2_[i];
			double p_s=prob(cov2_tmp, parameters);

			double term = 0;

			//a parameter linked to p p=exp(xb)/(1+exp(xb));
			if(f < (num_parameters_-1)) 
			{
				double factor = (1 - phi)*p_s*(1 - p_s)*cov2_[i][f];
				if (factor == 0) continue;

				for(int k = 0; k < met2_tmp; ++k)
					term += 1/((1 - phi)*p_s + phi*k);

				for(int k = 0; k < tot2_tmp- met2_tmp; ++k)
					term -= 1/((1 - phi)*(1 - p_s) + phi*k);

				deriv += term*factor;
			} 

			else 
			{ // the parameter linked to phi
				for(int k = 0; k < met2_tmp; ++k)
					term += (k - p_s)/((1 - phi)*p_s + phi*k);
      
				for(int k = 0; k < tot2_tmp - met2_tmp; ++k)
					term += (k - (1 - p_s))/((1 - phi)*(1 - p_s) + phi*k);

				for(int k = 0; k < tot2_tmp; ++k) 
				{
					term -= (k - 1)/(1 + phi*(k - 1));
				}
        
				deriv += term * phi * (1 - phi);
			}

		}
		//cout<<"Deriv"<<deriv<<endl;
		gsl_vector_set(output, f, deriv);

	}

}


void Regression::FullModel()
{
	// Incorporating treatment effect //;
	for(size_t i=0;i<cov1_.size();i++)
	{
		cov1_[i].insert(cov1_[i].begin(),1);
	}

	for(size_t i=0;i<cov2_.size();i++)
	{
		cov2_[i].insert(cov2_[i].begin(),0);
	}
	*((int*)(&num_parameters_)) =cov1_[0].size()+1;
};

bool Regression::low_coverage_check() 
{
  // Check if the site is covered by any of the samples //
	bool is_covered_in_treatment_samples = false; //sample 1;
	bool is_covered_in_control_samples = false;   //sample 2;
	
	for(size_t i=0; i<tot1_.size(); i++) 
	{
		if(tot1_[i]!=0) is_covered_in_treatment_samples =true;
	}

	for(size_t i=0; i<tot2_.size(); i++) 
	{
		if(tot2_[i]!=0) is_covered_in_control_samples =true;
	}
	return !is_covered_in_treatment_samples || !is_covered_in_control_samples;
}

bool Regression::extreme_coverage_check() 
{
	//Completely methylated or completely unmethylated.
  
	bool is_maximally_met = true;
	bool is_unmet = true;
  
	if( is_maximally_met | is_unmet)
	{
		for(size_t i=0; i<tot1_.size(); i++) 
		{
			if(tot1_[i] != met1_[i]) is_maximally_met = false;
			if(met1_[i] !=0 ) is_unmet =false;
		}
	}
	if( is_maximally_met | is_unmet)
	{
		for(size_t i=0; i<tot2_.size(); i++) 
		{
			if(tot2_[i] != met2_[i]) is_maximally_met = false;
			if(met2_[i] !=0 ) is_unmet =false;
		}
	}  
	return is_maximally_met || is_unmet;
}

