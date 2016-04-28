#include <iostream>
#include <string>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
//#include <armadillo>
#include <cmath>
#include <RcppArmadillo.h>

arma::vec MyDist(Rcpp::NumericMatrix Pos, Rcpp::NumericMatrix Z, Rcpp::NumericMatrix COV,Rcpp::NumericVector x_bar, Rcpp::NumericVector s, Rcpp::NumericMatrix G,Rcpp::NumericVector maxpos)
{
/* Make sure the Pos, Z and G has the same number of rows */
	int npos=Pos.nrow();
	int nz=Z.nrow();
	int ng=G.nrow();

	double X_bar=x_bar[0];
	double S=s[0];
	double Maxpos=maxpos[0];

	if((npos != nz) |(nz!=ng)) 	Rcpp::Rcout<<"Wrong Rows"<<std::endl;
	
	arma::mat pos = Rcpp::as<arma::mat>(Pos);
	arma::mat z = Rcpp::as<arma::mat>(Z);
	arma::mat cov = Rcpp::as<arma::mat>(COV);
	arma::mat g = Rcpp::as<arma::mat>(G);

	arma::mat dist_dif(npos,1,arma::fill::zeros);	

	int start1=0;
	int end1=0;
	bool run_forward=TRUE;
	bool run_backward=TRUE;
	double numer=1.0;
	double Deno=1.0;
	for(int i=0; i<npos; i++)
	{
		//Rcpp::Rcout<<i<<std::endl;
		dist_dif=abs(pos.col(0)-pos(i,0));
		start1=end1=i;
		run_forward=TRUE;
		run_backward=TRUE;

		while(run_forward | run_backward)
		{
			if(start1==0) run_backward=FALSE;
			if((end1==(npos-1))) run_forward=FALSE;

			if(run_backward)
			{
				if(dist_dif(start1-1,0)>=Maxpos) run_backward=FALSE;
				if(dist_dif(start1-1,0)<Maxpos) start1=start1-1;

			}
			if(run_forward)
			{
				if(dist_dif(end1+1,0)>=Maxpos) run_forward=FALSE;
				if(dist_dif(end1+1,0)<Maxpos) end1=end1+1;
			}
		}

		if(end1-start1==0)
		{
			numer=Z(i,0)-X_bar;
			Deno=S;
		}
		if(end1-start1>0)
		{
			arma::mat postmp = dist_dif.submat(start1, 0, end1, 0);
			arma::mat Ztmp = z.submat(start1, 0, end1, 0);
			arma::mat var(end1-start1+1,end1-start1+1,arma::fill::zeros);
			arma::mat dist_dif_tmp(var.n_rows,1,arma::fill::zeros);	
	
			for(int j=0; j<var.n_rows;j++)
			{
				dist_dif_tmp=abs(postmp.col(0)-postmp(j,0));
				for(int k=0; k<var.n_rows; k++)
				{
					var(j,k)=cov(dist_dif_tmp(k,0),0);
				//	Rcpp::Rcout<<var(j,k)<<"  ";
				}
				//Rcpp::Rcout<<std::endl;
			}
			
			arma::uvec q1 = find(postmp==0);
			arma::uvec q2 = find(postmp!=0);
			arma::mat var12=var.submat(q1,q2);
			arma::mat var22=var.submat(q2,q2);

			arma::vec post=postmp.col(0);
			arma::vec pos_sel=post.elem(q2);
			arma::vec z_st=Ztmp.col(0);
			arma::vec z_sel=z_st.elem(q2);


			std::vector<unsigned int> sel ;
			std::vector<int> sel_pos ;
			// Resize the matrix to exclude collinear //
			for(int j=0;j<pos_sel.size();j++)
			{				
				if(sel.size()==0) {sel.push_back(j);sel_pos.push_back(pos_sel(j));}
				if(sel.size()>0)
				{
 					if(std::find(sel_pos.begin(), sel_pos.end(), pos_sel(j)) == sel_pos.end())
					{
						sel_pos.push_back(pos_sel(j));
						sel.push_back(j);
					}
				}
			}
		
			arma::uvec SEL = arma::conv_to<arma::uvec>::from(sel);
			arma::mat var12_final=var12.cols(SEL);

			arma::mat var22_final=var22.submat(SEL,SEL);

			/*for(int k=0;k<var22_final.n_rows;k++)
			{
				for(int f=0;f<var22_final.n_cols;f++)
				Rcpp::Rcout<<var22_final(k,f)<<" ";
				Rcpp::Rcout<<std::endl;
			}*/

			arma::vec z_sel_final=z_sel.elem(SEL);

			arma::mat inv_var22_final=inv(var22_final) ;
			arma::mat Gtmp=var12_final * inv_var22_final * z_sel_final;
			arma::mat weight=var12_final * inv_var22_final;
			
		
			//arma::vec weight_vec(weight);
			double GG=Gtmp(0,0)+z(i,0);
			double sum_weight=arma::accu(weight)*1.0;
			double sum_weight_2=arma::accu(weight % weight)*1.0;
			double detmp=npos*(sum_weight_2+1)*1.0-(sum_weight+1)*(sum_weight+1)*1.0;
			detmp=detmp/(npos-1);
			detmp=sqrt(detmp);
			Deno=detmp*S;
			numer=GG-X_bar*(sum_weight+1);
			var.clear();
		}
		g(i,0)=numer/Deno*1.0;

	}
	
	return g;
}

RCPP_MODULE(MyDist)
{
	function("MyDist",&MyDist);
}
	

