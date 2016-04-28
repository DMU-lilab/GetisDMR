#include "GetisDMR_p_value.h"


#include <iostream>
#include <string>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <armadillo>
#include <cmath>
//#include <RcppArmadillo.h>

extern std::vector<std::string> SplitString(std::string line);

GetisDMR_P_Value::GetisDMR_P_Value(string DMRinputfile, string COVinputfile, string ALLinputfile)
{
	DMRinputfile_=DMRinputfile;
	COVinputfile_=COVinputfile;
	ALLinputfile_=ALLinputfile;
}

//load data //
void GetisDMR_P_Value::load_data()
{ 
	// Vector to store all the files //
	vector <string > parts;
	string tmp;

	std::ifstream dmrinputfile(DMRinputfile_.c_str(), std::ios::in);
	if(!dmrinputfile.is_open()) 
	{
		check_input=0;
		cout<<"DMR file can't be open! " <<DMRinputfile_ <<endl;
	}
	else
	{
		while(getline(dmrinputfile,tmp))
		{
			parts=SplitString(tmp);
			DMRINPUT_.push_back(parts);
		}
	}

	std::ifstream covinputfile(COVinputfile_.c_str(), std::ios::in);
	if(!covinputfile.is_open()) 
	{
		check_input=0;
		cout<<"COV file can't be open! " <<COVinputfile_ <<endl;
	}
	else
	{
		while(getline(covinputfile,tmp))
		{
			parts=SplitString(tmp);
			vector <double> parttmp;
			for(size_t i=0; i<parts.size();i++)
			{
				parttmp.push_back(stod(parts.at(i).c_str()));
			}
			COVINPUT_.push_back(parttmp);
		}
	}

	std::ifstream allinputfile(ALLinputfile_.c_str(), std::ios::in);
	if(!covinputfile.is_open()) 
	{
		check_input=0;
		cout<<"The file output from GetisDMR can't be open! " <<ALLinputfile_ <<endl;
	}
	else
	{
		getline(allinputfile,tmp); //The first line of the file is the header. 
		while(getline(allinputfile,tmp))
		{
			parts=SplitString(tmp);
			vector <double> parttmp;
			for(size_t i=0; i<parts.size();i++)
			{
				if(i==1) parttmp.push_back(stod(parts.at(i).c_str()));
				if(i==parts.size()-1) parttmp.push_back(stod(parts.at(i).c_str()));
			}
			ALLINPUT_.push_back(parttmp);
		}
	}
	
};

void GetisDMR_P_Value::print_load_data(bool PrintDMR, bool PrintCov, bool PrintAll)
{
	if(PrintDMR)
	{
		cout<<"DMR INPUT"<<endl;
		for(size_t i=0; i<DMRINPUT_.size();i++)
		{
			for(size_t j=0; j<DMRINPUT_[i].size();j++)
			cout<<DMRINPUT_[i][j]<<"\t";
			cout<<"\n";
		}
	}

	if(PrintCov)
	{
		cout<<"COV INPUT"<<endl;
		for(size_t i=0; i<COVINPUT_.size();i++)
		{
			for(size_t j=0; j<COVINPUT_[i].size();j++)
			cout<<COVINPUT_[i][j]<<"\t";
			cout<<"\n";
		};
	}
	if(PrintAll)
	{
		cout<<"The first 100 rows of the GetisDMR output which only includes location and g-score"<<endl;
		for(size_t i=0; i<100;i++)
		{
			for(size_t j=0; j<ALLINPUT_[i].size();j++)
			cout<<ALLINPUT_[i][j]<<"\t";
			cout<<"\n";
		};
	}
};



//arma::vec GetisAsympotic()
void GetisDMR_P_Value::GetisAsympotic()
{
	size_t npos=ALLINPUT_.size();
	size_t ncov=COVINPUT_.size();
	double Maxpos=COVINPUT_.at(ncov-1).at(0);

	arma::mat pos(npos,1,arma::fill::zeros);	
	arma::mat cov(ncov,1,arma::fill::zeros);
	arma::mat g(npos,1,arma::fill::zeros);	
	for(size_t i=0; i< npos; i++)
	{
		pos(i,0)=ALLINPUT_.at(i).at(0);
		g(i,0)=ALLINPUT_.at(i).at(1);
	}
	for(size_t i=0; i < ncov; i++)
	{
		cov(i,0)=COVINPUT_.at(i).at(1);
	}	

	arma::mat dist_dif(npos,1,arma::fill::zeros);	

	for(size_t dmr_no=0; dmr_no<DMRINPUT_.size(); dmr_no++)
//	for(size_t dmr_no=0; dmr_no<20; dmr_no++)
	{
	//	cout<<"dmr_no="<<dmr_no<<endl;
		arma::mat dist_start(npos,1,arma::fill::zeros);	
		arma::mat dist_end(npos,1,arma::fill::zeros);	
		dist_start=pos.col(0)-stod(DMRINPUT_.at(dmr_no).at(1));
		dist_end=pos.col(0)-stod(DMRINPUT_.at(dmr_no).at(2));
		bool run=1;
		bool check=0;
		size_t num=0;
		double g_sum=0;
		size_t SIZE_DMR=atoi(DMRINPUT_.at(dmr_no).at(2).c_str())-atoi(DMRINPUT_.at(dmr_no).at(1).c_str())+1;
		//cout<<"SIZE="<<DMRINPUT_.at(dmr_no).at(2)<<"2        "<<DMRINPUT_.at(dmr_no).at(1)<<SIZE_DMR+1000<<endl;
		//cout<<"SIZE="<<atoi(DMRINPUT_.at(dmr_no).at(2).c_str())<<"2"<<atoi(DMRINPUT_.at(dmr_no).at(1).c_str())<<SIZE_DMR+1000<<endl;
		arma::mat weight_tmp(SIZE_DMR+Maxpos*2,1,arma::fill::zeros);	
		arma::vec weight_indices(SIZE_DMR+Maxpos*2,arma::fill::zeros);	
		arma::mat Weight_DMR;
		size_t matnumcol=0;
		size_t matnumrow=SIZE_DMR+Maxpos*2;
		//cout<<"g_sum"<<endl;
		while(run)
		{
//if(check==1) cout<<"NUM="<<num<<"run="<<run<<endl;

			if((dist_start(num,0) >=0) & (dist_end(num,0)<=0 ))
			{
//			cout<<"NUM="<<num<<endl;
				check=1;
				size_t start1=0;
				size_t end1=0;
				bool run_forward=1;
				bool run_backward=1;
				//cout<<g_sum<<endl;
				g_sum=g_sum+g(num,0);
				dist_dif=abs(pos.col(0)-pos(num,0));

				start1=end1=num;
				run_forward=1;
				run_backward=1;

				while(run_forward | run_backward)
				{
//cout<<start1<<"  "<<end1<<endl;
					if(start1==0) run_backward=0;
					if((end1==(npos-1))) run_forward=0;

					if(run_backward)
					{
						if(dist_dif(start1-1,0)>=Maxpos) run_backward=0;
						if(dist_dif(start1-1,0)<Maxpos) start1=start1-1;
					}
					if(run_forward)
					{
						if(dist_dif(end1+1,0)>=Maxpos) run_forward=0;
						if(dist_dif(end1+1,0)<Maxpos) end1=end1+1;
					}
				}

				if(end1-start1==0)
				{
					// The loci does not have any neighbors within Maxpos //
				}

				if(end1-start1>0)
				{
					matnumcol=matnumcol+1;
					arma::mat postmp = dist_dif.submat(start1, 0, end1, 0);
					arma::mat PosTmp = pos.submat(start1, 0, end1, 0); // added for asymptoics//
//					arma::mat Gtmp = g.submat(start1, 0, end1, 0);
					arma::mat var(end1-start1+1,end1-start1+1,arma::fill::zeros);
					arma::mat dist_dif_tmp(var.n_rows,1,arma::fill::zeros);	
	
					for(size_t j=0; j<var.n_rows;j++)
					{
						dist_dif_tmp=abs(postmp.col(0)-postmp(j,0));
						for(size_t k=0; k<var.n_rows; k++)
						{
							var(j,k)=cov(dist_dif_tmp(k,0),0);
						}
					}

					arma::uvec q1 = find(postmp==0);
					arma::uvec q2 = find(postmp!=0);
					arma::mat var12=var.submat(q1,q2);
					arma::mat var22=var.submat(q2,q2);

					arma::vec post=postmp.col(0);
					arma::vec pos_sel=post.elem(q2);
					// added for asymptoics//
					arma::vec PosTmpvec=PosTmp.col(0);
					arma::vec PosTmp_sel=PosTmpvec.elem(q2);
					//

					std::vector<unsigned int> sel ;
					std::vector<int> sel_pos ;
					// added for asymptoics//
					std::vector<int> position_sel_weight;

					// Resize the matrix to exclude collinear //
					for(size_t j=0;j<pos_sel.size();j++)
					{				
						if(sel.size()==0) {sel.push_back(j);sel_pos.push_back(pos_sel(j));position_sel_weight.push_back(PosTmp_sel(j));}
						if(sel.size()>0)
						{
 							if(std::find(sel_pos.begin(), sel_pos.end(), pos_sel(j)) == sel_pos.end())
							{
								sel_pos.push_back(pos_sel(j));
								sel.push_back(j);
								position_sel_weight.push_back(PosTmp_sel(j));
							}
						}
					}
		
					arma::uvec SEL = arma::conv_to<arma::uvec>::from(sel);
					arma::mat var12_final=var12.cols(SEL);

					arma::mat var22_final=var22.submat(SEL,SEL);
					arma::mat inv_var22_final=inv(var22_final) ;
					arma::mat weight=var12_final * inv_var22_final;
					//cout<<weight.n_rows<<"    "<<weight.n_cols<<endl;

					//cout<<var12_final.n_rows<<"    "<<var12_final.n_cols<<endl;
					//cout<<inv_var22_final.n_rows<<"    "<<inv_var22_final.n_cols<<endl;

					weight_tmp.zeros();
					int aaa=pos(num,0)-atoi(DMRINPUT_.at(dmr_no).at(1).c_str())+Maxpos;
					weight_tmp(aaa)=1;
					weight_indices(aaa)=1;

					for(size_t ttt=0; ttt<position_sel_weight.size(); ttt++)
					{
						//cout<<weight_tmp.n_rows<<"    "<<weight.n_rows<<"    "<<position_sel_weight[ttt]-atoi(DMRINPUT_.at(dmr_no).at(1).c_str())+500<<endl;
						weight_tmp(position_sel_weight[ttt]-atoi(DMRINPUT_.at(dmr_no).at(1).c_str())+Maxpos,0)=weight(0,ttt);
						weight_indices(position_sel_weight[ttt]-atoi(DMRINPUT_.at(dmr_no).at(1).c_str())+Maxpos)=1;
					}


					Weight_DMR.reshape(matnumrow,matnumcol);

					Weight_DMR.col(matnumcol-1)=weight_tmp;

					var.clear();
				}

				//cout<<matnumrow<<" "<<matnumcol<<endl;
		
				//Calculate the weight //
			}
			else
			{
				//if(check==1) cout<<"NUM="<<num<<"   "<<check<<endl;
				if(check==1) run=0; // already covered the chunk 
			}
			
			//cout<<"NUM1"<<endl;			
			num=num+1;
			if(num==npos) {run=0; cout<<"ERROR: CHUNK "<<dmr_no+1<<" was not in the original data, please check"<<endl;}
		}
/*
		for(size_t iii=0; iii<matnumrow; iii++)
		{
			for(size_t jjj=0; jjj<matnumcol; jjj++)
			cout<<Weight_DMR(iii,jjj)<<"\t";
			cout<<"\n";
		}
*/		//			cout<<"NUM2"<<endl;
		arma::mat Weight_DMR_final=Weight_DMR.rows(find(weight_indices==1));
		/*cout<<"DDDDDDDDDDDDD"<<endl;
		for(size_t iii=0; iii<Weight_DMR_final.n_rows; iii++)
		{
			for(size_t jjj=0; jjj<Weight_DMR_final.n_cols; jjj++)
			cout<<Weight_DMR_final(iii,jjj)<<"\t";
			cout<<"\n";
		}*/
		arma::mat VCOV(Weight_DMR_final.n_cols,Weight_DMR_final.n_cols,arma::fill::zeros);
		for(size_t iii=0; iii<Weight_DMR_final.n_cols;iii++)
		for(size_t jjj=0; jjj<Weight_DMR_final.n_cols;jjj++)
		{
			arma::mat mat1tmp=Weight_DMR_final.col(iii);
			arma::mat mat2tmp=Weight_DMR_final.col(jjj);
			double sum_weight_1=arma::accu(mat1tmp)*1.0;
			double sum_weight_sq_1=arma::accu(mat1tmp % mat1tmp)*1.0;
			double sum_weight_2=arma::accu(mat2tmp)*1.0;
			double sum_weight_sq_2=arma::accu(mat2tmp % mat2tmp)*1.0;
			double sum_weight_12=arma::accu(mat1tmp % mat2tmp)*1.0;
			double numer=npos*sum_weight_12-sum_weight_1*sum_weight_2;
			double deno=sqrt(npos*sum_weight_sq_1-sum_weight_1*sum_weight_1)*sqrt(npos*sum_weight_sq_2-sum_weight_2*sum_weight_2);
			VCOV(iii,jjj)=numer/deno*1.0;
		}
			//cout<<"NUM3"<<endl;
	//	cout<<"NUMBER"<<Weight_DMR_final.n_cols<<"g_sum"<<g_sum<<endl;
		double sds=sqrt(accu(VCOV))/(Weight_DMR_final.n_cols-1);
		double numer=g_sum/(Weight_DMR_final.n_cols);
		double stats=numer/sds;
		boost::math::normal dist(0.0, 1.0);
		double p_value=2*(1-cdf(dist,abs(stats)));
		DMRINPUT_.at(dmr_no).push_back(to_string (numer));
		DMRINPUT_.at(dmr_no).push_back(to_string (sds));
		DMRINPUT_.at(dmr_no).push_back(to_string (stats));
		DMRINPUT_.at(dmr_no).push_back(to_string (p_value));
	//	for(size_t iii=0; iii<DMRINPUT_.at(dmr_no).size();iii++)
	//	cout<<DMRINPUT_.at(dmr_no).at(iii)<<"\t";
	//	cout<<"\n";
	//	cout<<2*(1-cdf(dist,1.975))<<"          "<<p_value<<endl;
		/*for(size_t iii=0; iii<VCOV.n_rows; iii++)
		{
			for(size_t jjj=0; jjj<VCOV.n_cols; jjj++)
			cout<<VCOV(iii,jjj)<<"\t";
			cout<<"\n";
		}*/;

	}
}

void GetisDMR_P_Value::output()
{
	for(size_t dmr_no=0; dmr_no<DMRINPUT_.size(); dmr_no++)
	{
		for(size_t iii=0; iii<DMRINPUT_.at(dmr_no).size();iii++)
		{
			if(iii<3) cout<<DMRINPUT_.at(dmr_no).at(iii)<<":";
			else cout<<DMRINPUT_.at(dmr_no).at(iii)<<"\t";
		}
		cout<<"\n";
	}
}

void GetisDMR_P_Value::output(string outdmr)
{
	ofstream myfile1;
	myfile1.open (outdmr.c_str());
	if(!myfile1.is_open()) cout<<"CAN'T open the output file"<<endl;
	else
	{
		myfile1<<"chr"<<":"<<"start"<<":"<<"end"<<":"<<"length"<<"\t";
		myfile1<<"stats"<<"\t"<<"sd"<<"\t"<<"test_stats"<<"\t"<<"p_value"<<"\t"<<endl;
		for(size_t dmr_no=0; dmr_no<DMRINPUT_.size(); dmr_no++)
		{
			for(size_t iii=0; iii<DMRINPUT_.at(dmr_no).size();iii++)
			{
				if(iii<3) myfile1<<DMRINPUT_.at(dmr_no).at(iii)<<":";
				else myfile1<<DMRINPUT_.at(dmr_no).at(iii)<<"\t";
			}
			myfile1<<"\n";
		}
	}
}

