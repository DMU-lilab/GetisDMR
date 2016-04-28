#include "multi_vs_multi.h"
//#include "one_vs_one.h"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <iostream>
#include <string>
#include <fstream>
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
#include "SplitString.h"
#include "compare2Prop.h"
#include "regression.h"
#include "BetaBinomialReg.h"
#include "Getis_Ord.h"
#include "FindDMR.h"
#include <Rcpp.h>
#include <RInside.h>
#include <gsl/gsl_matrix.h>
#include <armadillo>

extern std::vector<std::string> SplitString(std::string line);
extern string OutCovFile;
extern string CLINK_CPPFLAG;
extern string MydistSoftware;
extern string NewDistCppSoftware;
extern string FindDMRfile;
extern int MINTOT;
extern RInside RFun;

using namespace std;



// constructor //

MultivsMulti::MultivsMulti(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int, bool sorted_int, bool covinc_int)
{
	input1=input1_int;
	input2=input2_int;
	tolerance=tol_int;
	length=len_int;
	cutoff=cut_int;                                    
	outdir=outdir_int;
	sorted=sorted_int; 
	covinc=covinc_int;
}

MultivsMulti::MultivsMulti(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int, bool sorted_int, bool covinc_int, string covinput1_int, string covinput2_int)
{
	input1=input1_int;
	input2=input2_int;
	tolerance=tol_int;
	length=len_int;
	cutoff=cut_int;                                    
	outdir=outdir_int;
	sorted=sorted_int; 
	covinc=covinc_int;
	input1cov=covinput1_int;
	input2cov=covinput2_int;
}
// check the existance of the input files and out directory //

int MultivsMulti::check_input()
{
	int check_input=1; 


	string tmp;
	ifstream myfile1 (input1.c_str());
	int num1=0;
	if (myfile1.is_open())
	{
		while (getline (myfile1,tmp) )
		{
			num1=num1+1;
		}
		myfile1.close();
	}
	else cout << "Unable to open file:" <<input1<<endl; 

	ifstream myfile2 (input2.c_str());
	int num2=0;
	if (myfile2.is_open())
	{
		while ( getline (myfile2,tmp) )
		{
			num2=num2+1;
	
		}
		myfile2.close();
	}
	else cout << "Unable to open file:" <<input2<<endl; 


	if((num1<=1) | (num2<=1) )
	{
		cout<<"ERROR:You are trying to compare multiple sample vs multiple sample, but you give us more less than 2 samples for at least one experimental condition" <<endl;
		cout<<"No. for the first experiment:"<<num1<<endl;
		cout<<"No. for the second experiment:"<<num2<<endl;
		check_input=0;
	}
	else
	{
		string inputtmp;
		ifstream myfile1 (input1.c_str());
		while(getline (myfile1,inputtmp))
		{
			//open? 
			ifstream myfile1tmp (inputtmp.c_str());
			if (myfile1tmp.is_open())
			{
				input1file.push_back(inputtmp);
			}
			else
			{
				cout<<"FILE IN TREATMENT 1 Can't be open : " <<inputtmp <<endl;
				check_input=0;
			}
			myfile1tmp.close();
		}
		myfile1.close();

		ifstream myfile2 (input2.c_str());
		while(getline (myfile2,inputtmp))
		{
			//open? 
			ifstream myfile2tmp (inputtmp.c_str());
			if (myfile2tmp.is_open())
			{
				input2file.push_back(inputtmp);
			}
			else
			{
				cout<<"FILE IN TREATMENT 2 Can't be open : " <<inputtmp <<endl;
				check_input=0;
			}
			myfile2tmp.close();
		}
		myfile2.close();
	}

// Adding covariates //
	if(covinc)
	{

		int num1=0;
		ifstream myfile1 (input1cov.c_str());
		if (myfile1.is_open())
		{
			string tmp;
			while (getline (myfile1,tmp) )
			{
				if(!tmp.empty())
				{
					input1covfile=tmp;
					num1=num1+1;
				}
			}
			myfile1.close();
		}
		else cout << "Unable to open file:" <<input1cov<<endl; 

		int num2=0;
		ifstream myfile2 (input2cov.c_str());
		if (myfile2.is_open())
		{
			string tmp;	
			while (getline (myfile2,tmp) )
			{
				if(!tmp.empty())
				{	
					input2covfile=tmp;
					num2=num2+1;
				}
			}
			myfile2.close();
		}
		else cout << "Unable to open file:" <<input2cov<<endl; 

		if((num1!=1) | (num2!=1) )
		{
			cout<<"ERROR:For each treatment, you have more than one covariate or no covariates file" <<endl;
			cout<<"No. for the first experiment:"<<num1<<endl;
			cout<<"No. for the second experiment:"<<num2<<endl;
			check_input=0;
		}
		else
		{
			ifstream myfile1 (input1cov.c_str());
			getline (myfile1,input1covfile);
			myfile1.close();
			if(input1covfile.empty()) 
			{
				cout<<"ERROR: "<<input1cov<<" include empty line, please delete"<<endl;
				check_input=0;
			}			
			ifstream myfile2 (input2cov.c_str());
			getline (myfile2,input2covfile);
			myfile2.close();
			if(input2covfile.empty()) 
			{
				cout<<"ERROR: "<<input2cov<<" include empty line, please delete"<<endl;
				check_input=0;
			}
			//cout<<input1covfile<<input2covfile<<input1covfile.c_str()<<endl;
			ifstream myfile1tmp (input1covfile.c_str());
			if(!myfile1tmp.is_open()) 
			{
				cout << "Unable to open file:" <<input1covfile<<endl; 
				check_input=0;
			}
			myfile1tmp.close();
			
			ifstream myfile2tmp  (input2covfile.c_str());
			if(!myfile2tmp.is_open()) 
			{
				cout << "Unable to open file:" <<input2covfile<<endl; 
				check_input=0;
			}
			myfile2tmp.close();
		}

	}

	if ( access( outdir.c_str(), 0 ) == 0 )
	{
		struct stat status;
		stat( outdir.c_str(), &status );

		if (!( status.st_mode & S_IFDIR ))
		{
			cout << "ERROR:The outdir:"<<outdir<<" you entered is a file." << endl;
			check_input=0;
		}
	}
	else
	{
		cout << "ERROR:The outdir:"<<outdir<<" doesn't exist." << endl;
		check_input=0;
	}
	return check_input;

}
// Checking cutoff value, tolerance and length are realistic //
int MultivsMulti::check_cutoff()
{
	int check_cut=1;
	if(cutoff==0) {cout<<"ERROR: cutoff is exact zero"<<endl; check_cut=0;}
	if(cutoff<0) cutoff=abs(cutoff);
	if(tolerance > length) {cout<<"ERROR: tolerance is larger than considered length"<<endl; check_cut=0;}
	if(tolerance < 0) {cout<<"ERROR: tolerance is smaller than zero"<<endl; check_cut=0;}
	if(length < 0) {cout<<"ERROR: length is smaller than zero"<<endl; check_cut=0;}
	return check_cut;
}

// Print parameters //
void MultivsMulti::print()
{
	
	cout<<"Treatment 1 include :"<< "\n";
	for(size_t i=0; i<input1file.size(); i++) cout<<input1file[i]<<endl;
	cout<<"Treatment 2 include :"<< "\n";
	for(size_t i=0; i<input2file.size(); i++) cout<<input2file[i]<<endl;
	if(covinc)
	{
		cout<<"Cov 1 file : "<<input1covfile<<"\n";
		cout<<"Cov 2 file : "<<input2covfile<<"\n";
	}
	cout<<"outdir="<<outdir<<"\nlength="<<length<<"\ntolerance="<<tolerance<<"\ncutoff="<<cutoff<<endl;

}

//load data //
void MultivsMulti::load_data()
{ 
	// Vector to store all the files //
	vector <vector <string > > ALLDATA;
	size_t sample=input1file.size()+input2file.size();
	string tmp;
	for (size_t i = 0; i < input1file.size(); i++) 
	{
		vector <string> eachfile;
		std::ifstream f(input1file[i].c_str(), std::ios::in);
		while(getline(f,tmp))
		{
			eachfile.push_back(tmp);
		}
		ALLDATA.push_back(eachfile);
	}

	for (size_t i = 0; i < input2file.size(); i++) 
	{
		vector <string> eachfile;
		std::ifstream f(input2file[i].c_str(), std::ios::in);
		while(getline(f,tmp))
		{
			eachfile.push_back(tmp);
		}
		ALLDATA.push_back(eachfile);
	}

	vector <int> NoCpG;
	for(size_t i=0; i<ALLDATA.size(); i++) NoCpG.push_back(ALLDATA[i].size());
	vector <int> NoCurrent;
	for(size_t i=0; i<ALLDATA.size(); i++) NoCurrent.push_back(0);


	vector <bool> run;
	for (size_t i = 0; i < sample; i++) 	run.push_back(true);

	bool runany=true;
	while(runany)
	{
		for(size_t i=0; i<ALLDATA.size();i++)
		{
			if(run[i])
			{
 				string line1=ALLDATA[i][NoCurrent[i]];
				NoCurrent[i]=NoCurrent[i]+1;
				if(NoCurrent[i]>=NoCpG[i]) run[i]=0;

				std::vector<std::string> parts1;
				parts1=SplitString(line1);
				bool head1=parts1.at(0).find('#') != std::string::npos;
				if(!head1)
				{
					if(i==0) chromosome=parts1.at(0);
					double met=atoi(parts1.at(2).c_str());
					double tot=atoi(parts1.at(3).c_str());
					int pos=atoi(parts1.at(1).c_str());
					if((met<=tot) & (tot>=0) & (met >=0))
					{
						//////// if included in the map, then modify, otherwise insert into the map; 
						map<int, vector<int> >::iterator it;
						it=Tot.find(pos);
						if(it!=Tot.end())
						{
							vector <int> tot_int=it->second;
							tot_int[i]=tot;
							it->second=tot_int;
						}
						else
						{
 							vector <int> tot_int;
							for(size_t j=0; j<sample;j++)
							{
								if(j==i)tot_int.push_back(tot);
								else tot_int.push_back(0);
							}
       
							Tot.insert(pair<int, vector<int> >(pos,tot_int));
							//tot_int.clear();
						}
								
						it=Met.find(pos);
						if(it!=Met.end())
						{
							vector <int> met_int=it->second;
							met_int[i]=met;
							it->second=met_int;
						}
						else
						{
							vector <int> met_int;
							for(size_t j=0; j<sample;j++)
							{
								if(j==i)met_int.push_back(met);
								else met_int.push_back(0);
							}
       							Met.insert(pair<int, vector<int> >(pos,met_int));
							met_int.clear();
						} 
					}
				}
			}
		}
		bool runtmp=run[0];

		for(size_t j=1; j<run.size();j++)
		{
			runtmp = runtmp | run[j] ;
		}
		runany=runtmp;
	}

	if(!covinc)
	{
		vector<double> row;
		row.push_back(1);
		for (size_t i = 0; i < sample;i++)
		{
			COV.push_back(row);
		}

	}
	else
	{

		vector <string> covstr;
		vector <string> covstr1;
		vector <string> covstr2;
		string tmp;
		std::ifstream f1(input1covfile.c_str());
		arma::mat Mat;
		arma::mat Mat_tr;


		while(getline(f1,tmp))
		{
			covstr.push_back(tmp);
			covstr1.push_back(tmp);
		}

		std::ifstream f2(input2covfile.c_str());
		while(getline(f2,tmp))
		{
			covstr.push_back(tmp);
			covstr2.push_back(tmp);
		}

		if(covstr1.size()!=input1file.size())
		{
			cout<<"Under treatment condition 1, the number of individual files in covariate file does not match the number of individual in the CpG file"<<endl;
			COVERROR_=1;
		}

		if(covstr2.size()!=input2file.size())
		{
			cout<<"Under treatment condition 2, the number of individual files in covariate file does not match the number of individual in the CpG file"<<endl;
			COVERROR_=1;
		}

		if(!COVERROR_)
		{
			size_t Number_Cov=0;
			for(size_t i=0; i<covstr.size();i++)
			{
				string line1=covstr[i];
				std::vector<std::string> parts1;
				parts1=SplitString(line1);
				bool head1=parts1.at(0).find('#') != std::string::npos;
				if((!head1) & (!COVERROR_))
				{

					vector<double> row;
					row.push_back(1);
					if(Number_Cov==0) 
					{
						Number_Cov=parts1.size();
						Mat.set_size(covstr.size(), Number_Cov+1);
						Mat_tr.set_size(covstr.size(), Number_Cov+2);
  					}
					if(Number_Cov!=parts1.size()) 
					{
						cout<<"Number of Covariates does not match for each individual, please check"<<endl;
						COVERROR_=1;
					}
					else
					{
						Mat(i,0)=1;
						Mat_tr(i,0)=1;
						if(i<covstr1.size()) Mat_tr(i, 1)=1;
						else Mat_tr(i, 1)= 0;
						for(size_t j=0; j<parts1.size();j++)
						{
							row.push_back( stod(parts1.at(j)));
							Mat(i,j+1)=stod(parts1.at(j));
							Mat_tr(i,j+2)=stod(parts1.at(j));
						}
						COV.push_back(row);
					}
				}
			}


			double tmp1=abs(arma::det(Mat.t() * Mat));
			double tmp2=abs(arma::det(Mat_tr.t() * Mat_tr));

		
			bool a2=(tmp2<1.0e-7) & (tmp2>(-1.0e-7)) ;
			bool a1=(tmp1<1.0e-7) & (tmp1>(-1.0e-7));
			if(a1) 
			{
				cout<<"Mutli colinearity problem, please make sure you do not include intercept in your design matrix, intercept will be automatically added"<<endl;
				COVERROR_=1;
			}
			if(a2 & (!a1)) 
			{
				cout<<"Mutli colinearity problem, please make sure you do not include treatment effect in your design matrix, treatment effect will automatically added"<<endl;
				COVERROR_=1;
			}

		}

	}

	
};

bool MultivsMulti::Bool_Cov_Error()
{
	return COVERROR_;

}
void MultivsMulti::printdata()
{
	cout<<" TOT READS FOR EACH SAMPLE "<<endl;
	for(map<int, vector<int>  >::const_iterator it = Met.begin(); it != Met.end(); ++it)
	{
		cout << it->first << "\t";
		vector<int> tmp= it->second;
		for(size_t j=0; j<tmp.size(); j++)
		cout<<tmp[j]<<"\t";
		cout<<"\n";
	}
	cout<<" MET READS FOR EACH SAMPLE "<<endl;
	for(map<int, vector<int>  >::const_iterator it = Tot.begin(); it != Tot.end(); ++it)
	{
		cout << it->first << "\t";
		vector<int> tmp= it->second;
		for(size_t j=0; j<tmp.size(); j++)
		cout<<tmp[j]<<"\t";
		cout<<"\n";
	}
	

	cout<<" COVARIATES FOR EACH SAMPLE "<<endl;
	for(size_t j=0; j<COV.size(); j++)
	{
		for(size_t i=0; i<COV[j].size(); i++)
			cout<<COV[j][i]<<"\t";
		cout<<"\n";
	}

	
}


void MultivsMulti::compare()
{

	vector <int> pos_save;
	vector <double> z_save;

	vector < vector<double> > cov1int;
	vector<double> tmp;
//+input2file.size()
	for (size_t i=0;i<input1file.size();i++)
	{
		tmp=COV[i];
		cov1int.push_back(tmp);
	}
	
	vector < vector<double> > cov2int;
	for (size_t i=input1file.size();i<input1file.size()+input2file.size();i++)
	{
		tmp=COV[i];
		cov2int.push_back(tmp);
	}
	
	string outfile_exc=outdir+"/exclude.txt";
	ofstream myfile1;
	
	myfile1.open (outfile_exc.c_str());
	myfile1 <<"pos\t";
	string tmp1;
	string tmp2;
	for(size_t i=0; i<input1file.size();i++)
	{
		int tmp3=i+1;
		stringstream ss;
		ss << tmp3;
		string str = ss.str();
		tmp1="met1_"+str;
		tmp2="tot1_"+str;
		myfile1<<tmp1<<"\t"<<tmp2<<"\t";
	}
	for(size_t i=0; i<input2file.size();i++)
	{
		int tmp3=i+1;
		stringstream ss;
		ss << tmp3;
		string str = ss.str();
		tmp1="met2_"+str;
		tmp2="tot2_"+str;
		myfile1<<tmp1<<"\t"<<tmp2<<"\t";
	}
	myfile1<<"\n";
	for(map<int, vector<int>  >::const_iterator it = Tot.begin(); it != Tot.end(); ++it)
	{
		vector <int> tot=it->second;
		map<int, vector<int>  >::iterator itt2; 
		itt2=Met.find(it->first);
		if(itt2 == Met.end()) {cout<<"MET and TOT position does not match"<<endl;break;}
		else
		{
			vector <int> met=itt2->second;
		}
		vector <int>  met1int; 
		vector <int>  met2int; 
		vector <int>  tot1int; 
		vector <int>  tot2int; 

		for(size_t i=0; i<input1file.size(); i++)
		{
			met1int.push_back(itt2->second.at(i));
			tot1int.push_back(it->second.at(i));
		}
		for(size_t i=input1file.size(); i<input1file.size()+input2file.size(); i++)
		{
			met2int.push_back(itt2->second.at(i));
			tot2int.push_back(it->second.at(i));
		}

		Regression test;
		test.set_response(tot1int, met1int, tot2int, met2int,cov1int, cov2int);
		bool low_coverage=test.low_coverage_check();
		bool extreme_coverage=test.extreme_coverage_check();

		if(low_coverage)
		{
				myfile1 <<it->first<<"\t";
			for(size_t i=0; i<input1file.size(); i++)
			{
				myfile1 << met1int[i]<<"\t"<<tot1int[i]<<"\t";
			}
			for(size_t i=0; i<input2file.size(); i++)
			{
				myfile1 << met2int[i]<<"\t"<<tot2int[i]<<"\t";
			}
			myfile1<<"\n";
		}

		if(extreme_coverage)
		{
			my_result_multi tmp;
			tmp.chr=chromosome;
			tmp.pos=it->first;
			tmp.met1=met1int;
			tmp.tot1=tot1int;
			tmp.met2=met2int;
			tmp.tot2=tot2int;
			tmp.z_score=0;
			pos_save.insert(pos_save.end(),tmp.pos);
			z_save.insert(z_save.end(),0);
			mapResult.insert(pair<int, my_result_multi>(tmp.pos,tmp));
		}
		
		if((!low_coverage) & (!extreme_coverage))
		{
			vector <int> s1; 
			vector <int> s2;
			for(size_t i=0; i<tot1int.size();i++) if(tot1int[i]) s1.push_back(i);
			for(size_t i=0; i<tot2int.size();i++) if(tot2int[i]) s2.push_back(i);
			if(s1.size()==0) {cout<<"ERROR: treatment1 doesnot have coverage: Code is wrong"<<endl;break;}
			if(s2.size()==0) {cout<<"ERROR: treatment2 doesnot have coverage: Code is wrong"<<endl;break;}
			if((s1.size()==1) & (s2.size()==1)) 
			{
				// pair-wise comparison //
				int tmpmet1=met1int[s1[0]];
				int tmptot1=tot1int[s1[0]];
				int tmpmet2=met2int[s2[0]];
				int tmptot2=tot2int[s2[0]];
				Compare2Prop compare(tmpmet1,tmptot1,tmpmet2,tmptot2);
				double z_score=compare.compare2();
				my_result_multi tmp;
				tmp.chr=chromosome;
				tmp.pos=it->first;
				tmp.met1=met1int;
				tmp.tot1=tot1int;
				tmp.met2=met2int;
				tmp.tot2=tot2int;
				tmp.z_score=z_score;
				pos_save.insert(pos_save.end(),tmp.pos);
				z_save.insert(z_save.end(),tmp.z_score);
				mapResult.insert(pair<int, my_result_multi>(tmp.pos,tmp));				
			}
		
			if(!((s1.size()==1) & (s2.size()==1)))
			{
				BetaBinomialReg(test);
				double l_null=test.maximum_likelihood();
				test.FullModel();
			 	BetaBinomialReg(test);
				double l_alter=test.maximum_likelihood();
				//cout<<"FULL"<<l_alter<<"NULL"<<l_null<<endl;
				// chisq lrt //
				double log_lik_stat = -2*(l_null - l_alter);
				if(log_lik_stat<0) log_lik_stat=0;
				int met1tmp=0; int tot1tmp=0;
				for(size_t i=0; i<tot1int.size();i++) 
				{
					met1tmp=met1tmp+met1int[i];
					tot1tmp=tot1tmp+tot1int[i];
				}
				int met2tmp=0; int tot2tmp=0;
				for(size_t i=0; i<tot2int.size();i++) 
				{
					met2tmp=met2tmp+met2int[i];
					tot2tmp=tot2tmp+tot2int[i];
				}
//				double p1=met1tmp/tot1tmp*1.0;
//				double p2=met2tmp/tot2tmp*1.0;
				my_result_multi tmp;
				tmp.chr=chromosome;
				tmp.pos=it->first;
				tmp.met1=met1int;
				tmp.tot1=tot1int;
				tmp.met2=met2int;
				tmp.tot2=tot2int;
			//	int sign=1; if((p1-p2)<0) sign=-1;
				int sign=test.sign_treated_print();
				tmp.z_score=sqrt(log_lik_stat) * sign *1.0;
				pos_save.insert(pos_save.end(),tmp.pos);
				z_save.insert(z_save.end(),tmp.z_score);
				mapResult.insert(pair<int, my_result_multi>(tmp.pos,tmp));
			}

		}
	}
	myfile1.close();



	cout<<"Calculating Getis-ord....."<<endl;
	// Getis-ord calculation //
//	RInside R;
	Getis_Ord cal(pos_save,z_save,OutCovFile,CLINK_CPPFLAG,MydistSoftware,NewDistCppSoftware, RFun);
	cal.loadRFun();
	gi_final=cal.getis_ord();
	cal.clearR();

	vector<double> cut_range;
	cut_range.push_back(cutoff);

	// FindDMRs //
	cout<<"Finding DMRs....."<<endl;
	FindDMR finder(pos_save,gi_final,FindDMRfile,outdir,length,tolerance,chromosome,cut_range, RFun);
	finder.loadRFun();
	finder.FindDMR_implement();
	//finder.clearR();
}



void MultivsMulti::printFinalResult()
{
	cout<<"Output results for each CpG site....."<<endl;
	cout << "chr\t" <<"pos\t";
	string tmp1;
	string tmp2;
	for(size_t i=0; i<input1file.size();i++)
	{
		int tmp3=i+1;
		stringstream ss;
		ss << tmp3;
		string str = ss.str();
		tmp1="met1_"+str;
		tmp2="tot1_"+str;
		cout<<tmp1<<"\t"<<tmp2<<"\t";
	}
	for(size_t i=0; i<input2file.size();i++)
	{
		int tmp3=i+1;
		stringstream ss;
		ss << tmp3;
		string str = ss.str();
		tmp1="met2_"+str;
		tmp2="tot2_"+str;
		cout<<tmp1<<"\t"<<tmp2<<"\t";
	}
	cout<<"z_score\t"<<"g_score\t"<<"\n";
	int num=0;
	for(map<int, my_result_multi >::const_iterator it = mapResult.begin(); it != mapResult.end(); ++it)
	{
		cout << it->second.chr<< "\t" << it->second.pos<< "\t" ;
		for(size_t i=0; i<input1file.size();i++)
		{
			cout<<it->second.met1[i]<<"\t"<<it->second.tot1[i]<<"\t";
		}
		for(size_t i=0; i<input2file.size();i++)
		{
			cout<<it->second.met2[i]<<"\t"<<it->second.tot2[i]<<"\t";
		}


		cout<<it->second.z_score<<"\t"<<gi_final.at(num)<<"\n"; //
		num=num+1;
	}

};
void MultivsMulti::printFinalResult(string outdir)
{
	cout<<"Output results for each CpG site....."<<endl;

	string tmp=outdir+"/outputfinal.txt";
	ofstream myfile1;
	myfile1.open (tmp.c_str());
	myfile1 << "chr\t" <<"pos\t";
	string tmp1;
	string tmp2;
	for(size_t i=0; i<input1file.size();i++)
	{
		int tmp3=i+1;
		stringstream ss;
		ss << tmp3;
		string str = ss.str();
		tmp1="met1_"+str;
		tmp2="tot1_"+str;		
		myfile1<<tmp1<<"\t"<<tmp2<<"\t";
	}
	for(size_t i=0; i<input2file.size();i++)
	{
		int tmp3=i+1;
		stringstream ss;
		ss << tmp3;
		string str = ss.str();
		tmp1="met2_"+str;
		tmp2="tot2_"+str;
		myfile1<<tmp1<<"\t"<<tmp2<<"\t";
	}
	myfile1<<"z_score\t"<<"g_score\t"<<"\n";
	int num=0;
	for(map<int, my_result_multi >::const_iterator it = mapResult.begin(); it != mapResult.end(); ++it)
	{
		myfile1 << it->second.chr<< "\t" << it->second.pos<< "\t" ;
		for(size_t i=0; i<input1file.size();i++)
		{
			myfile1<<it->second.met1[i]<<"\t"<<it->second.tot1[i]<<"\t";
		}
		for(size_t i=0; i<input2file.size();i++)
		{
			myfile1<<it->second.met2[i]<<"\t"<<it->second.tot2[i]<<"\t";
		}


		myfile1<<it->second.z_score<<"\t"<<gi_final.at(num)<<"\n"; //
		num=num+1;
	}
	myfile1.close();
};
