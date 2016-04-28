#include "Getis_Ord.h"

extern string OutCovFile;
extern string CLINK_CPPFLAG;
extern string MydistSoftware;
extern string NewDistCppSoftware;


Getis_Ord::Getis_Ord(vector<int> pos_int,vector<double>z_score_int, string OutCovFile_int,string CLINK_CPPFLAGS_int,string MydistSoftware_int, string NewDistCppSoftware_int,RInside & mR ):R(mR)
{
	pos=pos_int;
	z_score=z_score_int;
	OutCovFile= OutCovFile_int;
	CLINK_CPPFLAG=CLINK_CPPFLAGS_int;
	MydistSoftware=MydistSoftware_int;
	NewDistCppSoftware=NewDistCppSoftware_int;
}

void Getis_Ord::loadRFun()
{
	R["ENVRCPPARMADILLO"] = CLINK_CPPFLAG;
	R["MydistSoftware"] = MydistSoftware;
	R["NewDistCppSoftware"] = NewDistCppSoftware;
	R.parseEvalQ("source(NewDistCppSoftware)");  
}

vector<double> Getis_Ord::getis_ord()
{

		SEXP gi;

	//	std::vector<int> re(pos, &pos[pos.size()]);
	//	R.assign(re, "re");	// assign STL vector to R's 're' variable
		
	//	std::vector<double> z(z_score, &z_score[z_score.size()]);
	//	R.assign(z, "z");	// assign STL vector to R's 'weightsvec' variable
	
		SEXP zz=Rcpp::wrap(z_score);
		R.assign(zz, "z");
		SEXP re=Rcpp::wrap(pos);
		R.assign(re,"re");

		R["OutCovFile"]=OutCovFile+"/Covarince.txt";
		string txt="GetisNew(re,z,OutCovFile)";

		if (R.parseEval(txt, gi))      
				throw std::runtime_error("R cannot evaluate '" + txt + "'");
 		vector<double> vector1=Rcpp::as<vector<double> >(gi);
		return vector1;
}


void Getis_Ord::clearR()
{
	R.parseEvalQ("rm(list=ls())");
}
