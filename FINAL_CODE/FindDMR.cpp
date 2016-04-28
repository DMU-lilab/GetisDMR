#include "FindDMR.h"

extern string OutCovFile;
extern string CLINK_CPPFLAG;
extern string MydistSoftware;
extern string NewDistCppSoftware;

// This function is mearly a wrap to call the R function FindDMRSoftware.r //
FindDMR::FindDMR(vector<int> pos_int,vector<double>g_score_int, string FindDMRfile_int,string OutCovFile_int,int length_int,int tolerance_int,string chrom_int, vector<double> cutrange_int,RInside & mR ):R(mR)
{
	pos=pos_int;
	g_score=g_score_int;
	OutCovFile= OutCovFile_int;
	FindDMRfile=FindDMRfile_int;	
	length=length_int;
	tolerance=tolerance_int;
	chrom=chrom_int;
	cutrange=cutrange_int;
}

void FindDMR::loadRFun()
{
	R["dmr_fun"] = FindDMRfile;
	R.parseEvalQ("source(dmr_fun)");  
}

void FindDMR::FindDMR_implement()
{

		SEXP gg=Rcpp::wrap(g_score);
		R.assign(gg, "Gi");
		string txt="G.i=list();G.i[[1]]=Gi;";
		R.parseEvalQ(txt);

		SEXP posi=Rcpp::wrap(pos);
		R.assign(posi,"posi");
	
		SEXP cut=Rcpp::wrap(cutrange);
		R.assign(cut,"range.di.high");

		R["cf.length"]=length;
		R["tolerance.length"]=	tolerance;
		R["chrom"]=chrom;

/*	
		SEXP len=length;
		R.assign(len,"cf.length");
		
		SEXP tol=tolerance;
		R.assign(tol,"tolerance.length");
		
		SEXP chrom=chr;
		R.assign(chrom,"chrom")
*/		



		
		if(strcmp(&OutCovFile[OutCovFile.length()-1], "/")) OutCovFile=OutCovFile+"/";

		R["outfile"]=OutCovFile;

//		txt="for (i in 1:length(range.di.high)) print(range.di.high[i]);";
//		R.parseEvalQ(txt);


		txt="findDMR(G.i,posi,range.di.high,range.di.low=NULL,cf.length=cf.length, tolerance.length=tolerance.length,outfile,chrom=chrom)";
		R.parseEvalQ(txt);
}



void FindDMR::clearR()
{
	R.parseEvalQ("rm(list=ls())");
}
