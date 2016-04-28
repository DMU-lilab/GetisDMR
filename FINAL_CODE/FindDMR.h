#ifndef FINDDMR_H
#define FINDDMR_H
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <algorithm> 
#include <boost/math/distributions/normal.hpp>
#include <Rcpp.h>
#include <RInside.h>

using namespace std;
class FindDMR
{
	private:
		vector<int> pos;
		vector<double>g_score;
		string FindDMRfile;
		string OutCovFile;
		int length;
		int tolerance;
		string chrom;
		vector<double> cutrange;
		RInside &R;              
	public:
		FindDMR(vector<int> pos_int,vector<double>g_score_int, string FindDMRfile_int,string OutCovFile_int,int length_int,int tolerance_int,string chrom_int, vector<double> cutrange_int,RInside & mR );
		void loadRFun();
		void FindDMR_implement();
		void clearR();
};

 
#endif
