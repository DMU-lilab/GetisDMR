#ifndef GETISORD_H
#define GETISORD_H
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
class Getis_Ord
{
	private:
		vector<int> pos;
		vector<double>z_score;
		vector<double>g_score;
		string OutCovFile;
		string CLINK_CPPFLAG;
		string MydistSoftware;
		string NewDistCppSoftware;
		RInside & R;              
	public:
		Getis_Ord(vector<int> pos,vector<double>z_score, string OutCovFile,string CLINK_CPPFLAGS_int,string MydistSoftware_int, string NewDistCppSoftware_int,RInside & mR);
		void loadRFun();
		vector<double> getis_ord();
		void clearR();
		//vector <double> getis_ord();

};

 
#endif
