#ifndef MULTIVSMULTI_H
#define MULTIVSMULTI_H

#include <vector>
#include <cmath>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

struct my_result_multi
{
	string chr;
	int pos;
	vector <int> met1;
	vector <int> tot1;
	vector <int> met2;
	vector <int> tot2;
	double z_score;
//	double getis_score;
};

class MultivsMulti
{
	private:
		string input1;
		string input2;
		vector <string> input1file;
		vector <string> input2file;
		map<int, vector<int> >Met;
		map<int, vector<int> >Tot;

		string input1cov;
		string input2cov;
		string input1covfile;
		string input2covfile;
		vector < vector<double> > COV;

		string chromosome;
		string outdir;
		bool sorted;
		int tolerance;
		int length;
		double cutoff;
		bool covinc;
		bool COVERROR_=0;

 		map<int, my_result_multi > mapResult; 
		vector <double> gi_final;
//		RInside R;

//map<int, vector<int> > Tot1; //pos, total reads for each sample in treatment1;
//map<int, vector<int> > Met1; //pos,   met reads for each sample in treatment1;
//map<int, vector<int> > Tot2; //pos, total reads for each sample in treatment2;
//map<int, vector<int> > Met2; //pos,   met reads for each sample in treatment2;

	public:

		MultivsMulti(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int,bool sorted_int, bool covinc_int);
		MultivsMulti(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int, bool sorted_int, bool covinc_int, string covinput1_int, string covinput2_int);
	    	int check_input();
		int check_cutoff();
		void print();
		void printdata();
		void compare();
		void load_data();
		void printFinalResult();
		void printFinalResult(string outdir);
		bool Bool_Cov_Error();
};
 
#endif

