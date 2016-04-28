#ifndef ONEVSONE_H
#define ONEVSONE_H
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <iostream>
#include <string>
#include <fstream>


using namespace std;
class OnevsOne
{
	private:
		string input1;
		string input2;
		string input1file;
		string input2file;
		string outdir;
		bool sorted;
		int tolerance;
		int length;
		double cutoff;
	public:
		OnevsOne(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int,bool sorted_int);
		void SetPara(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int,bool sorted_int);
	    	int check_input();
		int check_cutoff();
		void print();
		void printinput();
		void compare();
};
 
#endif

