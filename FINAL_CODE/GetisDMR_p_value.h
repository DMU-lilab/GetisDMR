#ifndef GETISDMR_P_VALUE_H
#define GETISDMR_P_VALUE_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 

using namespace std;

class GetisDMR_P_Value
{
	private: 
		string DMRinputfile_;
		string COVinputfile_;
		string ALLinputfile_;
		vector < vector<string> > DMRINPUT_;
		vector < vector<double> > COVINPUT_;
		vector < vector<double> > ALLINPUT_;
		bool check_input=1;
	public: 
		GetisDMR_P_Value(string DMRinputfile, string COVinputfile, string ALLinputfile);
		void load_data();
		void print_load_data(bool PrintDMR=1, bool PrintCov=1, bool PrintAll=1);
		void GetisAsympotic();
		void output();
		void output(string outdmr);
};

#endif 
