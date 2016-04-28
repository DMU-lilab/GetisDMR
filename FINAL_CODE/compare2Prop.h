#ifndef COMPARE2PROP_H
#define COMPARE2PROP_H
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <algorithm> 
#include <boost/math/distributions/normal.hpp>

using namespace std;
class Compare2Prop
{
	private:
		int met1;
		int met2;
		int tot1;
		int tot2;
	public:
		Compare2Prop(int met1_int, int tot1_int, int met2_int, int tot2_int);
		void SetPara(int met1_int, int tot1_int, int met2_int, int tot2_int);
		double compare2();
		double chisq();
		double fisher();


};

 
#endif
