#include <boost/algorithm/string.hpp>
#include <iostream>
#include <vector>
#include <string>
#include "SplitString.h"
using namespace std;
std::vector<std::string> SplitString(std::string line)
{
//std::string inputString("One!Two,Three sfds:Four   \tds");
	string delimiters("|,: \t /");
	vector<string> parts;
	boost::split(parts, line, boost::is_any_of(delimiters),boost::token_compress_on);
//	for(int i=0;i<parts.size();i++)
//	cout<<parts.at(i)<<":"<<endl;
	return parts;
}
