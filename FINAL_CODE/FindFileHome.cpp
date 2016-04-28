
#include <iostream>
#include <string>
#include <fstream>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"
#include <iostream>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>


namespace fs = boost::filesystem;
using namespace std;
string FindFileHome(string exp)
{
	string out="NOTFIND";

	char *homedir;

	if ((homedir = getenv("HOME")) == NULL) {
	    homedir = getpwuid(getuid())->pw_dir;
	}

	fs::path path=homedir;	

//	boost::filesystem::path path = boost::filesystem::current_path();
 	boost::filesystem::recursive_directory_iterator itr(path);
	bool run=1;
	while ((itr != boost::filesystem::recursive_directory_iterator()) & run)
	{
		string tmp=itr->path().string();
		if((tmp.find(exp) != std::string::npos))
		{
			string tmp1=tmp.substr(tmp.length()-exp.length(),tmp.length());
			if(tmp1==exp)
			{
				out=tmp;
				run=0;
			}
		}

        	++itr;
	}
	if(out=="NOTFIND") cout<<" Can't find "<<exp<<"anywhere from your home directory, please specify the location of your file"<<endl;
	return out;
}
