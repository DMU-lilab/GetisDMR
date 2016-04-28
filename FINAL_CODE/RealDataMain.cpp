#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <iostream>
#include <string>
#include <fstream>
#include "one_vs_one.h"
#include "SplitString.h"
#include "FindDMR.h"
#include <Rcpp.h>
#include <RInside.h>
#include "FindFileHome.h"
#include "multi_vs_multi.h"
#include "regression.h"
#include "BetaBinomialReg.h"
#include "GetisDMR_p_value.h"
#include <unistd.h>
using namespace std;



string OutCovFile = "";
string CLINK_CPPFLAG="-I /usr/lib/R/site-library/RcppArmadillo/include/";
string MydistSoftware="MydistSoftware.cpp";
string NewDistCppSoftware="NewDistCppSoftware.r";
string FindDMRfile="FindDMRSoftware.r";
int MINTOT=5;
RInside RFun;


//g++ GetisDMR_p_value.cpp multi_vs_multi.cpp BetaBinomialReg.cpp regression.cpp FindFileHome.cpp FindDMR.cpp SplitString.cpp compare2Prop.cpp Getis_Ord.cpp RealDataMain.cpp one_vs_one.cpp -I/usr/share/R/include -I/home/yalu/R/i686-pc-linux-gnu-library/3.0/Rcpp/include -I/home/yalu/R/i686-pc-linux-gnu-library/3.0/RInside/include -O3 -pipe -g -Wall  /usr/lib/i386-linux-gnu/libboost_program_options.so  -lboost_system  -lboost_filesystem -L/usr/lib/R/lib -lR  -lblas -llapack  -larmadillo -L/home/yalu/R/i686-pc-linux-gnu-library/3.0/RInside/lib -lRInside -Wl,-rpath,/home/yalu/R/i686-pc-linux-gnu-library/3.0/RInside/lib -lgsl -lgslcblas -lm -std=c++11  -o GetisDMR


namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace 
int main(int argc, char** argv) 
{ 

	char buffer[1000];
	char *answer = getcwd(buffer, sizeof(buffer));

	if (answer)
	{
 		OutCovFile = answer;
	}

	boost::filesystem::path full_path( boost::filesystem::current_path() );
	OutCovFile=full_path.string();
	//std::cout << "Current path is : " << OutCovFile << std::endl;

	int comparision_opt=0;
	string input1;
	string input2;
	string outdir;
	int tolerance=2;
	int length=6;
	double cutoff=1.645;
	bool sorted=1;
	string covinput1="NOINPUT";
	string covinput2="NOINPUT";
	bool asymptotic_run=0;
//	if(OutCovFile=="") cout<<"NO OUTPUT FILE SPECIFIED FOR COVARIANCE"<<endl;
	outdir=OutCovFile;
	try 
	{ 

		/** Define and parse the program options **/ 
		namespace po = boost::program_options; 
		po::options_description desc("Options"); 
		desc.add_options() 
			("help,h", "Print help messages") 
			("Comparison,c", po::value<int>(&comparision_opt)->required(),"1: 1 vs 1 (no biological replicates) \n2: with biological replicates but no covariates \n3: biological replicates and covariates") 
			("input1", po::value< string >(&input1)->required(),"The file including the locations of files under treatment condition 1")
			("input2", po::value< string >(&input2)->required(),"The file including the locations of files under treatment condition 2")
			("tolerance,t", po::value<int>(&tolerance), "tolerance")
			("length,l", po::value<int>(&length), "length")
			("cutoff", po::value<double>(&cutoff), "cutoff")
			("outfiledir,o",po::value< string >(&outdir), "The output folder: Default is the current folder")
			("sorted,s", po::value<bool>(&sorted), "The mtbr is sorted according to physical location")
		// ENVIROMENT FOR R //
			("CLINK_CPPFLAG",po::value< string >(&CLINK_CPPFLAG), "RcppArmadillo link")
			("MydistSoftware",po::value< string >(&MydistSoftware), "MydistSoftware directory")
			("NewDistCppSoftware",po::value< string >(&NewDistCppSoftware), "NewDistCppSoftware directory")
			("FindDMRfile",po::value< string >(&FindDMRfile), "FindDMRSoftware.r directory")
			("sorted",po::value< bool >(&sorted), "If the mtbr is sorted according to physical position, currently not support un-sorted")
			("MINTOT",po::value< int >(&MINTOT), "Minimum of the total number of reads required for each sample, default 5")
			("cov1",po::value< string >(&covinput1), "The covariates file of treatment condition 1")
			("cov2",po::value< string >(&covinput2), "The covariates file of treatment condition 2")
 			("P_value,p",po::value< bool >(&asymptotic_run), "Calculate p-value or not");

		po::variables_map vm; 
		try 
		{ 
			po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw

			/** --help option  **/ 
			if (vm.count("help")  ) 
			{ 
			std::cout << "Basic Command Line Parameters" << std::endl << desc << std::endl; 
			return SUCCESS; 
			} 
 
			po::notify(vm); // throws on error, so do after help in case 
		} 

		catch(po::error& e) 
		{ 
			std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
			std::cerr << desc << std::endl; 
			return ERROR_IN_COMMAND_LINE; 
		} 
 
		// application code here // 
		if(vm.count("tolerance"))
		tolerance=vm["tolerance"].as<int>();
		if(vm.count("length"))
		tolerance=vm["length"].as<int>();
		if(vm.count("cutoff"))
		tolerance=vm["cutoff"].as<double>();
		if(vm.count("outfiledir"))
		{
			outdir=vm["outfiledir"].as<string>();
			OutCovFile=outdir;
		}

		if(vm.count("Comparison"))
		comparision_opt=vm["Comparison"].as<int>();


		if(!vm.count("MydistSoftware"))
		{
			cout<<"Will search for MydistSoftware.cpp"<<endl;
			MydistSoftware=FindFileHome("/MydistSoftware.cpp");
			cout<<"This is what I found from your home direcotry: "<<MydistSoftware<<endl;
			if(MydistSoftware=="NOTFIND") return ERROR_IN_COMMAND_LINE; 
		}
		
		if(!vm.count("NewDistCppSoftware"))
		{
			cout<<"Will search for NewDistCppSoftware.r"<<endl;
			NewDistCppSoftware=FindFileHome("/NewDistCppSoftware.r");
			cout<<"This is what I found from your home direcotry: "<<NewDistCppSoftware<<endl;
			if(NewDistCppSoftware=="NOTFIND") return ERROR_IN_COMMAND_LINE; 
		}


	
		if(!vm.count("FindDMRfile"))
		{
			cout<<"Will search for FindDMRSoftware.r"<<endl;
			FindDMRfile=FindFileHome("/FindDMRSoftware.r");
			cout<<"This is what I found from your home direcotry: "<<FindDMRfile<<endl;
			if(FindDMRfile=="NOTFIND") return ERROR_IN_COMMAND_LINE; 
		}

		std::cout << "Current path is : " << OutCovFile << std::endl;
                if(OutCovFile=="") cout<<"NO OUTPUT FILE SPECIFIED FOR COVARIANCE"<<endl;



		if(!((comparision_opt==1) | (comparision_opt ==2) |(comparision_opt==3))) 
		{
			cout<<"comparison option is out of range, please check help "<<endl;
			return ERROR_IN_COMMAND_LINE; 
		}
		else
		{
			/** One to one comparision **/
			if(comparision_opt==1) 
			{
				cout<<"We are comparing: 1 vs 1 (no biological replicates)"<<endl;
				OnevsOne runonevsone(input1,input2, outdir, tolerance, length, cutoff,sorted);
				int check_input= runonevsone.check_input();
				int check_cutoff=runonevsone.check_cutoff();
				if((check_input==1) & (check_cutoff==1))
				{
					runonevsone.print();
					runonevsone.printinput();
					cout<<"FINISHING CHECKING ....."<<endl;
					cout<<"Start Analysis ......" <<endl;
					runonevsone.compare();
				
					if(asymptotic_run)
					{
						cout<<"Calculating p-value ......"<<endl;
						string outputfinal_asym=outdir+"/outputfinal.txt";
						std::string cutoff_str= std::to_string (cutoff);
						cutoff_str.erase ( cutoff_str.find_last_not_of('0') + 1, std::string::npos );
						string dmr1_asym=outdir+"/G1D"+cutoff_str.substr(0,4);
						string dmr2_asym=outdir+"/G1D-"+cutoff_str.substr(0,4);
						string cov_asym=OutCovFile+"/Covarince.txt";
						std::ifstream infile1_asym(dmr1_asym);
						std::ifstream infile2_asym(dmr2_asym);
						GetisDMR_P_Value asympotics1(dmr1_asym,cov_asym,outputfinal_asym);
						GetisDMR_P_Value asympotics2(dmr2_asym,cov_asym,outputfinal_asym);
						if(infile1_asym.good())
						{
							cout<<"Calculating p-value for "<<dmr1_asym<<endl;
							asympotics1.load_data();
							asympotics1.GetisAsympotic();
							asympotics1.output(dmr1_asym);
						}
						if(infile2_asym.good())
						{
							cout<<"Calculating p-value for "<<dmr2_asym<<endl;
							asympotics2.load_data();
							asympotics2.GetisAsympotic();
							asympotics2.output(dmr2_asym);
						}
					}

				}
			
			}
			else if(comparision_opt==2)
			{
				cout<<"We are comparing: Multi vs Multi (without covariates)"<<endl;
				bool covinc=0;

				MultivsMulti runmultivsmulti(input1,input2,outdir,tolerance, length, cutoff,sorted,covinc);
			
				int check_input= runmultivsmulti.check_input();
				int check_cutoff=runmultivsmulti.check_cutoff();
				
				if((check_input==1) & (check_cutoff==1))
				{
					runmultivsmulti.print();
					cout<<"FINISHING CHECKING ....."<<endl;
					cout<<"LOADING DATA"<<endl;
					runmultivsmulti.load_data();
					//runmultivsmulti.printdata();
					cout<<"Start Analysis ......" <<endl;
					runmultivsmulti.compare();
					runmultivsmulti.printFinalResult(outdir);
					//runmultivsmulti.printFinalResult();

					if(asymptotic_run)
					{
						cout<<"Calculating p-value ......"<<endl;
						string outputfinal_asym=outdir+"/outputfinal.txt";
						std::string cutoff_str= std::to_string (cutoff);
						cutoff_str.erase ( cutoff_str.find_last_not_of('0') + 1, std::string::npos );
						string dmr1_asym=outdir+"/G1D"+cutoff_str.substr(0,4);
						string dmr2_asym=outdir+"/G1D-"+cutoff_str.substr(0,4);
						string cov_asym=OutCovFile+"/Covarince.txt";
						std::ifstream infile1_asym(dmr1_asym);
						std::ifstream infile2_asym(dmr2_asym);
						GetisDMR_P_Value asympotics1(dmr1_asym,cov_asym,outputfinal_asym);
						GetisDMR_P_Value asympotics2(dmr2_asym,cov_asym,outputfinal_asym);
						if(infile1_asym.good())
						{
							cout<<"Calculating p-value for "<<dmr1_asym<<endl;
							asympotics1.load_data();
							asympotics1.GetisAsympotic();
							asympotics1.output(dmr1_asym);
						}
						if(infile2_asym.good())
						{
							cout<<"Calculating p-value for "<<dmr2_asym<<endl;
							asympotics2.load_data();
							asympotics2.GetisAsympotic();
							asympotics2.output(dmr2_asym);
						}
					}

				}
			}
			else
			{
				cout<<"with biological replicates and covariates"<<endl;
				bool covinc=1;
				if((covinput1=="NOINPUT") | (covinput2=="NOINPUT"))
				{
					if(covinput1=="NOINPUT") cout<<"Please provide the covariates input file for treatment condition 1"<<endl;
					if(covinput2=="NOINPUT") cout<<"Please provide the covariates input file for treatment condition 2"<<endl;
					return ERROR_IN_COMMAND_LINE; 
				}

				MultivsMulti runmultivsmulti(input1,input2,outdir,tolerance, length, cutoff,sorted,covinc,covinput1,covinput2);	
				int check_input= runmultivsmulti.check_input();
				int check_cutoff=runmultivsmulti.check_cutoff();

				if((check_input==1) & (check_cutoff==1))
				{
					runmultivsmulti.print();
					cout<<"FINISHING CHECKING ....."<<endl;
					cout<<"LOADING DATA"<<endl;
					runmultivsmulti.load_data();
					if(runmultivsmulti.Bool_Cov_Error()) 
					{
						cout<<"Something is wrong with the covariate files, please check!"<<endl;
						return ERROR_IN_COMMAND_LINE; 
					}
					if(!runmultivsmulti.Bool_Cov_Error()) 
					{
						//runmultivsmulti.printdata();
						cout<<"Start Analysis ......" <<endl;
						runmultivsmulti.compare();
						runmultivsmulti.printFinalResult(outdir);

						if(asymptotic_run)
						{
							cout<<"Calculating p-value ......"<<endl;
							string outputfinal_asym=outdir+"/outputfinal.txt";
							std::string cutoff_str= std::to_string (cutoff);
							cutoff_str.erase ( cutoff_str.find_last_not_of('0') + 1, std::string::npos );
							string dmr1_asym=outdir+"/G1D"+cutoff_str.substr(0,4);
							string dmr2_asym=outdir+"/G1D-"+cutoff_str.substr(0,4);
							string cov_asym=OutCovFile+"/Covarince.txt";
							std::ifstream infile1_asym(dmr1_asym);
							std::ifstream infile2_asym(dmr2_asym);
							GetisDMR_P_Value asympotics1(dmr1_asym,cov_asym,outputfinal_asym);
							GetisDMR_P_Value asympotics2(dmr2_asym,cov_asym,outputfinal_asym);
							if(infile1_asym.good())
							{
								cout<<"Calculating p-value for "<<dmr1_asym<<endl;
								asympotics1.load_data();
								asympotics1.GetisAsympotic();
								asympotics1.output(dmr1_asym);
							}
							if(infile2_asym.good())
							{
								cout<<"Calculating p-value for "<<dmr2_asym<<endl;
								asympotics2.load_data();
								asympotics2.GetisAsympotic();
								asympotics2.output(dmr2_asym);
							}
						}

					//runmultivsmulti.printFinalResult();*/;
					}
				}
			
			}
			
		}
	} 
	catch(std::exception& e) 
	{ 
		std::cerr << "Unhandled Exception reached the top of main: " << e.what() << ", application will now exit" << std::endl; 
		return ERROR_UNHANDLED_EXCEPTION; 
 
	} 


 	cout<<"Analysis is done!"<<endl;
	cout<<"Have a nice day :)"<<endl;

	
  return SUCCESS; 
 
} // main 
