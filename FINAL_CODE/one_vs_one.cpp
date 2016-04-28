#include "one_vs_one.h"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp" 
#include <iostream>
#include <string>
#include <fstream>
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
#include "SplitString.h"
#include "compare2Prop.h"
#include "Getis_Ord.h"
#include "FindDMR.h"
#include <Rcpp.h>
#include <RInside.h>

extern std::vector<std::string> SplitString(std::string line);
extern string OutCovFile;
extern string CLINK_CPPFLAG;
extern string MydistSoftware;
extern string NewDistCppSoftware;
extern string FindDMRfile;
extern int MINTOT;
extern RInside RFun;

using namespace std;
struct chr_met_tot
{
	double met;
	double tot;
	string chr;
};

struct my_result
{
	string chr;
	int pos;
	int met1;
	int tot1;
	int met2;
	int tot2;
	double z_score;
	double getis_score;
};


/* constructor */

OnevsOne::OnevsOne(string input1_int, string input2_int, string outdir_int, int tol_int, int len_int, double cut_int, bool sorted_int)
{
	SetPara(input1_int, input2_int, outdir_int, tol_int, len_int, cut_int,sorted_int);
}
void OnevsOne::SetPara(string input1_int, string input2_int,string outdir_int, int tol_int, int len_int, double cut_int, bool sorted_int)
{
	input1=input1_int;
	input2=input2_int;
	tolerance=tol_int;
	length=len_int;
	cutoff=cut_int;
	outdir=outdir_int;
	sorted=sorted_int;
}
/* check the existance of the input files and out directory */
int OnevsOne::check_input()
{
	int check_input=1;
	ifstream myfile1 (input1.c_str());
	int num1=0;
	if (myfile1.is_open())
	{
		while (getline (myfile1,input1file) )
		{
			num1=num1+1;
		}
		myfile1.close();
	}
	else cout << "Unable to open file:" <<input1<<endl; 

	ifstream myfile2 (input2.c_str());
	int num2=0;
	if (myfile2.is_open())
	{
		while ( getline (myfile2,input2file) )
		{
			num2=num2+1;
	
		}
		myfile2.close();
	}
	else cout << "Unable to open file:" <<input2<<endl; 


	if((num1!=1) | (num2!=1) )
	{
		cout<<"ERROR:You are trying to compare one sample vs another sample, but you give us more than 2 samples for at least one experimental condition" <<endl;
		cout<<"No. for the first experiment:"<<num1<<endl;
		cout<<"No. for the second experiment:"<<num2<<endl;
		check_input=0;
	}
	else
	{
		ifstream myfile1 (input1.c_str());
		getline (myfile1,input1file);
		myfile1.close();
		ifstream myfile2 (input2.c_str());
		getline (myfile2,input2file);
		myfile2.close();
	}


	if ( access( outdir.c_str(), 0 ) == 0 )
	{
		struct stat status;
		stat( outdir.c_str(), &status );

		if (!( status.st_mode & S_IFDIR ))
		{
			cout << "ERROR:The outdir:"<<outdir<<" you entered is a file." << endl;
			check_input=0;
		}
	}
	else
	{
		cout << "ERROR:The outdir:"<<outdir<<" doesn't exist." << endl;
		check_input=0;
	}
	return check_input;

}
/* Checking cutoff value, tolerance and length are realistic */
int OnevsOne::check_cutoff()
{
	int check_cut=1;
	if(cutoff==0) {cout<<"ERROR: cutoff is exact zero"<<endl; check_cut=0;}
	if(cutoff<0) cutoff=abs(cutoff);
	if(tolerance > length) {cout<<"ERROR: tolerance is larger than considered length"<<endl; check_cut=0;}
	if(tolerance < 0) {cout<<"ERROR: tolerance is smaller than zero"<<endl; check_cut=0;}
	if(length < 0) {cout<<"ERROR: length is smaller than zero"<<endl; check_cut=0;}
	return check_cut;
}

/* Print parameters */
void OnevsOne::print()
{
	cout<<"input 1="<<input1<<"\ninput 2="<<input2<<"\noutdir="<<outdir<<"\nlength="<<length<<"\ntolerance="<<tolerance<<"\ncutoff="<<cutoff<<endl;

}

/*Print the details for each input files to be compared */;
void OnevsOne::printinput()
{
	cout<<"input1file="<<input1file<<"\ninput2file="<<input2file<<endl;
}


/* Conduct one vs one comparison */
void OnevsOne::compare()
{

	int error=0;
	ifstream myfile1 (input1file.c_str());
	ifstream myfile2 (input2file.c_str());

	if ((!myfile1.is_open()) | (!myfile2.is_open()) )
	{
		if(!myfile1.is_open()) cout << "Unable to open file:" <<input1file<<endl; 
		if(!myfile2.is_open()) cout << "Unable to open file:" <<input2file<<endl; 
		error=1;
	}
	if(input1file==input2file) 
	{
		cout << "ERROR:The two input files are exactly the same. input1=" <<input1file<<",input2="<<input2file<<endl; 
		error=1;
	}

	if(error==0)
	{
		//Reading the data //; 
		if(!sorted) cout<< "ERROR: Current version can only handle the situation where the mtbr is sorted according to physical location, please sort"<<endl;
		else
		{

			bool run1=1;
			bool run2=1;
		
			map<int, my_result > mapResult; //chr,pos, met1,tot1,met2,tot2,z-score, getis-score
			map<int, chr_met_tot > record1; //pos, methylated, tot for file 1;
			map<int, chr_met_tot > record2; //pos, methylated, tot for file 1;
			
			vector<int> pos_save;
			vector<double> z_save;

			int num=0;
			bool head1=0; //0 no header, 1 header
			bool head2=0;
			while(run1 | run2)
			{
				num=num+1;
				string line1;
				getline(myfile1,line1);
				std::vector<std::string> parts1;
				if(myfile1.eof()) {run1=0; }//end1=1;}
				if(run1 & !line1.empty())
				{
					parts1=SplitString(line1);
					head1=parts1.at(0).find('#') != std::string::npos;
					if(!head1)
					{
					//	for(int i=0;i<parts1.size();i++) cout<<parts1.at(i)<<":"<<endl;
						double met=atoi(parts1.at(2).c_str());
						double tot=atoi(parts1.at(3).c_str());
						int chr=atoi(parts1.at(1).c_str());
						chr_met_tot tmp;
						tmp.met=met;
						tmp.tot=tot;
						tmp.chr=parts1.at(0);
						if((met<=tot) & (tot>0) & (met >=0))
						record1.insert(pair<int,  chr_met_tot> (chr, tmp));
					}
				}				

				string line2;
				getline(myfile2,line2);

				std::vector<std::string> parts2;
				if(myfile2.eof()) {run2=0; }//end2=1;}
				if(run2 & !line2.empty())
				{
					parts2=SplitString(line2);
					head2=parts2.at(0).find('#') != std::string::npos;
					if(!head2)
					{
						double met=atoi(parts2.at(2).c_str());
						double tot=atoi(parts2.at(3).c_str());
						int pos=atoi(parts2.at(1).c_str());
						chr_met_tot tmp;
						tmp.met=met;
						tmp.tot=tot;
						tmp.chr=parts2.at(0);
						if((met <= tot) & (tot>0) & (met >=0))
						record2.insert(pair<int,  chr_met_tot> (pos, tmp));
					}
				}

				// Will not check if they come from the same chromosome //
				// compare positons and conduct test for positions where both samples have data //
			
				if((record1.size()!=0) & (record2.size()!=0))
				{
					for(map<int, chr_met_tot >::iterator it=record1.begin();it!=record1.end();it++)
					{
						map<int, chr_met_tot >::iterator itt1;
						map<int, chr_met_tot >::iterator itt2;
						itt1=record2.find(it->first);
						if(itt1!=record2.end())
						{
							itt2=record1.find(itt1->first);
							my_result tmp;
							if((itt1->second.chr) != (itt2->second.chr)) tmp.chr="ERRORCHR";	
							else tmp.chr=itt1->second.chr;

							if((itt2->second.tot >= MINTOT) & (itt1->second.tot >= MINTOT))
							{
								tmp.pos=itt1->first;
								tmp.met1=itt2->second.met;
								tmp.tot1=itt2->second.tot;
								tmp.met2=itt1->second.met;
								tmp.tot2=itt1->second.tot;

								Compare2Prop compare(tmp.met1,tmp.tot1,tmp.met2,tmp.tot2);
								double z_score=compare.compare2();
							
								tmp.z_score=z_score;
								pos_save.insert(pos_save.end(),tmp.pos);
								z_save.insert(z_save.end(),z_score);
								mapResult.insert(pair<int, my_result>(tmp.pos,tmp));
								record1.erase (itt2);
								record2.erase (itt1);
							}						
						}
					}
				}

			}

			/* output positions existed only in file 1 or does not meet with the min total read requirement*/
 			ofstream myfile1;
			string tmp=outdir+"/exclude1.txt";
			myfile1.open (tmp.c_str());
			myfile1 << input1file <<": unmatched with met<tot and tot>0 and met>=0: \n";
			for(map<int, chr_met_tot >::const_iterator it = record1.begin(); it != record1.end(); ++it)
			{
				myfile1 << it->first << " " << it->second.met << " " << it->second.tot << "\n";
			}
			myfile1.close();

			/* output positions existed only in file 2 or does not meet with the min total read requirement*/
			tmp=outdir+"/exclude2.txt";
			myfile1.open (tmp.c_str());
			myfile1 << input2file <<": unmatched with met<tot and tot>0 and met>=0: \n";
			for(map<int, chr_met_tot >::const_iterator it = record2.begin(); it != record2.end(); ++it)
			{
				myfile1 << it->first << " " << it->second.met << " " << it->second.tot << "\n";
			}
			myfile1.close();

			cout<<"Calculating Getis-ord....."<<endl;
			// Getis-ord calculation //
//			RInside R;
			Getis_Ord cal(pos_save,z_save,OutCovFile,CLINK_CPPFLAG,MydistSoftware,NewDistCppSoftware, RFun);
			cal.loadRFun();
			vector <double> gi_final=cal.getis_ord();
			cal.clearR();

			cout<<"Output results for each CpG site....."<<endl;
			// Output results, including chr, pos, met1, tot1, met2, tot2, z_score, g_score //		
			tmp=outdir+"/outputfinal.txt";
			myfile1.open (tmp.c_str());
			myfile1 << "chr\t" <<"pos\t" <<"met1\t" <<"tot1\t" <<"met2\t" <<"tot2\t" <<"z_score\t"<<"g_score\t"<<"\n";
			int i=0;
			string chrom;
			if((mapResult.size())!=(gi_final.size())) myfile1 << "WRONG NOT MATCH" <<endl;
			else
			{
				for(map<int, my_result >::const_iterator it = mapResult.begin(); it != mapResult.end(); ++it)
				{
					myfile1 << it->second.chr<< "\t" << it->second.pos<< "\t" << it->second.met1<< "\t" << it->second.tot1<< "\t" << it->second.met2 << "\t" << it->second.tot2 << "\t"<<it->second.z_score<<"\t"<<gi_final.at(i)<<"\n";
					i=i+1;
					if(i==1) chrom=it->second.chr;
				}
			}


			vector<double> cut_range;
			cut_range.push_back(cutoff);

			// FindDMRs //
			cout<<"Finding DMRs....."<<endl;
			FindDMR finder(pos_save,gi_final,FindDMRfile,outdir,length,tolerance,chrom,cut_range, RFun);
			finder.loadRFun();
			finder.FindDMR_implement();
//			finder.clearR();
		}
	}

}
