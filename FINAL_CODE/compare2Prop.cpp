
#include "compare2Prop.h"
#include <algorithm> 
#include <boost/math/distributions/normal.hpp>
using namespace std;

double factorial(int n)
{
	double result=1;
	for(int i=1;i<=n;i++) result=result*i;
	return result;
};

Compare2Prop::Compare2Prop(int met1_int, int tot1_int, int met2_int, int tot2_int)
{
	SetPara( met1_int,  tot1_int,  met2_int,  tot2_int);
};
void Compare2Prop::SetPara(int met1_int, int tot1_int, int met2_int, int tot2_int)
{
	met1=met1_int;
	met2=met2_int;
	tot1=tot1_int;
	tot2=tot2_int;
};

/*
double Compare2Prop::compare2()
{
	double x1=double(met1);
	double x2=double(met2);
	double n1=double(tot1);
	double n2=double(tot2);
	bool correct=1;

	double est1=x1/n1*1.0;
	double est2=x2/n2*1.0;

	double YATES=0.0;
	if(correct) YATES=0.5;

	double DELTA=est1-est2;
	double signtmp=1;if(DELTA<0) signtmp=-1;
	double tmp=signtmp*(DELTA)/(1/n1+1/n2)*1.0;

	if(YATES>tmp) YATES=tmp;

	double p=0.00;
	p=(x1+x2)/(n1+n2)*1.0;
	double e11=n1*p;
	double e12=n1*(1-p);
	double e21=n2*p;
	double e22=n2*(1-p);
	double x11=x1;
	double x12=n1-x1;
	double x21=x2;
	double x22=n2-x2;
	double sign1=1; if(x11-e11<0) sign1=-1;
	double sign2=1; if(x12-e12<0) sign2=-1;
	double sign3=1; if(x21-e21<0) sign3=-1;
	double sign4=1; if(x22-e22<0) sign4=-1;
	//abs doesn't work for double for some reason, only work for numeric vector
	double c11=(sign1*(x11-e11)-YATES)*(sign1*(x11-e11)-YATES)/e11;
	double c12=(sign2*(x12-e12)-YATES)*(sign2*(x12-e12)-YATES)/e12;
	double c21=(sign3*(x21-e21)-YATES)*(sign3*(x21-e21)-YATES)/e21;
	double c22=(sign4*(x22-e22)-YATES)*(sign4*(x22-e22)-YATES)/e22;
	double stats=sqrt(c11+c12+c21+c22);
	if(DELTA<0) stats=(-1)*stats;
cout<<"FISHER"<<endl;
	cout<<stats<<endl;
	if(!((e11 >=5) & (e12 >=5) & (e21 >=5) & (e22 >=5)))  stats=fisher();
	return stats;
}*/;



double Compare2Prop::compare2()
{
	double x1=double(met1);
	double x2=double(met2);
	double n1=double(tot1);
	double n2=double(tot2);
	bool correct=1;
	double stats;


	if((n1==0) | (n2==0) ) stats=-999.0; //no coverage
	else
	{
		if((x1==0) & (x2==0)) stats=0.0;
		else
		{
			if((x1/n1==1) & (x2/n2==1)) stats=0.0;
			else
			{
				double est1=x1/n1*1.0;
				double est2=x2/n2*1.0;

				double YATES=0.0;
				if(correct) YATES=0.5;
			
				double DELTA=est1-est2;
				double signtmp=1;if(DELTA<0) signtmp=-1;
				double tmp=signtmp*(DELTA)/(1/n1+1/n2)*1.0;

				if(YATES>tmp) YATES=tmp;

				double p=0.00;
				p=(x1+x2)/(n1+n2)*1.0;
				double e11=n1*p;
				double e12=n1*(1-p);
				double e21=n2*p;
				double e22=n2*(1-p);
				double x11=x1;
				double x12=n1-x1;
				double x21=x2;
				double x22=n2-x2;
				double sign1=1; if(x11-e11<0) sign1=-1;
				double sign2=1; if(x12-e12<0) sign2=-1;
				double sign3=1; if(x21-e21<0) sign3=-1;
				double sign4=1; if(x22-e22<0) sign4=-1;
				//abs doesn't work for double for some reason, only work for numeric vector
				double c11=(sign1*(x11-e11)-YATES)*(sign1*(x11-e11)-YATES)/e11;
				double c12=(sign2*(x12-e12)-YATES)*(sign2*(x12-e12)-YATES)/e12;
				double c21=(sign3*(x21-e21)-YATES)*(sign3*(x21-e21)-YATES)/e21;
				double c22=(sign4*(x22-e22)-YATES)*(sign4*(x22-e22)-YATES)/e22;
				stats=sqrt(c11+c12+c21+c22);
				if(DELTA<0) stats=(-1)*stats;
			}
		}
	}
//	if(!((e11 >=5) & (e12 >=5) & (e21 >=5) & (e22 >=5)))  stats=fisher();
	return stats;
};


double Compare2Prop::chisq()
{
	double x1=double(met1);
	double x2=double(met2);
	double n1=double(tot1);
	double n2=double(tot2);
	bool correct=1;

	double est1=x1/n1*1.0;
	double est2=x2/n2*1.0;

	double YATES=0.0;
	if(correct) YATES=0.5;

	double DELTA=est1-est2;
	double signtmp=1;if(DELTA<0) signtmp=-1;
	double tmp=signtmp*(DELTA)/(1/n1+1/n2)*1.0;

	if(YATES>tmp) YATES=tmp;

	double p=0.00;
	p=(x1+x2)/(n1+n2)*1.0;
	double e11=n1*p;
	double e12=n1*(1-p);
	double e21=n2*p;
	double e22=n2*(1-p);
	double x11=x1;
	double x12=n1-x1;
	double x21=x2;
	double x22=n2-x2;
	double sign1=1; if(x11-e11<0) sign1=-1;
	double sign2=1; if(x12-e12<0) sign2=-1;
	double sign3=1; if(x21-e21<0) sign3=-1;
	double sign4=1; if(x22-e22<0) sign4=-1;
	//abs doesn't work for double for some reason, only work for numeric vector
	double c11=(sign1*(x11-e11)-YATES)*(sign1*(x11-e11)-YATES)/e11;
	double c12=(sign2*(x12-e12)-YATES)*(sign2*(x12-e12)-YATES)/e12;
	double c21=(sign3*(x21-e21)-YATES)*(sign3*(x21-e21)-YATES)/e21;
	double c22=(sign4*(x22-e22)-YATES)*(sign4*(x22-e22)-YATES)/e22;
	double stats=sqrt(c11+c12+c21+c22);
	if(DELTA<0) stats=(-1)*stats;
	if(!((e11 >=5) & (e12 >=5) & (e21 >=5) & (e22 >=5)))  stats=-99999;
	return stats;
}

double Compare2Prop::fisher()
{
	int met1_1=tot1-met1;
	int met2_1=tot2-met2;
	int n1_=met1+met2;
	int n2_=met1_1+met2_1;
	int n_1=tot1;
	int n_2=tot2;
	int n=tot1+tot2;
	int min_1=met1;
	if(min_1>met1_1) min_1=met1_1;
	if(min_1>met2) min_1=met2;
	if(min_1>met2_1) min_1=met2_1;

	double p=0.0;
	double deno=factorial(n)/factorial(n-n_1)/factorial(n_1);


	double a1;
	double a2;
	bool check=1;
	if((met1==min_1) & check)
	{
		for(int i=0; i<=min_1; i++)
		{
			int atmp=i;
			int ctmp=n_1-i;
			a1=factorial(n1_)/factorial(n1_-atmp)/factorial(atmp);
			a2=factorial(n2_)/factorial(n2_-ctmp)/factorial(ctmp);
			p=p+a1*a2/deno*1.0;
		}
		check=0 ;

	}
	else if ((met1_1==min_1) & check)
	{
		for(int i=0; i<=min_1; i++)
		{
			int atmp=i;
			int ctmp=n_1-i;
			a1=factorial(n2_)/factorial(n2_-atmp)/factorial(atmp);
			a2=factorial(n1_)/factorial(n1_-ctmp)/factorial(ctmp);
			p=p+a1*a2/deno*1.0;
		}
		check=0 ;

	}		
	else if ((met2==min_1) & check)
	{
		for(int i=0; i<=min_1; i++)
		{
			int atmp=i;
			int ctmp=n_2-i;
			a1=factorial(n2_)/factorial(n2_-atmp)/factorial(atmp);
			a2=factorial(n1_)/factorial(n1_-ctmp)/factorial(ctmp);
			p=p+a1*a2/deno*1.0;
		}
		check=0 ;

	}
	else
	{
		for(int i=0; i<=min_1; i++)
		{
			int atmp=i;
			int ctmp=n_2-i;
			a1=factorial(n2_)/factorial(n2_-atmp)/factorial(atmp);
			a2=factorial(n1_)/factorial(n1_-ctmp)/factorial(ctmp);
			p=p+a1*a2/deno*1.0;
		}
		check=0 ;
	}


	boost::math::normal dist(0.0, 1.0);
	if(p==1) p=0.99999;
	if(p==0) p=0.00001;
	double stats = quantile(dist, 1-p);
	return stats;

}



