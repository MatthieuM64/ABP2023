/*LIBRARY FOR RANDOM NUMBER GENERATOR*/

#ifndef DEF_RANDOM_MANGEAT_CPP
#define DEF_RANDOM_MANGEAT_CPP

#include <gsl/gsl_randist.h>

gsl_rng *GSL_r;

void init_gsl_ran()
{
	GSL_r=gsl_rng_alloc(gsl_rng_mt19937);
}

double ran()
{
	double r=gsl_rng_uniform(GSL_r);
	if (r!=0)
	{
		return r;
	}
	else
	{
		return ran();
	}
}

double gaussian()
{
	static bool eval=false;
	static double next;
	if (not eval)
	{
		double phi=2*M_PI*ran();
		double psi=ran();
		double rad=sqrt(-2*log(psi));
		eval=true;
		next=rad*sin(phi);
		return rad*cos(phi);
	}
	else
	{
		eval=false;
		return next;
	}
}

#endif


