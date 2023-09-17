/*C++ CODE - MANGEAT MATTHIEU - 2023 */
/*SEDIMENTATION OF INTERACTING ACTIVE BROWNIAN PARTICLES IN A BOX*/

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>

//Personal libraries.
#include "lib/random.cpp"
#include "lib/special_functions.cpp"

using namespace std;

void InformationRun();
string TimeRun(const string &c);
double square(const double &r);
double cube(const double &r);

//active Brownian particle.
class particle
{
	public:
	
	int index;
	double x,y,theta; //position+orientation.
	double a; //random radius.
	double vx,vy; //velocity.
	double Fpartx,Fparty; //particle-particle interaction.
	
	particle(const int &index_, const double &lsed, const double &sigma, const double &LX, const double &LY, const int &init);
	void updateTheta(const double &Ddt);
	void updatePos(const double &alpha, const double &F0, const double &LX, const double &LY, const double &dt);
};

//Initial condition.
particle::particle(const int &index_, const double &lsed, const double &sigma, const double &LX, const double &LY, const int &init)
{
	index=index_;
	a=0.5*sigma;
	vx=0, vy=0;
	Fpartx=0, Fparty=0;
	
	const double lx=LX-sigma, ly=LY-sigma;
	
	x=a+lx*ran();
	
	const double psi=ran();
	//Uniform initial condition.
	if (init==0)
	{
		y=a+ly*psi;
	}
	//Exponential initial condition (sedimentation profile for non-interacting ABPs, without vertical walls).
	else
	{
		y=a-lsed*log(1-psi+psi*exp(-ly/lsed));
	}
	
	theta=2*M_PI*ran();
}

//Update the orientation of the particle.
void particle::updateTheta(const double &Ddt)
{
	theta+= sqrt(2*Ddt)*gaussian();
	theta-=floor(theta/(2*M_PI))*2*M_PI;
}

//Update the position of the particle.
void particle::updatePos(const double &alpha, const double &F0, const double &LX, const double &LY, const double &dt)
{
	//Wall-particle interaction.
	static const double r0=pow(2,1./6)*a;
	double Fwallx=0, Fwally=0;
	if (x<r0)
	{
		const double invR2=square(a/x), invR6=cube(invR2);
		Fwallx+=48*(invR6/x)*(invR6-0.5);
	}
	else if (LX-x<r0)
	{
		const double invR2=square(a/(LX-x)), invR6=cube(invR2);
		Fwallx-=48*(invR6/(LX-x))*(invR6-0.5);
	}
	
	if (y<r0)
	{
		const double invR2=square(a/y), invR6=cube(invR2);
		Fwally+=48*(invR6/y)*(invR6-0.5);
	}
	else if (LY-y<r0)
	{
		const double invR2=square(a/(LY-y)), invR6=cube(invR2);
		Fwally-=48*(invR6/(LY-y))*(invR6-0.5);
	}
	
	//Velocity.
	vx=F0*Fpartx + Fwallx + cos(theta);
	vy=F0*Fparty + Fwally + sin(theta) - alpha;
	
	//Update of the position.
	x += vx*dt;
	y += vy*dt;
	
	if (x<0 or x>LX or y<0 or y>LY)
	{
		cerr << "Position outside the boundaries: index=" << index << " x=" << x << " y=" << y << endl;
	}
}

//Distance between two particles, divided by the sum of radii.
double distance2(const particle &ABP0, const particle &ABP1)
{
	return (square(ABP0.x-ABP1.x)+square(ABP0.y-ABP1.y))/square(ABP0.a+ABP1.a);
}

//Accept a new particle with hard-core criterion.
bool acceptParticle(const particle &ABP0, const vector<particle> &ABP)
{
	bool accept=true;
	for (int i=0; i<ABP.size(); i++)
	{
		if (distance2(ABP0,ABP[i])<1)
		{
			accept=false;
			break;
		}
	}
	return accept;
}

//Calculation of the particle-particle interaction.
void interparticleForce(const double &dmax, const int &lx, const int &ly, const int &Npart, vector<particle> &ABP)
{
	//Create a matrix with particle locations, box[i][j] regroups the particle indices with i*dmax<x<(i+1)*dmax and j*dmax<y<(j+1)*dmax.
	vector< vector< vector<int> > > box(lx,vector< vector<int> >(ly,vector<int>(0)));
	for (int i=0; i<Npart; i++)
	{
		box[int(ABP[i].x/dmax)][int(ABP[i].y/dmax)].push_back(i);
	}
	
	//Calculate the force F_i of each particles.
	for (int i=0; i<Npart; i++)
	{
		const int X0=int(ABP[i].x/dmax), Y0=int(ABP[i].y/dmax);
		//Take only the contribution for particles in neighbour boxes and with a distance smaller than 1.
		ABP[i].Fpartx=0., ABP[i].Fparty=0.;
		for (int X=max(0,X0-1);X<=min(X0+1,lx-1);X++)
		{
			for (int Y=max(0,Y0-1);Y<=min(Y0+1,ly-1);Y++)
			{
				const vector<int> neighbours=box[X][Y];
				for (int l=0; l<neighbours.size(); l++)
				{
					const int k=neighbours[l];
					if (k!=i and distance2(ABP[i],ABP[k])<1)
					{
						const double dx=ABP[i].x-ABP[k].x;
						const double dy=ABP[i].y-ABP[k].y;
						const double r=sqrt(square(dx)+square(dy));
						const double force=(ABP[i].a+ABP[k].a)/r-1.;
										
						ABP[i].Fpartx+=force*dx;
						ABP[i].Fparty+=force*dy;
					}
				}
			}
		}
	}
}

//Export the position of all particles.
void exportPosition(const double &Pe, const double &alpha, const double &F0, const double &LX, const double &LY, const int &init, const int &RAN, const double &texp, const int &Npart, const vector<particle> &ABP)
{
	int folder=system("mkdir -p data_ABP2d_int_position/");
	stringstream ss;
	ss << "./data_ABP2d_int_position/ABP2d_int_position_Pe=" << Pe << "_alpha=" << alpha << "_F0=" << F0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << texp << ".txt";
	string name=ss.str();
	
	//Write in the file Y:column X:lines.
	ofstream file(name.c_str(),ios::trunc);
	file.precision(6);
	
	for(int k=0;k<Npart;k++)
	{
		file << k+1 << " " << ABP[k].x << " " << ABP[k].y << " " << ABP[k].a << " " << ABP[k].theta << endl;
	}
	file.close();
}

//Update time-averaged steady-state quantities.
void updateSteadyState(vector< vector<int> > &RHO, vector< vector<double> > &PX, vector< vector<double> > &PY, vector< vector<double> > &JX, vector< vector<double> > &JY, const int &Npart, const vector<particle> &ABP)
{
	for(int k=0;k<Npart;k++)
	{
		const int X0=int(ABP[k].x), Y0=int(ABP[k].y);
		RHO[X0][Y0]++;
		PX[X0][Y0]+=cos(ABP[k].theta);
		PY[X0][Y0]+=sin(ABP[k].theta);
		JX[X0][Y0]+=ABP[k].vx;
		JY[X0][Y0]+=ABP[k].vy;
	}
}

//Export the density in a file.
void exportDensity(const vector< vector<int> > &RHO, const double &Pe, const double &alpha, const double &F0, const double &LX, const double &LY, const int &init, const int &RAN, const int &Nstat)
{
	int folderProba=system("mkdir -p data_ABP2d_int_density/");
	stringstream ssProba;
	ssProba << "./data_ABP2d_int_density/ABP2d_int_rho_Pe=" << Pe << "_alpha=" << alpha << "_F0=" << F0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nomProba=ssProba.str();
	
	//Write in the file Y:column X:lines.
	ofstream fileProba(nomProba.c_str(),ios::trunc);
	fileProba.precision(6);
	for (int y0=0; y0<LY; y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			fileProba << double(RHO[x0][y0])/Nstat << " ";
		}
		fileProba << endl;
	}
	fileProba.close();
}

//Export the polarization in a file.
void exportPolarization(const vector< vector<double> > &PX, const vector< vector<double> > &PY, const double &Pe, const double &alpha, const double &F0, const double &LX, const double &LY, const int &init, const int &RAN, const int &Nstat)
{
	int folderPolarization=system("mkdir -p data_ABP2d_int_polarization/");		
	stringstream ssPx, ssPy;
	ssPx << "./data_ABP2d_int_polarization/ABP2d_int_px_Pe=" << Pe << "_alpha=" << alpha << "_F0=" << F0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string namePx=ssPx.str();
	ssPy << "./data_ABP2d_int_polarization/ABP2d_int_py_Pe=" << Pe << "_alpha=" << alpha << "_F0=" << F0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string namePy=ssPy.str();
	
	//Write in the file Y:column X:lines.
	ofstream filePx(namePx.c_str(),ios::trunc);
	ofstream filePy(namePy.c_str(),ios::trunc);
	filePx.precision(6);
	filePy.precision(6);
	for (int y0=0; y0<LY; y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			filePx << PX[x0][y0]/Nstat << " ";
			filePy << PY[x0][y0]/Nstat << " ";
		}
		filePx << endl;
		filePy << endl;
	}
	filePx.close();
	filePy.close();
}

//Export the current in a file.
void exportCurrent(const vector< vector<double> > &JX, const vector< vector<double> > &JY, const double &Pe, const double &alpha, const double &F0, const double &LX, const double &LY, const int &init, const int &RAN, const int &Nstat)
{
	int folderCurrent=system("mkdir -p data_ABP2d_int_current/");		
	stringstream ssJx, ssJy;
	ssJx << "./data_ABP2d_int_current/ABP2d_int_jx_Pe=" << Pe << "_alpha=" << alpha << "_F0=" << F0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameJx=ssJx.str();
	ssJy << "./data_ABP2d_int_current/ABP2d_int_jy_Pe=" << Pe << "_alpha=" << alpha << "_F0=" << F0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameJy=ssJy.str();
	
	//Write in the file Y:column X:lines.
	ofstream fileJx(nameJx.c_str(),ios::trunc);
	ofstream fileJy(nameJy.c_str(),ios::trunc);
	fileJx.precision(6);
	fileJy.precision(6);
	for (int y0=0; y0<LY; y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			fileJx << JX[x0][y0]/Nstat << " ";
			fileJy << JY[x0][y0]/Nstat << " ";
		}
		fileJx << endl;
		fileJy << endl;
	}
	fileJx.close();
	fileJy.close();
}

//Read command line for physical and numerical parameters.
void ReadCommandLine(int argc, char** argv, double &dt, double &teq, double &tmax, int &Npart, double &Pe, double &alpha, double &F0, double &LX, double &LY, int &init, int &RAN)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-dt=" ))
		{
			dt=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-teq=" ))
		{
			teq=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-Npart=" ))
		{
			Npart=atoi(argv[i]+7);
		}
		else if (strstr(argv[i], "-Pe=" ))
		{
			Pe=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-alpha=" ))
		{
			alpha=atof(argv[i]+7);
		}
		else if (strstr(argv[i], "-F0=" ))
		{
			F0=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

//Main code.
int main(int argc, char *argv[])
{
	//Default parameters.
	//Physical parameters (Pe=Péclet number, alpha=ratio of velocities, F0=particle-particle force strenght, LX,LY=size of the box).
	double Pe=30, alpha=0.2, F0=100, LX=100, LY=400;
	//Numerical parameters (dt=integration time-step, teq=equilibration time, tmax=maximal time, Npart=number of particles, init=initial condition, RAN=index of RNG).
	//init=0: initial random height in the box, init=1: random height chosen with a sedimented profile such that p(y)=exp(-y/lsed).
	double dt=0.001, teq=100, tmax=1000;
	int Npart=5000, init=1, RAN=0;
	
	//Read the parameters taken in argument.
	ReadCommandLine(argc,argv,dt,teq,tmax,Npart,Pe,alpha,F0,LX,LY,init,RAN);

	//Start the random number generator.
	init_gsl_ran();
	cout << "GSL index = " << RAN << "\n";
	gsl_rng_set(GSL_r,RAN);
	
	//Maximal number of timesteps.
	const int Nsteps=int(tmax/dt)+1;
	const int NstepsEq=int(teq/dt)+1;
	
	//Molecular dynamics boxes.
	const double dmax=1.25; //check that this value is larger than the max distance needed and (LX/dmax, LY/dmax) are integers.
	const int lx=int(LX/dmax), ly=int(LY/dmax);
	if (fabs(LX/dmax-lx)>1e-8 or fabs(LY/dmax-ly)>1e-8)
	{
		cerr << "BAD MOLECULAR DYNAMICS BOXES: LX/dmax=" << LX/dmax << " LY/dmax=" << LY/dmax << " (must be integers)"  << endl;
		abort();
	}
	
	//Initialize all the ABPs with lsed given for non-interacting ABPs, without vertical walls.
	double Peg=Pe*alpha;
	double lsed=(1+0.5*Pe*Pe)/Peg;
	
	vector<particle> ABP;
	for(int k=0;k<Npart;k++)
	{
		
		const double sigma=0.8+0.4*ran(); //polydisperse particles.
		if (sigma>dmax)
		{
			cerr << "BAD MOLECULAR DYNAMICS BOXES: sigma/dmax=" << sigma/dmax << " (must be smaller than 1)" << endl;
			abort();
		}
		particle ABP0(k,lsed,sigma,LX,LY,init);
		while(not acceptParticle(ABP0,ABP))
		{
			ABP0=particle(k,lsed,sigma,LX,LY,init);
		}
		ABP.push_back(ABP0);
		//cout << k+1 << " " << ABP[k].x << " " << ABP[k].y << " " << ABP[k].a << endl;
	}
	
	//Steady-state quantities.
	vector< vector<int> > RHO(LX,vector<int>(LY,0.));
	vector< vector<double> > PX(LX,vector<double>(LY,0.)), PY(LX,vector<double>(LY,0.));
	vector< vector<double> > JX(LX,vector<double>(LY,0.)), JY(LX,vector<double>(LY,0.));
	int Nstat=0;
	
	//Time-evolution.
	for(int t=0;t<Nsteps+1;t++)
	{
		if (t%10000==0)
		{
			cout << "time=" << t*dt << running_time.TimeRun(" ") << endl;
		}
		
		//Export the position of particles (dynamics).
		if (RAN==0 and t%1000==0)
		{
			exportPosition(Pe,alpha,F0,LX,LY,init,RAN,t*dt,Npart,ABP);
		}
		
		//Update the steady-state quantities.
		if (t>NstepsEq and t%1000==0)
		{
			updateSteadyState(RHO,PX,PY,JX,JY,Npart,ABP);
			Nstat++;
			if (Nstat%100==0)
			{
				exportDensity(RHO,Pe,alpha,F0,LX,LY,init,RAN,Npart*Nstat);
				exportPolarization(PX,PY,Pe,alpha,F0,LX,LY,init,RAN,Npart*Nstat);
				exportCurrent(JX,JY,Pe,alpha,F0,LX,LY,init,RAN,Npart*Nstat);
				cout << "New exportation of data with Nstat=" << Nstat << " at time t=" << t*dt << endl;
			}
		}
		
		//Update the orientation of particles.
		for(int k=0;k<Npart;k++)
		{
			ABP[k].updateTheta(dt/Pe);
		}
		//Calculate the particle-particle interaction.
		interparticleForce(dmax,lx,ly,Npart,ABP);
		//Update the position of particles.
		for(int k=0;k<Npart;k++)
		{
			ABP[k].updatePos(alpha,F0,LX,LY,dt);
		}
	}
	
	return 0;
}
