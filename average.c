//#define _INFO
//#define ANDERS
/*********Standard libraries and header files.*********/
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits>
using namespace std;
#include <fftw.h>

/*********SU2 algebra and manipulations.***************/

#include "include/RpTimesSU2.h"

/*Lattice and global variables and definitions.********/

#define TRUE 1
#define FALSE 0

#include "include/macros.h"
#include "include/Mcinit.h"
#include "include/lattice.h"
#include "include/2hdm.h"
#include "include/initialFields.h"
#include "include/solveEOS.h"
#include "include/loadAndSaveData.h"
#include "include/codeTest.h"
#include "include/higgsWindingNumber.h"


struct Ran {
	unsigned long long int u,v,w;
	Ran(unsigned long long int j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline unsigned long long int int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline unsigned int uint32() { return (unsigned int)int64(); }
	inline int int32() { return (int)((unsigned int)int64()+numeric_limits<int>::min()); }
};

void getNcsSum(double *ncsSum){
	int i_save=0;int i_ncs=0;
	while (t<=t_stop) {
		
		solveEOS();
		calculateNcs();
		
		if(i_save%n_save==0){
			ncsSum[i_ncs]+=ncs;
			i_ncs++;
			//cout <<  t << "   " << ncs << "   " << calcPhi1sq() << "   " << calcPhi2sq() << "   " << higgs1WindingNumber() << "   " << higgs2WindingNumber() << "   " << calcKineticEnergyA() << "   " << calcKineticEnergyH() << "   " << calcPotentialEnergyA()  << "   " << calcPotentialEnergyH() << "   " << calcPotentialEnergyK()  << " " << calcGauss() <<endl;
		}
		i_save++;
		reposition();
		t+=deltat;
		
	}
}

void loadInitialParameters(int &_warmup){
	FILE *ifile;
	long _theseed;
	
	if((ifile = fopen(INITIALCFG,"r"))==NULL) 
    {    
		printf("Cannot open init.cfg.\n");
		exit(1);
    }
		
	fscanf(ifile,"%ld %i",
		   &_theseed,&_warmup);
	fclose(ifile);                           	
	
}

int main(int argc, char **argv)//argv[1] seed, argv[2] kappa
{
	if(argc!=4){
		cout << "Usage: ./average totalNumberAverage numberForSaving seed\n";exit(0);
	}
		
	int i_time, warmup,
	naverageSave=atoi(argv[2]),		//number for saving data
	naverage=atoi(argv[1]);			// number of initial field configurations(seeds)
	double *ncsKappa0,*ncsKappaMinus,*ncsKappaPlus,*ncsKappa0P,*ncsKappaMinusP,*ncsKappaPlusP;
	long seedForSeed=atoi(argv[3]),ranSeed;
	Ran ran(seedForSeed);
	ofstream out;
	string ifile;
	
	// Input parameters, once for all
	loadParameters(TRUE,atof(argv[2]));
	loadInitialParameters(warmup);
	n_save = (int)(t_save/deltat); deltat = t_save/n_save;i_time = (int)(t_stop/t_save)+1;
		
	// Initialized ncs
	ncsKappa0=new double[i_time];ncsKappaMinus=new double[i_time];ncsKappaPlus=new double[i_time];
	ncsKappa0P=new double[i_time];ncsKappaMinusP=new double[i_time];ncsKappaPlusP=new double[i_time];
	for(int i=0;i<i_time;i++){
		ncsKappa0[i]=0;ncsKappaMinus[i]=0;ncsKappaPlus[i]=0;
		ncsKappa0P[i]=0;ncsKappaMinusP[i]=0;ncsKappaPlusP[i]=0;
	}
	
	/* Time marching */
	
	running_time=clock()/CLOCKS_PER_SEC;
	t_stop+=deltat;
	
	for(int n=1;n<=naverage;n++){
		
		// Prepare initial fields
		do{
			ranSeed=ran.int32();
		}
		while(ranSeed==0);
		//cout << ranSeed << "\n";
		
		ifile=prepareInitialFields(seedForSeed,ranSeed, warmup);

		// kappa = 0
		reset(0);					// Reset first since loadLattice use the reset of prev, current and nextTime 
		loadLatticeFromFile(ifile);	// Load initial fields
		getNcsSum(ncsKappa0);	
		delete[] lattice;

		// kappa = 0
		reset(0);					// Reset first since loadLattice use the reset of prev, current and nextTime 
		loadLatticeFromFile(ifile);	// Load initial fields
		spatialInversion();
		getNcsSum(ncsKappa0P);		
		delete[] lattice;
		
		// kappa = 2.0
		reset(2.0);					// Reset first since loadLattice use the reset of prev, current and nextTime 
		loadLatticeFromFile(ifile);	// Load initial fields		
		getNcsSum(ncsKappaPlus);
		delete[] lattice;

		// kappa = 2.0
		reset(2.0);					// Reset first since loadLattice use the reset of prev, current and nextTime 
		loadLatticeFromFile(ifile);	// Load initial fields		
		spatialInversion();
		getNcsSum(ncsKappaPlusP);
		delete[] lattice;
		
		// kappa = -2.0
		reset(-2.0);				// Reset first since loadLattice use the reset of prev, current and nextTime 
		loadLatticeFromFile(ifile);	// Load initial fields		
		getNcsSum(ncsKappaMinus);
		delete[] lattice;
		
		// kappa = -2.0
		reset(-2.0);				// Reset first since loadLattice use the reset of prev, current and nextTime 
		loadLatticeFromFile(ifile);	// Load initial fields		
		spatialInversion();
		getNcsSum(ncsKappaMinusP);
		delete[] lattice;

		//Save the average Ncs
		if(n%naverageSave==0){
			ostringstream ostr;ostr << "results/N" << nx << "dt" << setfill('0') <<  setw(4) << int(1000.0*deltat) << "a" << setw(2) << int(100.0*avev) << "w" << warmup << "s" << argv[3] << "n" << n << ".dat";
			out.open(ostr.str().c_str(),ios::app);
			for(int j=0;j<i_time;j++) 
				out <<  t_save*j << "   " << ncsKappaPlus[j]/n\/*  << "   " << ncsKappaPlusP[j]/n\*/
				<< "   " << ncsKappa0[j]/n\/* << "   " << ncsKappa0P[j]/n\ 
				<< "   " << ncsKappaMinus[j]/n\/ << "   " << ncsKappaMinusP[j]/n\
				<<endl;
			out.close();
		}
		
		// Delete the temp file	
		system(string("rm "+ifile).c_str());

	}
	
	// Save computation time
	ostringstream ostr;
	ostr << "results/N" << nx << "dt" << setfill('0') <<  setw(4) << int(1000.0*deltat) << "a" << setw(2) << int(100.0*avev) << "w" << warmup << "s" << argv[3] << "RunningTime.dat";
	out.open(ostr.str().c_str(),ios::app);
	running_time=clock()/CLOCKS_PER_SEC - running_time;
	
	if(running_time>=60.0){
		if(running_time>=3600){
			out << "#Calculation is done and the running time is " << running_time/3600.0 << " hours.\n\n\n";
		}
		else {
			out << "#Calculation is done and the running time is " << running_time/60.0 << " minutes.\n\n\n";
		}
	}
	else {
		out << "#Calculation is done and the running time is " << running_time << " seconds.\n\n\n";
	}
		
	out.close();
	
	// Release meomory
	delete[] ncsKappa0,ncsKappaMinus,ncsKappaPlus; 
	
	return 0;


}
