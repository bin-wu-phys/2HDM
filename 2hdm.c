//#define _INFO
//#define ANDERS
//#define _C
//#define _P

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

int main(int argc, char **argv)//argv[1] seed, argv[2] kappa
{
	int i_save=0;ofstream out;
	ostringstream ostr;
	
	if(argc!=3){
		cout << "Usage: ./2hdm seed kappa\n";exit(0);
	}
	/* Time labels. */

	prevTime = 0;currentTime = 1;nextTime = 2;
	
	/* Setup lattice and initial conditions. Load in lattice */

	t=0;ncs=0;
	loadParameters(TRUE,atof(argv[2]));
	//loadLattice(atoi(argv[1]));//if N is small
	loadLattice(atoi(argv[1]));
	ostr << "N" << nx << "dt" << setfill('0') <<  setw(4) << int(1000.0*deltat) << "a" << setw(2) << int(100.0*avev);
	if(kappa>=0)
		ostr << "k" << setw(6)<< int(1.0e6*kappa);
	else
		ostr << "km" << setw(6)<< int(-1.0e6*kappa);
	ostr << "w" << warmup;
	if(theseed>=0)
		ostr << "s" << theseed;
	else
		ostr << "sm" << -theseed;	
	
	ncsdot=0;
	n_save = (int)(t_save/deltat); deltat = t_save/n_save;
	int n_save0=n_save;
	
	#ifdef _INFO 
	cout <<"\nResults are going to be saved in " << ostr.str().c_str() << "...\n" << endl;
	#endif
/*
	out.open(ostr.str().c_str(),ios::app);
	out << "#Loading initial fields with the following parameters:\n";
	out << "#N = " << nx << ", a_t/a = " << deltat << ", a*v = " << avev << ", kappa = " << kappa << ", warmup = " << warmup << ", seed = " << theseed << endl;
	out.close();
*/	/* Time marching */
	
	running_time=clock()/CLOCKS_PER_SEC;
	//t_stop=0;
	t_stop+=deltat;
#ifdef _P
	spatialInversion();
#ifdef _C
	chargeConjugate();
	ostr << "CP";
#else
	ostr << "P";
#endif
#else
#ifdef _C
	chargeConjugate();
	ostr << "C";	
#endif
#endif
	out.open(string(ostr.str()+".dat").c_str(),ios::app);
	double aInvGeV=avev/vev;aInvGeV*=aInvGeV;
	out << "# parameters: " << gg << "   " << am11sq/aInvGeV << "   " << am22sq/aInvGeV << "   " << am12sqR/aInvGeV  << "   " << am12sqI/aInvGeV << "   " << lambda1 << "   " << lambda2 << "   " << lambda3 << "   " << lambda4 << "   " << lambda5R << "   " << lambda5I<<"\n";
	out.close();
	//Turn off gauge field
	//turnOffGaugeField();
	//amplifyFieldAmplidues(20.0);
	while (t<=t_stop) {
		//If turn off gauge fields
		//phisAtNextTime();
		// If turn off Higgs field
		//linksAtNextTimeOffHiggses();//
		solveEOS();
		calculateNcs();
		
		if(i_save%n_save==0){
			out.open(string(ostr.str()+".dat").c_str(),ios::app);			
			out <<  t << "   " << ncs << "   "
			<< higgs1WindingNumber() << "   " 
			<< higgs2WindingNumber() << "   " 
			<< calcPhi1sq() << "   "  
			<< calcPhi2sq() << "   " << calcTheta() << "   "
			<< calcKineticEnergy() << "   " << calcPotentialEnergy() << "   " 
			<< calcKineticEnergyA() << "   "
			<< calcPotentialEnergyA()  << "   "
			<< calcKineticEnergyH() << "   " << calcPotentialEnergyH() << "   " << calcPotentialEnergyK()  << " " 
			<< calcGauss() 
			<<endl;
			out.close();
			calcPhi1sqHisgram(string(ostr.str()+"Phi1sq.dat"));
			calcPhi2sqHisgram(string(ostr.str()+"Phi2sq.dat"));
			calcThetaHisgram(string(ostr.str()+"Theta.dat"));
			calcImPhi1aPhi2Hisgram(string(ostr.str()+"ImPhi1aPhi2.dat"));
		}

		i_save++;
		reposition();
		t+=deltat;
		
	}
	//outputPhi2();
	out.open(string(ostr.str()+".dat").c_str(),ios::app);
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
	
	delete[] lattice;
 
	return 0;

}
