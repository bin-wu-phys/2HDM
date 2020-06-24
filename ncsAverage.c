#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#define NUM 300

int main(){
	ifstream inName,in;
	ofstream out;
	string fname;
	double t, ncs, hwind1, hwind2, phi1, phi2, theta, energyk, energyp, energyak, energyap, energyhk, energyhp,energypk, err;
	double ncsAverage[NUM+1],ncsErr[NUM+1],phi1sq[NUM+1],phi2sq[NUM+1],hw1[NUM+1],hw2[NUM+1];unsigned int n=0,i;
	
	for(i=0;i<=NUM;i++){
		ncsAverage[i]=0;ncsErr[i]=0;phi1sq[i]=0;phi2sq[i]=0;hw1[i]=0;hw2[i]=0;
	}
	
	inName.open("tmp");
	while((inName >> fname)){
		in.open(string("results/final/CPV/max/para02/k2/"+fname).c_str());
		cout << fname << "\n";
		t=0;n++;i=0;
		do{
			in >> t >> ncs >> hwind1 >> hwind2 >> phi1 >> phi2 >>  theta >> energyk >> energyp >> energyak >> energyap >> energyhk >> energyhp >> energypk >> err;
			ncsAverage[i]+=ncs;ncsErr[i]+=(ncs*ncs);phi1sq[i]+=phi1;phi2sq[i]+=phi2;hw1[i]+=hwind1;hw2[i]+=hwind2;i++;
		}while(t<NUM);

		in.close();
	}
	
	out.open("ncsAverage.dat");
	t=0;
	for(i=0;i<=NUM;i++){
		ncsAverage[i]/=n;
		out << t << "   " << ncsAverage[i] << "   " << sqrt(ncsErr[i]/n-ncsAverage[i]*ncsAverage[i]) << "   "  << phi1sq[i]/n  << "   " << phi2sq[i]/n << "   "  << hw1[i]/n << "   "  << hw2[i]/n << "\n";
		t+=1.0;
	}
	
	return 0;
}