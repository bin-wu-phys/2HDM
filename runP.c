#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
using namespace std;


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

#define BEGNUM 0

int main(int argc, char **argv){// 
	ifstream in;
	int i=0+BEGNUM,n;
	unsigned long seedForSeed=atoi(argv[1]);
	long ranSeed;
	double kappa=atof(argv[2]);
	Ran ran(seedForSeed);

	for(int j=0;j<BEGNUM;j++){
		do{
			ranSeed=ran.int32();
		}
		while(ranSeed==0);
	}

	do{
		ostringstream num;
		num << "num" <<seedForSeed;
		if(kappa>0)
			num << "k" << setw(6)<< int(1.0e6*kappa);
		else if(kappa==0)
			num << "k0";
		else
			num << "km" << setw(6)<< int(-1.0e6*kappa);
		num<<".wb";
		system(string("ps aux|grep -c ./2hdmp > "+ num.str()).c_str());
		in.open(num.str().c_str());
		in >> n;
		in.close();
		while (n<27) {
			ostringstream ostr;
			do{
				ranSeed=ran.int32();
			}
			while(ranSeed==0);
			ostr << "./2hdmp " << ranSeed << " " << kappa << "&";
			cout << ostr.str() << "\n";
			system(ostr.str().c_str());
			n++;i++;
		}			
		sleep(300);
	}while(i<25+BEGNUM);
}
