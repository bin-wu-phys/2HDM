#include <iostream>
#include <fstream>
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


int main(int argc, char **argv){// 
	int i=0;long ranSeed;
	unsigned long seedForSeed=atoi(argv[1]);
	Ran ran(seedForSeed);

	do{
			do{
				ranSeed=ran.int32();
			}
			while(ranSeed==0);
			if(ranSeed < 0) ranSeed=-ranSeed;
			cout << ranSeed << "\n";
			i++;
	}while(i<=100);
}
