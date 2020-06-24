/*
 *  RpTimesSU2.h
 *  
 *
 *  Created by Bin Wu on 12/22/11.
 *  Copyright 2011 Universitaet Bielefeld. All rights reserved.
 *
 */

#include <complex>

class RpTimesSU2{
public:
	double e[4];
	RpTimesSU2(){};
	RpTimesSU2(double a0, double a1, double a2, double a3);
	
	RpTimesSU2(const RpTimesSU2& u);
		
	void operator+=(RpTimesSU2 u);
	
	void operator-=(RpTimesSU2 u);
	
	void setUnitary();
	
	RpTimesSU2 UMinusUa();
	
	void setZero(){for(int i=0;i<4;i++)e[i]=0;};
	
	double getNorm();
	
	RpTimesSU2 chargeConjugate(){
		return RpTimesSU2(e[0],-e[1],e[2],-e[3]);
	}
	
};


// Implementation of RpTimesSU2 member methods

RpTimesSU2::RpTimesSU2(const RpTimesSU2& u){
	for(int i=0;i<4;i++) e[i]=u.e[i];
}

RpTimesSU2::RpTimesSU2(double a0, double a1, double a2, double a3){
	e[0]=a0;e[1]=a1;e[2]=a2;e[3]=a3;
}

void RpTimesSU2::operator+=(RpTimesSU2 u){
	for(int i=0;i<4;i++)e[i]+=u.e[i];
}

void RpTimesSU2::operator-=(RpTimesSU2 u){
	for(int i=0;i<4;i++)e[i]-=u.e[i];
}

void RpTimesSU2::setUnitary(){
	e[0]=sqrt(1.0-e[1]*e[1]-e[2]*e[2]-e[3]*e[3]);
}

RpTimesSU2 RpTimesSU2::UMinusUa(){// U - U^\dagger
	return RpTimesSU2(0,2.0*e[1],2.0*e[2],2.0*e[3]);
}

/* Dagger(conjugate Transpose) of RpTimesSU2 */

RpTimesSU2 operator~(const RpTimesSU2& u){
	return RpTimesSU2(u.e[0],-u.e[1],-u.e[2],-u.e[3]);
}

RpTimesSU2 operator+(const RpTimesSU2& u1, const RpTimesSU2& u2){
	return RpTimesSU2(u1.e[0]+u2.e[0],u1.e[1]+u2.e[1],u1.e[2]+u2.e[2],u1.e[3]+u2.e[3]);
}

RpTimesSU2 operator-(const RpTimesSU2& u1, const RpTimesSU2& u2){
	return RpTimesSU2(u1.e[0]-u2.e[0],u1.e[1]-u2.e[1],u1.e[2]-u2.e[2],u1.e[3]-u2.e[3]);
}

RpTimesSU2 operator*(const RpTimesSU2& u1, const RpTimesSU2& u2){
	return RpTimesSU2(u1.e[0]*u2.e[0]-u1.e[1]*u2.e[1]-u1.e[2]*u2.e[2]-u1.e[3]*u2.e[3],\
					  u1.e[0]*u2.e[1]+u1.e[1]*u2.e[0]-u1.e[2]*u2.e[3]+u1.e[3]*u2.e[2], \
					  u1.e[0]*u2.e[2]+u1.e[2]*u2.e[0]-u1.e[3]*u2.e[1]+u1.e[1]*u2.e[3], \
					  u1.e[0]*u2.e[3]+u1.e[3]*u2.e[0]-u1.e[1]*u2.e[2]+u1.e[2]*u2.e[1]);
}

RpTimesSU2 operator*(const double& c1, const RpTimesSU2& u){
	return RpTimesSU2(c1*u.e[0],c1*u.e[1],c1*u.e[2],c1*u.e[3]);
}

RpTimesSU2 operator/(const RpTimesSU2& u, const double& c1){
	return RpTimesSU2(u.e[0]/c1,u.e[1]/c1,u.e[2]/c1,u.e[3]/c1);
}

double Tr(const RpTimesSU2& u){
	return 2.0*u.e[0];
}

double Tr(const RpTimesSU2& a,const RpTimesSU2& b){//Tr(a*b)
	return 2.0*(a.e[0]*b.e[0] - a.e[1]*b.e[1] - a.e[2]*b.e[2] - a.e[3]*b.e[3]);
}

double Trsq(const RpTimesSU2& u){// Tr(~u*u)
	return 2.0*(u.e[0]*u.e[0]+u.e[1]*u.e[1]+u.e[2]*u.e[2]+u.e[3]*u.e[3]);
}

double TrAdaggerB(const RpTimesSU2& a, const RpTimesSU2& b){// Tr(~a*b)
	return 2.0*(a.e[0]*b.e[0] + a.e[1]*b.e[1] + a.e[2]*b.e[2] + a.e[3]*b.e[3]);
}

ostream& operator<<(ostream& out, const RpTimesSU2& u){
	return out << "( " << u.e[0] << ", " << u.e[1] << ", " << u.e[2] << ", " << u.e[3] << " )";
}

double max(const RpTimesSU2& u){
	double maxdb=0;
	for(int a=1;a<=3;a++)
		if(fabs(u.e[a])>maxdb) maxdb=u.e[a];
	return maxdb;
}

double cmp(const RpTimesSU2& u1,const RpTimesSU2& u2){
	double maxdb=0,diff;
	for(int a=1;a<=3;a++){
		if(u2.e[a]!=0){
			diff=fabs(1.0-u1.e[a]/u2.e[a]);
		}
		else {
			diff = fabs(u1.e[a]);
		}
		if(diff>maxdb) maxdb=diff;
	}
	return maxdb;
}

#ifdef ANDERS
double relativeDifference(const RpTimesSU2& u1,const RpTimesSU2& u2){// Anders' check
	double err,df;
	RpTimesSU2 sum,diff;
	sum=u1+u2;diff=u1-u2;
	err=sum.e[1]*sum.e[1]+sum.e[2]*sum.e[2]+sum.e[3]*sum.e[3];
	df=diff.e[1]*diff.e[1]+diff.e[2]*diff.e[2]+diff.e[3]*diff.e[3];
	if(err!=0)
		df=sqrt(df/err);
	else
		df=sqrt(df);
	return df;
}
#endif

double RpTimesSU2::getNorm(){
	return sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]+e[3]*e[3]);
}

RpTimesSU2 operator*(const complex<double>& c1, const RpTimesSU2& u){
	return RpTimesSU2(c1.real()*u.e[0]-c1.imag()*u.e[3],c1.real()*u.e[1]-c1.imag()*u.e[2],c1.real()*u.e[2]+c1.imag()*u.e[1],c1.real()*u.e[3]+c1.imag()*u.e[0]);
}
