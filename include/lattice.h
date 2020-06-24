/******************************  lattice.h ******************************
 *																		*
 *  Header file for lattice simulation									*
 *	A copy of anders' code with redundant variables omitted				*
 *																		*
 ************************************************************************/


 
#define PI 3.14159265358979323846
#define XUP 0
#define YUP 1
#define ZUP 2
#define ZDOWN 3
#define YDOWN 4
#define XDOWN 5
#define OPP_DIR(dir)    (5-(dir))       // Opposite direction
#define NUMTSLICES	3					// numbers of time slices
#define NUMSDIR		3					// numbers of spatial directions
#define TEMPTIME    3

/* Pauli matrices times i */

RpTimesSU2 iSigma[]={RpTimesSU2(0,1.0,0,0), RpTimesSU2(0,0,1.0,0), RpTimesSU2(0,0,0,1.0)};

/* The lattice is an array of sites.  */

class site {
public:
	struct site *buddy[6];					/* 6 spatial directions */
	RpTimesSU2 link[NUMTSLICES][NUMSDIR];	/* 3 time slices, 3 spatial directions */
	RpTimesSU2 efield[NUMTSLICES+1][NUMSDIR];	/* 3 time slices, 3 spatial directions, U_0n*/
	RpTimesSU2 U0nBar[NUMSDIR];				// \bar{U}_0n at currentTime
	RpTimesSU2 Umn[NUMTSLICES][NUMSDIR];	// Umn[0]=U12, Umn[1]=U23, Umn[2]=U31	
	RpTimesSU2 UmnBar[NUMTSLICES][NUMSDIR];	// Umn[0]=U12, Umn[1]=U23, Umn[2]=U31	
	RpTimesSU2 phi1[NUMTSLICES];			/* 3 time slices + 1 for solving implicit EOS*/
	RpTimesSU2 phiUni;						// used for calculating higgs winding number
	RpTimesSU2 phi2[NUMTSLICES];			/* 3 time slices + 1 for solving implicit EOS*/
	RpTimesSU2 phi2Uni;
	//RpTimesSU2 gaugeUnitary;		// phi_2 = gaugeUnitary*(0,phi2.norm);

	double ncs_loc;
	
	RpTimesSU2 phi1iSigma3[NUMTSLICES];		/* Phi_1 times i sigma^3 */
	RpTimesSU2 phi2iSigma3[NUMTSLICES];		/* Phi_2 times i sigma^3 */
	double TrPhi1sq[NUMTSLICES];			/* Tr(Phi_1 Phi_1^dagger) */
	double TrPhi2sq[NUMTSLICES];			/* Tr(Phi_2 Phi_2^dagger) */
	double TrPhi2Phi1a[NUMTSLICES];			/* Tr(Phi_2 Phi_1^dagger)=Tr(Phi_1 Phi_2^dagger)*/
	double TrPhi2iSigma3Phi1a[NUMTSLICES];	/* Tr(Phi_2 isigma^3 Phi_1^dagger)=- Tr(Phi_1 isigma^3 Phi_2^dagger)*/
	
	// Quantities used in solving EOS to same computation time
	RpTimesSU2 Tr231aUmnBar[NUMTSLICES][NUMSDIR];	// TrPhi2iSigma3Phi1a[NUMTSLICES]*UmnBar[NUMTSLICES][NUMSDIR];
	RpTimesSU2 Tr231aU0nBar[NUMSDIR];				// TrPhi2iSigma3Phi1a[currentTime]*\bar{U}_0n at currentTime
	RpTimesSU2 dU0iUjk1[NUMTSLICES][NUMSDIR];		// Used in linksAtNextTime()
	RpTimesSU2 dU0iUjk2[NUMTSLICES][NUMSDIR];		// Used in linksAtNextTime()

	RpTimesSU2 prod1[6];							// Used in linksAtNextTime()
	RpTimesSU2 prod2[6];							// Used in linksAtNextTime()
	RpTimesSU2 prod3[6];							// Used in linksAtNextTime()
	RpTimesSU2 prod4[6];							// Used in linksAtNextTime()
	
	RpTimesSU2 getUmn(int iTime,int m, int n);
	RpTimesSU2 getUmnBar(int iTime,int m, int n);
	RpTimesSU2 getTr231aUmnBar(int iTime,int m, int n);
	RpTimesSU2 getprod1(int m, int n);
	RpTimesSU2 getprod2(int m, int n);
	RpTimesSU2 getprod3(int m, int n);
	RpTimesSU2 getprod4(int m, int n);
	
	void setUmn(int iTime,int m, int n, RpTimesSU2 UmniTime);
	void setUmnBar(int iTime,int m, int n, RpTimesSU2 UmniTime);
	void setprod1(int m, int n, RpTimesSU2 UmniTime);	
	void setprod2(int m, int n, RpTimesSU2 UmniTime);	
	void setprod3(int m, int n, RpTimesSU2 UmniTime);	
	void setprod4(int m, int n, RpTimesSU2 UmniTime);	
	void updatePhiExpression(int iTime);
	void updateTr231aUmnBar(int iTime);
	void updateTr231aU0nBar(int iTime);
};

RpTimesSU2 site::getprod1(int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=prod1[0];
	else if(m==XUP&&n==ZUP)		
		retUmn=prod1[1];
	else if(m==YUP&&n==XUP)
		retUmn=prod1[2];
	else if(m==YUP&&n==ZUP)
		retUmn=prod1[3];	
	else if(m==ZUP&&n==XUP)
		retUmn=prod1[4];
	else if(m==ZUP&&n==YUP)
		retUmn=prod1[5];
	return retUmn;
}

void site::setprod1(int m, int n, RpTimesSU2 UmniTime){
	if(m==XUP&&n==YUP)
		prod1[0]=UmniTime;
	else if(m==XUP&&n==ZUP)		
		prod1[1]=UmniTime;
	else if(m==YUP&&n==XUP)
		prod1[2]=UmniTime;
	else if(m==YUP&&n==ZUP)
		prod1[3]=UmniTime;	
	else if(m==ZUP&&n==XUP)
		prod1[4]=UmniTime;
	else if(m==ZUP&&n==YUP)
		prod1[5]=UmniTime;
}

RpTimesSU2 site::getprod2(int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=prod2[0];
	else if(m==XUP&&n==ZUP)		
		retUmn=prod2[1];
	else if(m==YUP&&n==XUP)
		retUmn=prod2[2];
	else if(m==YUP&&n==ZUP)
		retUmn=prod2[3];	
	else if(m==ZUP&&n==XUP)
		retUmn=prod2[4];
	else if(m==ZUP&&n==YUP)
		retUmn=prod2[5];
	return retUmn;
}

void site::setprod2(int m, int n, RpTimesSU2 UmniTime){
	if(m==XUP&&n==YUP)
		prod2[0]=UmniTime;
	else if(m==XUP&&n==ZUP)		
		prod2[1]=UmniTime;
	else if(m==YUP&&n==XUP)
		prod2[2]=UmniTime;
	else if(m==YUP&&n==ZUP)
		prod2[3]=UmniTime;	
	else if(m==ZUP&&n==XUP)
		prod2[4]=UmniTime;
	else if(m==ZUP&&n==YUP)
		prod2[5]=UmniTime;
}


RpTimesSU2 site::getprod3(int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=prod3[0];
	else if(m==XUP&&n==ZUP)		
		retUmn=prod3[1];
	else if(m==YUP&&n==XUP)
		retUmn=prod3[2];
	else if(m==YUP&&n==ZUP)
		retUmn=prod3[3];	
	else if(m==ZUP&&n==XUP)
		retUmn=prod3[4];
	else if(m==ZUP&&n==YUP)
		retUmn=prod3[5];
	return retUmn;
}

void site::setprod3(int m, int n, RpTimesSU2 UmniTime){
	if(m==XUP&&n==YUP)
		prod3[0]=UmniTime;
	else if(m==XUP&&n==ZUP)		
		prod3[1]=UmniTime;
	else if(m==YUP&&n==XUP)
		prod3[2]=UmniTime;
	else if(m==YUP&&n==ZUP)
		prod3[3]=UmniTime;	
	else if(m==ZUP&&n==XUP)
		prod3[4]=UmniTime;
	else if(m==ZUP&&n==YUP)
		prod3[5]=UmniTime;
}


RpTimesSU2 site::getprod4(int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=prod4[0];
	else if(m==XUP&&n==ZUP)		
		retUmn=prod4[1];
	else if(m==YUP&&n==XUP)
		retUmn=prod4[2];
	else if(m==YUP&&n==ZUP)
		retUmn=prod4[3];	
	else if(m==ZUP&&n==XUP)
		retUmn=prod4[4];
	else if(m==ZUP&&n==YUP)
		retUmn=prod4[5];
	return retUmn;
}

void site::setprod4(int m, int n, RpTimesSU2 UmniTime){
	if(m==XUP&&n==YUP)
		prod4[0]=UmniTime;
	else if(m==XUP&&n==ZUP)		
		prod4[1]=UmniTime;
	else if(m==YUP&&n==XUP)
		prod4[2]=UmniTime;
	else if(m==YUP&&n==ZUP)
		prod4[3]=UmniTime;	
	else if(m==ZUP&&n==XUP)
		prod4[4]=UmniTime;
	else if(m==ZUP&&n==YUP)
		prod4[5]=UmniTime;
}

RpTimesSU2 site::getUmn(int iTime,int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=Umn[iTime][0];
	else if(m==YUP&&n==XUP)
		retUmn=~Umn[iTime][0];
	else if(m==YUP&&n==ZUP)
		retUmn=Umn[iTime][1];
	else if(m==ZUP&&n==YUP)
		retUmn=~Umn[iTime][1];	
	else if(m==ZUP&&n==XUP)
		retUmn=Umn[iTime][2];
	else if(m==XUP&&n==ZUP)
		retUmn=~Umn[iTime][2];
	return retUmn;
}

RpTimesSU2 site::getUmnBar(int iTime,int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=UmnBar[iTime][0];
	else if(m==YUP&&n==XUP)
		retUmn=~UmnBar[iTime][0];
	else if(m==YUP&&n==ZUP)
		retUmn=UmnBar[iTime][1];
	else if(m==ZUP&&n==YUP)
		retUmn=~UmnBar[iTime][1];	
	else if(m==ZUP&&n==XUP)
		retUmn=UmnBar[iTime][2];
	else if(m==XUP&&n==ZUP)
		retUmn=~UmnBar[iTime][2];
	return retUmn;
}

RpTimesSU2 site::getTr231aUmnBar(int iTime,int m, int n){
	RpTimesSU2 retUmn;
	if(m==XUP&&n==YUP)
		retUmn=Tr231aUmnBar[iTime][0];
	else if(m==YUP&&n==XUP)
		retUmn=~Tr231aUmnBar[iTime][0];
	else if(m==YUP&&n==ZUP)
		retUmn=Tr231aUmnBar[iTime][1];
	else if(m==ZUP&&n==YUP)
		retUmn=~Tr231aUmnBar[iTime][1];	
	else if(m==ZUP&&n==XUP)
		retUmn=Tr231aUmnBar[iTime][2];
	else if(m==XUP&&n==ZUP)
		retUmn=~Tr231aUmnBar[iTime][2];
	return retUmn;
}

void site::setUmn(int iTime,int m, int n, RpTimesSU2 UmniTime){
	if(m==XUP&&n==YUP)
		Umn[iTime][0]=UmniTime;
	else if(m==YUP&&n==ZUP)
		Umn[iTime][1]=UmniTime;
	else if(m==ZUP&&n==XUP)
		Umn[iTime][2]=UmniTime;
}

void site::setUmnBar(int iTime,int m, int n, RpTimesSU2 UmniTime){
	if(m==XUP&&n==YUP)
		UmnBar[iTime][0]=UmniTime;
	else if(m==YUP&&n==ZUP)
		UmnBar[iTime][1]=UmniTime;
	else if(m==ZUP&&n==XUP)
		UmnBar[iTime][2]=UmniTime;
}

void site::updatePhiExpression(int iTime){
	phi1iSigma3[iTime] = phi1[iTime]*iSigma[2];
	phi2iSigma3[iTime] = phi2[iTime]*iSigma[2];
	TrPhi1sq[iTime] = Trsq(phi1[iTime]);
	TrPhi2sq[iTime] = Trsq(phi2[iTime]);
	TrPhi2Phi1a[iTime] = TrAdaggerB(phi1[iTime],phi2[iTime]);
	TrPhi2iSigma3Phi1a[iTime] = TrAdaggerB(phi1[iTime],phi2iSigma3[iTime]);
}

void site::updateTr231aUmnBar(int iTime){
	Tr231aUmnBar[iTime][0]=TrPhi2iSigma3Phi1a[iTime]*UmnBar[iTime][0];
	Tr231aUmnBar[iTime][1]=TrPhi2iSigma3Phi1a[iTime]*UmnBar[iTime][1];
	Tr231aUmnBar[iTime][2]=TrPhi2iSigma3Phi1a[iTime]*UmnBar[iTime][2];
}

void site::updateTr231aU0nBar(int iTime){
	Tr231aU0nBar[0]=TrPhi2iSigma3Phi1a[iTime]*(U0nBar[0].UMinusUa());
	Tr231aU0nBar[1]=TrPhi2iSigma3Phi1a[iTime]*(U0nBar[1].UMinusUa());
	Tr231aU0nBar[2]=TrPhi2iSigma3Phi1a[iTime]*(U0nBar[2].UMinusUa());
}
