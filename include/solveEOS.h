#include "tensors.h"
#define ERROR 1.0e-6
#define NTRIAL 1000

// Calculate ncs_loc

void calculateNcsLoc(){
	for(int x=0;x<volume1;x++){
		lattice[x].ncs_loc=( Tr( lattice[x].U0nBar[XUP].UMinusUa(), lattice[x].getUmnBar(currentTime,YUP, ZUP))\
							+ Tr( lattice[x].U0nBar[YUP].UMinusUa(),lattice[x].getUmnBar(currentTime,ZUP, XUP) )\
							+ Tr( lattice[x].U0nBar[ZUP].UMinusUa(),lattice[x].getUmnBar(currentTime,XUP, YUP) ) );
	}
}

/* Calculate phi1 and phi2 at the next time step */

void phisAtNextTime(){
	int i,dir;
	RpTimesSU2 ddphi;
	
	//if(kappa!=0)
	calculateNcsLoc();
	/* phi1 and phi2 at the next time step */	
	
	for (i=0; i<volume1; i++) {
		ddphi.setZero();
		for (dir=XUP; dir<=ZUP; dir++) {
			ddphi=ddphi+lattice[i].link[currentTime][dir]*(lattice[i].buddy[dir]->phi1[currentTime])\
			+(~lattice[i].buddy[OPP_DIR(dir)]->link[currentTime][dir])*lattice[i].buddy[OPP_DIR(dir)]->phi1[currentTime];
		}
		lattice[i].phi1[nextTime]=  deltatsq*ddphi - lattice[i].phi1[prevTime]\
		+ (2.0 - deltatsq*( 0.5*(lattice[i].TrPhi1sq[currentTime] - am11sq) + 6.0) - beta31*lattice[i].TrPhi2sq[currentTime])*lattice[i].phi1[currentTime]\
		+ (beta121R - beta21a1*lattice[i].TrPhi2Phi1a[currentTime] - beta51I*lattice[i].TrPhi2iSigma3Phi1a[currentTime])*lattice[i].phi2[currentTime]\
		+ (beta121I + beta231a1*lattice[i].TrPhi2iSigma3Phi1a[currentTime] - beta51I*lattice[i].TrPhi2Phi1a[currentTime] + betak1*lattice[i].ncs_loc)*lattice[i].phi2iSigma3[currentTime];
		
		ddphi.setZero();
		for (dir=XUP; dir<=ZUP; dir++) {
			ddphi=ddphi+lattice[i].link[currentTime][dir]*(lattice[i].buddy[dir]->phi2[currentTime])\
			+(~lattice[i].buddy[OPP_DIR(dir)]->link[currentTime][dir])*lattice[i].buddy[OPP_DIR(dir)]->phi2[currentTime];
		}
		lattice[i].phi2[nextTime]= deltatsq*ddphi - lattice[i].phi2[prevTime]\
		+ (2.0 - deltatsq*( 0.5*(lattice[i].TrPhi2sq[currentTime] - am22sq) + 6.0) - beta32*lattice[i].TrPhi1sq[currentTime])*lattice[i].phi2[currentTime]\
		+ (beta122R - beta21a2*lattice[i].TrPhi2Phi1a[currentTime] - beta52I*lattice[i].TrPhi2iSigma3Phi1a[currentTime])*lattice[i].phi1[currentTime]\
		+ ( beta52I*lattice[i].TrPhi2Phi1a[currentTime] - beta122I - beta231a2*lattice[i].TrPhi2iSigma3Phi1a[currentTime] -betak2*lattice[i].ncs_loc)*lattice[i].phi1iSigma3[currentTime];
		
		lattice[i].updatePhiExpression(nextTime);
		lattice[i].updateTr231aUmnBar(nextTime);
	}

}

/* Calculate links at the next time step */
void linksAtNextTimeKappa0(){//kappa = 0
	for (int x=0; x<volume1; x++ ) {
		for (int j=XUP; j<=ZUP; j++) {
			lattice[x].efield[TEMPTIME][j]=lattice[x].efield[prevTime][j];
			for (int n=XUP; n<=ZUP; n++) {
				if (n!=j) {
					lattice[x].efield[TEMPTIME][j]-=deltatsq*(lattice[x].getUmn(currentTime,j,n)\
					-(~lattice[x].buddy[OPP_DIR(n)]->link[currentTime][n])\
					*lattice[x].buddy[OPP_DIR(n)]->getUmn(currentTime,j,n)\
					*lattice[x].buddy[OPP_DIR(n)]->link[currentTime][n]);
					//( spatialPlaquette(j,n,x,currentTime) - spatialBackwardPlaquette(j,n,x,currentTime) ) ));
				}
			}
			lattice[x].efield[TEMPTIME][j]-=betaE1*(lattice[x].link[currentTime][j]*lattice[x].buddy[j]->phi1[currentTime]*(~lattice[x].phi1[currentTime]));
			lattice[x].efield[TEMPTIME][j]-=betaE2*(lattice[x].link[currentTime][j]*lattice[x].buddy[j]->phi2[currentTime]*(~lattice[x].phi2[currentTime]));
			lattice[x].efield[TEMPTIME][j].setUnitary();
			
			lattice[x].efield[currentTime][j]=lattice[x].efield[TEMPTIME][j];
			lattice[x].link[nextTime][j]=lattice[x].efield[currentTime][j]*lattice[x].link[currentTime][j];
		}
	}
	updateUmn(nextTime);updateUmnBar(nextTime);updateU0nBar();updateTr231aU0nBar();
}

void setProds(){
	for (int x=0; x<volume1; x++ ) {
		for (int n=XUP; n<=ZUP; n++) {
			for (int k=XUP; k<=ZUP; k++) {
				if(k!=n){
					lattice[x].setprod1(n,k,(~lattice[x].link[currentTime][n])*lattice[x].getUmn(currentTime,n,k));
					lattice[x].setprod2(k,n,(~(lattice[x].link[currentTime][k]*lattice[x].buddy[k]->link[currentTime][n])));
					lattice[x].setprod3(k,n,lattice[x].buddy[OPP_DIR(k)]->getUmn(currentTime,n,k)*lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k]);
					lattice[x].setprod4(k,n,lattice[x].getUmn(currentTime,n,k)-(~lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k])*lattice[x].buddy[OPP_DIR(k)]->getUmn(currentTime,k,n)*lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k]);
				}
			}
		}
	}
}

void linksAtNextTime(){
	
	//get trial solution (nextTime)
	for (int x=0; x<volume1; x++ ) {
		for (int n=XUP; n<=ZUP; n++) {
			lattice[x].efield[nextTime][n]=lattice[x].efield[TEMPTIME][n];
			if(betaEk!=0){//
				RpTimesSU2 rhon(0,0,0,0);
				for (int i=XUP; i<=ZUP; i++) {
					if(i!=n){
						for (int k=i+1; k<=ZUP; k++) {
							if((k!=n)){
								//cout << n << "   " << i << "   " << k << "   " << eabc(n,i,k) << "\n";
								lattice[x].dU0iUjk1[currentTime][n]=( (lattice[x].getTr231aUmnBar(currentTime,i,k)+lattice[x].getTr231aUmnBar(nextTime,i,k) ).UMinusUa());
								lattice[x].dU0iUjk2[currentTime][n]=lattice[x].link[nextTime][n]*(( lattice[x].buddy[n]->getTr231aUmnBar(currentTime,i,k)+lattice[x].buddy[n]->getTr231aUmnBar(nextTime,i,k)).UMinusUa())*(~lattice[x].link[currentTime][n]);
								if(eabc(n,i,k)==1.0){
									rhon+=(lattice[x].dU0iUjk1[currentTime][n]*lattice[x].efield[currentTime][n]-lattice[x].efield[prevTime][n]*lattice[x].dU0iUjk1[prevTime][n]\
										   \
										   +lattice[x].dU0iUjk2[currentTime][n]-lattice[x].dU0iUjk2[prevTime][n]
										   \
										   \
										   +lattice[x].getprod4(k,n)*(lattice[x].Tr231aU0nBar[i])\
										   \
										   -lattice[x].getprod4(i,n)*(lattice[x].Tr231aU0nBar[k])\
										   \
										   +lattice[x].link[currentTime][n]*(\
																			 //
											(lattice[x].buddy[n]->Tr231aU0nBar[i])*lattice[x].getprod1(n,k)\
											\
											+lattice[x].buddy[n]->link[currentTime][k]\
											*(lattice[x].buddy[n]->buddy[k]->Tr231aU0nBar[i])*lattice[x].getprod2(k,n)\
											\
											-(lattice[x].buddy[n]->Tr231aU0nBar[k])*lattice[x].getprod1(n,i)\
											\
											-lattice[x].buddy[n]->link[currentTime][i]\
											*(lattice[x].buddy[n]->buddy[i]->Tr231aU0nBar[k])*lattice[x].getprod2(i,n)\
											//
											)\
										   \
										   +lattice[x].link[currentTime][k]*(lattice[x].buddy[k]->Tr231aU0nBar[i])*lattice[x].getprod1(k,n)\
										   \
										   -lattice[x].link[currentTime][i]*(lattice[x].buddy[i]->Tr231aU0nBar[k])*lattice[x].getprod1(i,n)\
										   \
										   -(~lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k])*(\
											//
											(lattice[x].buddy[OPP_DIR(k)]->Tr231aU0nBar[i])*lattice[x].getprod3(k,n)\
											\
											+lattice[x].buddy[OPP_DIR(k)]->link[currentTime][n]*(/**/lattice[x].buddy[n]->buddy[OPP_DIR(k)]->link[currentTime][k]\
											*(lattice[x].buddy[n]->Tr231aU0nBar[i])\
											\
											+(lattice[x].buddy[n]->buddy[OPP_DIR(k)]->Tr231aU0nBar[i])\
											*lattice[x].buddy[n]->buddy[OPP_DIR(k)]->link[currentTime][k]/**/)*(~lattice[x].link[currentTime][n])\
											//
											)\
										   \
										   +(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])*(\
											//
											(lattice[x].buddy[OPP_DIR(i)]->Tr231aU0nBar[k])*lattice[x].getprod3(i,n)\
											\
											+lattice[x].buddy[OPP_DIR(i)]->link[currentTime][n]*(/**/lattice[x].buddy[n]->buddy[OPP_DIR(i)]->link[currentTime][i]\
											*(lattice[x].buddy[n]->Tr231aU0nBar[k])\
											\
											+(lattice[x].buddy[n]->buddy[OPP_DIR(i)]->Tr231aU0nBar[k])\
											*lattice[x].buddy[n]->buddy[OPP_DIR(i)]->link[currentTime][i]/**/)*(~lattice[x].link[currentTime][n])\
											//
											)\
										   );
								}
								else if(eabc(n,i,k)==-1.0){
									rhon-=(lattice[x].dU0iUjk1[currentTime][n]*lattice[x].efield[currentTime][n]-lattice[x].efield[prevTime][n]*lattice[x].dU0iUjk1[prevTime][n]\
										   \
										   +lattice[x].dU0iUjk2[currentTime][n]-lattice[x].dU0iUjk2[prevTime][n]
										   \
										   \
										   +lattice[x].getprod4(k,n)*(lattice[x].Tr231aU0nBar[i])\
										   \
										   -lattice[x].getprod4(i,n)*(lattice[x].Tr231aU0nBar[k])\
										   \
										   +lattice[x].link[currentTime][n]*(\
																			 //
																			 (lattice[x].buddy[n]->Tr231aU0nBar[i])*lattice[x].getprod1(n,k)\
																			 \
																			 +lattice[x].buddy[n]->link[currentTime][k]\
																			 *(lattice[x].buddy[n]->buddy[k]->Tr231aU0nBar[i])*lattice[x].getprod2(k,n)\
																			 \
																			 -(lattice[x].buddy[n]->Tr231aU0nBar[k])*lattice[x].getprod1(n,i)\
																			 \
																			 -lattice[x].buddy[n]->link[currentTime][i]\
																			 *(lattice[x].buddy[n]->buddy[i]->Tr231aU0nBar[k])*lattice[x].getprod2(i,n)\
																			 //
																			 )\
										   \
										   +lattice[x].link[currentTime][k]*(lattice[x].buddy[k]->Tr231aU0nBar[i])*lattice[x].getprod1(k,n)\
										   \
										   -lattice[x].link[currentTime][i]*(lattice[x].buddy[i]->Tr231aU0nBar[k])*lattice[x].getprod1(i,n)\
										   \
										   -(~lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k])*(\
																								   //
																								   (lattice[x].buddy[OPP_DIR(k)]->Tr231aU0nBar[i])*lattice[x].getprod3(k,n)\
																								   \
																								   +lattice[x].buddy[OPP_DIR(k)]->link[currentTime][n]*(/**/lattice[x].buddy[n]->buddy[OPP_DIR(k)]->link[currentTime][k]\
																																						*(lattice[x].buddy[n]->Tr231aU0nBar[i])\
																																						\
																																						+(lattice[x].buddy[n]->buddy[OPP_DIR(k)]->Tr231aU0nBar[i])\
																																						*lattice[x].buddy[n]->buddy[OPP_DIR(k)]->link[currentTime][k]/**/)*(~lattice[x].link[currentTime][n])\
																								   //
																								   )\
										   \
										   +(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])*(\
																								   //
																								   (lattice[x].buddy[OPP_DIR(i)]->Tr231aU0nBar[k])*lattice[x].getprod3(i,n)\
																								   \
																								   +lattice[x].buddy[OPP_DIR(i)]->link[currentTime][n]*(/**/lattice[x].buddy[n]->buddy[OPP_DIR(i)]->link[currentTime][i]\
																																						*(lattice[x].buddy[n]->Tr231aU0nBar[k])\
																																						\
																																						+(lattice[x].buddy[n]->buddy[OPP_DIR(i)]->Tr231aU0nBar[k])\
																																						*lattice[x].buddy[n]->buddy[OPP_DIR(i)]->link[currentTime][i]/**/)*(~lattice[x].link[currentTime][n])\
																								   //
																								   )\
										   );
								}
								 
								 /*For non-symmetrized plaquette
								rhon+=4.0*eabc(n,i,k)*( lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].getUmnBar(currentTime,i,k)*lattice[x].link[nextTime][n]*(~lattice[x].link[currentTime][n])\
												   \
												   -lattice[x].TrPhi2iSigma3Phi1a[prevTime]*lattice[x].link[currentTime][n]*(~lattice[x].link[prevTime][n])*lattice[x].getUmnBar(prevTime,i,k)\
												   \
												   \
												   +lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].getUmn(currentTime,n,k)*(lattice[x].U0nBar[i].UMinusUa())\
												   \
													-lattice[x].buddy[OPP_DIR(k)]->TrPhi2iSigma3Phi1a[currentTime]*(~lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k])\
													*(lattice[x].buddy[OPP_DIR(k)]->U0nBar[i].UMinusUa())*lattice[x].buddy[OPP_DIR(k)]->getUmn(currentTime,n,k)*lattice[x].buddy[OPP_DIR(k)]->link[currentTime][k]\
													\
												   );*/
							}
						}//k loop
					}
				}//i loop
				lattice[x].efield[nextTime][n]+=betaEk*rhon;lattice[x].efield[nextTime][n].setUnitary();
			}//endif
		}//n loop
	}
	
	//put trial solution (nextTime) back to currentTime
	for(int x=0;x<volume1;x++){
			for (int j=XUP; j<=ZUP; j++) {	
				//lattice[x].efield[currentTime][j]=lattice[x].efield[nextTime][j];
				lattice[x].link[nextTime][j]=lattice[x].efield[nextTime][j]*lattice[x].link[currentTime][j];//also links to nextTime	
			}
	}
	updateUmn(nextTime);updateUmnBar(nextTime);updateU0nBar();updateTr231aU0nBar();
}

// Solve implicit EOS

double getMaxError(){
	double maxErr=0,diff;
	for(int x=0;x<volume1;x++){
		for(int j=XUP;j<=ZUP;j++){
			diff=cmp(lattice[x].efield[nextTime][j],lattice[x].efield[currentTime][j]);
			if(diff>maxErr) maxErr=diff;
			lattice[x].efield[currentTime][j]=lattice[x].efield[nextTime][j];
		}
	}
	return maxErr;
}

/* Check Gauss constraint */

double calcGauss(){
	double maxErr=0,diff;
	RpTimesSU2 dE,dE1,dEk;
	int x;
	for(x=0;x<volume1;x++)
		//for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++)
    {
		//x=i*ny*nz+j*nz+k;
		dE.setZero();
		for (int n=XUP; n<=ZUP; n++) {
			dE+=lattice[x].efield[currentTime][n]\
			- (~lattice[x].buddy[OPP_DIR(n)]->link[currentTime][n])\
			*lattice[x].buddy[OPP_DIR(n)]->efield[currentTime][n]\
			*lattice[x].buddy[OPP_DIR(n)]->link[currentTime][n];
		}
		
		dE1=(betaGauss1*lattice[x].phi1[currentTime]*(~lattice[x].phi1[nextTime])\
			 +betaGauss2*lattice[x].phi2[currentTime]*(~lattice[x].phi2[nextTime]));
		if(kappa!=0){
			dEk.setZero();
			for(int i=XUP;i<=ZUP;i++){
				for(int j=XUP;j<=ZUP;j++){
					if(j!=i){
						for(int k=XUP;k<=ZUP;k++){
							if((k!=i)&&(k!=j)){
								
								dEk+=eabc(i,j,k)*(lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].efield[currentTime][i]*lattice[x].getUmnBar(currentTime,j,k)\
												  \
												  -lattice[x].TrPhi2iSigma3Phi1a[nextTime]*(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])\
												  *lattice[x].buddy[OPP_DIR(i)]->link[nextTime][i]*lattice[x].getUmnBar(nextTime,j,k)\
												  \
												  +lattice[x].TrPhi2iSigma3Phi1a[nextTime]*lattice[x].getUmnBar(nextTime,j,k)*lattice[x].efield[currentTime][i]\
												  \
												  -lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].getUmnBar(currentTime,j,k)\
												  *(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])*lattice[x].buddy[OPP_DIR(i)]->link[nextTime][i]\
												  \
												  +lattice[x].link[nextTime][i]*( lattice[x].buddy[i]->TrPhi2iSigma3Phi1a[currentTime]*lattice[x].buddy[i]->getUmnBar(currentTime,j,k)\
																				 +lattice[x].buddy[i]->TrPhi2iSigma3Phi1a[nextTime]*lattice[x].buddy[i]->getUmnBar(nextTime,j,k))*(~lattice[x].link[currentTime][i])\
												  \
												  -(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])*( lattice[x].buddy[OPP_DIR(i)]->TrPhi2iSigma3Phi1a[currentTime]*lattice[x].buddy[OPP_DIR(i)]->getUmnBar(currentTime,j,k)\
																										  +lattice[x].buddy[OPP_DIR(i)]->TrPhi2iSigma3Phi1a[nextTime]*lattice[x].buddy[OPP_DIR(i)]->getUmnBar(nextTime,j,k))*lattice[x].buddy[OPP_DIR(i)]->link[nextTime][i]\
												  );
								
								/*For non-symmetrized plaquette
								dEk+=4.0*eabc(i,j,k)*(lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].efield[currentTime][i]*lattice[x].getUmnBar(currentTime,j,k)\
												\
												-lattice[x].buddy[OPP_DIR(i)]->TrPhi2iSigma3Phi1a[currentTime]*(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])\
												*lattice[x].buddy[OPP_DIR(i)]->getUmnBar(currentTime,j,k)*lattice[x].buddy[OPP_DIR(i)]->link[nextTime][i]\
												  );*/
							}
						}
					}
				}
			}
			//cout << "\n";
			dEk=betaEk*dEk;
		}
		else {
			dEk.setZero();
		}
		
#ifdef  ANDERS
		maxErr+=relativeDifference(dE,dE1+dEk);
#else
		diff=cmp(dE,dE1+dEk);
		//if(diff>1.0)
		//	cout << t << ": " << i << ", " << j << ", " << k << " : " << diff << "\n";
		if(fabs(diff)>maxErr){
			maxErr=diff;
			//cout << i << ", " << j << ", " << k << ": " << dE << "/" << dE1 << " = ( " << dE.e[1]/dE1.e[1] << ", " << dE.e[2]/dE1.e[2] << ", " << dE.e[3]/dE1.e[3] << " )" << endl;
		}
#endif
	}
#ifdef  ANDERS
	maxErr/=volume;
#endif
	
	return maxErr;
}

void solveEOS(){
	int nIteration=0;
	double maxerr,gauss;
	
	// EField(U_{0j}) at currentTime
	linksAtNextTimeKappa0();

	phisAtNextTime();
 
	if(kappa!=0){
		setProds();
		do{
			linksAtNextTime();
			phisAtNextTime();
			
			if(nIteration>=NTRIAL){
				cout << "Not convergent...\n";exit(0);
			}
			nIteration++;
			//cout << nIteration << "\n";
			maxerr=getMaxError();
			//gauss=calcGauss();
			//cout << t << "   " << maxerr << "   " << gauss << "\n";
		}
		while((maxerr>ERROR));
	}
		
}

/* Calculate the Chern-Simons number */

void calculateNcs(){//20120217
	
	ncsdotPrev=ncsdot;ncsdot=0;
	
	for (int x=0; x<volume1; x++) {
		ncsdot-=lattice[x].ncs_loc;
		
		/*For non-symmetrized plaquette
		ncsdot-=(Tr( temporalPlaquetteBar(XUP,x).UMinusUa(),spatialPlaquetteBar(YUP, ZUP, x, currentTime) )\
				 + Tr( temporalPlaquetteBar(YUP,x).UMinusUa(),spatialPlaquetteBar(ZUP, XUP, x, currentTime) )\
				 + Tr( temporalPlaquetteBar(ZUP,x).UMinusUa(),spatialPlaquetteBar(XUP, YUP, x, currentTime) ));*/
	}
	ncs+=(ncsdotPrev + ncsdot)/(16.0*PI*PI);
}

/*********** Repositions for next iteration *********/

void reposition()
{
	int interm;
	
	/* Shift timelabel cyclicly. */
	
	interm = prevTime;
	prevTime = currentTime;
	currentTime = nextTime;
	nextTime = interm;
	
}

/* Calculate Phi_1^2 and Phi_2^2 */

#define NHH 40

double calcPhi1sq(){
	double phisq=0;
	for (int x=0; x<volume1; x++) {
		phisq+=lattice[x].TrPhi1sq[currentTime];
	}
	return phisq/(av1*av1*volume*lambda1);
}


double calcPhi1sqHisgram(string fname){
	unsigned int hisgram[NHH]; 
	double dh=4.0/NHH,h,pref=1.0/(av1*av1*lambda1),phisq;
		
	// Initialize all numbers into zero
	for(int i=0;i<NHH;i++) hisgram[i]=0;
	
	// count
	for(int x=0;x<volume1;x++){
		h=0.5*dh;
		phisq=pref*lattice[x].TrPhi1sq[currentTime];
		for(int i=0;i<NHH;i++){
			if( (phisq>= (h-0.5*dh))&&(phisq< (h+0.5*dh))){ 
				hisgram[i]+=1;break;
			}
			
			h+=dh;
		}
	}
	
	// output
	ofstream out;
	out.open(fname.c_str(),ios::app);
	h=0.5*dh;//int count=0;
	out << "# t =" << t << "\n";
	for(int i=0;i<NHH;i++){
		out << h << "   " << hisgram[i]  << "\n";
		//count+=hisgram[i];
		h+=dh;
	}
	//cout << count << endl;
	out << "\n\n";
	
}

double calcPhi2sq(){
	double phisq=0;
	for (int x=0; x<volume1; x++) {
		phisq+=lattice[x].TrPhi2sq[currentTime];
	}
	return phisq/(av2*av2*volume*lambda2);
}

double calcPhi2sqHisgram(string fname){
	unsigned int hisgram[NHH]; 
	double dh=4.0/NHH,h,pref=1.0/(av2*av2*lambda2),phisq;
	
	// Initialize all numbers into zero
	for(int i=0;i<NHH;i++) hisgram[i]=0;
	
	// count
	for(int x=0;x<volume1;x++){
		h=0.5*dh;
		phisq=pref*lattice[x].TrPhi2sq[currentTime];
		for(int i=0;i<NHH;i++){
			if( (phisq>= (h-0.5*dh))&&(phisq< (h+0.5*dh))){ 
				hisgram[i]+=1;break;
			}
			
			h+=dh;
		}
	}
	
	// output
	ofstream out;
	out.open(fname.c_str(),ios::app);
	h=0.5*dh;
	out << "# t =" << t << "\n";
	for(int i=0;i<NHH;i++){
		out << h << "   " << hisgram[i]  << "\n";
		h+=dh;
	}
	out << "\n\n";
	
}

double calcTheta(){// Relative phase between the two Higgses
	double theta=0,thetatmp;
	for(int x=0;x<volume1;x++){
		thetatmp=asin(lattice[x].TrPhi2iSigma3Phi1a[currentTime]/sqrt(lattice[x].TrPhi1sq[currentTime]*lattice[x].TrPhi2sq[currentTime]));
		theta+=(thetatmp);
	}
	return theta/(theta12*volume1);
}

double calcThetaHisgram(string fname){
	unsigned int hisgram[NHH]; 
	double dh=2.0*PI/NHH,h,thetatmp;
	
	// Initialize all numbers into zero
	for(int i=0;i<NHH;i++) hisgram[i]=0;
	
	// count
	for(int x=0;x<volume1;x++){
		h=-PI+0.5*dh;
		thetatmp=asin(lattice[x].TrPhi2iSigma3Phi1a[currentTime]/sqrt(lattice[x].TrPhi1sq[currentTime]*lattice[x].TrPhi2sq[currentTime]));
		for(int i=0;i<NHH;i++){
			if( (thetatmp>= (h-0.5*dh))&&(thetatmp< (h+0.5*dh))){ 
				hisgram[i]+=1;break;
			}
			h+=dh;
		}
	}
	
	// output
	ofstream out;
	out.open(fname.c_str(),ios::app);
	h=-PI+0.5*dh;
	out << "# t =" << t << "\n";
	for(int i=0;i<NHH;i++){
		out << h << "   " << hisgram[i]  << "\n";
		h+=dh;
	}
	out << "\n\n";
}

double calcImPhi1aPhi2Hisgram(string fname){
	unsigned int hisgram[NHH]; 
	double dh=4.0/NHH,h,thetatmp;
	
	// Initialize all numbers into zero
	for(int i=0;i<NHH;i++) hisgram[i]=0;
	
	// count
	for(int x=0;x<volume1;x++){
		h=-2.0+0.5*dh;
		thetatmp=0.5*lattice[x].TrPhi2iSigma3Phi1a[currentTime]/(sqrt(lambda1*lambda2)*av1*av2);
		for(int i=0;i<NHH;i++){
			if( (thetatmp>= (h-0.5*dh))&&(thetatmp< (h+0.5*dh))){ 
				hisgram[i]+=1;break;
			}
			h+=dh;
		}
	}
	
	// output
	ofstream out;
	out.open(fname.c_str(),ios::app);
	h=-2.0+0.5*dh;
	out << "# t =" << t << "\n";
	for(int i=0;i<NHH;i++){
		out << h << "   " << hisgram[i]  << "\n";
		h+=dh;
	}
	out << "\n\n";
}
