/*
 *  higgsWindingNumber.h
 *  
 *
 *  Created by Bin Wu on 1/12/12.
 *  Copyright 2012 Universitaet Bielefeld. All rights reserved.
 *
 */

double gaugeWindingNumber(){
	double gwn=0;
	
	for(int x=0;x<volume1;x++){
		gwn+=Tr(lattice[x].link[currentTime][XUP],(lattice[x].link[currentTime][YUP]*lattice[x].link[currentTime][ZUP] -lattice[x].link[currentTime][ZUP]*lattice[x].link[currentTime][YUP] ));
		//Tr((lattice[x].link[currentTime][XUP]-(~lattice[x].link[currentTime][XUP]))\
				*((lattice[x].link[currentTime][YUP]-(~lattice[x].link[currentTime][YUP]))*(lattice[x].link[currentTime][ZUP]-(~lattice[x].link[currentTime][ZUP]))\
				  -(lattice[x].link[currentTime][ZUP]-(~lattice[x].link[currentTime][ZUP]))*(lattice[x].link[currentTime][YUP]-(~lattice[x].link[currentTime][YUP]))));
	}
	
	return gwn/(8.0*PI*PI);
}


double gaugeWindingNumberAF(){
	double gwn=0;
	
	for(int x=0;x<volume1;x++){
		gwn+=Tr(lattice[x].link[currentTime][XUP],lattice[x].getUmn(currentTime,YUP,ZUP).UMinusUa())\
				+Tr(lattice[x].link[currentTime][YUP],lattice[x].getUmn(currentTime,ZUP,XUP).UMinusUa())\
					+Tr(lattice[x].link[currentTime][ZUP],lattice[x].getUmn(currentTime,XUP,YUP).UMinusUa());
		//Tr((lattice[x].link[currentTime][XUP]-(~lattice[x].link[currentTime][XUP]))\
		*((lattice[x].link[currentTime][YUP]-(~lattice[x].link[currentTime][YUP]))*(lattice[x].link[currentTime][ZUP]-(~lattice[x].link[currentTime][ZUP]))\
		-(lattice[x].link[currentTime][ZUP]-(~lattice[x].link[currentTime][ZUP]))*(lattice[x].link[currentTime][YUP]-(~lattice[x].link[currentTime][YUP]))));
	}
	
	return -gwn/(16.0*PI*PI);
}

double higgs1WindingNumber(){
	RpTimesSU2 dPhidiPhia[NUMSDIR];
	double hwn=0;
	
	for(int x=0;x<volume1;x++)lattice[x].phiUni=lattice[x].phi1[currentTime]/lattice[x].phi1[currentTime].getNorm();
	
	for(int x=0;x<volume1;x++){
		for(int j=XUP;j<=ZUP;j++){
			dPhidiPhia[j]= 0.5*(lattice[x].buddy[j]->phiUni - lattice[x].buddy[OPP_DIR(j)]->phiUni)*(~lattice[x].phiUni);
		}
		for(int i=XUP;i<=ZUP;i++){
			for(int j=XUP;j<=ZUP;j++){
				if(j!=i){
					for(int k=XUP;k<=ZUP;k++){
						if(k!=i&&k!=j)
							hwn-=eabc(i,j,k)*Tr(dPhidiPhia[i],dPhidiPhia[j]*dPhidiPhia[k]);
					}
				}
			}
		}
	}
	return hwn/(24.0*PI*PI);
}

double higgs2WindingNumber(){
	RpTimesSU2 dPhidiPhia[NUMSDIR];
	double hwn=0;
	
	for(int x=0;x<volume1;x++)lattice[x].phiUni=lattice[x].phi2[currentTime]/lattice[x].phi2[currentTime].getNorm();
	
	for(int x=0;x<volume1;x++){
		for(int j=XUP;j<=ZUP;j++){
			dPhidiPhia[j]= 0.5*(lattice[x].buddy[j]->phiUni - lattice[x].buddy[OPP_DIR(j)]->phiUni)*(~lattice[x].phiUni);
		}
		for(int i=XUP;i<=ZUP;i++){
			for(int j=XUP;j<=ZUP;j++){
				if(j!=i){
					for(int k=XUP;k<=ZUP;k++){
						if(k!=i&&k!=j)
							hwn-=eabc(i,j,k)*Tr(dPhidiPhia[i],dPhidiPhia[j]*dPhidiPhia[k]);
					}
				}
			}
		}
	}
	return hwn/(24.0*PI*PI);
}


void outputNH1(){
	RpTimesSU2 dPhidiPhia[NUMSDIR];
	double afm = 0.1975*avev/vev;//a in fm
	
	for(int x=0;x<volume1;x++)lattice[x].phiUni=lattice[x].phi1[currentTime]/lattice[x].phi1[currentTime].getNorm();
	
	for(int z=0;z<nz;z++){
		cout << "# z = " << afm*z << "\n";
		for(int x=0;x<nx;x++){
			for(int y=0;y<ny;y++){
				double hwn=0;
				int i = x*ny*nz+y*nz+z;
				for(int j=XUP;j<=ZUP;j++){
					dPhidiPhia[j]= 0.5*(lattice[i].buddy[j]->phiUni - lattice[i].buddy[OPP_DIR(j)]->phiUni)*(~lattice[i].phiUni);
				}
				for(int l=XUP;l<=ZUP;l++){
					for(int m=XUP;m<=ZUP;m++){
						if(m!=l){
							for(int n=XUP;n<=ZUP;n++){
								if(n!=l&&n!=m)
									hwn-=eabc(l,m,n)*Tr(dPhidiPhia[l],dPhidiPhia[m]*dPhidiPhia[n]);
							}
						}
					}
				}
				cout << afm*x << "   " << afm*y << "   " << hwn/(24.0*PI*PI) << "\n";
			}
		}
		cout << "\n\n";
	}
}

void outputPhi2(){
	ostringstream ostr;ofstream out;
	double afm = 0.1975*avev/vev;//a in fm
	
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
	
#ifdef _P
#ifdef _C
	ostr << "CPH.dat";
#else
	ostr << "PH.dat";
#endif
#else
#ifdef _C
	ostr << "CH.dat";	
#else
	ostr << "H.dat";
#endif
#endif
	out.open(ostr.str().c_str());
	
	for(int x=0;x<volume1;x++)lattice[x].phiUni=lattice[x].phi2[currentTime]/lattice[x].phi2[currentTime].getNorm();
	
	for(int z=0;z<nz;z++){
		out << "# z = " << afm*z << "\n";
		for(int x=0;x<nx;x++){
			for(int y=0;y<ny;y++){
				int i = x*ny*nz+y*nz+z;
				out << afm*x << "   " << afm*y << "   " << lattice[i].phiUni.e[0] << "   " << lattice[i].phiUni.e[1] << "   "<< lattice[i].phiUni.e[2] << "   "<< lattice[i].phiUni.e[3] << "\n";
			}
		}
		out << "\n\n";
	}
}

void outputNH2(){
	ostringstream ostr;ofstream out;
	RpTimesSU2 dPhidiPhia[NUMSDIR];
	double afm = 0.1975*avev/vev;//a in fm
	double hwntot=0;
	
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

#ifdef _P
#ifdef _C
	ostr << "CPH.dat";
#else
	ostr << "PH.dat";
#endif
#else
#ifdef _C
	ostr << "CH.dat";	
#else
	ostr << "H.dat";
#endif
#endif
	out.open(ostr.str().c_str());
	
	for(int x=0;x<volume1;x++)lattice[x].phiUni=lattice[x].phi2[currentTime]/lattice[x].phi2[currentTime].getNorm();
	
	for(int z=0;z<nz;z++){
		out << "# z = " << afm*z << "\n";
		for(int x=0;x<nx;x++){
			for(int y=0;y<ny;y++){
				double hwn=0;
				int i = x*ny*nz+y*nz+z;
				for(int j=XUP;j<=ZUP;j++){
					dPhidiPhia[j]= 0.5*(lattice[i].buddy[j]->phiUni - lattice[i].buddy[OPP_DIR(j)]->phiUni)*(~lattice[i].phiUni);
				}
				for(int l=XUP;l<=ZUP;l++){
					for(int m=XUP;m<=ZUP;m++){
						if(m!=l){
							for(int n=XUP;n<=ZUP;n++){
								if(n!=l&&n!=m)
									hwn-=eabc(l,m,n)*Tr(dPhidiPhia[l],dPhidiPhia[m]*dPhidiPhia[n]);
							}
						}
					}
				}
				hwn=hwn/(24.0*PI*PI);hwntot+=hwn;
				out << afm*x << "   " << afm*y << "   " << hwn << "\n";
			}
		}
		out << "\n\n";
	}
	out << "NH2 = " << hwntot << "\n";
}

/*
#define NHH2 40

void hisgramh2(string fname){
	unsigned int hisgramH2[4][NHH2]; double dh=2.0/NHH2,h;
	
	// Find SU(2) matrix
	for(int x=0;x<volume1;x++)lattice[x].phiUni=lattice[x].phi1[currentTime]/lattice[x].phi1[currentTime].getNorm();
	
	// Initialize all numbers into zero
	for(int n=0;n<NHH2;n++){
		hisgramH2[0][n]=0;hisgramH2[1][n]=0;hisgramH2[2][n]=0;hisgramH2[3][n]=0;
	}
	
	// count
	for(int x=0;x<volume1;x++){
		h=-1.0+0.5*dh;
		for(int i=0;i<NHH2;i++){
			for(int a=0;a<4;a++){
				if( (lattice[x].phiUni.e[a]>= (h-0.5*dh))&&(lattice[x].phiUni.e[a]< (h+0.5*dh))){
					hisgramH2[a][i]+=1;
				}
			}
			//cout << lattice[currentTime][x] << " " << tInterval << "   " << (tInterval-0.1*a) << "   " << (tInterval+0.1*a) << "   " << numInt[i] << "\n";
			h+=dh;
		}
	}
	
	// output
	ofstream out;
	out.open(fname.c_str(),ios::app);
	h=-1.0+0.5*dh;
	out << "# t =" << t << "\n";
	for(int i=0;i<NHH2;i++){
		out << h << "   " << hisgramH2[0][i]  << "   " << hisgramH2[1][i]  << "   " << hisgramH2[2][i]  << "   " << hisgramH2[3][i]  << "\n";
		h+=dh;
	}
	out << "\n\n";
}
*/