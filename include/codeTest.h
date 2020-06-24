// Code Testing Funcitons

void amplifyFieldAmplidues(double prefactor){
	double pref2=prefactor*prefactor;
	for(int x=0;x<volume1;x++){
		lattice[x].phi1[prevTime]=prefactor*lattice[x].phi1[prevTime];
		lattice[x].phi1[currentTime]=prefactor*lattice[x].phi1[currentTime];
		lattice[x].phi2[prevTime]=prefactor*lattice[x].phi2[prevTime];
		lattice[x].phi2[currentTime]=prefactor*lattice[x].phi2[currentTime];
		
		lattice[x].updatePhiExpression(prevTime);
		lattice[x].updatePhiExpression(currentTime);
		
		for(int dir=XUP;dir<=ZUP;dir++){
			lattice[x].efield[prevTime][dir]=pref2*lattice[x].efield[prevTime][dir];
			lattice[x].efield[prevTime][dir].setUnitary();lattice[x].efield[currentTime][dir]=lattice[x].efield[prevTime][dir];//Random guess for efield at currentTime = at prevTime
			lattice[x].link[currentTime][dir]=lattice[x].efield[prevTime][dir]*lattice[x].link[prevTime][dir];
		}
		
	}
	updateUmn(currentTime);updateUmnBar(currentTime);
	updateUmn(prevTime);updateUmnBar(prevTime);}

// Spatial inversion initial fields

void spatialInversion(){
	
 	int i,j,k; 
	int x,y,z; 
	site *myself;	
	
	//Active point of view
   	
	// Allocate space for lattice
	
	site *latticeP = new site[volume1];
	if(latticeP==NULL)
	{
		printf("No room for lattice\n");
		exit(1);
	}
	
	//cout << "latticeP is allocated!" << endl;
	
	// Find neighbors
	
	for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++)
	{
	    myself=latticeP+i*ny*nz+j*nz+k;
	    myself->buddy[ZUP]=latticeP+i*ny*nz+j*nz+(k+1)%nz;
	    myself->buddy[YUP]=latticeP+i*ny*nz+((j+1)%ny)*nz+k;
	    myself->buddy[XUP]=latticeP+((i+1)%nx)*ny*nz+j*nz+k;
	    (myself->buddy[ZUP])->buddy[ZDOWN]=myself;
	    (myself->buddy[YUP])->buddy[YDOWN]=myself;
	    (myself->buddy[XUP])->buddy[XDOWN]=myself;
	}
	//cout << "Buddy is done!" << endl;	
	
	for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++){
	    x = i*ny*nz+j*nz+k;
	    y = ((nx-i)%nx)*ny*nz+((ny-j)%ny)*nz+(nz-k)%nz;
	    //cout << "( " << i << ", " << j << ", " << k << " ) vs ( " << (nx-i)%nx << ", " << (ny-j)%ny << ", " << (nz-k)%nz << ")" << endl;
	    for(int tInt=prevTime;tInt<=nextTime;tInt++){
			latticeP[x].phi1[tInt]=lattice[y].phi1[tInt];latticeP[x].phi2[tInt]=lattice[y].phi2[tInt];
			latticeP[x].phi1iSigma3[tInt]=lattice[y].phi1iSigma3[tInt];latticeP[x].phi2iSigma3[tInt]=lattice[y].phi2iSigma3[tInt];
			latticeP[x].TrPhi1sq[tInt]=lattice[y].TrPhi1sq[tInt];latticeP[x].TrPhi2sq[tInt]=lattice[y].TrPhi2sq[tInt];
			latticeP[x].TrPhi2Phi1a[tInt]=lattice[y].TrPhi2Phi1a[tInt];latticeP[x].TrPhi2iSigma3Phi1a[tInt]=lattice[y].TrPhi2iSigma3Phi1a[tInt];
	    }			
		for(z=XUP;z<=ZUP;z++){
			latticeP[x].link[prevTime][z]=~(lattice[y].buddy[OPP_DIR(z)]->link[prevTime][z]);
			latticeP[x].efield[prevTime][z]=~(lattice[y].buddy[OPP_DIR(z)]->efield[prevTime][z]);
			latticeP[x].link[currentTime][z]=latticeP[x].efield[prevTime][z]*latticeP[x].link[prevTime][z];
		}
	}
	delete[] lattice;
	lattice = latticeP;
/*
	
	// Passive point of view	
	
	// Reverse the coordinate
	for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++)
	{
	    myself=lattice+i*ny*nz+j*nz+k;
	    myself->buddy[ZDOWN]=lattice+i*ny*nz+j*nz+(k+1)%nz;
	    myself->buddy[YDOWN]=lattice+i*ny*nz+((j+1)%ny)*nz+k;
	    myself->buddy[XDOWN]=lattice+((i+1)%nx)*ny*nz+j*nz+k;
	    (myself->buddy[ZDOWN])->buddy[ZUP]=myself;
	    (myself->buddy[YDOWN])->buddy[YUP]=myself;
	    (myself->buddy[XDOWN])->buddy[XUP]=myself;
	}

	// A_i -> -A_i
	for(x=0;x<volume1;x++){
		for(z=XUP;z<=ZUP;z++){
			lattice[x].link[prevTime][z]=~lattice[x].link[prevTime][z];
			lattice[x].efield[prevTime][z]=~lattice[x].efield[prevTime][z];
			lattice[x].link[currentTime][z]=lattice[x].efield[prevTime][z]*lattice[x].link[prevTime][z];
		}
	}
*/	
	updateUmn(currentTime);updateUmnBar(currentTime);
	updateUmn(prevTime);updateUmnBar(prevTime);	
}

void chargeConjugate(){
	for(int x=0;x<volume1;x++){
		for(int tInt=prevTime;tInt<=nextTime;tInt++){
			lattice[x].phi1[tInt]=lattice[x].phi1[tInt].chargeConjugate();lattice[x].phi2[tInt]=lattice[x].phi2[tInt].chargeConjugate();
		}
		for(int z=XUP;z<=ZUP;z++){
			lattice[x].link[prevTime][z]=lattice[x].link[prevTime][z].chargeConjugate();
			lattice[x].efield[prevTime][z]=lattice[x].efield[prevTime][z].chargeConjugate();
			lattice[x].link[currentTime][z]=lattice[x].efield[prevTime][z]*lattice[x].link[prevTime][z];
		}
		lattice[x].updatePhiExpression(prevTime);
		lattice[x].updatePhiExpression(currentTime);
	}
	updateUmn(currentTime);updateUmnBar(currentTime);
	updateUmn(prevTime);updateUmnBar(prevTime);
}

// Evolution of Higgs Fields without gauge fields

void turnOffGaugeField(){
	for (int x=0; x<volume1; x++ ) {
		for (int j=XUP; j<=ZUP; j++) {
			lattice[x].link[prevTime][j]=RpTimesSU2(1,0,0,0);lattice[x].link[currentTime][j]=RpTimesSU2(1,0,0,0);lattice[x].link[nextTime][j]=RpTimesSU2(1,0,0,0);
			lattice[x].efield[prevTime][j]=RpTimesSU2(1,0,0,0);lattice[x].efield[currentTime][j]=RpTimesSU2(1,0,0,0);lattice[x].efield[nextTime][j]=RpTimesSU2(1,0,0,0);
		}
	}
}

// Evolution of Gauge field without Higgses

void turnOffHiggses(){
	double	pref1=10.0, pref2=10.0, pref3=10.0,
	f1[volume1], f2[volume1], f3[volume1];
	ifstream in;
	in.open("include/initialGauge.txt");
	int x;

	for(x=0;x<volume1;x++){
		in >> f1[x];
		cout << "f1[" << x << "] = " << "   " << f1[x] << endl;
	}
	for(x=0;x<volume1;x++){
		in >> f2[x];
		//cout << "f2[" << x << "] = " << "   " << f2[x] << endl;
	}
	for(x=0;x<volume1;x++){
		in >> f3[x];
		//cout << "f3[" << x << "] = " << "   " << f3[x] << endl;
	}
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				x=i*ny*nz+j*nz+k;
				int xp=((i+1)%nx)*ny*nz+j*nz+k;
				int yp=i*ny*nz+((j+1)%ny)*nz+k;
				int zp=i*ny*nz+j*nz+((k+1)%nz);
				lattice[x].efield[prevTime][0]=RpTimesSU2(0,0.5,0,0);lattice[x].efield[prevTime][0].setUnitary();
				lattice[x].efield[prevTime][1]=RpTimesSU2(0,0.5,0,0);lattice[x].efield[prevTime][1].setUnitary();
				lattice[x].efield[prevTime][2]=RpTimesSU2(0,0.5,0,0);lattice[x].efield[prevTime][2].setUnitary();
				//cout << ((1.0-0.5*(Tr(lattice[x].efield[prevTime][XUP])))+2.0-0.5*(Tr(lattice[x].efield[prevTime][YUP])+Tr(lattice[x].efield[prevTime][ZUP]))) << endl;
				for(int n=XUP;n<=ZUP;n++){ 
					lattice[x].link[currentTime][n]=lattice[x].efield[prevTime][n]*lattice[x].link[prevTime][n];
				}
				
				lattice[x].phi1[prevTime]=RpTimesSU2(0,0,0,0);lattice[x].phi1[currentTime]=RpTimesSU2(0,0,0,0);
				lattice[x].phi2[prevTime]=RpTimesSU2(0,0,0,0);lattice[x].phi2[currentTime]=RpTimesSU2(0,0,0,0);
				
				lattice[x].updatePhiExpression(currentTime);
			}
		}
	}
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				x=i*ny*nz+j*nz+k;	
		cout << i << "   " << j << "   " << k << "   " << (lattice[x].efield[prevTime][XUP].e[1]+lattice[x].efield[prevTime][YUP].e[1]+lattice[x].efield[prevTime][ZUP].e[1])/\
		(lattice[x].buddy[OPP_DIR(XUP)]->efield[prevTime][XUP].e[1]+lattice[x].buddy[OPP_DIR(YUP)]->efield[prevTime][YUP].e[1]\
		 +lattice[x].buddy[OPP_DIR(ZUP)]->efield[prevTime][ZUP].e[1]) << endl;
			}
		}
	}
	updateUmn(prevTime);
	updateUmn(currentTime);
	
}

/* Calculate links at the next time step */
void linksAtNextTimeOffHiggses(){//kappa = 0
	for (int x=0; x<volume1; x++ ) {
		for (int j=XUP; j<=ZUP; j++) {
			lattice[x].efield[currentTime][j]=lattice[x].efield[prevTime][j];
			for (int n=XUP; n<=ZUP; n++) {
				if (n!=j) {
					lattice[x].efield[currentTime][j]-=deltatsq*(lattice[x].getUmn(currentTime,j,n)\
															  -(~lattice[x].buddy[OPP_DIR(n)]->link[currentTime][n])\
															  *lattice[x].buddy[OPP_DIR(n)]->getUmn(currentTime,j,n)\
															  *lattice[x].buddy[OPP_DIR(n)]->link[currentTime][n]);
				}
			}
			lattice[x].efield[currentTime][j].setUnitary();
			lattice[x].link[nextTime][j]=lattice[x].efield[currentTime][j]*lattice[x].link[currentTime][j];
		}
	}
	updateUmn(nextTime);
}

// Kinetic Energy of Gauge Fields

double calcKineticEnergyA(){
	double ke=0;
	for(int x=0;x<volume1;x++){
		ke+=betagt*(3.0-0.5*(Tr(lattice[x].efield[currentTime][XUP])+Tr(lattice[x].efield[currentTime][YUP])+Tr(lattice[x].efield[currentTime][ZUP])));
	}
	return ke*vev/(avev*deltat);//in GeV
}

// Kinetic Energy of Higgs Fields

double calcKineticEnergyH(){
	double ke=0;
	for(int x=0;x<volume1;x++){
		ke+=betah1t*Trsq(lattice[x].phi1[nextTime]-lattice[x].phi1[currentTime])+betah2t*Trsq(lattice[x].phi2[nextTime]-lattice[x].phi2[currentTime]);
	}
	return ke*vev/(avev*deltat);//in GeV
}

// Total Kinetic Energy

double calcKineticEnergy(){
	double ke=0;
	for(int x=0;x<volume1;x++){
		ke+=betagt*(3.0-0.5*(Tr(lattice[x].efield[currentTime][XUP])+Tr(lattice[x].efield[currentTime][YUP])+Tr(lattice[x].efield[currentTime][ZUP])))\
		+betah1t*Trsq(lattice[x].phi1[nextTime]-lattice[x].phi1[currentTime])+betah2t*Trsq(lattice[x].phi2[nextTime]-lattice[x].phi2[currentTime]);
	}
	return ke*vev/(avev*deltat);//in GeV
}

double TrUmn(int m, int n, int x, int tInt){
	return TrAdaggerB(lattice[x].link[tInt][n]*lattice[x].buddy[n]->link[tInt][m],lattice[x].link[tInt][m]*lattice[x].buddy[m]->link[tInt][n]);
}

// Potential of energy from covariant spatial derivatives

double calcPotentialEnergyK(){
	double pe=0;
	
	for(int x=0;x<volume1;x++){
		for(int n=XUP;n<=ZUP;n++){
			pe+=betah1x*Trsq(lattice[x].link[currentTime][n]*lattice[x].buddy[n]->phi1[currentTime] - lattice[x].phi1[currentTime])\
			+betah2x*Trsq(lattice[x].link[currentTime][n]*lattice[x].buddy[n]->phi2[currentTime] - lattice[x].phi2[currentTime]);
		}
	}
	return pe*vev/(avev*deltat);//in GeV
}

// Potential energy of gauge fields

double calcPotentialEnergyA(){
	double pe=0;
	
	for(int x=0;x<volume1;x++){
		pe+=betagx*(3.0-0.5*(TrUmn(XUP,YUP,x,currentTime)+TrUmn(YUP,ZUP,x,currentTime)+TrUmn(ZUP,XUP,x,currentTime)));
	}
	return pe*vev/(avev*deltat);//in GeV
}

// Potential energy of higgs fields

double calcPotentialEnergyH(){
	double pe=0;
	
	for(int x=0;x<volume1;x++){
		pe+=beta1r*(lattice[x].TrPhi1sq[currentTime]-am11sq)*(lattice[x].TrPhi1sq[currentTime]-am11sq)\
		+beta2r*(lattice[x].TrPhi2sq[currentTime]-am22sq)*(lattice[x].TrPhi2sq[currentTime]-am22sq)\
		+beta3*lattice[x].TrPhi1sq[currentTime]*lattice[x].TrPhi2sq[currentTime]\
		+(beta4+beta5R)*lattice[x].TrPhi2Phi1a[currentTime]*lattice[x].TrPhi2Phi1a[currentTime]\
		+(beta4-beta5R)*lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].TrPhi2iSigma3Phi1a[currentTime]\
		+beta5I*lattice[x].TrPhi2Phi1a[currentTime]*lattice[x].TrPhi2iSigma3Phi1a[currentTime]\
		-beta12R*lattice[x].TrPhi2Phi1a[currentTime]-beta12I*lattice[x].TrPhi2iSigma3Phi1a[currentTime];
	}
	return pe*vev/(avev*deltat);//in GeV
}

// Total potential energy

double calcPotentialEnergy(){
	double pe=0;
	
	for(int x=0;x<volume1;x++){
		for(int n=XUP;n<=ZUP;n++){
			pe+=betah1x*Trsq(lattice[x].link[currentTime][n]*lattice[x].buddy[n]->phi1[currentTime] - lattice[x].phi1[currentTime])\
			+betah2x*Trsq(lattice[x].link[currentTime][n]*lattice[x].buddy[n]->phi2[currentTime] - lattice[x].phi2[currentTime]);
		}
		pe+=betagx*(3.0-0.5*(TrUmn(XUP,YUP,x,currentTime)+TrUmn(YUP,ZUP,x,currentTime)+TrUmn(ZUP,XUP,x,currentTime)))\
		+beta1r*(lattice[x].TrPhi1sq[currentTime]-am11sq)*(lattice[x].TrPhi1sq[currentTime]-am11sq)\
		+beta2r*(lattice[x].TrPhi2sq[currentTime]-am22sq)*(lattice[x].TrPhi2sq[currentTime]-am22sq)\
		+beta3*lattice[x].TrPhi1sq[currentTime]*lattice[x].TrPhi2sq[currentTime]\
		+(beta4+beta5R)*lattice[x].TrPhi2Phi1a[currentTime]*lattice[x].TrPhi2Phi1a[currentTime]\
		+(beta4-beta5R)*lattice[x].TrPhi2iSigma3Phi1a[currentTime]*lattice[x].TrPhi2iSigma3Phi1a[currentTime]\
		+beta5I*lattice[x].TrPhi2Phi1a[currentTime]*lattice[x].TrPhi2iSigma3Phi1a[currentTime]\
		-beta12R*lattice[x].TrPhi2Phi1a[currentTime]-beta12I*lattice[x].TrPhi2iSigma3Phi1a[currentTime];
	}
	return pe*vev/(avev*deltat);//in GeV
}

