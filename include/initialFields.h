// U_{ij,x}

RpTimesSU2 spatialPlaquette(int diri, int dirj, int xInt, int timeInt){
	return lattice[xInt].link[timeInt][diri]*lattice[xInt].buddy[diri]->link[timeInt][dirj]\
	*(~lattice[xInt].buddy[dirj]->link[timeInt][diri])*(~lattice[xInt].link[timeInt][dirj]);
}

/* Calculate symmetrized U_{i j, x} */

RpTimesSU2 spatialPlaquetteBar(int i, int j, int x, int timeInt){
	return 0.25*(spatialPlaquette(i,j,x,timeInt)\
				 +lattice[x].link[timeInt][j]*(~lattice[x].buddy[OPP_DIR(i)]->buddy[j]->link[timeInt][i])\
				 *(~lattice[x].buddy[OPP_DIR(i)]->link[timeInt][j])*lattice[x].buddy[OPP_DIR(i)]->link[timeInt][i]\
				 +(~lattice[x].buddy[OPP_DIR(j)]->link[timeInt][j])*lattice[x].buddy[OPP_DIR(j)]->link[timeInt][i]\
				 *lattice[x].buddy[OPP_DIR(j)]->buddy[i]->link[timeInt][j]*(~lattice[x].link[timeInt][i])\
				 +(~lattice[x].buddy[OPP_DIR(i)]->link[timeInt][i])*(~lattice[x].buddy[OPP_DIR(i)]->buddy[OPP_DIR(j)]->link[timeInt][j])\
				 *lattice[x].buddy[OPP_DIR(i)]->buddy[OPP_DIR(j)]->link[timeInt][i]*lattice[x].buddy[OPP_DIR(j)]->link[timeInt][j]);
}

RpTimesSU2 spatialBackwardPlaquette(int diri, int dirj, int xInt, int timeInt){
	return (~lattice[xInt].buddy[OPP_DIR(dirj)]->link[timeInt][dirj])*(lattice[xInt].buddy[OPP_DIR(dirj)]->link[timeInt][diri])\
	*(lattice[xInt].buddy[diri]->buddy[OPP_DIR(dirj)]->link[timeInt][dirj])*(~lattice[xInt].link[timeInt][diri]);
	
	//lattice[xInt].link[timeInt][dirj]*(~lattice[xInt].buddy[OPP_DIR(diri)]->buddy[dirj]->link[timeInt][diri])\
	*(~lattice[xInt].buddy[OPP_DIR(diri)]->link[timeInt][dirj])*(lattice[xInt].buddy[OPP_DIR(diri)]->link[timeInt][diri]);
}

RpTimesSU2 temporalPlaquetteBar(int i, int x){
	return 0.25*( lattice[x].link[nextTime][i]*(~lattice[x].link[currentTime][i])\
				 +lattice[x].link[currentTime][i]*(~lattice[x].link[prevTime][i])\
				 +(~lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i])*lattice[x].buddy[OPP_DIR(i)]->link[nextTime][i]\
				 +(~lattice[x].buddy[OPP_DIR(i)]->link[prevTime][i])*lattice[x].buddy[OPP_DIR(i)]->link[currentTime][i]);
}

void setUmn(int m, int n, int x, int iTime){
	lattice[x].setUmn(iTime,m,n,spatialPlaquette(m, n, x, iTime));
}

void updateUmn(int iTime){
	for(int x=0;x<volume1;x++){
		setUmn(XUP, YUP, x, iTime);setUmn(YUP, ZUP, x, iTime);setUmn(ZUP, XUP, x, iTime);
	}
}

void setUmnBar(int m, int n, int x, int iTime){
	lattice[x].setUmnBar(iTime,m,n,spatialPlaquetteBar(m, n, x, iTime));
	
	//For non-symmetrized plaquette
	//lattice[x].setUmnBar(iTime,m,n,spatialPlaquette(m, n, x, iTime));
}

void updateUmnBar(int iTime){
	for(int x=0;x<volume1;x++){
		setUmnBar(XUP, YUP, x, iTime);setUmnBar(YUP, ZUP, x, iTime);setUmnBar(ZUP, XUP, x, iTime);
	}
}

void updateTr231aUmnBar(int iTime){
	for(int x=0;x<volume1;x++){
		lattice[x].updateTr231aUmnBar(iTime);
	}
}

void updateU0nBar(){
	for(int x=0;x<volume1;x++){
		for(int n=XUP;n<=ZUP;n++){
			lattice[x].U0nBar[n]=temporalPlaquetteBar(n, x);
			
			//For non-symmetrized plaquette
			//lattice[x].U0nBar[n]=lattice[x].efield[currentTime][n];
		}
	}
}
void updatedU0iUjk(){
	for(int x=0;x<volume1;x++){
		for(int n=XUP;n<=ZUP;n++){
			for (int i=XUP; i<=ZUP; i++) {
				if(i!=n){
					for (int k=i+1; k<=ZUP; k++) {
						if((k!=n)){
							lattice[x].dU0iUjk1[prevTime][n]=((lattice[x].getTr231aUmnBar(currentTime,i,k)+lattice[x].getTr231aUmnBar(prevTime,i,k) ).UMinusUa());
							lattice[x].dU0iUjk2[prevTime][n]=lattice[x].link[currentTime][n]*(( lattice[x].buddy[n]->getTr231aUmnBar(currentTime,i,k)+lattice[x].buddy[n]->getTr231aUmnBar(prevTime,i,k)).UMinusUa())*(~lattice[x].link[prevTime][n]);
						}
					}
				}
			}
		}
	}
}

void updateTr231aU0nBar(){
	for(int x=0;x<volume1;x++){
		lattice[x].updateTr231aU0nBar(currentTime);
	}
}
void reset(double kappadb){
	prevTime = 0;currentTime = 1;nextTime = 2;	t=0;ncs=0;ncsdot=0;
	kappa=kappadb;betaEk=0.5*gg*kappa*deltat/(av1*av2*sqrt(lambda1*lambda2));
	betak1=kappa*deltat*sqrt(lambda1/lambda2)/(av1*av2);betak2=kappa*deltat*sqrt(lambda2/lambda1)/(av1*av2);
}

//update parameters used in solving EOS
void updateLatticePara(){
	
	/* Convert parameters to lattice parameters and back again. */
	deltatsq=deltat*deltat;
	
	betagt=4.0/(gg*deltat); betagx=4.0*deltat/gg; 
	betah1t=0.5/(lambda1*deltat); betah1x=0.5*deltat/lambda1; 
	betah2t=0.5/(lambda2*deltat); betah2x=0.5*deltat/lambda2; 
	beta1r=0.25*betah1x; beta2r=0.25*betah2x; 
	beta3=0.25*lambda3*deltat/(lambda1*lambda2); beta4=0.25*lambda4*deltat/(lambda1*lambda2);
	beta5R=0.25*lambda5R*deltat/(lambda1*lambda2);beta5I=0.5*lambda5I*deltat/(lambda1*lambda2);
	beta12R=0.5*am12sqR*deltat/sqrt(lambda1*lambda2);beta12I=0.5*am12sqI*deltat/sqrt(lambda1*lambda2);
	
	// parameters for Higgs time marching
	beta31=0.5*lambda3*deltatsq/lambda2;beta32=0.5*lambda3*deltatsq/lambda1;
	beta21a1=0.5*deltatsq*(lambda4+lambda5R)/lambda2;beta21a2=0.5*deltatsq*(lambda4+lambda5R)/lambda1;
	beta231a1=0.5*deltatsq*(lambda5R-lambda4)/lambda2;beta231a2=0.5*deltatsq*(lambda5R-lambda4)/lambda1;
	betak1=kappa*deltat*sqrt(lambda1/lambda2)/(av1*av2);betak2=kappa*deltat*sqrt(lambda2/lambda1)/(av1*av2);
	
	//beta41=0.5*lambda4*deltatsq/lambda2;beta42=0.5*lambda4*deltatsq/lambda1;
	//beta51R=0.5*lambda5R*deltatsq/lambda2;beta52R=0.5*lambda5R*deltatsq/lambda1;
	beta51I=0.5*lambda5I*deltatsq/lambda2;beta52I=0.5*lambda5I*deltatsq/lambda1;
	beta121R=0.5*am12sqR*deltatsq*sqrt(lambda1/lambda2);beta122R=0.5*am12sqR*deltatsq*sqrt(lambda2/lambda1);
	beta121I=0.5*am12sqI*deltatsq*sqrt(lambda1/lambda2);beta122I=0.5*am12sqI*deltatsq*sqrt(lambda2/lambda1);
	
	betaE1=0.5*gg*deltatsq/lambda1;betaE2=0.5*gg*deltatsq/lambda2;betaEk=0.125*gg*kappa*deltat/(av1*av2*sqrt(lambda1*lambda2));	//20120214// parameters for E-field time marching, additional 1/2 because of the different definition of E-field
	
	betaGauss1=0.5*gg/lambda1;betaGauss2=0.5*gg/lambda2;			// parameters for checking Gauss constraints
	
}

/***********The inputfiles are opened and read*********/

void loadParameters(bool arg=FALSE,double arckappa=0)
{
	FILE *ifile;
	double aInvGeV;
	
	if((ifile = fopen(LATTICEPARAMETERS,"r"))==NULL) 
    {    
		printf("Cannot open ");printf(LATTICEPARAMETERS);printf(".\n");
		exit(1);
    }
	fscanf(ifile,"%lf %lf %lf",&deltat, &avev, &vev);fclose(ifile);aInvGeV=avev/vev;
	
#ifdef _INFO
	cout <<"\nParameters loaded: a_t/a = " << deltat << ", a*v = " << avev << "\n"; 
#endif
	
	if((ifile = fopen(PARAMETERS,"r"))==NULL)  
    {    
		printf("Cannot open ");printf(PARAMETERS);printf(".\n");
		exit(1);
    }
	fscanf(ifile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &gg,&am11sq,&am22sq, &am12sqR, &am12sqI,&lambda1,&lambda2, &lambda3, &lambda4, &lambda5R, &lambda5I, &av1, &av2, &theta12, &kappa);
	fclose(ifile); 
	
	if(arg==TRUE) kappa=arckappa;
	
	av1*=aInvGeV;av2*=aInvGeV;aInvGeV*=aInvGeV;am11sq*=aInvGeV; am22sq*=aInvGeV; am12sqR*=aInvGeV; am12sqI*=aInvGeV;
	
	if((ifile = fopen(TWOHDM,"r"))==NULL)  
    {    
		printf("Cannot open ");printf(TWOHDM);printf(".\n");
		exit(1);
    }
	fscanf(ifile,"%lf %lf",&t_stop,&t_save);
	fclose(ifile); 
	updateLatticePara();
}

/***********Make the lattice and nearest neighbours.*********/

void makeLattice()
{
	int i,j,k; 
	int x,y,z; 
	site *myself;	
	
	/* Find neighbors. */
	
	for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++)
    {
		myself=lattice+i*ny*nz+j*nz+k;
		myself->buddy[ZUP]=lattice+i*ny*nz+j*nz+(k+1)%nz;
		myself->buddy[YUP]=lattice+i*ny*nz+((j+1)%ny)*nz+k;
		myself->buddy[XUP]=lattice+((i+1)%nx)*ny*nz+j*nz+k;
		(myself->buddy[ZUP])->buddy[ZDOWN]=myself;
		(myself->buddy[YUP])->buddy[YDOWN]=myself;
		(myself->buddy[XUP])->buddy[XDOWN]=myself;
    }
	
	return;
}

/*********** Load an old lattice configuration.*********/

string prepareInitialFields(long seedForSeed,long lseed, int warmup){
	
	initialFieldsMC *iField=new initialFieldsMC(gg,  am11sq,  am22sq, am12sqR, am12sqI, lambda1, lambda2,lseed, warmup, deltat, avev, vev);
	ostringstream ostr;
	
	ostr << "para/N" << iField->getN() << "dt" << setfill('0') <<  setw(4) << int(1000.0*deltat) << "a" << setw(2) << int(100.0*avev) << "w" << warmup << "ss" << seedForSeed;
	if(lseed>0)
		ostr << "s" << lseed;
	else
		ostr << "sm" << -lseed;
	
	iField->generate();iField->dump_config(ostr.str());
	
	delete iField;
	
	return ostr.str();
}

void loadLatticeFromFile(string fname=INITIALFIELDS)
{
	int x,y,z,dir,i,j;
	FILE* conffile;
	
	// Higgs basis transformation
	RpTimesSU2 phi1o,phi2o;
	
	double am12sq=sqrt(am12sqR*am12sqR+am12sqI*am12sqI),norm,u1;
	
	complex<double> c(am12sqR/am12sq,am12sqI/am12sq);
	
	complex<double> u11,u12,u21,u22;
	
	u1=am22sq-am11sq+sqrt(4.0*am12sq*am12sq+(am11sq-am22sq)*(am11sq-am22sq));
	norm=1.0/sqrt( 4.0*am12sq*am12sq+u1*u1);
	u11=norm*u1*c;u21=-2.0*am12sq*norm;
	
	u1=am11sq-am22sq+sqrt(4.0*am12sq*am12sq+(am11sq-am22sq)*(am11sq-am22sq));
	norm=1.0/sqrt( 4*am12sq*am12sq+u1*u1);
	u12=norm*u1;u22=2.0*am12sq*norm*conj(c);
	
	//Load from file
	if((conffile = fopen(fname.c_str(),"rb"))==NULL) 
	{
		printf("cannot open file");printf(INITIALFIELDS);printf(".\n");
		exit(1);
	}

	if(fread(&deltat,sizeof deltat,1,conffile)!=1){printf("Write error14.\n");}
	if(fread(&avev,sizeof avev,1,conffile)!=1){printf("Write error14.\n");}
	if(fread(&warmup,sizeof warmup,1,conffile)!=1){printf("Write error14.\n");}
	if(fread(&theseed,sizeof theseed,1,conffile)!=1){printf("Write error14.\n");}
	
	if(fread(&nx, sizeof nx,1,conffile)!=1){printf("Read error.\n");}
	if(fread(&ny, sizeof ny,1,conffile)!=1){printf("Read error.\n");}
	if(fread(&nz, sizeof nz,1,conffile)!=1){printf("Read error.\n");}
	
#ifdef _INFO
	cout << "\nLoading initial fields with the following parameters:\n";
	cout << "N = " << nx << ", a_t/a = " << deltat << ", a*v = " << avev << ", kappa = " << kappa << ", warmup = " << warmup << ", seed = " << theseed << endl;
#endif
	
	volume1 = nx * ny * nz;
	volume = (double)(nx * ny * nz);
	
	/* Allocate space for lattice. */
	
	lattice = new site[volume1];
	if(lattice==NULL)
    {
		printf("No room for lattice\n");
		exit(1);
    }
	
	makeLattice();
	
	for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
		
		i = x*ny*nz+y*nz+z;
		
		/*Links at the initial time prevTime*/
		
		for(dir=XUP;dir<=ZUP;dir++)
			for(j=0;j<4;j++)
				if(fread(&lattice[i].link[prevTime][dir].e[j],sizeof lattice[i].link[prevTime][dir].e[j],1,conffile)!=1){printf("Read error.\n");}
		
		/*Electric Fields at prevTime*/
		
		for(dir=XUP;dir<=ZUP;dir++){
			for(j=1;j<4;j++){
				if(fread(&lattice[i].efield[prevTime][dir].e[j],sizeof lattice[i].efield[prevTime][dir].e[j],1,conffile)!=1){printf("Read error.\n");}
				lattice[i].efield[prevTime][dir].e[j]=-0.5*lattice[i].efield[prevTime][dir].e[j];//one half of the electric field in Anders' code
			}
			lattice[i].efield[prevTime][dir].setUnitary();lattice[i].efield[currentTime][dir]=lattice[i].efield[prevTime][dir];//Random guess for efield at currentTime = at prevTime
			lattice[i].link[currentTime][dir]=lattice[i].efield[prevTime][dir]*lattice[i].link[prevTime][dir];
		}
		
		/*Higgs Fields at the initial time prevTime*/
		
		if(fread(&lattice[i].phi1[prevTime].e[0],sizeof lattice[i].phi1[prevTime].e[0],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi1[prevTime].e[3],sizeof lattice[i].phi1[prevTime].e[3],1,conffile)!=1){printf("Read error.\n");}
		
		if(fread(&lattice[i].phi1[prevTime].e[2],sizeof lattice[i].phi1[prevTime].e[2],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi1[prevTime].e[1],sizeof lattice[i].phi1[prevTime].e[1],1,conffile)!=1){printf("Read error.\n");}		
		lattice[i].phi1[prevTime].e[2]=-lattice[i].phi1[prevTime].e[2];
		
		if(fread(&lattice[i].phi2[prevTime].e[0],sizeof lattice[i].phi2[prevTime].e[0],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi2[prevTime].e[3],sizeof lattice[i].phi2[prevTime].e[3],1,conffile)!=1){printf("Read error.\n");}
		
		if(fread(&lattice[i].phi2[prevTime].e[2],sizeof lattice[i].phi2[prevTime].e[2],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi2[prevTime].e[1],sizeof lattice[i].phi2[prevTime].e[1],1,conffile)!=1){printf("Read error.\n");}		
		lattice[i].phi2[prevTime].e[2]=-lattice[i].phi2[prevTime].e[2];
		
		
		/*Higgs Fields at the next time currentTime*/
		
		if(fread(&lattice[i].phi1[currentTime].e[0],sizeof lattice[i].phi1[currentTime].e[0],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi1[currentTime].e[3],sizeof lattice[i].phi1[currentTime].e[3],1,conffile)!=1){printf("Read error.\n");}
		
		
		if(fread(&lattice[i].phi1[currentTime].e[2],sizeof lattice[i].phi1[currentTime].e[2],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi1[currentTime].e[1],sizeof lattice[i].phi1[currentTime].e[1],1,conffile)!=1){printf("Read error.\n");}		
		lattice[i].phi1[currentTime].e[2]=-lattice[i].phi1[currentTime].e[2];lattice[i].phi1[currentTime]=(lattice[i].phi1[prevTime])+deltat*lattice[i].phi1[currentTime]; 
		
		if(fread(&lattice[i].phi2[currentTime].e[0],sizeof lattice[i].phi2[currentTime].e[0],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi2[currentTime].e[3],sizeof lattice[i].phi2[currentTime].e[3],1,conffile)!=1){printf("Read error.\n");}
		
		if(fread(&lattice[i].phi2[currentTime].e[2],sizeof lattice[i].phi2[currentTime].e[2],1,conffile)!=1){printf("Read error.\n");}
		if(fread(&lattice[i].phi2[currentTime].e[1],sizeof lattice[i].phi2[currentTime].e[1],1,conffile)!=1){printf("Read error.\n");}		
		lattice[i].phi2[currentTime].e[2]=-lattice[i].phi2[currentTime].e[2];lattice[i].phi2[currentTime]=(lattice[i].phi2[prevTime])+deltat*lattice[i].phi2[currentTime];  
		
		// Scale Higgs fields
		lattice[i].phi1[prevTime]=(0.5*sqrt(2.0)*lattice[i].phi1[prevTime]);lattice[i].phi2[prevTime]=(0.5*sqrt(2.0)*lattice[i].phi2[prevTime]);
		phi1o=u11*lattice[i].phi1[prevTime]/sqrt(lambda1)+u12*lattice[i].phi2[prevTime]/sqrt(lambda2);phi2o=u21*lattice[i].phi1[prevTime]/sqrt(lambda1)+u22*lattice[i].phi2[prevTime]/sqrt(lambda2);
		lattice[i].phi1[prevTime]=sqrt(lambda1)*phi1o;lattice[i].phi2[prevTime]=sqrt(lambda2)*phi2o;
		
		lattice[i].phi1[currentTime]=(0.5*sqrt(2.0)*lattice[i].phi1[currentTime]);lattice[i].phi2[currentTime]=(0.5*sqrt(2.0)*lattice[i].phi2[currentTime]);
		phi1o=u11*lattice[i].phi1[currentTime]/sqrt(lambda1)+u12*lattice[i].phi2[currentTime]/sqrt(lambda2);phi2o=u21*lattice[i].phi1[currentTime]/sqrt(lambda1)+u22*lattice[i].phi2[currentTime]/sqrt(lambda2);
		lattice[i].phi1[currentTime]=sqrt(lambda1)*phi1o;lattice[i].phi2[currentTime]=sqrt(lambda2)*phi2o;
				
		/* phi1 i sigma^3, phi2 i sigma^3 etc. */
		lattice[i].updatePhiExpression(prevTime);
		lattice[i].updatePhiExpression(currentTime);
	}
	
	updateUmn(currentTime);updateUmnBar(currentTime);
	updateUmn(prevTime);updateUmnBar(prevTime);
	updateTr231aUmnBar(prevTime);updateTr231aUmnBar(currentTime);updatedU0iUjk();
	
	fclose(conffile);
		
	return;
}



/*********** Load an old lattice configuration.*********/

void loadLattice(long lseed)
{
	int x,y,z,i,j;
	
	// Higgs basis transformation
	RpTimesSU2 phi1o,phi2o;
	
	double am12sq=sqrt(am12sqR*am12sqR+am12sqI*am12sqI),norm,u1;
	
	complex<double> c(am12sqR/am12sq,am12sqI/am12sq);
		
	complex<double> u11,u12,u21,u22;
	
	u1=am22sq-am11sq+sqrt(4.0*am12sq*am12sq+(am11sq-am22sq)*(am11sq-am22sq));
	norm=1.0/sqrt( 4.0*am12sq*am12sq+u1*u1);
	u11=norm*u1*c;u21=-2.0*am12sq*norm;
	
	u1=am11sq-am22sq+sqrt(4.0*am12sq*am12sq+(am11sq-am22sq)*(am11sq-am22sq));
	norm=1.0/sqrt( 4*am12sq*am12sq+u1*u1);
	u12=norm*u1;u22=2.0*am12sq*norm*conj(c);
	
	// Generate initial fields
	initialFieldsMC *iField=new initialFieldsMC(lseed);
	
	iField->generate();
	
	deltat=iField->deltat;avev=iField->avev;warmup=iField->warmup;theseed=iField->theseed;
	
	nx=N;ny=N;nz=N;
	
#ifdef _INFO
	cout << "\nLoading initial fields with the following parameters:\n";
	cout << "N = " << nx << ", a_t/a = " << deltat << ", a*v = " << avev << ", kappa = " << kappa << ", warmup = " << warmup << ", seed = " << theseed << endl;
#endif
	
	volume1 = nx * ny * nz;
	volume = (double)(nx * ny * nz);
	
	/* Allocate space for lattice. */
	
	lattice = new site[volume1];
	if(lattice==NULL)
    {
		printf("No room for lattice\n");
		exit(1);
    }
	
	makeLattice();
	
	for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
		
		i = x*ny*nz+y*nz+z;
		
		/*Links at the initial time prevTime*/
		
		for(int dir=XUP;dir<=ZUP;dir++)
			lattice[i].link[prevTime][dir]=RpTimesSU2(1,0,0,0);
		
		/*Electric Fields at prevTime*/
		
		lattice[i].efield[prevTime][XUP]=RpTimesSU2(0,-0.5*iField->Elec1_1[x][y][z],-0.5*iField->Elec1_2[x][y][z],-0.5*iField->Elec1_3[x][y][z]);//one half of the electric field in Anders' code
		lattice[i].efield[prevTime][XUP].setUnitary();lattice[i].efield[currentTime][XUP]=lattice[i].efield[prevTime][XUP];//Random guess for efield at currentTime = at prevTime
		lattice[i].link[currentTime][XUP]=lattice[i].efield[prevTime][XUP]*lattice[i].link[prevTime][XUP];
		
		lattice[i].efield[prevTime][YUP]=RpTimesSU2(0,-0.5*iField->Elec2_1[x][y][z],-0.5*iField->Elec2_2[x][y][z],-0.5*iField->Elec2_3[x][y][z]);//one half of the electric field in Anders' code
		lattice[i].efield[prevTime][YUP].setUnitary();lattice[i].efield[currentTime][YUP]=lattice[i].efield[prevTime][YUP];//Random guess for efield at currentTime = at prevTime
		lattice[i].link[currentTime][YUP]=lattice[i].efield[prevTime][YUP]*lattice[i].link[prevTime][YUP];		
		
		lattice[i].efield[prevTime][ZUP]=RpTimesSU2(0,-0.5*iField->Elec3_1[x][y][z],-0.5*iField->Elec3_2[x][y][z],-0.5*iField->Elec3_3[x][y][z]);//one half of the electric field in Anders' code
		lattice[i].efield[prevTime][ZUP].setUnitary();lattice[i].efield[currentTime][ZUP]=lattice[i].efield[prevTime][ZUP];//Random guess for efield at currentTime = at prevTime
		lattice[i].link[currentTime][ZUP]=lattice[i].efield[prevTime][ZUP]*lattice[i].link[prevTime][ZUP];		
		/*Higgs Fields at the initial time prevTime*/
		
		lattice[i].phi1[prevTime].e[0]=iField->out1[x][y][z].re;
		lattice[i].phi1[prevTime].e[3]=iField->out2[x][y][z].re;
		
		lattice[i].phi1[prevTime].e[2]=-iField->out3[x][y][z].re;
		lattice[i].phi1[prevTime].e[1]=iField->out4[x][y][z].re;		
		
		lattice[i].phi2[prevTime].e[0]=iField->out21[x][y][z].re;
		lattice[i].phi2[prevTime].e[3]=iField->out22[x][y][z].re;
		
		lattice[i].phi2[prevTime].e[2]=-iField->out23[x][y][z].re;
		lattice[i].phi2[prevTime].e[1]=iField->out24[x][y][z].re;			
		
		/*Higgs Fields at the next time currentTime*/
		
		lattice[i].phi1[currentTime].e[0]=iField->out5[x][y][z].re;
		lattice[i].phi1[currentTime].e[3]=iField->out6[x][y][z].re;
		
		
		lattice[i].phi1[currentTime].e[2]=-iField->out7[x][y][z].re;
		lattice[i].phi1[currentTime].e[1]=iField->out8[x][y][z].re;	
		lattice[i].phi1[currentTime]=(lattice[i].phi1[prevTime])+deltat*lattice[i].phi1[currentTime]; 
		
		lattice[i].phi2[currentTime].e[0]=iField->out25[x][y][z].re;
		lattice[i].phi2[currentTime].e[3]=iField->out26[x][y][z].re;
		
		lattice[i].phi2[currentTime].e[2]=-iField->out27[x][y][z].re;
		lattice[i].phi2[currentTime].e[1]=iField->out28[x][y][z].re;
		lattice[i].phi2[currentTime]=(lattice[i].phi2[prevTime])+deltat*lattice[i].phi2[currentTime]; 		// Scale Higgs fields
				
		lattice[i].phi1[prevTime]=(0.5*sqrt(2.0)*lattice[i].phi1[prevTime]);lattice[i].phi2[prevTime]=(0.5*sqrt(2.0)*lattice[i].phi2[prevTime]);
		phi1o=u11*lattice[i].phi1[prevTime]/sqrt(lambda1)+u12*lattice[i].phi2[prevTime]/sqrt(lambda2);phi2o=u21*lattice[i].phi1[prevTime]/sqrt(lambda1)+u22*lattice[i].phi2[prevTime]/sqrt(lambda2);
		lattice[i].phi1[prevTime]=sqrt(lambda1)*phi1o;lattice[i].phi2[prevTime]=sqrt(lambda2)*phi2o;

		lattice[i].phi1[currentTime]=(0.5*sqrt(2.0)*lattice[i].phi1[currentTime]);lattice[i].phi2[currentTime]=(0.5*sqrt(2.0)*lattice[i].phi2[currentTime]);
		phi1o=u11*lattice[i].phi1[currentTime]/sqrt(lambda1)+u12*lattice[i].phi2[currentTime]/sqrt(lambda2);phi2o=u21*lattice[i].phi1[currentTime]/sqrt(lambda1)+u22*lattice[i].phi2[currentTime]/sqrt(lambda2);
		lattice[i].phi1[currentTime]=sqrt(lambda1)*phi1o;lattice[i].phi2[currentTime]=sqrt(lambda2)*phi2o;
		
		/* phi1 i sigma^3, phi2 i sigma^3 etc. */
		lattice[i].updatePhiExpression(prevTime);
		lattice[i].updatePhiExpression(currentTime);
	}
	
	delete iField;

	updateUmn(currentTime);updateUmnBar(currentTime);
	updateUmn(prevTime);updateUmnBar(prevTime);
	updateTr231aUmnBar(prevTime);updateTr231aUmnBar(currentTime);updatedU0iUjk();

}
