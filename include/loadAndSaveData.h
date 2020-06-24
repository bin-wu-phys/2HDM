

/* Save Lattice at t */
#define ERRSAVELATTICEATT "Writing error in saveLatticeAtT\n"

void saveLatticeAtT(){
	FILE* ofile;ofstream out;
	
	if((ofile = fopen(LATTICEATT,"wb"))==NULL) 
    {
		printf("cannot open file ");printf(LATTICEATT);printf(".\n");
		exit(1);
    }
	
	/* Saving global scalars. */
	FWRITE(ofile, nx,ERRSAVELATTICEATT);
	FWRITE(ofile, ny,ERRSAVELATTICEATT);
	FWRITE(ofile, nz,ERRSAVELATTICEATT);
	FWRITE(ofile, volume1,ERRSAVELATTICEATT);
	FWRITE(ofile, n_save,ERRSAVELATTICEATT);
	FWRITE(ofile, t,ERRSAVELATTICEATT);
	FWRITE(ofile, t_stop,ERRSAVELATTICEATT);
	FWRITE(ofile, t_save,ERRSAVELATTICEATT);
	FWRITE(ofile, ncs_current,ERRSAVELATTICEATT);
	FWRITE(ofile, nextTime, ERRSAVELATTICEATT);
	FWRITE(ofile, currentTime, ERRSAVELATTICEATT);
	FWRITE(ofile, prevTime, ERRSAVELATTICEATT);
	FWRITE(ofile, volume, ERRSAVELATTICEATT);
	FWRITE(ofile, running_time,ERRSAVELATTICEATT);
	
	
	/* Parameters in the 2HDM */
	
	// parameters for solving higgs fields at next time step
	FWRITE(ofile, beta31,ERRSAVELATTICEATT);
	FWRITE(ofile, beta32,ERRSAVELATTICEATT);
	FWRITE(ofile, beta21a1,ERRSAVELATTICEATT);
	FWRITE(ofile, beta21a2,ERRSAVELATTICEATT);
	FWRITE(ofile, beta231a1,ERRSAVELATTICEATT);
	FWRITE(ofile, beta231a2,ERRSAVELATTICEATT);
	FWRITE(ofile, beta51I,ERRSAVELATTICEATT);
	FWRITE(ofile, beta52I,ERRSAVELATTICEATT);
	FWRITE(ofile, beta121R,ERRSAVELATTICEATT);
	FWRITE(ofile, beta121I,ERRSAVELATTICEATT);
	FWRITE(ofile, beta122R,ERRSAVELATTICEATT);
	FWRITE(ofile, beta122I,ERRSAVELATTICEATT);
	
	//parameter for solving E at current time
	FWRITE(ofile, betaE1,ERRSAVELATTICEATT);
	FWRITE(ofile, betaE2,ERRSAVELATTICEATT);
	
	//parameters for checking Gauss constraint
	FWRITE(ofile, betaGauss1,ERRSAVELATTICEATT);
	FWRITE(ofile, betaGauss2,ERRSAVELATTICEATT);
	
	FWRITE(ofile, gg,ERRSAVELATTICEATT);
	FWRITE(ofile, lambda1,ERRSAVELATTICEATT);
	FWRITE(ofile, lambda2,ERRSAVELATTICEATT);
	FWRITE(ofile, lambda3,ERRSAVELATTICEATT);
	FWRITE(ofile, lambda4,ERRSAVELATTICEATT);
	FWRITE(ofile, lambda5R,ERRSAVELATTICEATT);
	FWRITE(ofile, lambda5I,ERRSAVELATTICEATT);
	FWRITE(ofile, am11sq,ERRSAVELATTICEATT);
	FWRITE(ofile, am22sq,ERRSAVELATTICEATT);
	FWRITE(ofile, am12sqR,ERRSAVELATTICEATT);
	FWRITE(ofile, am12sqI,ERRSAVELATTICEATT);
	FWRITE(ofile, av1,ERRSAVELATTICEATT);
	FWRITE(ofile, av2,ERRSAVELATTICEATT);
	
	FWRITE(ofile, deltatsq,ERRSAVELATTICEATT);
	FWRITE(ofile, deltat,ERRSAVELATTICEATT);
	FWRITE(ofile, avev,ERRSAVELATTICEATT);
	FWRITE(ofile, vev,ERRSAVELATTICEATT);
	
	FWRITE(ofile, ncs,ERRSAVELATTICEATT);
	
	if(fwrite(iSigma,sizeof iSigma[0],3,ofile)!=3)printf(ERRSAVELATTICEATT);	
	if(fwrite(lattice,sizeof lattice[0],volume1,ofile)!=volume1)printf(ERRSAVELATTICEATT);
	
	fflush(ofile);
	fclose(ofile);
	
	out.open("N16dt001.dat",ios::app);
	out <<  t << "   " << ncs << "   " << calcPhi1sq() << "   " << calcPhi2sq() << " " << calcGauss() <<endl;
}


/* Load Lattice at t */
#define ERRREADLATTICEATT "Reading error in loadLatticeAtT\n"

void loadLatticeAtT(){
	FILE* ofile;
	
	if((ofile = fopen(LATTICEATT,"rb"))==NULL) 
    {
		printf("cannot open file ");printf(LATTICEATT);printf(".\n");
		exit(1);
    }
	
	/* Saving global scalars. */
	FREAD(ofile, nx,ERRREADLATTICEATT);
	FREAD(ofile, ny,ERRREADLATTICEATT);
	FREAD(ofile, nz,ERRREADLATTICEATT);
	FREAD(ofile, volume1,ERRREADLATTICEATT);
	FREAD(ofile, n_save,ERRREADLATTICEATT);
	
	lattice = new site[volume1];
	if(lattice==NULL)
    {
		printf("No room for lattice\n");
		exit(1);
    }	
	
	FREAD(ofile, t,ERRREADLATTICEATT);
	FREAD(ofile, t_stop,ERRREADLATTICEATT);
	FREAD(ofile, t_save,ERRREADLATTICEATT);
	FREAD(ofile, ncs_current,ERRREADLATTICEATT);
	FREAD(ofile, nextTime, ERRREADLATTICEATT);
	FREAD(ofile, currentTime, ERRREADLATTICEATT);
	FREAD(ofile, prevTime, ERRREADLATTICEATT);
	FREAD(ofile, volume, ERRREADLATTICEATT);
	FREAD(ofile, running_time,ERRREADLATTICEATT);
	
	
	/* Parameters in the 2HDM */
	
	// parameters for solving higgs fields at next time step
	FREAD(ofile, beta31,ERRSAVELATTICEATT);
	FREAD(ofile, beta32,ERRSAVELATTICEATT);
	FREAD(ofile, beta21a1,ERRSAVELATTICEATT);
	FREAD(ofile, beta21a2,ERRSAVELATTICEATT);
	FREAD(ofile, beta231a1,ERRSAVELATTICEATT);
	FREAD(ofile, beta231a2,ERRSAVELATTICEATT);
	FREAD(ofile, beta51I,ERRSAVELATTICEATT);
	FREAD(ofile, beta52I,ERRSAVELATTICEATT);
	FREAD(ofile, beta121R,ERRSAVELATTICEATT);
	FREAD(ofile, beta121I,ERRSAVELATTICEATT);
	FREAD(ofile, beta122R,ERRSAVELATTICEATT);
	FREAD(ofile, beta122I,ERRSAVELATTICEATT);
	
	//parameter for solving E at current time
	FREAD(ofile, betaE1,ERRSAVELATTICEATT);
	FREAD(ofile, betaE2,ERRSAVELATTICEATT);
	
	//parameters for checking Gauss constraint
	FREAD(ofile, betaGauss1,ERRSAVELATTICEATT);
	FREAD(ofile, betaGauss2,ERRSAVELATTICEATT);
	
	FREAD(ofile, gg,ERRREADLATTICEATT);
	FREAD(ofile, lambda1,ERRREADLATTICEATT);
	FREAD(ofile, lambda2,ERRREADLATTICEATT);
	FREAD(ofile, lambda3,ERRREADLATTICEATT);
	FREAD(ofile, lambda4,ERRREADLATTICEATT);
	FREAD(ofile, lambda5R,ERRREADLATTICEATT);
	FREAD(ofile, lambda5I,ERRREADLATTICEATT);
	FREAD(ofile, am11sq,ERRREADLATTICEATT);
	FREAD(ofile, am22sq,ERRREADLATTICEATT);
	FREAD(ofile, am12sqR,ERRREADLATTICEATT);
	FREAD(ofile, am12sqI,ERRREADLATTICEATT);
	FREAD(ofile, av1,ERRREADLATTICEATT);
	FREAD(ofile, av2,ERRREADLATTICEATT);
	
	FREAD(ofile, deltatsq,ERRREADLATTICEATT);
	FREAD(ofile, deltat,ERRREADLATTICEATT);
	FREAD(ofile, avev,ERRREADLATTICEATT);
	FREAD(ofile, vev,ERRREADLATTICEATT);
	
	FREAD(ofile, ncs,ERRREADLATTICEATT);
	
	if(fread(iSigma,sizeof iSigma[0],3,ofile)!=3)printf(ERRREADLATTICEATT);
	if(fread(lattice,sizeof lattice[0],volume1,ofile)!=volume1)printf(ERRREADLATTICEATT);
	
	fclose(ofile);
	
	makeLattice();
}
