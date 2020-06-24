//#define _INFO_INIT
/*********Standard libraries and header files.*********/

#define INITIALCFG "para/init.cfg"
/*********Global variables.****************************/

#define N	NX              
#define N2 (N/2)
#define factor 1.0                                 /* Coefficient on choice of da, db...*/
#define PI 3.14159265358979323846
/*************Random numbers.**************************/
double mersenne();
void seed_mersenne(long seed);
#define dran() mersenne()

class initialFieldsMC{
public:
	
/*********Files.***************************************/

FILE *conffile, *pfile, *ifile, *sfile;                    /* in- and output files */

/*********Arrays.**************************************/

fftw_complex in1[N][N][N],in2[N][N][N];                 /* Mode arrays */
fftw_complex in3[N][N][N],in4[N][N][N];
fftw_complex in5[N][N][N],in6[N][N][N];
fftw_complex in7[N][N][N],in8[N][N][N];
fftw_complex in21[N][N][N],in22[N][N][N];                 /* Mode arrays */
fftw_complex in23[N][N][N],in24[N][N][N];
fftw_complex in25[N][N][N],in26[N][N][N];
fftw_complex in27[N][N][N],in28[N][N][N];
fftw_complex out1[N][N][N],out2[N][N][N];
fftw_complex out3[N][N][N],out4[N][N][N];
fftw_complex out5[N][N][N],out6[N][N][N];
fftw_complex out7[N][N][N],out8[N][N][N];
fftw_complex out21[N][N][N],out22[N][N][N];
fftw_complex out23[N][N][N],out24[N][N][N];
fftw_complex out25[N][N][N],out26[N][N][N];
fftw_complex out27[N][N][N],out28[N][N][N];
fftw_complex rho1[N][N][N],rho2[N][N][N],rho3[N][N][N]; /* Charge density array, x and k space */
fftw_complex chi1[N][N][N],chi2[N][N][N],chi3[N][N][N]; /* Electric potential array, x and k */

double a[N][N][N][8];                                   /* Mode coefficients, real, phi */
double b[N][N][N][8];                                   /* Mode coefficients, imag, phi */
double c[N][N][N][8];                                   /* Mode coefficients, real, pi */
double d[N][N][N][8];                                   /* Mode coefficients, imag, pi */
double om_nk[N][N][N][2];                               /* omega_k [0], nk [1] */

double Elec1_1[N][N][N], Elec2_1[N][N][N], Elec3_1[N][N][N]; /* Efield arrays */
double Elec1_2[N][N][N], Elec2_2[N][N][N], Elec3_2[N][N][N]; /* Elec1_2 is flavour 2, lorentz 1*/
double Elec1_3[N][N][N], Elec2_3[N][N][N], Elec3_3[N][N][N];

fftwnd_plan planen;                                      /* FFTW plan variable */

double f1, f2, f3, f4, a_sq, f_sq, G, S;				 /* Global action contributions */
int hits, accept;                                 /* Counting variables */
double gg, betagt, betah1t, betah2t, lambda1, lambda2, deltat, am11sq, am22sq, am12sqR, am12sqI, avev;
long theseed; 
int warmup;

int higgs;												/*wb: higgs = 0 for higgs doublet 1 and higgs = 4 for higgs doublet 2*/

/**********Off we go...**********************************/

	initialFieldsMC(){// Initialize from configuration files
#ifdef _DEBUG
		cout << "Parameters inputing...\n";
#endif
		input();
#ifdef _DEBUG
		cout << "Parameters inputed\n";
#endif
	}
	
	initialFieldsMC(long lseed){// Initialize from configuration files except theseed
#ifdef _DEBUG
		cout << "Parameters inputing...\n";
#endif
		input();
		theseed=lseed;
#ifdef _DEBUG
		cout << "Parameters inputed\n";
#endif
	}
	
	initialFieldsMC(double _gg, double _am11sq, double _am22sq,double _am12sqR, double _am12sqI,double _lambda1,double _lambda2,// Parameters in the 2hdm
					long _theseed,int _warmup,// Parameters for intial fields
					double _deltat,double _avev,double vev){// Parameters for lattice
		
		input(_gg,  _am11sq,  _am22sq,_am12sqR,_am12sqI, _lambda1, _lambda2,_theseed, _warmup, _deltat, _avev, vev);		
	}
	
	int getN(){
		return N;
	}
	
	void generate(){
		int steps;
#ifdef _DEBUG
		cout << "Starting generating initial fields...\n";
#endif
		
		seed_mersenne(theseed);
		
		masssqDiagonal();
		
		/*Higgs 1*/
		higgs = 0;
		setup();
		
		if(am11sq>0){
			hits = accept = 0;
		
			for(steps=0;steps<warmup;steps++)
			{
				metropolis_loop();
				get_C();
			}
		
			metropolis_loop();
			get_C();
			check_charge();
		}
		
		/*Higgs 1*/
		higgs = 4;
		setup();
		
		if(am22sq>0){
		
			hits = accept = 0;
		
			for(steps=0;steps<warmup;steps++)
			{
				metropolis_loop();
				get_C();
			}
		
			metropolis_loop();
			get_C();
			check_charge();
		}
		
		/*Save initial higgs fields*/
		prepare_output();
		find_rho();
#ifdef _DEBUG
		cout << "Generating initial fields is done!\n";
#endif
		
	}

	void masssqDiagonal(){
		double am11sqtmp,am22sqtmp,am12sqsq,sqrtm;
		
		am12sqsq=am12sqR*am12sqR+am12sqI*am12sqI;
		sqrtm=sqrt( 4.0*am12sqsq+(am11sq-am22sq)*(am11sq-am22sq));
		am11sqtmp=am11sq;am22sqtmp=am22sq;
		
		am11sq=0.5*(am11sqtmp+am22sqtmp-sqrtm);am22sq=0.5*(am11sqtmp+am22sqtmp+sqrtm);
#ifdef	_INFO_INIT
		cout << am11sq << "   "  << am11sqtmp << "   " << am22sq << "   " << am22sqtmp << endl;
#endif
	}
	
/***********The inputfiles are opened and read*********/
	
void input()
{	
	double vev, aInvGeV;
	if((pfile = fopen(PARAMETERS,"r"))==NULL)  
    {    
		printf("Cannot open ");printf(PARAMETERS);printf(".\n");
		exit(1);
    }
	
	if((ifile = fopen(INITIALCFG,"r"))==NULL) 
    {    
		printf("Cannot open init.cfg.\n");
		exit(1);
    }
	
	if((sfile = fopen(LATTICEPARAMETERS,"r"))==NULL) 
    {    
		printf("Cannot open ");printf(LATTICEPARAMETERS);printf(".\n");
		exit(1);
    }
	
	fscanf(pfile,"%lf %lf %lf %lf %lf %lf %lf",
		   &gg,&am11sq,&am22sq,&am12sqR,&am12sqI,&lambda1,&lambda2);
	fclose(pfile);                        
	
	fscanf(ifile,"%ld %i",
		   &theseed,&warmup);
	fclose(ifile);                           
	
	fscanf(sfile,"%lf %lf %lf",
		   &deltat,&avev,&vev);
	fclose(sfile);                           
	aInvGeV=avev/vev;aInvGeV*=aInvGeV;
	am11sq*=(aInvGeV);am22sq*=aInvGeV;
	am12sqR*=aInvGeV;am12sqI*=aInvGeV;

#ifdef _INFO_INIT
	printf("\n");
	printf(" The input parameters are:\n");
	printf(" g squared		  = %f\n",gg);
	printf(" at/a			  = %f\n",deltat);
	printf(" a^2*m11^2        = %f\n",am11sq);
	printf(" a^2*m22^2        = %f\n",am22sq);
	printf(" a^2*m12^2        = %f+%f i\n",am12sqR,am12sqI);
	printf(" lambda_1         = %f\n",lambda1);
	printf(" lambda_2         = %f\n",lambda2);
	
	printf("\n");
#endif
	
	betagt = 4.0/(gg*deltat);
	betah1t = 0.25*gg*betagt/lambda1;
	betah2t = 0.25*gg*betagt/lambda2;

#ifdef _INFO_INIT
	
	printf(" The lattice beta coefficients:\n");
	printf(" beta_gt = %f\n",betagt);
	printf(" beta_h1t = %f\n",betah1t);
	printf(" beta_h2t = %f\n",betah2t);
	printf("\n");
#endif	
	return;
}

	
	void input(double _gg, double _am11sq, double _am22sq, double _am12sqR, double _am12sqI,//dimensionless
			   double _lambda1,double _lambda2,long _theseed,int _warmup,double _deltat,double _avev,double vev)
	{	
		double aInvGeV;
		
		gg=_gg;am11sq=_am11sq;am22sq=_am22sq;am12sqR=_am12sqR;am12sqI=_am12sqI;lambda1=_lambda1;lambda2=_lambda2;
		
		theseed=_theseed;warmup=_warmup;
		
		deltat=_deltat;avev=_avev;
		
/*		aInvGeV=avev/vev;aInvGeV*=aInvGeV;
		am11sq*=aInvGeV;am22sq*=aInvGeV;
*/		
		betagt = 4.0/(gg*deltat);
		betah1t = 0.25*gg*betagt/lambda1;
		betah2t = 0.25*gg*betagt/lambda2;
		
	}
	

/********Initializes arrays and nk, omega*******************/

void setup()
{
	double pref, musq, pref3, pref4, prefag, amHsq;
	int i,k, l, m, mink, minl, minm;
	
	(higgs==0)?(amHsq=am11sq):(amHsq=am22sq);
	
	musq=amHsq*0.5+6.0;
	pref=2.0*PI/N;
	S = 0.0;
	pref3 = 1.0;
	prefag = sqrt(1.0);
	
	f1 = f2 = f3 = f4 = 0.0;
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {
		mink=(N-k)%N;
		minl=(N-l)%N;
		minm=(N-m)%N;
		
		if(((amHsq*0.5-6.0+2.0*(cos(pref*k)+cos(pref*l)+cos(pref*m)))>0.0)&&(amHsq>0))
		{
			om_nk[k][l][m][0]=sqrt(musq-2.0*(cos(pref*k)+cos(pref*l)+cos(pref*m)));
			om_nk[k][l][m][1]=0.5;
#ifdef _INFO_INIT			
			printf("%i %i %i %f %f\n",k,l,m,om_nk[k][l][m][0],om_nk[k][l][m][1]);	
#endif			
			for(i=0;i<4;i++){
				a[k][l][m][i+higgs]=prefag*(0.5-dran())*sqrt(om_nk[k][l][m][1]);
				a[mink][minl][minm][i+higgs]=a[k][l][m][i+higgs];
			}
			for(i=0;i<4;i++){
				b[k][l][m][i+higgs]=prefag*(0.5-dran())*sqrt(om_nk[k][l][m][1]);
				b[mink][minl][minm][i+higgs]=-b[k][l][m][i+higgs];
			}			
			for(i=0;i<4;i++){
				c[k][l][m][i+higgs]=prefag*(0.5-dran())*sqrt(om_nk[k][l][m][1]);
				c[mink][minl][minm][i+higgs]=c[k][l][m][i+higgs];
			}			
			for(i=0;i<4;i++){
				d[k][l][m][i+higgs]=prefag*(0.5-dran())*sqrt(om_nk[k][l][m][1]);
				d[mink][minl][minm][i+higgs]=-d[k][l][m][i+higgs];
			}
		}
		else
		{
			om_nk[k][l][m][0]=10.0;
			om_nk[k][l][m][1]=0.00000000000001;
			
			for(i=0;i<4;i++){
				a[k][l][m][i+higgs]=0.0;			
				b[k][l][m][i+higgs]=0.0;			
				c[k][l][m][i+higgs]=0.0;			
				d[k][l][m][i+higgs]=0.0;
			}
		}
    }
	
	/*Get initial action contributions, and set corner b's and d's to zero ****/
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {
		if((k==0)&&(l==0)&&(m==0))
		{
			a_sq = a[0][0][0][0+higgs]*a[0][0][0][0+higgs]
			+a[0][0][0][1+higgs]*a[0][0][0][1+higgs]
			+a[0][0][0][2+higgs]*a[0][0][0][2+higgs]
			+a[0][0][0][3+higgs]*a[0][0][0][3+higgs];
			
			b[0][0][0][0+higgs]=b[0][0][0][1+higgs]=b[0][0][0][2+higgs]=b[0][0][0][3+higgs]=0.0;
			d[0][0][0][0+higgs]=d[0][0][0][1+higgs]=d[0][0][0][2+higgs]=d[0][0][0][3+higgs]=0.0;
			
		}	
		else
		{ 
			if(((k==N2)||(k==0))&&((l==N2)||(l==0))&&((m==N2)||(m==0)))
			{
				b[k][l][m][0+higgs]=b[k][l][m][1+higgs]=b[k][l][m][2+higgs]=b[k][l][m][3+higgs]=0.0;
				d[k][l][m][0+higgs]=d[k][l][m][1+higgs]=d[k][l][m][2+higgs]=d[k][l][m][3+higgs]=0.0; 
				pref3 = 1.0;
			}
			else
			{
				pref3 = 0.5;
			}
			
			f1+=pref3*(a[k][l][m][0+higgs]*c[k][l][m][0+higgs]+a[k][l][m][1+higgs]*c[k][l][m][1+higgs]
					   +a[k][l][m][2+higgs]*c[k][l][m][2+higgs]+a[k][l][m][3+higgs]*c[k][l][m][3+higgs]
					   +b[k][l][m][0+higgs]*d[k][l][m][0+higgs]+b[k][l][m][1+higgs]*d[k][l][m][1+higgs]
					   +b[k][l][m][2+higgs]*d[k][l][m][2+higgs]+b[k][l][m][3+higgs]*d[k][l][m][3+higgs]);
			f2+=pref3*(a[k][l][m][2+higgs]*c[k][l][m][1+higgs]-a[k][l][m][1+higgs]*c[k][l][m][2+higgs]
					   +a[k][l][m][0+higgs]*c[k][l][m][3+higgs]-a[k][l][m][3+higgs]*c[k][l][m][0+higgs]
					   +b[k][l][m][2+higgs]*d[k][l][m][1+higgs]-b[k][l][m][1+higgs]*d[k][l][m][2+higgs]
					   +b[k][l][m][0+higgs]*d[k][l][m][3+higgs]-b[k][l][m][3+higgs]*d[k][l][m][0+higgs]);
			f3+=pref3*(a[k][l][m][2+higgs]*c[k][l][m][0+higgs]-a[k][l][m][0+higgs]*c[k][l][m][2+higgs]
					   +a[k][l][m][3+higgs]*c[k][l][m][1+higgs]-a[k][l][m][1+higgs]*c[k][l][m][3+higgs]
					   +b[k][l][m][2+higgs]*d[k][l][m][0+higgs]-b[k][l][m][0+higgs]*d[k][l][m][2+higgs]
					   +b[k][l][m][3+higgs]*d[k][l][m][1+higgs]-b[k][l][m][1+higgs]*d[k][l][m][3+higgs]);
			f4+=pref3*(a[k][l][m][0+higgs]*c[k][l][m][1+higgs]-a[k][l][m][1+higgs]*c[k][l][m][0+higgs]
					   +a[k][l][m][3+higgs]*c[k][l][m][2+higgs]-a[k][l][m][2+higgs]*c[k][l][m][3+higgs]
					   +b[k][l][m][0+higgs]*d[k][l][m][1+higgs]-b[k][l][m][1+higgs]*d[k][l][m][0+higgs]
					   +b[k][l][m][3+higgs]*d[k][l][m][2+higgs]-b[k][l][m][2+higgs]*d[k][l][m][3+higgs]);
		}
    }  
	
	f_sq = f1*f1 + f2*f2 + f3*f3 + f4*f4;
	G = f_sq/(a_sq*2.0*om_nk[0][0][0][1]);
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {
		
		if((k==0)&&(l==0)&&(m==0)){pref4=0.0;}else{pref4=1.0;}
		
		if(((k==N2)||(k==0))&&((l==N2)||(l==0))&&((m==N2)||(m==0)))
	    {pref3 = 1.0;}else{pref3 = 0.5;}
		
		S+=(1.0/(2.0*om_nk[k][l][m][1]))*( a[k][l][m][0+higgs]*a[k][l][m][0+higgs]
										  +b[k][l][m][0+higgs]*b[k][l][m][0+higgs]
										  +pref4*c[k][l][m][0+higgs]*c[k][l][m][0+higgs]
										  +d[k][l][m][0+higgs]*d[k][l][m][0+higgs]
										  
										  +a[k][l][m][1+higgs]*a[k][l][m][1+higgs]
										  +b[k][l][m][1+higgs]*b[k][l][m][1+higgs]
										  +pref4*c[k][l][m][1+higgs]*c[k][l][m][1+higgs]
										  +d[k][l][m][1+higgs]*d[k][l][m][1+higgs]
										  
										  +a[k][l][m][2+higgs]*a[k][l][m][2+higgs]
										  +b[k][l][m][2+higgs]*b[k][l][m][2+higgs]
										  +pref4*c[k][l][m][2+higgs]*c[k][l][m][2+higgs]
										  +d[k][l][m][2+higgs]*d[k][l][m][2+higgs]
										  
										  +a[k][l][m][3+higgs]*a[k][l][m][3+higgs]
										  +b[k][l][m][3+higgs]*b[k][l][m][3+higgs]
										  +pref4*c[k][l][m][3+higgs]*c[k][l][m][3+higgs]
										  +d[k][l][m][3+higgs]*d[k][l][m][3+higgs]
										  );
    }
	
	S+=G;
	
	return;
}

/* Get the c[0][0][0] from constraint */

void get_C()
{
	c[0][0][0][0+higgs]=-(a[0][0][0][0+higgs]*f1-a[0][0][0][3+higgs]*f2+a[0][0][0][2+higgs]*f3-a[0][0][0][1+higgs]*f4)/a_sq;
	c[0][0][0][1+higgs]=-(a[0][0][0][1+higgs]*f1+a[0][0][0][2+higgs]*f2+a[0][0][0][3+higgs]*f3+a[0][0][0][0+higgs]*f4)/a_sq;
	c[0][0][0][2+higgs]=-(a[0][0][0][2+higgs]*f1-a[0][0][0][1+higgs]*f2-a[0][0][0][0+higgs]*f3+a[0][0][0][3+higgs]*f4)/a_sq;
	c[0][0][0][3+higgs]=-(a[0][0][0][3+higgs]*f1+a[0][0][0][0+higgs]*f2-a[0][0][0][1+higgs]*f3-a[0][0][0][2+higgs]*f4)/a_sq;
	
	return;
}


/* Check that charge is conserved in a,b,c,d */

void check_charge()
{
	double the_charge_1, the_charge_2, the_charge_3, pref3;
	int k,l,m;
	
	the_charge_1 = the_charge_2 = the_charge_3 = 0.0;
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {
		
		if(((k==N2)||(k==0))&&((l==N2)||(l==0))&&((m==N2)||(m==0)))
	    {pref3 = 1.0;}else{pref3 = 0.5;}
		
		the_charge_1+=pref3*(c[k][l][m][1+higgs]*a[k][l][m][2+higgs]-c[k][l][m][2+higgs]*a[k][l][m][1+higgs]+
							 d[k][l][m][1+higgs]*b[k][l][m][2+higgs]-d[k][l][m][2+higgs]*b[k][l][m][1+higgs]+
							 c[k][l][m][3+higgs]*a[k][l][m][0+higgs]-c[k][l][m][0+higgs]*a[k][l][m][3+higgs]+
							 d[k][l][m][3+higgs]*b[k][l][m][0+higgs]-d[k][l][m][0+higgs]*b[k][l][m][3+higgs]);
		the_charge_2+=pref3*(c[k][l][m][0+higgs]*a[k][l][m][2+higgs]-c[k][l][m][2+higgs]*a[k][l][m][0+higgs]+
							 d[k][l][m][0+higgs]*b[k][l][m][2+higgs]-d[k][l][m][2+higgs]*b[k][l][m][0+higgs]+
							 c[k][l][m][1+higgs]*a[k][l][m][3+higgs]-c[k][l][m][3+higgs]*a[k][l][m][1+higgs]+
							 d[k][l][m][1+higgs]*b[k][l][m][3+higgs]-d[k][l][m][3+higgs]*b[k][l][m][1+higgs]);
		the_charge_3+=pref3*(c[k][l][m][1+higgs]*a[k][l][m][0+higgs]-c[k][l][m][0+higgs]*a[k][l][m][1+higgs]+
							 d[k][l][m][1+higgs]*b[k][l][m][0+higgs]-d[k][l][m][0+higgs]*b[k][l][m][1+higgs]+
							 c[k][l][m][2+higgs]*a[k][l][m][3+higgs]-c[k][l][m][3+higgs]*a[k][l][m][2+higgs]+
							 d[k][l][m][2+higgs]*b[k][l][m][3+higgs]-d[k][l][m][3+higgs]*b[k][l][m][2+higgs]);
    }
	
	/*printf("%.15e %.15e %.15e\n",the_charge_1,the_charge_2,the_charge_3); */
	
	return;
}


/*Choose a random site and field. Check if it is a corner */

void metropolis_loop()
{
	int i, k, l, m, field;
	double amHsq;
	
	(higgs==0)?(amHsq=am11sq):(amHsq=am22sq);	
	for(i=0;i<(N*N*N*8);i++)
    {
		do{k = (int)(N*dran());}while(k==N);
		do{l = (int)(N*dran());}while(l==N);
		do{m = (int)(N*dran());}while(m==N);
		
		hits++;
		
		if(om_nk[k][l][m][0]<sqrt(amHsq))
		{
			if(((k==0)||(k==N2))&&((l==0)||(l==N2))&&((m==0)||(m==N2)))
			{
				do{field = (int)(8*dran());}while(field==8);
				hit_corner(k,l,m,field);
			}
			else
			{
				do{field = (int)(16*dran());}while(field==16);
				hit(k,l,m,field);
			}
		}
    }
	
	return;
}


/* Choose hit subroutine given field value */

void hit(int x, int y, int z, int deg)
{
	if(deg==0){hit_a_1(x,y,z,1.0);}
	if(deg==1){hit_a_2(x,y,z,1.0);}
	if(deg==2){hit_a_3(x,y,z,1.0);}
	if(deg==3){hit_a_4(x,y,z,1.0);}
	if(deg==4){hit_b_1(x,y,z,1.0);}
	if(deg==5){hit_b_2(x,y,z,1.0);}
	if(deg==6){hit_b_3(x,y,z,1.0);}
	if(deg==7){hit_b_4(x,y,z,1.0);}
	if(deg==8){hit_c_1(x,y,z,1.0);}
	if(deg==9){hit_c_2(x,y,z,1.0);}
	if(deg==10){hit_c_3(x,y,z,1.0);}
	if(deg==11){hit_c_4(x,y,z,1.0);}
	if(deg==12){hit_d_1(x,y,z,1.0);}
	if(deg==13){hit_d_2(x,y,z,1.0);}
	if(deg==14){hit_d_3(x,y,z,1.0);}
	if(deg==15){hit_d_4(x,y,z,1.0);}
	
	return;
}

/* If chosen site is a corner, only a and c are allowed to be updated*/

void hit_corner(int x, int y, int z, int deg)
{
	
	if((x==0)&&(y==0)&&(z==0))
    {
		hit_a_zero((deg%4)+higgs);
    }
	else
    {
		if(deg==0){hit_a_1(x,y,z,0.5);}
		if(deg==1){hit_a_2(x,y,z,0.5);}
		if(deg==2){hit_a_3(x,y,z,0.5);}
		if(deg==3){hit_a_4(x,y,z,0.5);}
		if(deg==4){hit_c_1(x,y,z,0.5);}
		if(deg==5){hit_c_2(x,y,z,0.5);}
		if(deg==6){hit_c_3(x,y,z,0.5);}
		if(deg==7){hit_c_4(x,y,z,0.5);}
    } 
	
	return;
}

/* If chosen site is 0,0,0, only a is allowed to be updated, and in a special way */

void hit_a_zero(int degr)
{
	double da, da_s, df1, df2, df3, df4, ds, da_sq, df_sq, a_sqs;
	
	da = factor*(0.5-dran())*sqrt(om_nk[0][0][0][1]);
	da_s = (((2.0*a[0][0][0][degr]+da)*da)+a_sq);
	da_sq = 1.0/da_s-1.0/a_sq;
	df1 = 0.0;
	df2 = 0.0;
	df3 = 0.0;
	df4 = 0.0;
	df_sq = 0.0;
	ds = ((2.0*a[0][0][0][degr]+da)*da)/(2.0*om_nk[0][0][0][1])
	+f_sq*da_sq/(2.0*om_nk[0][0][0][1])+2.0*log(da_s/a_sq);
	
	if(dran()<exp(-ds))
    {
		a[0][0][0][degr]+=da;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=(2.0*a[0][0][0][degr]-da)*da;
		
		accept++;
    }
	
	return;
}

/* Update a, field 1 */

void hit_a_1(int kx, int ky, int kz, double prefs)
{
	double da, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	da = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= c[kx][ky][kz][0+higgs]*da;
	df2= c[kx][ky][kz][3+higgs]*da;
	df3= -c[kx][ky][kz][2+higgs]*da;
	df4= c[kx][ky][kz][1+higgs]*da;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*a[kx][ky][kz][0+higgs]+da)*da)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		a[kx][ky][kz][0+higgs]+=da;
		a[minkx][minky][minkz][0+higgs]+=da;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update a, field 2 */

void hit_a_2(int kx, int ky, int kz, double prefs)
{
	double da, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	da = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= c[kx][ky][kz][1+higgs]*da;
	df2= -c[kx][ky][kz][2+higgs]*da;
	df3= -c[kx][ky][kz][3+higgs]*da;
	df4= -c[kx][ky][kz][0+higgs]*da;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*a[kx][ky][kz][1+higgs]+da)*da)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		a[kx][ky][kz][1+higgs]+=da;
		a[minkx][minky][minkz][1+higgs]+=da;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update a, field 3 */

void hit_a_3(int kx, int ky, int kz, double prefs)
{
	double da, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	da = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= c[kx][ky][kz][2+higgs]*da;
	df2= c[kx][ky][kz][1+higgs]*da;
	df3= c[kx][ky][kz][0+higgs]*da;
	df4= -c[kx][ky][kz][3+higgs]*da;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*a[kx][ky][kz][2+higgs]+da)*da)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		a[kx][ky][kz][2+higgs]+=da;
		a[minkx][minky][minkz][2+higgs]+=da;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update a, field 4 */

void hit_a_4(int kx, int ky, int kz, double prefs)
{
	double da, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	da = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= c[kx][ky][kz][3+higgs]*da;
	df2= -c[kx][ky][kz][0+higgs]*da;
	df3= c[kx][ky][kz][1+higgs]*da;
	df4= c[kx][ky][kz][2+higgs]*da;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*a[kx][ky][kz][3+higgs]+da)*da)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		a[kx][ky][kz][3+higgs]+=da;
		a[minkx][minky][minkz][3+higgs]+=da;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update b, field 1 */

void hit_b_1(int kx, int ky, int kz, double prefs)
{
	double db, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	db = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= d[kx][ky][kz][0+higgs]*db;
	df2= d[kx][ky][kz][3+higgs]*db;
	df3= -d[kx][ky][kz][2+higgs]*db;
	df4= d[kx][ky][kz][1+higgs]*db;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*b[kx][ky][kz][0+higgs]+db)*db)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		b[kx][ky][kz][0+higgs]+=db;
		b[minkx][minky][minkz][0+higgs]-=db;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update b, field 21 */

void hit_b_2(int kx, int ky, int kz, double prefs)
{
	double db, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	db = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= d[kx][ky][kz][1+higgs]*db;
	df2= -d[kx][ky][kz][2+higgs]*db;
	df3= -d[kx][ky][kz][3+higgs]*db;
	df4= -d[kx][ky][kz][0+higgs]*db;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*b[kx][ky][kz][1+higgs]+db)*db)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		b[kx][ky][kz][1+higgs]+=db;
		b[minkx][minky][minkz][1+higgs]-=db;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update b, field 3 */

void hit_b_3(int kx, int ky, int kz, double prefs)
{
	double db, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	db = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= d[kx][ky][kz][2+higgs]*db;
	df2= d[kx][ky][kz][1+higgs]*db;
	df3= d[kx][ky][kz][0+higgs]*db;
	df4= -d[kx][ky][kz][3+higgs]*db;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*b[kx][ky][kz][2+higgs]+db)*db)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		b[kx][ky][kz][2+higgs]+=db;
		b[minkx][minky][minkz][2+higgs]-=db;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update b, field 4 */

void hit_b_4(int kx, int ky, int kz, double prefs)
{
	double db, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	db = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= d[kx][ky][kz][3+higgs]*db;
	df2= -d[kx][ky][kz][0+higgs]*db;
	df3= d[kx][ky][kz][1+higgs]*db;
	df4= d[kx][ky][kz][2+higgs]*db;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*b[kx][ky][kz][3+higgs]+db)*db)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		b[kx][ky][kz][3+higgs]+=db;
		b[minkx][minky][minkz][3+higgs]-=db;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update c, field 1 */

void hit_c_1(int kx, int ky, int kz, double prefs)
{
	double dc, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dc = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= a[kx][ky][kz][0+higgs]*dc;
	df2= -a[kx][ky][kz][3+higgs]*dc;
	df3= a[kx][ky][kz][2+higgs]*dc;
	df4= -a[kx][ky][kz][1+higgs]*dc;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*c[kx][ky][kz][0+higgs]+dc)*dc)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		c[kx][ky][kz][0+higgs]+=dc;
		c[minkx][minky][minkz][0+higgs]+=dc;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update c, field 2 */

void hit_c_2(int kx, int ky, int kz, double prefs)
{
	double dc, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dc = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= a[kx][ky][kz][1+higgs]*dc;
	df2= a[kx][ky][kz][2+higgs]*dc;
	df3= a[kx][ky][kz][3+higgs]*dc;
	df4= a[kx][ky][kz][0+higgs]*dc;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*c[kx][ky][kz][1+higgs]+dc)*dc)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		c[kx][ky][kz][1+higgs]+=dc;
		c[minkx][minky][minkz][1+higgs]+=dc;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update c, field 3 */

void hit_c_3(int kx, int ky, int kz, double prefs)
{
	double dc, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dc = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= a[kx][ky][kz][2+higgs]*dc;
	df2= -a[kx][ky][kz][1+higgs]*dc;
	df3= -a[kx][ky][kz][0+higgs]*dc;
	df4= a[kx][ky][kz][3+higgs]*dc;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*c[kx][ky][kz][2+higgs]+dc)*dc)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		c[kx][ky][kz][2+higgs]+=dc;
		c[minkx][minky][minkz][2+higgs]+=dc;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update c, field 4 */

void hit_c_4(int kx, int ky, int kz, double prefs)
{
	double dc, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dc = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= a[kx][ky][kz][3+higgs]*dc;
	df2= a[kx][ky][kz][0+higgs]*dc;
	df3= -a[kx][ky][kz][1+higgs]*dc;
	df4= -a[kx][ky][kz][2+higgs]*dc;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*c[kx][ky][kz][3+higgs]+dc)*dc)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		c[kx][ky][kz][3+higgs]+=dc;
		c[minkx][minky][minkz][3+higgs]+=dc;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update d, field 1 */

void hit_d_1(int kx, int ky, int kz, double prefs)
{
	double dd, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dd = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= b[kx][ky][kz][0+higgs]*dd;
	df2= -b[kx][ky][kz][3+higgs]*dd;
	df3= b[kx][ky][kz][2+higgs]*dd;
	df4= -b[kx][ky][kz][1+higgs]*dd;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*d[kx][ky][kz][0+higgs]+dd)*dd)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		d[kx][ky][kz][0+higgs]+=dd;
		d[minkx][minky][minkz][0+higgs]-=dd;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update d, field 2 */

void hit_d_2(int kx, int ky, int kz, double prefs)
{
	double dd, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dd = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= b[kx][ky][kz][1+higgs]*dd;
	df2= b[kx][ky][kz][2+higgs]*dd;
	df3= b[kx][ky][kz][3+higgs]*dd;
	df4= b[kx][ky][kz][0+higgs]*dd;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*d[kx][ky][kz][1+higgs]+dd)*dd)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		d[kx][ky][kz][1+higgs]+=dd;
		d[minkx][minky][minkz][1+higgs]-=dd;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update d, field 3 */

void hit_d_3(int kx, int ky, int kz, double prefs)
{
	double dd, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dd = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= b[kx][ky][kz][2+higgs]*dd;
	df2= -b[kx][ky][kz][1+higgs]*dd;
	df3= -b[kx][ky][kz][0+higgs]*dd;
	df4= b[kx][ky][kz][3+higgs]*dd;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*d[kx][ky][kz][2+higgs]+dd)*dd)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds))
    {
		d[kx][ky][kz][2+higgs]+=dd;
		d[minkx][minky][minkz][2+higgs]-=dd;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	return;
}

/* Update d, field 4 */

void hit_d_4(int kx, int ky, int kz, double prefs)
{
	double dd, df1, df2, df3, df4, ds, da_sq, df_sq;
	int minkx, minky, minkz;
	
	
	dd = factor*(0.5-dran())*sqrt(om_nk[kx][ky][kz][1]);
	df1= b[kx][ky][kz][3+higgs]*dd;
	df2= b[kx][ky][kz][0+higgs]*dd;
	df3= -b[kx][ky][kz][1+higgs]*dd;
	df4= -b[kx][ky][kz][2+higgs]*dd;
	df_sq = ((2.0*f1+df1)*df1+(2.0*f2+df2)*df2+(2.0*f3+df3)*df3+(2.0*f4+df4)*df4);
	ds = ((2.0*d[kx][ky][kz][3+higgs]+dd)*dd)/(2.0*om_nk[kx][ky][kz][1])
	+df_sq/(2.0*om_nk[0][0][0][1]*a_sq);
	da_sq = 0.0;
	
	minkx = (N-kx)%N; minky = (N-ky)%N; minkz = (N-kz)%N;
	
	if(dran()<exp(-ds)) 
    {
		d[kx][ky][kz][3+higgs]+=dd;
		d[minkx][minky][minkz][3+higgs]-=dd;
		S+=ds;
		f1+=df1;
		f2+=df2;
		f3+=df3;
		f4+=df4;
		f_sq+=df_sq;
		a_sq+=da_sq;
		
		accept++;
    }
	
	return;
}



/*Get the phi_k and pi_k from a,b,c,d and Fourier tranform to phi_x and pi_x*/

void prepare_output()
{
	int k, l, m;
	double pref1, pref2, prefactor1, prefactor2;
	int mink, minl, minm;
	int x,y,z;
	int nx, ny, nz;
	
	nx=ny=nz=N;
	prefactor1 = sqrt(lambda1)/(sqrt((double)(nx*ny*nz)));
	prefactor2 = sqrt(lambda2)/(sqrt((double)(nx*ny*nz)));
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {
		if(((k==0)||(k==N2))&&((l==0)||(l==N2))&&((m==0)||(m==N2)))
		{
			pref1 = 1.0/sqrt(om_nk[k][l][m][0]);
			pref2 = sqrt(om_nk[k][l][m][0]);
		}
		else
		{
			pref1 = 1.0/sqrt(2.0*om_nk[k][l][m][0]);
			pref2 = sqrt(om_nk[k][l][m][0]/2.0);
		}
		
		in1[k][l][m].re = pref1*a[k][l][m][0];
		in1[k][l][m].im = pref1*b[k][l][m][0];
		in2[k][l][m].re = pref1*a[k][l][m][1];
		in2[k][l][m].im = pref1*b[k][l][m][1];
		in3[k][l][m].re = pref1*a[k][l][m][2];
		in3[k][l][m].im = pref1*b[k][l][m][2];
		in4[k][l][m].re = pref1*a[k][l][m][3];
		in4[k][l][m].im = pref1*b[k][l][m][3];
		
		in5[k][l][m].re = pref2*c[k][l][m][0];
		in5[k][l][m].im = pref2*d[k][l][m][0];
		in6[k][l][m].re = pref2*c[k][l][m][1];
		in6[k][l][m].im = pref2*d[k][l][m][1];
		in7[k][l][m].re = pref2*c[k][l][m][2];
		in7[k][l][m].im = pref2*d[k][l][m][2];
		in8[k][l][m].re = pref2*c[k][l][m][3];
		in8[k][l][m].im = pref2*d[k][l][m][3];
		
		in21[k][l][m].re = pref1*a[k][l][m][4];
		in21[k][l][m].im = pref1*b[k][l][m][4];
		in22[k][l][m].re = pref1*a[k][l][m][5];
		in22[k][l][m].im = pref1*b[k][l][m][5];
		in23[k][l][m].re = pref1*a[k][l][m][6];
		in23[k][l][m].im = pref1*b[k][l][m][6];
		in24[k][l][m].re = pref1*a[k][l][m][7];
		in24[k][l][m].im = pref1*b[k][l][m][7];
		
		in25[k][l][m].re = pref2*c[k][l][m][4];
		in25[k][l][m].im = pref2*d[k][l][m][4];
		in26[k][l][m].re = pref2*c[k][l][m][5];
		in26[k][l][m].im = pref2*d[k][l][m][5];
		in27[k][l][m].re = pref2*c[k][l][m][6];
		in27[k][l][m].im = pref2*d[k][l][m][6];
		in28[k][l][m].re = pref2*c[k][l][m][7];
		in28[k][l][m].im = pref2*d[k][l][m][7];
		
    }
	
	planen = fftw3d_create_plan(nx,ny,nz,FFTW_BACKWARD,FFTW_MEASURE);
	
	fftwnd_one(planen, &in1[0][0][0], &out1[0][0][0]);
	fftwnd_one(planen, &in2[0][0][0], &out2[0][0][0]);
	fftwnd_one(planen, &in3[0][0][0], &out3[0][0][0]);
	fftwnd_one(planen, &in4[0][0][0], &out4[0][0][0]);
	fftwnd_one(planen, &in5[0][0][0], &out5[0][0][0]);
	fftwnd_one(planen, &in6[0][0][0], &out6[0][0][0]);
	fftwnd_one(planen, &in7[0][0][0], &out7[0][0][0]);
	fftwnd_one(planen, &in8[0][0][0], &out8[0][0][0]);
	
	fftwnd_one(planen, &in21[0][0][0], &out21[0][0][0]);
	fftwnd_one(planen, &in22[0][0][0], &out22[0][0][0]);
	fftwnd_one(planen, &in23[0][0][0], &out23[0][0][0]);
	fftwnd_one(planen, &in24[0][0][0], &out24[0][0][0]);
	fftwnd_one(planen, &in25[0][0][0], &out25[0][0][0]);
	fftwnd_one(planen, &in26[0][0][0], &out26[0][0][0]);
	fftwnd_one(planen, &in27[0][0][0], &out27[0][0][0]);
	fftwnd_one(planen, &in28[0][0][0], &out28[0][0][0]);
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {  
		z = m;
		y = l;
		x = k;
		
		out1[x][y][z].re=prefactor1*out1[x][y][z].re;
		out1[x][y][z].im=prefactor1*out1[x][y][z].im;
		out2[x][y][z].re=prefactor1*out2[x][y][z].re;
		out2[x][y][z].im=prefactor1*out2[x][y][z].im;
		out3[x][y][z].re=prefactor1*out3[x][y][z].re;
		out3[x][y][z].im=prefactor1*out3[x][y][z].im;
		out4[x][y][z].re=prefactor1*out4[x][y][z].re;
		out4[x][y][z].im=prefactor1*out4[x][y][z].im;
		out5[x][y][z].re=prefactor1*out5[x][y][z].re;
		out5[x][y][z].im=prefactor1*out5[x][y][z].im;
		out6[x][y][z].re=prefactor1*out6[x][y][z].re;
		out6[x][y][z].im=prefactor1*out6[x][y][z].im;
		out7[x][y][z].re=prefactor1*out7[x][y][z].re;
		out7[x][y][z].im=prefactor1*out7[x][y][z].im;
		out8[x][y][z].re=prefactor1*out8[x][y][z].re;
		out8[x][y][z].im=prefactor1*out8[x][y][z].im;
		
		out21[x][y][z].re=prefactor2*out21[x][y][z].re;
		out21[x][y][z].im=prefactor2*out21[x][y][z].im;
		out22[x][y][z].re=prefactor2*out22[x][y][z].re;
		out22[x][y][z].im=prefactor2*out22[x][y][z].im;
		out23[x][y][z].re=prefactor2*out23[x][y][z].re;
		out23[x][y][z].im=prefactor2*out23[x][y][z].im;
		out24[x][y][z].re=prefactor2*out24[x][y][z].re;
		out24[x][y][z].im=prefactor2*out24[x][y][z].im;
		out25[x][y][z].re=prefactor2*out25[x][y][z].re;
		out25[x][y][z].im=prefactor2*out25[x][y][z].im;
		out26[x][y][z].re=prefactor2*out26[x][y][z].re;
		out26[x][y][z].im=prefactor2*out26[x][y][z].im;
		out27[x][y][z].re=prefactor2*out27[x][y][z].re;
		out27[x][y][z].im=prefactor2*out27[x][y][z].im;
		out28[x][y][z].re=prefactor2*out28[x][y][z].re;
		out28[x][y][z].im=prefactor2*out28[x][y][z].im;
#ifdef _INFO_INIT
		cout << out1[x][y][z].re << "   " << out2[x][y][z].re << "   " << out3[x][y][z].re << "   " << out4[x][y][z].re << endl;
#endif
    }
	
	fftwnd_destroy_plan(planen);
	
	return;
}


/*Solve for E-fields by En_x = d_n \chi_x, \chi_k=-rho_k/k^2. Requires som Fourier tranforms*/

void find_rho()
{
	int k,l,m;
	double prefactor,prefactor1,prefactor2,rho_1,rho_2,rho_3;
	double Elec_1, Elec_2, Elec_3;
	int xup, xdn, yup, ydn, zup, zdn;
	int nx, ny, nz;
	double k_squared;
	
	nx = ny = nz = N;
	prefactor1 = 2.0*deltat*betah1t/betagt;
	prefactor2 = 2.0*deltat*betah2t/betagt;
	
	rho_1 = rho_2 = rho_3 = 0.0;
	Elec_1 = Elec_2 = Elec_3 = 0.0;
	
	/* Get charge density with proper coefficients */
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {  
		rho1[k][l][m].re=prefactor1*(out2[k][l][m].re*out7[k][l][m].re
									 -out3[k][l][m].re*out6[k][l][m].re
									 +out4[k][l][m].re*out5[k][l][m].re
									 -out1[k][l][m].re*out8[k][l][m].re)+
		prefactor2*(out22[k][l][m].re*out27[k][l][m].re
					-out23[k][l][m].re*out26[k][l][m].re
					+out24[k][l][m].re*out25[k][l][m].re
					-out21[k][l][m].re*out28[k][l][m].re);
		
		rho1[k][l][m].im=0.0;
		
		rho2[k][l][m].re=prefactor1*(out1[k][l][m].re*out7[k][l][m].re
									 -out3[k][l][m].re*out5[k][l][m].re
									 +out2[k][l][m].re*out8[k][l][m].re
									 -out4[k][l][m].re*out6[k][l][m].re)+
		prefactor2*(out21[k][l][m].re*out27[k][l][m].re
					-out23[k][l][m].re*out25[k][l][m].re
					+out22[k][l][m].re*out28[k][l][m].re
					-out24[k][l][m].re*out26[k][l][m].re);
		
		rho2[k][l][m].im=0.0;
		
		rho3[k][l][m].re=prefactor1*(out2[k][l][m].re*out5[k][l][m].re
									 -out1[k][l][m].re*out6[k][l][m].re
									 +out3[k][l][m].re*out8[k][l][m].re
									 -out4[k][l][m].re*out7[k][l][m].re)+
		prefactor2*(out22[k][l][m].re*out25[k][l][m].re
					-out21[k][l][m].re*out26[k][l][m].re
					+out23[k][l][m].re*out28[k][l][m].re
					-out24[k][l][m].re*out27[k][l][m].re);
		
		rho3[k][l][m].im=0.0;
    }
	
	/* Check that total charge is zero */
	
	/*  for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
	 {
	 rho_1+=rho1[k][l][m].re;
	 rho_2+=rho2[k][l][m].re;
	 rho_3+=rho3[k][l][m].re;
	 }
	 printf("%.15e %.15e %.15e\n",rho_1,rho_2,rho_3);*/
	
	/* Get rho_k */
	
	planen = fftw3d_create_plan(nx,ny,nz,FFTW_FORWARD,FFTW_MEASURE|FFTW_IN_PLACE);
	
	fftwnd_one(planen, &rho1[0][0][0], NULL);
	fftwnd_one(planen, &rho2[0][0][0], NULL);
	fftwnd_one(planen, &rho3[0][0][0], NULL);
	
	fftwnd_destroy_plan(planen);
	
	prefactor = 1.0/sqrt((double)(N*N*N));
	
	/* Normalize and get chi_k */
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {  
		if((k==0)&&(l==0)&&(m==0))
		{
		}
		else
		{
			rho1[k][l][m].re=prefactor*rho1[k][l][m].re;
			rho1[k][l][m].im=prefactor*rho1[k][l][m].im;
			rho2[k][l][m].re=prefactor*rho2[k][l][m].re;
			rho2[k][l][m].im=prefactor*rho2[k][l][m].im;
			rho3[k][l][m].re=prefactor*rho3[k][l][m].re;
			rho3[k][l][m].im=prefactor*rho3[k][l][m].im;
			
			k_squared = 6.0-2.0*cos(2.0*PI*k/N)-2.0*cos(2.0*PI*l/N)-2.0*cos(2.0*PI*m/N);
			
			chi1[k][l][m].re=-rho1[k][l][m].re/k_squared;
			chi2[k][l][m].re=-rho2[k][l][m].re/k_squared;
			chi3[k][l][m].re=-rho3[k][l][m].re/k_squared;
			chi1[k][l][m].im=-rho1[k][l][m].im/k_squared;
			chi2[k][l][m].im=-rho2[k][l][m].im/k_squared;
			chi3[k][l][m].im=-rho3[k][l][m].im/k_squared;
		}
    }
	
	planen = fftw3d_create_plan(nx,ny,nz,FFTW_BACKWARD,FFTW_MEASURE|FFTW_IN_PLACE);
	
	/* Get chi_x */
	
	fftwnd_one(planen, &chi1[0][0][0], NULL);
	fftwnd_one(planen, &chi2[0][0][0], NULL);
	fftwnd_one(planen, &chi3[0][0][0], NULL);
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {  
		chi1[k][l][m].re=prefactor*chi1[k][l][m].re;
		chi1[k][l][m].im=prefactor*chi1[k][l][m].im;
		chi2[k][l][m].re=prefactor*chi2[k][l][m].re;
		chi2[k][l][m].im=prefactor*chi2[k][l][m].im;
		chi3[k][l][m].re=prefactor*chi3[k][l][m].re;
		chi3[k][l][m].im=prefactor*chi3[k][l][m].im;
    }
	
	fftwnd_destroy_plan(planen);
	
	/* And finally E_x */
	
	for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
    {  
		xup=(k+1)%N;;
		yup=(l+1)%N;;
		zup=(m+1)%N;
		
		Elec1_1[k][l][m]=-(chi1[xup][l][m].re-chi1[k][l][m].re);
		Elec2_1[k][l][m]=-(chi1[k][yup][m].re-chi1[k][l][m].re);
		Elec3_1[k][l][m]=-(chi1[k][l][zup].re-chi1[k][l][m].re);
		Elec1_2[k][l][m]=-(chi2[xup][l][m].re-chi2[k][l][m].re);
		Elec2_2[k][l][m]=-(chi2[k][yup][m].re-chi2[k][l][m].re);
		Elec3_2[k][l][m]=-(chi2[k][l][zup].re-chi2[k][l][m].re);
		Elec1_3[k][l][m]=-(chi3[xup][l][m].re-chi3[k][l][m].re);
		Elec2_3[k][l][m]=-(chi3[k][yup][m].re-chi3[k][l][m].re);
		Elec3_3[k][l][m]=-(chi3[k][l][zup].re-chi3[k][l][m].re);
		
    }
	
	/*Check that charge is zero in E-fields */
	
    for(k=0;k<N;k++)for(l=0;l<N;l++)for(m=0;m<N;m++)
	{  
		xdn=(N+k-1)%N;
		ydn=(N+l-1)%N;
		zdn=(N+m-1)%N;
		Elec_1 +=Elec1_1[k][l][m]-Elec1_1[xdn][l][m]
		+Elec2_1[k][l][m]-Elec2_1[k][ydn][m]
		+Elec3_1[k][l][m]-Elec3_1[k][l][zdn];
		Elec_2 +=Elec1_2[k][l][m]-Elec1_2[xdn][l][m]
		+Elec2_2[k][l][m]-Elec2_2[k][ydn][m]
		+Elec3_2[k][l][m]-Elec3_2[k][l][zdn];
		Elec_3 +=Elec1_3[k][l][m]-Elec1_3[xdn][l][m]
		+Elec2_3[k][l][m]-Elec2_3[k][ydn][m]
		+Elec3_3[k][l][m]-Elec3_3[k][l][zdn];
	}
#ifdef _INFO_INIT
    printf("%.15e %.15e %.15e\n",Elec_1,Elec_2,Elec_3); 
#endif 
	return;
}


void dump_config(string fname=INITIALFIELDS)
{
	
	double value1, value2, value3;
	int x,y,z,nx,ny,nz;
		
	nx = ny = nz = N;
	
	if((conffile = fopen(fname.c_str(),"wb"))==NULL) 
    {
		printf("cannot open file conf_su2.\n");
		exit(1);
    }
	
	/* Parameters characterize initial fields */
	
	if(fwrite(&deltat,sizeof deltat,1,conffile)!=1){printf("Write error14.\n");}
	if(fwrite(&avev,sizeof avev,1,conffile)!=1){printf("Write error14.\n");}
	if(fwrite(&warmup,sizeof warmup,1,conffile)!=1){printf("Write error14.\n");}
	if(fwrite(&theseed,sizeof theseed,1,conffile)!=1){printf("Write error14.\n");}
	
	if(fwrite(&nx,sizeof nx,1,conffile)!=1){printf("Write error14.\n");}
	if(fwrite(&ny,sizeof ny,1,conffile)!=1){printf("Write error14.\n");}
	if(fwrite(&nz,sizeof nz,1,conffile)!=1){printf("Write error14.\n");}
	
	for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
		
		/* Links are one */
		
		value1 = 1.0; value2 = 0.0;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		/* Efields as just calculated */
		
		value1 = Elec1_1[x][y][z]; value2 = Elec1_2[x][y][z]; value3 = Elec1_3[x][y][z];   
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value3,sizeof value3,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 = Elec2_1[x][y][z]; value2 = Elec2_2[x][y][z]; value3 = Elec2_3[x][y][z];   
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value3,sizeof value3,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 = Elec3_1[x][y][z]; value2 = Elec3_2[x][y][z]; value3 = Elec3_3[x][y][z];   
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value3,sizeof value3,1,conffile)!=1){printf("Write error14.\n");}
		
		/* Higgs fields */
		
		value1 =    out1[x][y][z].re; value2 =    out2[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 =    out3[x][y][z].re; value2 =    out4[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 =    out21[x][y][z].re; value2 =    out22[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 =    out23[x][y][z].re; value2 =    out24[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		/* Higgs momentum */
		
		value1 =    out5[x][y][z].re; value2 =    out6[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 =    out7[x][y][z].re; value2 =    out8[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}	 
		
		value1 =    out25[x][y][z].re; value2 =    out26[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}
		
		value1 =    out27[x][y][z].re; value2 =    out28[x][y][z].re;
		if(fwrite(&value1,sizeof value1,1,conffile)!=1){printf("Write error14.\n");}
		if(fwrite(&value2,sizeof value2,1,conffile)!=1){printf("Write error14.\n");}	 
		
	}
	
	fflush(conffile);
	fclose(conffile);
	
	return;
}

};