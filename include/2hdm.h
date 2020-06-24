/******************************  2hdm.h *********************************
*																		*
*  Header file for 2hdm simulation										*
*  Version 1, Dec 18 2011												*
*																		*
************************************************************************/

#define TWOHDM			"para/2hdm.cfg"
#define LATTICEATT		"latticeAtT"
#define FWRITE(fname, var,err) if(fwrite(&var,sizeof var,1,fname)!=1)printf(err);
#define FREAD(fname, var,err) if(fread(&var,sizeof var,1,fname)!=1)printf(err);

/* The following are global scalars. */
int nx,ny,nz;	/* Lattice dimensions. */
int volume1;	/* Volume of lattice = nx*ny*nz. */
double t, t_stop, t_save, ncs_current;
int n_save;
int nextTime, currentTime, prevTime;
double volume;
double running_time;

int warmup;
long theseed;								/* parameters in initial fields */

/* Parameters in the 2HDM */

double betah1t, betah1x, betah2t, betah2x, betagt, betagx, beta1r, beta2r, beta3, beta4, beta5R, beta5I, beta12R, beta12I;
double	beta31, beta32, beta21a1, beta21a2, beta231a1, beta231a2, beta51I, beta52I, beta121R, beta121I, beta122R, beta122I,betak1, betak2,// parameters for solving higgs fields at next time step
		betaE1, betaE2,betaEk,//parameter for solving E at current time
		betaGauss1, betaGauss2;//parameters for checking Gauss constraint
double	gg, lambda1, lambda2, lambda3, lambda4, lambda5R, lambda5I, am11sq, am22sq, am12sqR, am12sqI, av1, av2, theta12, kappa;
double	deltatsq, deltat, avev, vev;
double	ncs, ncsdotPrev, ncsdot;

/* The lattice is a single global variable. */

site *lattice;


