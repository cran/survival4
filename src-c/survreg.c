/* @(#)survreg.c	4.8 03/30/93  */
/*
** Fit one of several censored data distributions
**
** Input
**      maxiter - max # of iterations allowed
**      n       - number of subjects
**      nvar    - number of variables in the x matrix
**      y       - matrix of start, stop, event
**      ny      - # columns of y.  =3 if there is interval censored data
**      event   - 1=exact, 0= right censored, 2=left censored, 3=interval
**      covar   - covariates for patient i
**                   Note that S sends this in column major order
**      offset  - offset vector (usually 0)
**      beta    - initial values for the parameters
**      np      - number of "extra" parameters in the distribution.  Always
**                  >=1.
**      pars    - An array of length npar*2.  The first column contains
**                  initial parameter estimates, the second is =1 if this
**                  parameter is fixed, i.e. not part of the iteration.
**      eps     - tolerance for convergence.  Iteration continues until the
**                  relative change in the deviance is <= eps.
**      dist    -  1=extreme value, 2=logistic, 3=gaussian, 4=cauchy
**
**  Output
**      beta    - the final coef vector
**      maxiter - the number of iterations consumed
**      imat    - the information matrix
**      loglik  - the final deviance (log lik + constant)
**      flag    - success flag  0 =ok
**                              -1= did not converge
**      deriv   - 3 or 6 column array: terms of the deviance, first deriv
**                  wrt eta, and second deriv wrt eta.  If there are any
**                  interval censored obs, then the next 3 cols contain, for
**                  those observations only, d/sigma, d/sigma^2, and
**                  d/(sigma eta).
**
**  Work arrays
**      newbeta(nvar)- always contains the "next iteration"
**      u(nvar)      - first deriv of the loglik
**      savediag(nvar)- for rnewton
**    these are all passed as a single vector, and then broken out.
**
**  Trickery: in order to use rnewton without having to pass all the data
**    arrays through it and out the other side, I declare many things as
**    static external.  This means that all the routines in this text file
**    can see them, but no others (they might conflict with S internals!).
**    If you break this file in two, they won't work.
**
**  For the calculations involving sigma and an interval censored datum, I
**    can't calculate the sigma derivatives simply from the eta derivs.
**
**  To add a new distribution:
**              add a new "static void" declaration
**              add it to the "switch(*dist)" list, (2 places)
**              add the new subroutine to the bottom of the code, see
**                      logist_d as an example
*/
#include <math.h>
#include <stdio.h>
#define  PI     3.141592653589793
#define  SPI    2.506628274631001     /* sqrt(2*pi) */
#define  ROOT_2 1.414213562373095

static int debug =0;       /* normally set to zero */

double **dmatrix();
static void exvalue_d();
static void logistic_d();
static void gauss_d();
static void cauchy_d();
static void (*sreg_gg)();
void sreg_g();
void sreg_deriv();

static int    nvar, np;
static double **covar;
static double *time2, *time, *status;
static double *g, *dg, *ddg;
static double *dsig, *ddsig, *dsg;
static double *offset;
static double *parms, *pfixed;

void survreg(maxiter, nx, nvarx, y, ny, covar2, offset2,
		 beta, npx, parmsx, u, imatx, loglik, flag, eps,
		 deriv, dist)
long    *nx,
	*maxiter,
	*nvarx,
	*ny,
	*npx,
	*flag,
	*dist;
double  *covar2,
	*offset2,
	*parmsx,
	*deriv,
	*imatx;
double  beta[],
	*u,
	*y,
	loglik[2],
	*eps;
    {
    register int n, nvar2, i;
    int maxiter2;
    double *newbeta,
	   *savediag;
    double **imat;

    maxiter2 = *maxiter;
    n = *nx;
    nvar = *nvarx;
    np   = *npx;
    offset = offset2;
    parms  = parmsx;
    pfixed = parmsx + np;

    nvar2 = nvar;
    for (i=0; i<np; i++) if (pfixed[i]!=1) nvar2++;
    imat = dmatrix(imatx, nvar2, nvar2);
    covar = dmatrix(covar2, n, nvar);
    newbeta = u+nvar2;
    savediag= newbeta + nvar2;

    g = deriv;
    dg= g+n;
    ddg = dg +n;
    dsig= ddg +n;
    ddsig = dsig +n;
    dsg = ddsig +n;

    if (*ny==2) {
	time=y;
	status = y+n;
	}
    else {
	time=y;
	time2 = time + n;
	status = time2 +n;
	}

    switch(*dist) {
	case 1: sreg_gg = exvalue_d;  break;
	case 2: sreg_gg = logistic_d; break;
	case 3: sreg_gg = gauss_d;    break;
	case 4: sreg_gg = cauchy_d;   break;
	}

    *flag = rnewton(&maxiter2, n, nvar2, beta, u, imat, loglik, *eps,
		    sreg_g,  sreg_deriv, newbeta, savediag, debug);
    *maxiter = maxiter2;
    }


void sreg_g(n, nvar2, beta, loglik)
int n, nvar2;
double beta[];
double *loglik;
    {
    register int person, i,j;
    double  eta,
	    sigma;
    double  z, zu,
	    temp;
    double  sig2;
    static double  funs[4], ufun[4];

    /*
    ** calculate the first and second derivative wrt eta
    */
    if (pfixed[0]==1) sigma = exp(parms[0]);
    else sigma = exp(beta[nvar]);
    i= nvar;
    for (j=0; j<np; j++) if (pfixed[j]!=1) parms[j] = beta[i++];

    *loglik =0;
    sig2 = 1/(sigma*sigma);
    for (person=0; person<n; person++) {
	eta =0;
	for (i=0; i<nvar; i++) eta += beta[i] * covar[i][person];
	eta += offset[person];
	z = (time[person] - eta) /sigma;

	j = status[person];       /*convert to integer */
	switch(j) {
	    case 1:                             /* exact */
		(*sreg_gg)(z, funs,1);
		g[person] = log(funs[1])  - log(sigma);
		temp = funs[2]/sigma;
		dg[person] = -temp;
		ddg[person]= funs[3]*sig2 - temp*temp;
		break;
	    case 0:                             /* right censored */
		(*sreg_gg)(z, funs,2);
		g[person] = log(funs[1]);
		temp = funs[2]/(funs[1]*sigma);
		dg[person] = temp;
		ddg[person]= -funs[3]*sig2/funs[1] - temp*temp;
		break;
	    case 2:                             /* left censored */
		(*sreg_gg)(z, funs,2);
		g[person] = log(funs[0]);
		temp = funs[2]/(funs[0]*sigma);
		dg[person]= -temp;
		ddg[person]= funs[3]*sig2/funs[0] - temp*temp;
		break;
	    case 3:                             /* interval censored */
		zu = (time2[person] - eta)/sigma;  /*upper endpoint */
		(*sreg_gg)(z, funs, 2);
		(*sreg_gg)(zu,ufun ,2);
		if (z>0)  temp = funs[1] - ufun[1]; /*stop roundoff in tails*/
		else      temp = ufun[0] - funs[0];
		g[person] = log(temp);
		dg[person]  = -(ufun[2] - funs[2])/(temp*sigma);
		ddg[person] = (ufun[3] - funs[3])*sig2/temp -
						     dg[person]*dg[person];
		if (pfixed[0]==0) {
		    dsig[person] = (z*funs[2] - zu*ufun[2])/temp;
		    ddsig[person]= ((zu*zu*ufun[3] - z*z*funs[3])
				      +(zu*ufun[2] - z*funs[2])) /temp -
				    dsig[person]*dsig[person];
		    dsg[person] = ((zu*ufun[3] - z*funs[3])
				   + (ufun[2] - funs[2])) /(temp*sigma)  -
				      dsig[person]*dg[person];
		    }
		break;
	    }
	*loglik += g[person];
	}
    }

/*
**  This routine works for all the routines that have ony a scale
**   parameter as the "extra".
*/
void sreg_deriv(n, nvar2, beta, u, imat)
int     n;
int     nvar2;
double  beta[],
	u[],
	**imat;
    {
    register int i,j, person;
    double eta,
	   sz,
	   sigma;

    for (i=0; i<nvar2; i++) {
	u[i]=0;
	for (j=0; j<nvar2; j++)
	    imat[i][j]=0;
	}

    if (pfixed[0]==1) sigma = parms[0];
    else              sigma = exp(beta[nvar]);
    i= nvar;
    for (j=0; j<np; j++) if (pfixed[j]!=1) parms[j] = beta[i++];

    for (person=0; person<n; person++) {
	for (i=0; i<nvar; i++) {
	    u[i] += dg[person] * covar[i][person];
	    for (j=0; j<=i; j++) {
		imat[j][i] -= covar[i][person] *covar[j][person] *ddg[person];
		}
	    }

	if (pfixed[0]!=1) {   /* need derivative wrt log sigma */
	    eta =0;
	    for (i=0; i<nvar; i++) eta += beta[i] *covar[i][person];
	    eta += offset[person];
	    sz = (time[person] - eta);   /* sigma * z */

	    if (status[person]==3) {/* interval censored */
		u[nvar] += dsig[person];
		for (i=0; i<nvar; i++) {
		    imat[i][nvar] -= dsg[person] * covar[i][person];
		    }
		imat[nvar][nvar] -=  ddsig[person];
		imat[nvar][nvar] +=  dsig[person]*dsig[person];
		}
	    else {
		if (status[person]==1) u[nvar] += sz*dg[person] -1;
		else                   u[nvar] += sz*dg[person];
		for (i=0; i<nvar; i++) {
		    imat[i][nvar] -= (sz*ddg[person] - dg[person]) *
						covar[i][person];
		    }
		imat[nvar][nvar] -= sz*sz*ddg[person] - sz*dg[person];
		}
	    }
	}
    }

void survreg_g(nx, y, ny, eta, parms, deriv, ncol, dist)
long    *nx,
	*ny,
	*ncol,
	*dist;
double  eta[],
	*deriv;
double  parms[],
	*y;
    {
    register int n;
    register int person, j;
    double  sigma;
    double  z, zu,
	    temp;
    double  sig2;
    static double  funs[4], ufun[4];

    n = *nx;
    g = deriv;
    dg= g+n;
    ddg = dg +n;
    dsig= ddg +n;
    ddsig = dsig +n;
    dsg = ddsig +n;

    if (*ny==2) {
	time=y;
	status = y+n;
	}
    else {
	time=y;
	time2 = time + n;
	status = time2 +n;
	}

    switch(*dist) {
	case 1: sreg_gg = exvalue_d;  break;
	case 2: sreg_gg = logistic_d; break;
	case 3: sreg_gg = gauss_d;    break;
	case 4: sreg_gg = cauchy_d;   break;
	}

    /*
    ** calculate the first and second derivative wrt eta
    */
    sigma = exp(parms[0]);
    sig2 = 1/(sigma*sigma);
    for (person=0; person<n; person++) {
	z = (time[person] - eta[person]) /sigma;

	j = status[person];       /*convert to integer */
	switch(j) {
	    case 1:                             /* exact */
		(*sreg_gg)(z, funs,1);
		g[person] = log(funs[1])  - log(sigma);
		temp = -funs[2]/sigma;
		dg[person] = temp;
		ddg[person]= funs[3]*sig2 - temp*temp;
		if (*ncol==6) {
		    dsig[person] = temp*sigma*z -1;
		    dsg[person] =  sigma*z*ddg[person] - temp;
		    ddsig[person]= sigma*z* dsg[person];
		    }
		break;
	    case 0:                             /* right censored */
		(*sreg_gg)(z, funs,2);
		g[person] = log(funs[1]);
		temp = funs[2]/(funs[1]*sigma);
		dg[person] = temp;
		ddg[person]= -funs[3]*sig2/funs[1] - temp*temp;
		if (*ncol==6) {
		    dsig[person] = temp*sigma*z;
		    dsg[person] =  sigma*z*ddg[person] - temp;
		    ddsig[person]= sigma*z* dsg[person];
		    }
		break;
	    case 2:                             /* left censored */
		(*sreg_gg)(z, funs,2);
		g[person] = log(funs[0]);
		temp = funs[2]/(funs[0]*sigma);
		dg[person]= -temp;
		ddg[person]= funs[3]*sig2/funs[0] - temp*temp;
		if (*ncol==6) {
		    dsig[person] = temp*sigma*z;
		    dsg[person] =  sigma*z*ddg[person] - temp;
		    ddsig[person]= sigma*z* dsg[person];
		    }
		break;
	    case 3:                             /* interval censored */
		zu = (time2[person] - eta[person])/sigma;  /*upper endpoint */
		(*sreg_gg)(z, funs, 2);
		(*sreg_gg)(zu,ufun ,2);
		if (z>0)  temp = funs[1] - ufun[1]; /*stop roundoff in tails*/
		else      temp = ufun[0] - funs[0];
		g[person] = log(temp);
		dg[person]  = -(ufun[2] - funs[2])/(temp*sigma);
		ddg[person] = (ufun[3] - funs[3])*sig2/temp -
						     dg[person]*dg[person];
		if (*ncol==6) {
		    dsig[person] = (z*funs[2] - zu*ufun[2])/temp;
		    ddsig[person]= ((zu*zu*ufun[3] - z*z*funs[3])
				      +(zu*ufun[2] - z*funs[2])) /temp -
				    dsig[person]*dsig[person];
		    dsg[person] = ((zu*ufun[3] - z*funs[3])
				   + (ufun[2] - funs[2])) /(temp*sigma)  -
				      dsig[person]*dg[person];
		    }
		break;
	    }
	}
    }

/*
**  Case      ans[0]    ans[1]       ans[2]     ans[3]
**   1                    f          f'/f        f''/ f
**   2          F        1-F         f           f'
**
**  We do both F and 1-F to avoid the error in (1-F) for F near 1
*/

static void logistic_d(z, ans, j)
double z, ans[4];
int j;
    {
    double w, temp;
    int    sign, ii;

    /*
    ** The symmetry of the logistic allows me to be careful, and never take
    **  exp(large number).  This routine should be very accurate.
    */
    if (z>0)  {
	w = exp(-z);
	sign = -1;
	ii=0;
	}
    else {
	w = exp(z);
	sign = 1;
	ii=1;
	}
    temp = 1+w;
    switch(j) {
	case 1:  ans[1] = w/(temp*temp);
		 ans[2] = sign*(1-w)/temp;
		 ans[3] = (w*w -4*w +1)/(temp*temp);
		 break;
	case 2:  ans[1-ii] = w/temp;
		 ans[ii]   = 1/temp;
		 ans[2] = w/(temp*temp);
		 ans[3] = sign*ans[2]*(1-w)/temp;
		 break;
	}
    }

static void gauss_d(z, ans, j)
double z, ans[4];
int j;
    {
    double f;

    f = exp(-z*z/2) /SPI;
    switch(j) {
	case 1: ans[1] =f;
		ans[2] = -z;
		ans[3] = z*z -1;
		break;
	case 2: if (z>0) {
		    ans[0] = (1 + erf(z/ROOT_2))/2;
		    ans[1] =  erfc(z/ROOT_2) /2;
		    }
		else {
		    ans[1] = (1 + erf(-z/ROOT_2))/2;
		    ans[0] =  erfc(-z/ROOT_2) /2;
		    }
		ans[2] = f;
		ans[3] = -z*f;
		break;
	}
    }

/*
** In the Gaussian and logistic cases, I could avoid numeric disaster by only
**   evaluating exp(x) for x<0.  By symmetry, I could get what I need for
**   x >0.  The extreme value dist is asymmetric, and I don't yet see the
**   numeric tricks that I need.
** Probobly, a Taylor series will need to be used for large z.
*/

static void exvalue_d(z, ans, j)
double z, ans[4];
int j;
    {
    double temp;
    double w;
    w = exp(z);
    temp = exp(-w);
    switch(j) {
	case 1:  ans[1] = w*temp;
		 ans[2] = 1-w;
		 ans[3] = w*(w-3) +1;
		 break;
	case 2:  ans[0] = 1-temp;
		 ans[1] = temp;
		 ans[2] = w*temp;
		 ans[3] = w*temp*(1-w);
		 break;
	}
    }

static void cauchy_d(z, ans, j)
double z, ans[4];
int j;
    {
    double temp;

    temp = 1/(1 + z*z);
    switch(j) {
	case 1:  ans[1] = temp/PI;
		 ans[2] = -2*z*temp;
		 ans[3] = (6*z*z -2) * temp * temp;
		 break;
	case 2:  ans[0] = 0.5 + atan(z)/PI;
		 ans[1] = 1 - ans[0];
		 ans[2] = temp/PI;
		 ans[3] = -2*z*temp*temp/PI;
		 break;
	}
    }

