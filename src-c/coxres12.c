/* SCCS  @(#)coxres12.c	4.4  01/12/93      */
/*
** Does only the  integral(xbar(t) dM_i(t)) part of the score residuals
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of time and status values
**      strata  =1 for the last obs of each strata
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      method  ==1 for efron method
**
** Output
**      x is replaced with the appropriate value
**
** Scratch
**      scratch,  from which a and a2 are carved
*/
#include <stdio.h>
extern double **dmatrix();

void coxres12(nx, nvarx, y, covar2, strata, score, method, scratch)
long    nx[1],
	nvarx[1],
	*method,
	strata[];
double  y[],
	*covar2,
	*scratch,
	score[];

    {
    register int i,j, k;
    register double temp;
    int n, nvar;
    int deaths;
    double *time, *status;
    double *a, *a2;
    double denom, efron_wt;
    double **covar;
    double hazard;

    n = *nx;
    nvar  = *nvarx;
    time = y;
    status = y+n;
    a = scratch;
    a2 = a+nvar;
    /*
    **  Set up the ragged array
    */
    covar=  dmatrix(covar2, n, nvar);

    efron_wt=0;
    deaths=0;
    for (i=0; i<nvar; i++) a2[i] =0;
    strata[n-1] =1;  /*failsafe */
    for (i=n-1; i >=0; i--) {
	if (strata[i]==1) {
	    denom =0;
	    for (j=0; j<nvar; j++) a[j] =0;
	    }

	denom += score[i];
	if (status[i]==1) {
	    deaths++;
	    efron_wt += score[i];
	    for (j=0; j<nvar; j++) a2[j] += score[i]*covar[j][i];
	    }
	for (j=0; j<nvar; j++) {
	    a[j] += score[i] * covar[j][i];
	    covar[j][i] =0;
	    }

	if (deaths>0 && (i==0 || strata[i-1]==1 || time[i]!=time[i-1])){
	    /* last obs of a set of tied death times */
	    hazard = deaths/ denom;
	    temp = *method * (deaths-1) / 2.0;
	    for (j=0; j<nvar; j++)    /*create the weighted mean in a2 */
		a2[j] = (a[j]-temp*a2[j]) / (denom - temp*efron_wt);

	    for (k=i; k<n; k++) {
		if (time[k]==time[i]) {
		    for (j=0; j<nvar; j++)
			covar[j][k] = a2[j]*(status[k] - score[k]*hazard);
		    }
		else {
		    for (j=0; j<nvar; j++) covar[j][k] -= a2[j]*score[k]*hazard;
		    }
		}
	    efron_wt =0;
	    deaths =0;
	    for (j=0; j<nvar; j++)  a2[j] =0;
	    }
	}
    }
