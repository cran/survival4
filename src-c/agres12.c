/* SCCS: @(#)agres12.c	1.5 01/12/93
/*
** Do the score residuals
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of start, stop, and event
**      strata  =1 for the last obs of each strata
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      method  ==1 for efron approx
**
** Output
**      resid2  matrix of score residuals, same "shape" of matrix as covar2
**
** Scratch
**      a       vector of length 3*nvar
*/
#include <stdio.h>
extern double **dmatrix();

void agres12(nx, nvarx, y, covar2, strata, score,
		method, resid2, a)
long    nx[1],
	nvarx[1],
	method[1],
	strata[];
double  y[],
	*covar2,
	score[],
	*a,
	*resid2;

    {
    register int i,k;
    register double temp;
    int n, nvar;
    int person;
    double denom, time;
    double *a2, *mean;
    double efron_wt;
    double hazard;
    int    deaths;
    double *start, *stop, *event;
    double **covar,
	   **resid;

    n = *nx;
    nvar  = *nvarx;
    start =y;
    stop  = y+n;
    event = y+(n+n);
    /*
    **  Set up the ragged arrays
    */
    covar=  dmatrix(covar2, n, nvar);
    resid = dmatrix(resid2, n, nvar);
    a2  = a+nvar;
    mean= a2 + nvar;

    for (person=0; person<n; ) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean over the risk set, also hazard at this time
	    */
	    denom =0;
	    efron_wt =0;
	    deaths =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		}
	    time = stop[person];
	    for (k=person; k<n; k++) {
		if (start[k] < time) {
		    denom += score[k];
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + score[k]*covar[i][k];
			}
		     if (stop[k]==time && event[k]==1) {
			deaths++;
			efron_wt += score[k];
			for (i=0; i<nvar; i++)
			    a2[i] = a2[i] + score[k]*covar[i][k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /* compute the weighted mean */
	    temp = *method * (deaths -1)/ 2.0;
	    for (i=0; i<nvar; i++)
		mean[i] = (a[i] - temp*a2[i]) / (denom - temp*efron_wt);

	    hazard = deaths/denom;

	    /* add it in for everyone in the risk set*/
	    for (k=person; k<n; k++) {
		if (start[k] < time) {
		    for (i=0; i<nvar; i++)
			resid[i][k] -= (covar[i][k] -mean[i])*score[k]*hazard;
		    if (stop[k]==time) {
			person++;
			if (event[k]==1)
			    for (i=0; i<nvar; i++)
				resid[i][k] += (covar[i][k] -mean[i]);
			}
		    }
		if (strata[k]==1) break;
		}
	    }
	}
    }
