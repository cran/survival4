/* SCCS @(#)coxfit_null.c	4.6 06/17/93  */
/*
** Special case: fit the "Null" model.  All that is needed are the loglik
**     and the residual  -- 90% of the work is the residual
**  the input parameters are
**
**       nused        :number of people
**       method       : ==1 for efron
**       time(n)      :time of event or censoring for person i
**       status(n)    :status for the ith person    1=dead , 0=censored
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       score(n)    :the risk score
**       weights(n)  :case weights
**
**  returned parameters
**       loglik       : the log-likelihood for the data
**       resid        : the martingale residual for each subject
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>

void coxfit_null( nusedx, method, time, status, score, weights, strata,
		  loglik, resid)

long    *nusedx,
	*method,
	strata[],
	status[];
double  score[],
	weights[],
	time[];
double  loglik[],  /* returned values */
	resid[];
{
    register int i,j;
    double deaths;
    int n;
    double  denom;
    double  hazard;
    double  e_denom, temp, downwt;
    int     lastone;
    double  meanwt;

    n = *nusedx;
    /*
    ** pass 1- resid will contain the risk sum
    */
    strata[n-1] =1;  /* just in case */
    for (i=n-1; i>=0; i--) {
	if (strata[i] == 1) denom = 0;
	denom += score[i] * weights[i];
	if (i==0 || strata[i-1]==1 || time[i-1]!=time[i])
	    resid[i] = denom;
	else resid[i] =0;
	}

    /*
    ** pass2, now do all the work
    */
    deaths=0;
    e_denom=0;
    hazard =0;
    meanwt =0;
    lastone = 0;
    loglik[0] =0;
    for (i= 0; i<n; i++) {
	if (resid[i]!=0) denom = resid[i];
	resid[i] = status[i];
	if (status[i]==1) {
	    deaths++;
	    e_denom += score[i] * weights[i];
	    meanwt += weights[i];
	    loglik[0] += weights[i] *log(score[i]);
	    }
	if (strata[i]==1 ||  time[i+1]!=time[i]) {
	    /*last subject of a set of tied times */
	    if (deaths<2 || *method==0) {
		hazard += meanwt/denom;
		loglik[0] -= meanwt * log(denom);
		for (j=lastone; j<=i; j++) resid[j] -= score[j]*hazard;
		}
	    else {
		temp = hazard;
		meanwt /= deaths;
		for (j=0; j<deaths; j++) {
		    downwt = j /deaths;
		    hazard +=  meanwt/(denom - e_denom* downwt);
		    temp   +=  meanwt*(1-downwt)/(denom - e_denom* downwt);
		    loglik[0] -=  meanwt *log(denom - e_denom* downwt);
		    }
		for (j=lastone; j<=i; j++) {
		    if (status[j]==0) resid[j] = -score[j]*hazard;
		    else  resid[j] -=  score[j]* temp;
		    }
		}
	    lastone =i +1;
	    deaths =0;
	    e_denom =0;
	    meanwt =0;
	    }
	}

    /* finish any "trailing" censored obs */
    for (j=lastone; j<n; j++)  resid[j] -= score[j]*hazard;
    }
