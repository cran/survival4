/* SCCS @(#)agsurv2.c	4.7 05/23/94  */
/*
** Fit the survival curve, the special case of an Anderson-Gill style data
**   This program differs from survfit in several key ways:
**       Only returns data at the event times, not at censoring times
**       Fewer work arrays, but it is slower.
**
**   This is similar to survfit in that a complete curve is produced for
**        each strata.  If there are multiple 'subjects' in the newdata
**        list, then a matrix of survival curves is produced, with one
**        column for each input vector.
**
**  Input
**    n=# of subjects
**    nvar - number of vars in xmat -- will be zero if se is not desired
**             (the calling routine knows that xmat is only needed for the
**              correct second term of the se)
**    y - 3 column matrix containing strart, stop, event
**    score[n] - vector of weights
**    strata[n] - ==1 at the last obs of each strata
**    xmat   = data matrix that generated the Cox fit
**    varcov   = variance matrix of the coefs
**    nsurv    = the method 1=Kalbfleisch/Prentice  2= Tsiatis/Breslow
**                              3= Tsiatis, Efron approx
**
**    ncurve  = # of curves to produce
**    newx(nvar, ncurve) =  new subject x matrix
**    newrisk(ncurve)  = weights for the new subjects
**
** Output
**    surv  - the survival - of length ncurve*nsurv
**    varh  - the variance of the hazard function
**    nsurv - returned, number of survival time points
**    y[1,] - contains the survival times
**    y[2,] - the number of subjects at risk at that survival time
**    y[3,]  - the number of events at that time
**    strata[0]= # of strata, strata[1:n]= last obs strata 1,2, etc
**
**  Work
**    d[3*nvar]
**
**  Input must be sorted by (event before censor) within stop time within strata,
*/
#include <math.h>
double **dmatrix();

void agsurv2(sn, snvar, y, score, strata, surv, varh,
		  xmat, varcov, snsurv, d,
		  sncurve, newx, newrisk)
long *sn, *snvar;
long *snsurv, *sncurve;
long strata[];
double score[], y[], xmat[];
double varh[];
double surv[];
double d[], varcov[];
double newx[], newrisk[];
{
    register int i,j,k,l;
    double hazard, varhaz;
    double *start, *stop, *event;
    int n, nvar;
    int nsurv, method;
    int kk, psave;
    int deaths;
    double *a, *a2;
    int ncurve;
    double **covar,
	   **imat,
	   **covar2;
    int column,
	nrisk,
	nstrat,
	nsave,
	person;
    double time,
	   weight,
	   e_denom,
	   denom;
    double crisk,
	   guess, inc,
	   sumt,
	   km;
    double temp,
	   downwt,
	   d2;

    n = *sn;  nvar = *snvar;
    ncurve = *sncurve;
    method = *snsurv;
    start =y;
    stop  = y+n;
    event = y+n+n;
    a = d+nvar;
    a2 = a+nvar;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(xmat, n, nvar);
    imat = dmatrix(varcov,  nvar, nvar);
    covar2 = dmatrix(newx, ncurve, nvar);
    nsurv =0;
    nstrat =0;

    for (column=0; column<ncurve; column++) {
	crisk = newrisk[column];
	hazard  =0;
	varhaz  =0;
	km =1;
	for (i=0; i<nvar; i++) d[i] =0;
	nsave = nsurv;
	for (person=0; person<n;) {
	    if (event[person]==0) person++;
	    else {
		/*
		** compute the mean and denominator over the risk set
		*/
		denom =0;
		e_denom=0;
		for(i=0; i<nvar; i++){
		    a[i] =0;
		    a2[i]=0;
		    }
		time = stop[person];
		nrisk =0;
		deaths=0;
		for (k=person; k<n; k++) {
		    if (start[k] < time) {
			nrisk++;
			weight = score[k]/crisk;
			denom += weight;
			for (i=0; i<nvar; i++) {
			    a[i] += weight*(covar[i][k]- covar2[i][column]);
			    }
			 }
		    if (stop[k]==time && event[k]==1) {
			deaths++;
			e_denom += weight;
			for (i=0; i<nvar; i++) {
			    a2[i] += weight*(covar[i][k]- covar2[i][column]);
			    }
			}
		    if (strata[k]==1) break;
		    }

		/*
		** Add results all events at this time point
		*/
		psave = person;  /* for KM case below */
		temp =0;
		for (k=person; k<n && stop[k]==time; k++) {
		    if (event[k]==1) {
			kk =k ;      /*save for km case */
			if (method==3) downwt = temp++/deaths;
			else           downwt =0;
			d2 = (denom - downwt*e_denom);
			hazard += 1/d2;
			varhaz += 1/(d2*d2);
			for (i=0; i<nvar; i++)
			    d[i] += (a[i]- downwt*a2[i])/ (d2*d2);
			}
		    person++;
		    if (strata[k]==1) break;
		    }

		if (method==1) {
		    /*
		    ** kalbfleisch estimator is harder;
		    */
		    if (deaths ==1) {
			km *= pow(1- score[kk]/(crisk*denom), crisk/score[kk]);
			}
		    else {           /*find the zero of an equation */
			guess = .5;
			inc = .25;
			for (l=0; l<35; l++) { /* bisect it to death */
			    sumt =0;
			    for (k=psave; k<person; k++) {
				if (event[k] ==1) {
				    temp = score[k]/crisk;
				    sumt +=  temp/(1-pow(guess, temp));
				    }
				}
			    if (sumt < denom)  guess += inc;
				 else          guess -= inc;
			    inc = inc/2;
			    }
			km *= guess;
			}
		    surv[nsurv] = km;;
		    }
		else surv[nsurv] = exp(-hazard);

		temp =0;
		for (i=0; i<nvar; i++)
		    for (j=0; j< nvar; j++)
			temp += d[i]*d[j]*imat[i][j];
		varh[nsurv] = varhaz + temp;
		if (column==(ncurve-1)) {
		    /* on the last pass, I can overwrite old data with new */
		    i = nsurv - nsave;
		    start[i] = time;
		    stop[i] = nrisk;
		    event[i]= deaths;
		    }
		nsurv++;
		}

	    if (strata[person-1]==1) {
		if (column==(ncurve-1)) {
		    nstrat++;
		    strata[nstrat]= nsurv-nsave;
		    }
		km=1;
		hazard  =0;
		varhaz  =0;
		}
	    }
	}
    *snsurv = nsurv/ ncurve;
    strata[0] = nstrat;
    }



