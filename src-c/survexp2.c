/*SCCS @(#)survexp2.c	4.3  08/10/92
** Calcultate expected survival for a single column of p's
**
** Input
**       ntime: number of output time points
**       times: the output time points, in days
**       dim:  the dimensions of the hazard table: nage, nsex, nyear
**       ages:  the values of the ages, in days, increasing.  ages[1]==0.
**       year :  the first date that this 'year' should be used, in days since 1/1/60
**       hazard: the hazard rate, per day: c style dimension is (nyear,2,nage)
**       n : number of subjects
**       entry: entry date for each subject, in days since 1/1/60
**       birth: birth date for each subject, in days since 1/1/60
**       sex:   1=male, 2=female
**         (The program can handle any value >0 here, up to second dim of haz)
**       nsurv: second dimension of surv.  Will be either 1 or n.  In the
**                latter case there is a separate set of expecteds per person
**       special: if !=0, the special case is run, where times has one value
**                  per subject, and I want each subject's expected at his t.
**
** Output
**        surv: the expected survival for the cohort, at the specified times.
*/
#include <math.h>
void survexp2(ntime, times, dim, ages, year, hazard, n,
			   entry, birth, sex, work, nsurv, special, surv)
long    *ntime, *nsurv,  *n;
long    *special;
long    times[];
long    year[];
long    ages[];
long    dim[];
long    entry[], birth[], sex[];
float   hazard[];
double  surv[];
double  work[];
    {
    int    i,j,k;
    int    si;          /* offset, for my own subscripting */
    int    iy;          /* current year, or row of the table */
    int    isurv;       /* indexing for the 'surv' matrix */
    int    dtime, stime;
    int    entage;
    int    nage, nsex, nyear;
    double temp;

    nage = dim[0];  nsex=dim[1]; nyear=dim[2];
    isurv=0;
    for (i=0; i< *ntime * *nsurv; i++) surv[i] =0;

    /*
    ** From here on down, it is essentially two separate programs.  One for
    **   for the special case, and one for the ordinary case.
    */
    if (*special) goto spcase;
    for (i=0; i< *n; i++) {
	for (j=0; j< *ntime; j++) work[j]=0;

	si = (sex[i]-1) * nage;
	iy =1;
	entage = entry[i] - birth[i];
	stime=0 ;               /*time since entry, for the surv intervals */
	dtime=0 ;               /*time since entry, for the age pointer    */
	/*
	** Logic - The interval being mapped into is always (stime, time[k]).
	**   The cumulative hazard for this interval, or some portion of this
	**   interval, will be added to work[k], and stime will be moved
	**   ahead.  Dtime marks how far ahead I can go before getting another
	**   hazard rate value from the hazard table.  At that point in time,
	**   the year is checked as well, and incremented if necessary.
	*/
	k=0;
	for(j=0 ; j< (nage-1); j++) {
	    dtime =  ages[j+1] - entage;  /*time since entry */
	    if (dtime <0) continue;
	    while(iy< nyear  &&  (birth[i]+ages[j]) >= year[iy]) {
		iy++;
		si += nsex* nage;
		}
	    for (; k< *ntime; k++) {
		if (times[k] > dtime) {
		    work[k] += hazard[j+si]*(dtime - stime);
		    stime=dtime;
		    break;
		    }
		else {
		    work[k] += hazard[j+si]*(times[k]-stime);
		    stime=times[k];
		    }
		}
	    }


	/*
	** finish it off-  times that are greater than the last dtime
	*/
	while( k < *ntime)  {
	    work[k] += hazard[j+si]*(times[k]-stime);
	    stime = times[k];
	    k++;
	    }

	/* Now add it in to the total survival   */
	temp =0 ;
	for (k=0; k< *ntime; k++) {
	    temp -= work[k];
	    surv[k+isurv] += exp(temp);   /* exp(-cumulative hazard) */
	    }
	if (*nsurv >1) isurv += *ntime;
	}

    if (*nsurv ==1)
	for (k=0; k< *ntime; k++)
	    surv[k] = surv[k]/ *n;
    return;

spcase:
    for (i=0; i< *n; i++) {
	si = (sex[i]-1) * nage;
	iy =1;
	entage = entry[i] - birth[i];
	stime=0 ;               /*time since entry, for the surv intervals */
	dtime=0 ;               /*time since entry, for the age pointer    */
	surv[i] =0;
	for(j=0 ; j< (nage-1); j++) {
	    dtime =  ages[j+1] - entage;  /*time since entry */
	    if (dtime <0) continue;
	    while(iy< nyear  &&  (birth[i]+ages[j]) >= year[iy]) {
		iy++;
		si += nsex* nage;
		}
	    if (times[i] > dtime) {
		surv[i] += hazard[j+si]*(dtime - stime);
		stime=dtime;
		}
	    else {
		surv[i] += hazard[j+si]*(times[i]-stime);
		break;
		}
	    }
	if (times[i] > dtime)  surv[i] +=  hazard[j+si]*(times[i]-stime);

	surv[i] = exp(-surv[i]);
	}
    }
