/* SCCS @(#)rnewton.c	4.7 9/20/92 */
/*
** Ridge stabilized Newton iteration
**
** Input
**      maxiter     max # of iterations
**      n           number of obs (not used, just passed along)
**      nvar        number of parameters
**      beta        initial values of the parameters
**      eps         convergence criteria.  Iteration continues until the
**                     percent change in loglik is <= eps
**      dolk        pointer to a function that will return the loglik
**      doimat      pointer to a function that will compute the first
**                     derivative vector and -1*second deriv
**      debug       if =1, print out lots of messages as I go
**
** Output
**      beta        the coef vector
**      u           first derivative vector
**      imat        the information matrix at beta
**      loglik(2)   loglik at beta=start and beta=final
**      maxiter     actual # of iterations used
**
**  Return value
**      0 = ok
**      1 to nvar = non SPD matrix detected on variable __.
**              Because of the ridge stabilization, this only can happen if
**                      maxiter =0.
**      -1= failed to converge
**
**  Work arrays
**      newbeta(nvar)   always contains the "next iteration" of beta
**      savediag(nvar)  the diagonal of imat
**
**  Calls functions:  cholesky2, chsolve2
*/
#include <math.h>
#include <stdio.h>
#define POWER 2      /*how fast to increase or decrease the ridge parameter*/

int rnewton(maxiter, n, nvar, beta, u, imat, loglik, eps,
		    dolk, doimat, newbeta, savediag, debug)
int     *maxiter,
	n,
	nvar,
	debug;
double  u[],
	beta[],
	loglik[2],
	newbeta[],
	savediag[],
	eps;
double  **imat;

void    (*dolk)(),
	(*doimat)();
    {
    register int i,j;
    int     ierr, iter;
    double  newlk;
    int     halving;    /*are we doing step halving at the moment? */
    int     numhalf;    /* number of step halvings done so far  */
    int     levenberg;  /*is ridge stabilization turned on? */
    double  tau;        /*amount of ridge */

    /*
    ** do the initial iteration step
    */
    (*dolk)(n, nvar, beta, &newlk);
    loglik[0] = newlk;
    loglik[1] = newlk;

    if (debug>0) {
	fprintf(stderr, "\nCall to rnewton: LL=%f, coef= ", newlk);
	for (i=0; i<nvar; i++) fprintf(stderr, "%f ", beta[i]);
	fprintf(stderr, "\n");
	fflush(stderr);
	}

    (*doimat)(n, nvar, beta, u, imat);
    for (i=0; i<nvar; i++) savediag[i] = imat[i][i];
    tau =1;

    /*
    ** Update the betas and test for convergence
    */
    ierr = cholesky2(imat, nvar);
    if (ierr <nvar) {
	if (*maxiter==0) {   /* Bail out */
	    for (i=0; i<nvar; i++) {
		imat[i][i] = savediag[i];
		for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
		}
	    return(ierr);
	    }
	/*
	** Turn on ridge stabilization
	*/
	levenberg =1;
	do {
	    for (i=0; i<nvar; i++)
		   imat[i][i] = savediag[i] + fabs(savediag[i]) *tau;
	    ierr = cholesky2(imat, nvar);
	    tau *= POWER;
	    }  while (ierr <nvar);
	tau /= POWER;
	}
    else levenberg=0;

    chsolve2(imat, nvar, u);  /*replaced by u * inverse(imat) */
    for (i=0; i<nvar; i++)
	newbeta[i] = beta[i] + u[i];

    /*
    **  Never, never complain about convergence on the first step.  Thay way,
    **   if someone has to they can force one iter at a time.
    */
    if (*maxiter==0) {
	for (i=0; i<nvar; i++) {
	    imat[i][i] = savediag[i];
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	    }
	return(0);  /* and leave the old beta in peace */
	}

    /*
    ** here is the main loop
    */
    numhalf =0;
    for (iter=1; iter<=*maxiter; iter++) {
	halving=0;
	(*dolk)(n, nvar, newbeta, &newlk);
	if (debug>0) {
	    fprintf(stderr, "Iter %2d: LL=%f, tau=%f, coef= ",
				    iter, newlk, tau*levenberg);
	    for (i=0; i<nvar; i++) fprintf(stderr, "%f ", newbeta[i]);
	    fprintf(stderr, "\n");
	    fflush(stderr);
	    }

	if ((eps + newlk) <= loglik[1]) do {   /* step halving */
	    halving++;
	    if (numhalf++ > 20) {levenberg=1; tau *=POWER;}
	    if (halving > 20) {    /*stop an infinite loop */
		for (i=0; i<nvar; i++) newbeta[i] = beta[i];
		levenberg=1;
		tau *= POWER;
		goto retry;  /*this ain't working -- try ridging! */
		}
	    for (i=0; i<nvar; i++)
		newbeta[i] = (newbeta[i] + beta[i]) /2;
	    (*dolk)(n, nvar, newbeta, &newlk);

	    if (debug>0) {
		fprintf(stderr, "Halving! LL=%f, tau=%f, coef= ",
					 newlk, tau*levenberg);
		for (i=0; i<nvar; i++) fprintf(stderr, "%f ", newbeta[i]);
		fprintf(stderr, "\n");
		fflush(stderr);
		}
	    } while ( (eps + newlk) <= loglik[1]) ;
	tau /= 2*POWER;

	(*doimat)(n, nvar, newbeta, u, imat);
	for (i=0; i<nvar; i++)  savediag[i] = imat[i][i];

	/* Check out the information matrix */
	if (!levenberg) {
	    ierr = cholesky2(imat,nvar);
	    if (ierr <nvar) levenberg=1;
	    }
retry:  if (levenberg) {
	    do {
		for (i=0; i<nvar; i++)
		    imat[i][i] = savediag[i] + fabs(savediag[i])* tau;
		ierr = cholesky2(imat, nvar);
		tau *= POWER;
		if (debug>0 && ierr !=0)
			fprintf(stderr, "\tBump tau to %f\n", tau);
		}  while (ierr <nvar);
	    tau /= POWER;
	    }

	/* Am I done ? */
	if (halving==0 && fabs(1-(loglik[1]/newlk))<=eps) {
	    loglik[1] = newlk;
	    for (i=0; i<nvar; i++) {  /* restore imat */
		imat[i][i] = savediag[i];
		for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
		beta[i] = newbeta[i];
		}
	    *maxiter = iter;
	    return(0);
	    }

	/* Not done --- set up for another pass */
	loglik[1] = newlk;
	chsolve2(imat, nvar, u);
	for (i=0; i<nvar; i++) {
	    beta[i] = newbeta[i];
	    newbeta[i] += u[i];
	    }
	}

    /*
    **  Fell out the bottom-- didn't make it in the requested iterations
    */
    loglik[1] = newlk;
    for (i=0; i<nvar; i++) {  /* restore imat */
	imat[i][i] = savediag[i];
	for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	beta[i] = newbeta[i];
	}
    *maxiter = iter-1;
    return(-1);
    }
