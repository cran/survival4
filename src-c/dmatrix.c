#include <S.h>
/* SCCS @(#)dmatrix.c	4.4 08/25/92  */
/*
** set up ragged arrays, with #of columns and #of rows
*/
double **dmatrix(array, ncol, nrow)
double  *array;
int ncol, nrow;
    {
    register int i;
    register double **pointer;

    pointer = (double **) S_alloc(nrow, sizeof(double *));
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }
