/* SCCS @(#)chinv2.c	2.5 12/28/94  */
/*
** matrix inversion, given the cholesky decomposition
**
** input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**
** returned: the upper triangle will contain the inverse
**            below the diagonal will be junk
**
**  Terry Therneau
*/

chinv2(matrix ,n)
int  n;
double **matrix;
     {
     register double temp;
     register int i,j,k;

     /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
     for (i=0; i<n; i++){
	  if (matrix[i][i] >0) {
	      matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
	      for (j= (i+1); j<n; j++) {
		   matrix[j][i] = -matrix[j][i];
		   for (k=0; k<i; k++)     /*sweep operator */
			matrix[j][k] += matrix[j][i]*matrix[i][k];
		   }
	      }
	  }

     /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     */
     for (i=0; i<n; i++) {
	  if (matrix[i][i]==0) {  /* singular row */
		for (j=0; j<i; j++) matrix[j][i]=0;
		for (j=i; j<n; j++) matrix[i][j]=0;
		}
	  else {
	      for (j=(i+1); j<n; j++) {
		   temp = matrix[j][i]*matrix[j][j];
		   if (j!=i) matrix[i][j] = temp;
		   for (k=i; k<j; k++)
			matrix[i][k] += temp*matrix[j][k];
		   }
	      }
	  }
     }
