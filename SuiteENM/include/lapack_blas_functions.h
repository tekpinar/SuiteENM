//Purpose: To keep all functions of Newton-Raphson minimization
//         procedures I developed during my PhD in one place.
//Author : Mustafa Tekpinar
//Date   : 15 December 2011
//History: 
//-makePosDefANDsolveLapack_v2() has been added to the functions. 23 January 2012.

#ifndef  LAPACK_BLAS_FUNCTIONS_H
#define  LAPACK_BLAS_FUNCTIONS_H
double ddot (int N, double *X, int INCX, double *Y, int INCY);
int    dposv (char UPLO, int N, int NRHS, double *A, int LDA, double *B, int LDB);
void   makePosDefANDsolveLapack(int N, double *hessFortran, double **hessC, double *delta_x, double *grad, double beta);
void   makePosDefANDsolveLapack_v2(int N, double *hessFortran, double **hessC, double *delta_x, double *grad, double beta, \
				   int printDetails);
#endif
