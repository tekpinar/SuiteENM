//Purpose: lapack_blas_functions.c keeps all functions of Newton-Raphson minimization
//         procedures I developed during my PhD in one place.
//Author : Mustafa Tekpinar
//Date   : 15 December 2011
/* Copyright (C) 2012  Mustafa Tekpinar

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; 
version 2.1 of the License.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
02110-1301  USA

Additionally, you can obtain a copy of the licence from 
http://www.gnu.org/licenses/lgpl-2.1.html

*/

//Email: tekpinar@buffalo.edu

//History: 
//-makePosDefANDsolveLapack_v2() has been added to the functions. 23 January 2012.

#include <stdio.h>
//Wrap BLAS functions by Jeffrey Hafner
/* Double Vector Dot Product */
double ddot (int N, double *X, int INCX, double *Y, int INCY) 
{
  extern double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);
  return ddot_(&N, X, &INCX, Y, &INCY);
}

/* Ax=B solution of a symmetric positive definite matrix A by Cholesky*/
int dposv (char UPLO, int N, int NRHS, double *A, int LDA, double *B, int LDB)
{
  extern void dposv_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
  int INFO;
  dposv_ (&UPLO, &N, &NRHS, A, &LDA, B, &LDB, &INFO);
  return INFO;
}

void makePosDefANDsolveLapack(int N, double *hessFortran, double **hessC, double *delta_x, double *grad, double beta)
{
  unsigned i=0, j=0;
  //Make way for direct call to Lapack routine=================================//
  for(i=0; i<3*N*3*N; i++)
    hessFortran[i]=0.0;
  
  for(i=0; i<3*N; i++)
    for(j=0; j<=i; j++)
      hessFortran[3*N*i+j]=hessC[i][j];
  
  int info;
  int posDefCounter=0;
  //int dposv (char UPLO, int N, int NRHS, double *A, int LDA, double *B, int LDB)
  while(1)
    {
      for(i=0; i<3*N; i++)
	delta_x[i]=grad[i];
      info=dposv( 'U',  (3*N),  1, hessFortran,  (3*N), delta_x,  (3*N));
      if( info > 0 )
	{
	  //	  fprintf(stdout, "The leading minor of order %i is not positive definite;\n", info );
	  //	  fprintf(stdout, "the solution could not be computed.\n" );
	  //	  fprintf(stdout, "Warning: Not pos def!\n" );
	  for(i=0; i<3*N*3*N; i++)
	    hessFortran[i]=0.0;

	  for(i=0; i<3*N; i++)
	    for(j=0; j<=i; j++)
	      {
		if(i==j)
		  hessFortran[3*N*i+j]= (hessC[i][j] + beta);
		else
		  hessFortran[3*N*i+j]= hessC[i][j];
	      }
	  //	  beta=beta*10;

	  posDefCounter=+1;
	  if (posDefCounter==100)
	    {
	      fprintf(stderr, "ERROR: Can not make it pos def\n");
	      break;
	    }
	  beta=beta*(posDefCounter+1);
	  //?? exit(1);
	}
      
      if (info==0)
	{
	  //	  fprintf(stdout, "Finally, it is pos def now!\n");
	  break;
	}
    }
}


void makePosDefANDsolveLapack_v2(int N, double *hessFortran, double **hessC, double *delta_x, double *grad, double beta, int printDetails)
{

  //Purpose: The purpose is pretty clear: To convert a non positive definite matrix to a positive definite matrix
  //         by adding a small amount of 'double beta' to diagonal elements. After making it pos def, we solve
  //         Ax=B by using Cholesky factorization. The difference between previous version and this is in their way
  //         they increase 'double beta' if beta is not sufficient to make it positive definite!
  unsigned i=0, j=0;
  //Make way for direct call to Lapack routine=================================//
  for(i=0; i<3*N*3*N; i++)
    hessFortran[i]=0.0;
  
  for(i=0; i<3*N; i++)
    for(j=0; j<=i; j++)
      hessFortran[3*N*i+j]=hessC[i][j];
  
  int info;
  int posDefCounter=0;
  //int dposv (char UPLO, int N, int NRHS, double *A, int LDA, double *B, int LDB)
  while(1)
    {
      for(i=0; i<3*N; i++)
	delta_x[i]=grad[i];
      info=dposv( 'U',  (3*N),  1, hessFortran,  (3*N), delta_x,  (3*N));
      if( info > 0 )
	{
	  //	  fprintf(stdout, "The leading minor of order %i is not positive definite;\n", info );
	  //	  fprintf(stdout, "the solution could not be computed.\n" );
	  if(printDetails)
	    fprintf(stdout, "Warning: Not pos def!\n" );
	  for(i=0; i<3*N*3*N; i++)
	    hessFortran[i]=0.0;

	  for(i=0; i<3*N; i++)
	    for(j=0; j<=i; j++)
	      {
		if(i==j)
		  hessFortran[3*N*i+j]= (hessC[i][j] + beta);
		else
		  hessFortran[3*N*i+j]= hessC[i][j];
	      }
	  beta=beta*10;
	  posDefCounter=+1;
	  if (posDefCounter==100)
	    {
	      fprintf(stderr, "ERROR: Can not make it pos def\n");
	      break;
	    }

	  //?? exit(1);
	}
      
      if (info==0)
	{
	  if(printDetails)
	    fprintf(stdout, "Finally, it is pos def now!\n");
	  break;
	}
    }
}
