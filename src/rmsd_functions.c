//Purpose: rmsd_functions.c contains a bunch of functions that calculate rmsd and 
//         superimpose proteins.  

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
//Personal header file 
#include <structures.h>
#include <defines.h>

extern void u3b_(double *W, double *X, double *Y, int *N,int *MODE, double *RMS, double *U, double *T,int  *IER);

void superpose_v0(int N/*Number of mapped residues*/, CAcoord *atom1, CAcoord *atom2, int *map1, FILE *FID_debug)
{
  int i=0, j=0;
  char buffer[256];
  char *linePtr;
  double *x = (double*)calloc(3*N, sizeof(double));
  if(x==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for x array");
      exit(EXIT_FAILURE);
    }
  double *y = (double*)calloc(3*N, sizeof(double));
  if(y==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for y array");
      exit(EXIT_FAILURE);
    }
  FILE *CAconf1=fopen("CAcoordConf1.txt", "r");
  if(CAconf1==NULL)
    {
      fprintf(stderr, "CAcoordConf1.txt can not be opened\n");
      exit(EXIT_FAILURE);
    }
  
  FILE *CAconf2=fopen("CAcoordConf2.txt", "r");
  if(CAconf2==NULL)
    {
      fprintf(stderr, "CAcoordConf2.txt can not be opened\n");
      exit(EXIT_FAILURE);
    }

  for(i=0;i<N;i++)
    {
      linePtr=fgets(buffer,sizeof(buffer), CAconf1); 
      if(linePtr==NULL)
	{
	  fprintf(stderr, "fgets can not read lines of CA coordinates -1 \n");
	  exit(EXIT_FAILURE);
	}
      sscanf(buffer, "%lf\t%lf\t%lf\n",&x[3*i], &x[3*i+1] , &x[3*i+2]);
    }
  fclose(CAconf1);
  int unmapped=0;    

  for(i=0;i<N;i++)
    {
      linePtr=fgets(buffer,sizeof(buffer), CAconf2); 
      if(linePtr==NULL)
	{
	  fprintf(stderr, "fgets can not read lines of CA coordinates -2\n");
	  exit(EXIT_FAILURE);
	}
      sscanf(buffer, "%lf\t%lf\t%lf\n",&y[3*i], &y[3*i+1] , &y[3*i+2]);

      if( (y[3*i] ==INVALID_CRD) && (y[3*i+1] ==INVALID_CRD) && (y[3*i+2] ==INVALID_CRD))
	{
	  map1[i]=(-2);
	  fprintf(stdout, "This is unmapped%d\n", i);
	  unmapped++;
	}
      else
	{
	  map1[i]=i;
	}
    }
  fclose(CAconf2);
  fprintf(stdout, "Number of unmapped=%d\n", unmapped);
  int mappedNumber=N-unmapped;
  fprintf(stdout, "Number of all=%d\n", N);
  fprintf(stdout, "Number of mapped=%d\n", mappedNumber);
  double *a = (double*)calloc(3*(mappedNumber), sizeof(double));
  if(a==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for a array");
      exit(EXIT_FAILURE);
    }
  double *b = (double*)calloc(3*(mappedNumber), sizeof(double));
  if(b==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for b array");
      exit(EXIT_FAILURE);
    }
  
  int l=0;
  
  for (i=0;i<N;++i)
    {
      if(map1[i]!=(-2))
	{
	  a[3*l]  = x[3*i];
	  a[3*l+1]= x[3*i+1];
	  a[3*l+2]= x[3*i+2];
	  
	  b[3*l]  = y[3*i];
	  b[3*l+1]= y[3*i+1];
	  b[3*l+2]= y[3*i+2];
	  l++;
	}
    }

  /*Call fortran subroutine for superposition: u3b.f*/
  int ier=0, mode=1;
  double sum_wt=0, rms=999.0, w[mappedNumber], u1[9], t1[3];  
  
  for (i=0;i<mappedNumber;++i) 
    {
      w[i]=1.0;
      sum_wt+=w[i];
    }
  
  
  u3b_(w, b, a, &mappedNumber, &mode, &rms, u1, t1, &ier);
  rms=sqrt(rms/sum_wt);
  printf("RMSD:%lf\n", rms);
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
  
  /* for (j=0;j<9;++j) */
  /*   { */
  /*     //      fprintf(FID_debug, "u1[%d]=%lf\t",j, u1[j]); */
  /*     if((j+1)%3==0) */
  /* 	fprintf(FID_debug,"\n"); */
  /*   } */
  
  double U[3][3]; 
  for (i=0;i<3;++i)
    {
      for (j=0;j<3;++j)
	{
	  U[i][j]=u1[3*j+i];
	  //	  fprintf(FID_debug, "U[%d][%d]=%lf\t", i, j, U[i][j]);
	}
      //      fprintf(FID_debug, "\n");
    }
  
  double *rotatedB = (double*)calloc(3*mappedNumber, sizeof(double));
  if(rotatedB==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for rotatedB array");
      exit(EXIT_FAILURE);
    }
  
  for(i=0;i<mappedNumber;i++)
    for(j=0;j<3;j++)
      {
	rotatedB[3*i+j]=((U[j][0])*(b[3*i]) +  (U[j][1])*(b[3*i+1])+ (U[j][2])*( b[3*i+2]));
	rotatedB[3*i+j]=(t1[j]+ rotatedB[3*i+j]);
      }
  
  for (i=0;i<3;++i)
    fprintf(stdout, "T[%d]=%lf\n", i, t1[i]);
  
  FILE * f = fopen ("CAcoordConf2.txt", "w");
  if(f==NULL)
    {
      fprintf(stderr, "CAcoordConf2.txt can not be created\n");
      exit(EXIT_FAILURE);
    }
  l=0;
  FILE *residueInfo1=fopen("map1.pdb", "w");
  assert(residueInfo1!=NULL);
  FILE *residueInfo2=fopen("map2.pdb", "w");
  assert(residueInfo2!=NULL);
  
  for(i=0;i<N;i++)
    {
      if(map1[i]==-2)
      	{
      	  fprintf(f, "%.3lf\t%.3lf\t%.3lf\n", INVALID_CRD, INVALID_CRD, INVALID_CRD);
	  
	  fprintf(residueInfo1,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n" , (i+1), " CA ", \
		  atom1[i].residname,  atom1[i].chain, atom1[i].resNo, atom1[i].x,  atom1[i].y, atom1[i].z, 1.00, 0.00);
	  
	  fprintf(residueInfo2,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " CA ", \
		  atom2[i].residname,  atom2[i].chain, atom2[i].resNo, INVALID_CRD,  INVALID_CRD, INVALID_CRD, 1.00, 0.00); 
	  
      	}
      else
      	{
	  fprintf(f, "%.3lf\t%.3lf\t%.3lf\n", rotatedB[3*l], rotatedB[3*l+1], rotatedB[3*l+2]);
	  fprintf(residueInfo1,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n" , (i+1), " CA ", \
		  atom1[i].residname,  atom1[i].chain, atom1[i].resNo, atom1[i].x,  atom1[i].y, atom1[i].z, 1.00, 0.00); 
	  fprintf(residueInfo2,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " CA ", \
		  atom2[i].residname,  atom2[i].chain, atom2[i].resNo, rotatedB[3*l],  rotatedB[3*l+1], rotatedB[3*l+2], 1.00, 0.00); 

      	  l++;
      	}
    }
  
  fclose (f);

  fclose(residueInfo1);
  fclose(residueInfo2);
  free(rotatedB);
  free(x);
  free(y);
  free(a);
  free(b);
}

double *superpose_v1(int N/*Number of mapped residues*/, CAcoord *atomUpdated, double *resultVec, FILE *FID_debug)
{
  /*Purpose: This function superposes updated coordinates to 2nd conformation and */
  /*         return superposed coordinates as an 1D array.                        */
  int i=0, j=0;
  char buffer[256];
  char *linePtr;

  double *x = (double*)calloc(3*N, sizeof(double));
  if(x==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for x array");
      exit(EXIT_FAILURE);
    }
  double *y = (double*)calloc(3*N, sizeof(double));
  if(y==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for y array");
      exit(EXIT_FAILURE);
    }
  
  FILE *CAconf1=fopen("CAcoordConf1.txt", "r");
  if(CAconf1==NULL)
    {
      fprintf(stderr, "CAcoordConf1.txt can not be opened\n");
      exit(EXIT_FAILURE);
    }
  
  for(i=0;i<N;i++)
    {
      x[3*i]   = atomUpdated[i].x;
      x[3*i+1] = atomUpdated[i].y;
      x[3*i+2] = atomUpdated[i].z;
    }
  
  for(i=0;i<N;i++)
    {
      linePtr=fgets(buffer,sizeof(buffer), CAconf1); 
      if(linePtr==NULL)
	{
	  fprintf(stderr, "fgets can not read lines of CA coordinates\n");
	  exit(EXIT_FAILURE);
	}
      sscanf(buffer, "%lf\t%lf\t%lf\n",&y[3*i], &y[3*i+1] , &y[3*i+2]);
    }
  fclose(CAconf1);
  
  /*Call fortran subroutine for superposition: u3b.f*/
  int ier=0, mode=1;
  double sum_wt=0, rms=999.0, w[N], u1[9], t1[3];  
  
  for (i=0;i<N;++i) 
    {
      w[i]=1.0;
      sum_wt+=w[i];
    }
  
  
  u3b_(w, x, y, &N, &mode, &rms, u1, t1, &ier);
  rms=sqrt(rms/sum_wt);
  printf("RMSD:%lf\n", rms);
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
  
  /* for (j=0;j<9;++j) */
  /*   { */
  /*     fprintf(FID_debug, "u1[%d]=%lf\t",j, u1[j]); */
  /*     if((j+1)%3==0) */
  /* 	fprintf(stdout, "\n"); */
  /*   } */
  
  double U[3][3];
  for (i=0;i<3;++i)
    {
      for (j=0;j<3;++j)
	{
	  U[i][j]=u1[3*j+i];
	  //	  fprintf(FID_debug, "U[%d][%d]=%lf\t", i, j, U[i][j]);
	}
      //      fprintf(FID_debug, "\n");
    }
  
  for(i=0;i<N;i++)
    {
      for(j=0;j<3;j++)
	{
	  resultVec[3*i+j]=((U[j][0])*(x[3*i]) +  (U[j][1])*(x[3*i+1])+ (U[j][2])*( x[3*i+2]));
	  resultVec[3*i+j]=(t1[j]+ resultVec[3*i+j]);
	}
    }
  
  /* for (i=0;i<3;++i) */
  /*   fprintf(stdout, "T[%d]=%lf\n", i, t1[i]); */
  
  /*   for(i=0;i<N;i++) */
  /*     { */
  /*       for(j=0;j<3;j++) */
  /* 	   fprintf(stdout, "%.6lf\t", resultVec[3*i+j]); */
  
  /*       fprintf(stdout, "\n"); */
  /*     } */
  
  free(x);
  free(y);
  return resultVec;
}

double superpose(int N/*Number of mapped residues*/, CAcoord *atom1, CAcoord *atom2, int *map1, double *w/*Weight for RMSD*/)
{
  //  FILE *FID=fopen("");
  int i=0, j=0;
  char buffer[256];
  char *linePtr;
  double *x = (double*)calloc(3*N, sizeof(double));
  if(x==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for x array");
      exit(EXIT_FAILURE);
    }
  double *y = (double*)calloc(3*N, sizeof(double));
  if(y==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for y array");
      exit(EXIT_FAILURE);
    }
  FILE *CAconf1=fopen("CAcoordConf1.txt", "r");
  if(CAconf1==NULL)
    {
      fprintf(stderr, "CAcoordConf1.txt can not be opened\n");
      exit(EXIT_FAILURE);
    }
  
  FILE *CAconf2=fopen("CAcoordConf2.txt", "r");
  if(CAconf2==NULL)
    {
      fprintf(stderr, "CAcoordConf2.txt can not be opened\n");
      exit(EXIT_FAILURE);
    }

  for(i=0;i<N;i++)
    {
      linePtr=fgets(buffer,sizeof(buffer), CAconf1); 
      //      fprintf(stdout, "i=%d\n", i);
      if(linePtr==NULL)
	{
	  fprintf(stderr, "Hangi fgets can not read lines of CA coordinates\n");
	  exit(EXIT_FAILURE);
	}
      sscanf(buffer, "%lf\t%lf\t%lf\n",&x[3*i], &x[3*i+1] , &x[3*i+2]);
    }
  fclose(CAconf1);
  int unmapped=0;    

  for(i=0;i<N;i++)
    {
      linePtr=fgets(buffer,sizeof(buffer), CAconf2); 
      if(linePtr==NULL)
	{
	  fprintf(stderr, "fgets can not read lines of CA coordinates\n");
	  exit(EXIT_FAILURE);
	}
      sscanf(buffer, "%lf\t%lf\t%lf\n",&y[3*i], &y[3*i+1] , &y[3*i+2]);

      if( (y[3*i] ==INVALID_CRD) && (y[3*i+1] ==INVALID_CRD) && (y[3*i+2] ==INVALID_CRD))
	{
	  map1[i]=(-2);
	  fprintf(stdout, "This is unmapped%d\n", i);
	  unmapped++;
	}
      else
	{
	  map1[i]=i;
	}
    }
  fclose(CAconf2);
  fprintf(stdout, "Number of unmapped=%d\n", unmapped);
  int mappedNumber=N-unmapped;
  fprintf(stdout, "Number of all=%d\n", N);
  fprintf(stdout, "Number of mapped=%d\n", mappedNumber);
  double *a = (double*)calloc(3*(mappedNumber), sizeof(double));
  if(a==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for a array");
      exit(EXIT_FAILURE);
    }
  double *b = (double*)calloc(3*(mappedNumber), sizeof(double));
  if(b==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for b array");
      exit(EXIT_FAILURE);
    }
  
  int l=0;
  
  for (i=0;i<N;++i)
    {
      if(map1[i]!=(-2))
	{
	  a[3*l]  = x[3*i];
	  a[3*l+1]= x[3*i+1];
	  a[3*l+2]= x[3*i+2];
	  
	  b[3*l]  = y[3*i];
	  b[3*l+1]= y[3*i+1];
	  b[3*l+2]= y[3*i+2];
	  l++;
	}
    }

  /*Call fortran subroutine for superposition: u3b.f*/
  int ier=0, mode=1;
  double sum_wt=0, rms=999.0, u1[9], t1[3]; // w[mappedNumber], 
  
  for (i=0;i<mappedNumber;++i) 
    {
      //      w[i]=1.0;
      sum_wt+=w[i];
    }
  
  
  u3b_(w, b, a, &mappedNumber, &mode, &rms, u1, t1, &ier);
  rms=sqrt(rms/sum_wt);
  printf("RMSD:%lf\n", rms);
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
  
/*   for (j=0;j<9;++j) */
/*     { */
/*       fprintf(FID_debug, "u1[%d]=%lf\t",j, u1[j]); */
/*       if((j+1)%3==0) */
/* 	fprintf(FID_debug,"\n"); */
/*     } */
  
  double U[3][3]; 
  for (i=0;i<3;++i)
    {
      for (j=0;j<3;++j)
	{
	  U[i][j]=u1[3*j+i];
	  //	  fprintf(FID_debug, "U[%d][%d]=%lf\t", i, j, U[i][j]);
	}
      //      fprintf(FID_debug, "\n");
    }
  
  double *rotatedB = (double*)calloc(3*mappedNumber, sizeof(double));
  if(rotatedB==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for rotatedB array");
      exit(EXIT_FAILURE);
    }
  
  for(i=0;i<mappedNumber;i++)
    for(j=0;j<3;j++)
      {
	rotatedB[3*i+j]=((U[j][0])*(b[3*i]) +  (U[j][1])*(b[3*i+1])+ (U[j][2])*( b[3*i+2]));
	rotatedB[3*i+j]=(t1[j]+ rotatedB[3*i+j]);
      }
  
  for (i=0;i<3;++i)
    fprintf(stdout, "T[%d]=%lf\n", i, t1[i]);
  
  FILE * f = fopen ("CAcoordConf2.txt", "w");
  if(f==NULL)
    {
      fprintf(stderr, "CAcoordConf2.txt can not be created\n");
      exit(EXIT_FAILURE);
    }
  l=0;
  FILE *residueInfo1=fopen("map1.pdb", "w");
  assert(residueInfo1!=NULL);
  FILE *residueInfo2=fopen("map2.pdb", "w");
  assert(residueInfo2!=NULL);
  
  for(i=0;i<N;i++)
    {
      if(map1[i]==-2)
      	{
      	  fprintf(f, "%.3lf\t%.3lf\t%.3lf\n", INVALID_CRD, INVALID_CRD, INVALID_CRD);
	  
	  fprintf(residueInfo1,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n" , (i+1), " CA ", \
		  atom1[i].residname,  atom1[i].chain, atom1[i].resNo, atom1[i].x,  atom1[i].y, atom1[i].z, 1.00, 0.00);
	  
	  fprintf(residueInfo2,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " CA ", \
		  atom2[i].residname,  atom2[i].chain, atom2[i].resNo, INVALID_CRD,  INVALID_CRD, INVALID_CRD, 1.00, 0.00); 
	  
      	}
      else
      	{
	  fprintf(f, "%.3lf\t%.3lf\t%.3lf\n", rotatedB[3*l], rotatedB[3*l+1], rotatedB[3*l+2]);
	  fprintf(residueInfo1,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n" , (i+1), " CA ", \
		  atom1[i].residname,  atom1[i].chain, atom1[i].resNo, atom1[i].x,  atom1[i].y, atom1[i].z, 1.00, 0.00); 
	  fprintf(residueInfo2,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " CA ", \
		  atom2[i].residname,  atom2[i].chain, atom2[i].resNo, rotatedB[3*l],  rotatedB[3*l+1], rotatedB[3*l+2], 1.00, 0.00); 

      	  l++;
      	}
    }
  
  fclose (f);

  fclose(residueInfo1);
  fclose(residueInfo2);
  free(rotatedB);
  free(x);
  free(y);
  free(a);
  free(b);

  return rms;
}
double calculateRMS_v1(int N/*System size*/, double *CAcoordSet1, CAcoord *atomUpdated, FILE *FID_debug)
{
  
  double *CAcoordSet2 = (double*) calloc (3*N,sizeof(double));
  if (CAcoordSet2==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }
  int i=0;
  
  for(i=0;i<N;i++)
    {
      CAcoordSet2[3*i]   = atomUpdated[i].x;
      CAcoordSet2[3*i+1] = atomUpdated[i].y;
      CAcoordSet2[3*i+2] = atomUpdated[i].z;
    }
  
  //--Call Fortran routine-------------------    
  int k=0, ier=0, mode=0;
  double rms=0, sum_wt=0, w[N], u1[9], t1[3];
  for (k=0;k<N;++k)
    {
      w[k]=1.0;
      sum_wt+=w[k];
    }
  
  u3b_(w, CAcoordSet1, CAcoordSet2, &N, &mode, &rms, u1, t1, &ier);
  
  rms=sqrt(rms/sum_wt);
  fprintf(FID_debug, "RMSD:%lf\n", rms);
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
  //--End of Fortran routine call                                                                                                                                                  
  free(CAcoordSet2);
  
  return rms;
}

double calculateRMS_v2(int N/*System size*/, CAcoord *atomSet1, CAcoord *atomSet2, double *w/*Weight for each atom*/)
{
  int i=0;  
  double *CAcoordSet1 = (double*)calloc (3*N,sizeof(double));
  if (CAcoordSet1==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }
  
  double *CAcoordSet2 = (double*)calloc (3*N,sizeof(double));
  if (CAcoordSet2==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }

  
  for(i=0;i<N;i++)
    {
      CAcoordSet1[3*i]   = atomSet1[i].x;
      CAcoordSet1[3*i+1] = atomSet1[i].y;
      CAcoordSet1[3*i+2] = atomSet1[i].z;

      CAcoordSet2[3*i]   = atomSet2[i].x;
      CAcoordSet2[3*i+1] = atomSet2[i].y;
      CAcoordSet2[3*i+2] = atomSet2[i].z;
    }
  
  //--Call Fortran routine-------------------    
  int    ier=0, mode=0;
  double rms=0, sum_wt=0, u1[9], t1[3];
  for(i=0;i<N;i++)
    {
      //      w[i]=1.0;
      sum_wt+=w[i];
    }
  
  u3b_(w, CAcoordSet1, CAcoordSet2, &N, &mode, &rms, u1, t1, &ier);
  
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    {
    fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
    exit(EXIT_FAILURE);
    }
  //--End of Fortran routine call 

  rms=sqrt(rms/sum_wt);
  //  fprintf(FID_debug, "RMSD:%lf\n", rms);

  free(CAcoordSet1);
  free(CAcoordSet2);
  
  return rms;
}





double *superpose_v2(int N/*Number of mapped residues*/, CAcoord *atomUpdated, double *resultVec, double *w)
{
  /*Purpose: This function superposes updated coordinates to 2nd conformation and */
  /*         return superposed coordinates as an 1D array.                        */
  int i=0, j=0;
  char buffer[256];
  char *linePtr;

  double *x = (double*)calloc(3*N, sizeof(double));
  if(x==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for x array");
      exit(EXIT_FAILURE);
    }
  double *y = (double*)calloc(3*N, sizeof(double));
  if(y==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for y array");
      exit(EXIT_FAILURE);
    }
  
  FILE *CAconf2=fopen("CAcoordConf2.txt", "r");
  if(CAconf2==NULL)
    {
      fprintf(stderr, "CAcoordConf2.txt can not be opened\n");
      exit(EXIT_FAILURE);
    }
  
  for(i=0;i<N;i++)
    {
      x[3*i]   = atomUpdated[i].x;
      x[3*i+1] = atomUpdated[i].y;
      x[3*i+2] = atomUpdated[i].z;
    }
  
  for(i=0;i<N;i++)
    {
      linePtr=fgets(buffer,sizeof(buffer), CAconf2); 
      if(linePtr==NULL)
	{
	  fprintf(stderr, "fgets can not read lines of CA coordinates\n");
	  exit(EXIT_FAILURE);
	}
      sscanf(buffer, "%lf\t%lf\t%lf\n",&y[3*i], &y[3*i+1] , &y[3*i+2]);
    }
  fclose(CAconf2);
  
  /*Call fortran subroutine for superposition: u3b.f*/
  int ier=0, mode=1;
  double sum_wt=0, rms=999.0,  u1[9], t1[3];  //w[N],
  
  for (i=0;i<N;++i) 
    {
      //      w[i]=1.0;
      sum_wt+=w[i];
    }
  
  
  u3b_(w, x, y, &N, &mode, &rms, u1, t1, &ier);
  rms=sqrt(rms/sum_wt);
  fprintf(stdout, "RMSD:%lf\n", rms);
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
  
/*   for (j=0;j<9;++j) */
/*     { */
/*       fprintf(FID_debug, "u1[%d]=%lf\t",j, u1[j]); */
/*       if((j+1)%3==0) */
/* 	fprintf(stdout, "\n"); */
/*     } */
  
  double U[3][3];
  for (i=0;i<3;++i)
    {
      for (j=0;j<3;++j)
	{
	  U[i][j]=u1[3*j+i];
	  //	  fprintf(FID_debug, "U[%d][%d]=%lf\t", i, j, U[i][j]);
	}
      //      fprintf(FID_debug, "\n");
    }
  
  for(i=0;i<N;i++)
    {
      for(j=0;j<3;j++)
	{
	  resultVec[3*i+j]=((U[j][0])*(x[3*i]) +  (U[j][1])*(x[3*i+1])+ (U[j][2])*( x[3*i+2]));
	  resultVec[3*i+j]=(t1[j]+ resultVec[3*i+j]);
	}
    }
  
  for (i=0;i<3;++i)
    fprintf(stdout, "T[%d]=%lf\n", i, t1[i]);
  
  /*   for(i=0;i<N;i++) */
  /*     { */
  /*       for(j=0;j<3;j++) */
  /* 	   fprintf(stdout, "%.6lf\t", resultVec[3*i+j]); */
  
  /*       fprintf(stdout, "\n"); */
  /*     } */
  
  free(x);
  free(y);
  return resultVec;
}

double superimpose(int N/*System size*/, CAcoord *atomSet1, CAcoord *atomSet2, double *w/*Weight for each atom*/, int mode)
{
  int i=0;  
  double *CAcoordSet1 = (double*)calloc (3*N,sizeof(double));
  if (CAcoordSet1==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }
  double *CAcoordSet2 = (double*)calloc (3*N,sizeof(double));
  if (CAcoordSet2==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }
  
  for(i=0;i<N;i++)
    {
      CAcoordSet1[3*i]   = atomSet1[i].x;
      CAcoordSet1[3*i+1] = atomSet1[i].y;
      CAcoordSet1[3*i+2] = atomSet1[i].z;
      
      CAcoordSet2[3*i]   = atomSet2[i].x;
      CAcoordSet2[3*i+1] = atomSet2[i].y;
      CAcoordSet2[3*i+2] = atomSet2[i].z;
    }
  
  //--Call Fortran routine-------------------    
  int    ier=0, j=0;
  double rms=0, sum_wt=0, u1[9], t1[3];
  for(i=0;i<N;i++)
    {
      sum_wt+=w[i];
    }
  
  u3b_(w, CAcoordSet1, CAcoordSet2, &N, &mode, &rms, u1, t1, &ier);
  
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    {
      fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
      exit(EXIT_FAILURE);
    }
  //--End of Fortran routine call 
  if(mode==1)
    {
      double U[3][3];
      for (i=0;i<3;++i)
	{
	  for (j=0;j<3;++j)
	    {
	      U[i][j]=u1[3*j+i];
	      //	  fprintf(FID_debug, "U[%d][%d]=%lf\t", i, j, U[i][j]);
	    }
	  //  fprintf(FID_debug, "\n");
	}
      
      for(i=0;i<N;i++)
	{
	  //After this point, you dont need second coordinate set. So, use it as temp storage!
	  for(j=0;j<3;j++)
	    {
	      CAcoordSet2[3*i+j]=((U[j][0])*(CAcoordSet1[3*i]) +  (U[j][1])*(CAcoordSet1[3*i+1])+ (U[j][2])*( CAcoordSet1[3*i+2]));
	      CAcoordSet2[3*i+j]=(t1[j]+ CAcoordSet2[3*i+j]);
	    }
	}
      
      for(i=0;i<N;i++)
	{
	  atomSet1[i].x=CAcoordSet2[3*i];
	  atomSet1[i].y=CAcoordSet2[3*i+1];
	  atomSet1[i].z=CAcoordSet2[3*i+2];
	}

      if(0)
	for (i=0;i<3;++i)
	  fprintf(stdout, "T[%d]=%lf\n", i, t1[i]);
    }
  rms=sqrt(rms/sum_wt);
  //  fprintf(FID_debug, "RMSD:%lf\n", rms);

  free(CAcoordSet1);
  free(CAcoordSet2);
  
  return rms;
}
double superimpose_allatoms(int N/*# of all atoms*/, pdb_v23 *info1, pdb_v23 *info2, double *w/*Weight for each atom*/, int mode)
{

  //Purpose and whats new: 
  //This function calculates rmsd over all atoms! In previous versions, it was generally just over CA atoms. 
  int i=0;  
  double *allAtomSet1 = (double*)calloc (3*N,sizeof(double));
  if (allAtomSet1==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }
  double *allAtomSet2 = (double*)calloc (3*N,sizeof(double));
  if (allAtomSet2==NULL) 
    {
      fprintf(stderr, "No memory allocation fo RMSD calculation\n");
      exit (EXIT_FAILURE);
    }
  
  for(i=0;i<N;i++)
    {
      allAtomSet1[3*i]   = info1[i].x;
      allAtomSet1[3*i+1] = info1[i].y;
      allAtomSet1[3*i+2] = info1[i].z;
      
      allAtomSet2[3*i]   = info2[i].x;
      allAtomSet2[3*i+1] = info2[i].y;
      allAtomSet2[3*i+2] = info2[i].z;
    }
  
  //--Call Fortran routine-------------------    
  int    ier=0, j=0;
  double rms=0, sum_wt=0, u1[9], t1[3];
  for(i=0;i<N;i++)
    {
      sum_wt+=w[i];
    }
  
  u3b_(w, allAtomSet1, allAtomSet2, &N, &mode, &rms, u1, t1, &ier);
  
  if (ier==-1)
    fprintf(stderr, "\nError: SUPERPOSITION IS NOT UNIQUE!\n");
  
  else if (ier==-2)
    {
      fprintf(stderr, "\nError: NO RESULT OBTAINED!\n");
      exit(EXIT_FAILURE);
    }
  //--End of Fortran routine call 
  if(mode==1)
    {
      double U[3][3];
      for (i=0;i<3;++i)
	{
	  for (j=0;j<3;++j)
	    {
	      U[i][j]=u1[3*j+i];
	      //	  fprintf(FID_debug, "U[%d][%d]=%lf\t", i, j, U[i][j]);
	    }
	  //  fprintf(FID_debug, "\n");
	}
      
      for(i=0;i<N;i++)
	{
	  //After this point, you dont need second coordinate set. So, use it as temp storage!
	  for(j=0;j<3;j++)
	    {
	      allAtomSet2[3*i+j]=((U[j][0])*(allAtomSet1[3*i]) +  (U[j][1])*(allAtomSet1[3*i+1])+ (U[j][2])*( allAtomSet1[3*i+2]));
	      allAtomSet2[3*i+j]=(t1[j]+ allAtomSet2[3*i+j]);
	    }
	}
      
      for(i=0;i<N;i++)
	{
	  info1[i].x=allAtomSet2[3*i];
	  info1[i].y=allAtomSet2[3*i+1];
	  info1[i].z=allAtomSet2[3*i+2];
	}

      if(0)
	for (i=0;i<3;++i)
	  fprintf(stdout, "T[%d]=%lf\n", i, t1[i]);
    }
  rms=sqrt(rms/sum_wt);
  //  fprintf(FID_debug, "RMSD:%lf\n", rms);

  free(allAtomSet1);
  free(allAtomSet2);
  
  return rms;
}
