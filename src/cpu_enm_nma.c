//Purpose: This program will calculate normal modes of proteins by 
//         using elastic network model. For the first time, I will 
//         use a Lennard-Jones (6-12) type of potential for normal mode
//         analysis.
//Date   : July 24, 2013.

/* Copyright (C) 2013  Mustafa Tekpinar

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
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>

//Personal header files
//#include <time_functions.h>
#define KBT 2.49434192
#define M_PI 3.14159265358979323846

typedef struct
{
  /*Purpose: To keep all 'ATOM' records data in a pdb file according to PDB FORMAT Version 2.3*/
  int    serial;                         //Atom serial number
  char   name[5];                        //Atom name
  char   altLoc;                         //Alternate location indicator
  char   resName[4];                     //Residue name: 
  char   chainID;                        //Chain identifier
  int    resSeq;                         //Residue sequence number
  char   iCode;                          //Code for insertion of residues
  double x;                              //Orthogonal coordinates for X in Angstroms
  double y;                              //Orthogonal coordinates for Y in Angstroms
  double z;                              //Orthogonal coordinates for Z in Angstroms
  double occupancy;                      //Occupancy 
  double tempFactor;                     //Temperature factor
  char   element[2];                     //Element symbol, right-justified
  char   charge[2];                      //Charge on the atom. I dont know why it is char?????
  int    selTag;
}pdb_v23;

typedef struct
{
  /*Purpose: This structure keeps CA information of a residue*/
  int    resNo;
  int    selTag; //This tag will be used to select a subset of CA coordinates. 
                 //Example: CA of a domain in protein. 
  char   residname[4];
  char   chain;
  double x;
  double y;
  double z;
  double occ;
  double beta;
  double mass;
}CAcoord;

int  *scanpdb(int *atomCount, char *inFileName)
{
/*Purpose: To determine the total number of atoms and residues for memory allocation*/
/*In future I may use it to keep protein chain and segment infos*/

  int i=0;
  int numberofResids=0;
  int numberofChains=1; /*Remember that if there is just one chain, you will not see any */
                        /*TER signal to count chain numbers.                             */
  int numberofHelices=0;
  int numberofSheets=0;
  char line[100];
  
  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }

  while(1)
    {
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{

	  //	  if(((strncmp(" CA", line+12, 3))==0))
	  if(((strncmp(" CA ", line+12, 4))==0 ) && ((line[16]==' ') || (line[16]=='A')  ) )
	    numberofResids+=1;
	  
	  i+=1; 
	}
      if((strncmp("TER   ", line, 6))==0)
	{
	  numberofChains+=1;
	}

      if((strncmp("HELIX ", line, 6))==0)
	{
	  numberofHelices+=1;
	}

      if((strncmp("SHEET ", line, 6))==0)
	{
	  numberofSheets+=1;
	}
    }
  atomCount[0]=i;               //atomCount[0] is the total number of atoms
  atomCount[1]=numberofResids;  //atomCount[1] is the total number of residues
  atomCount[2]=numberofChains;  //atomCount[2] is the total number of chains
  atomCount[3]=numberofHelices; //atomCount[3] is the total number of helices
  atomCount[4]=numberofSheets;  //atomCount[4] is the total number of sheets

  if((atomCount[0]==0) || (atomCount[1]==0))
    {
      fprintf(stderr, "Error: Can not read atoms and residues\n");
      exit(EXIT_FAILURE);
    }
  else
    {
      fprintf(stdout, "Total number of atoms:%d\n", atomCount[0]);
      fprintf(stdout, "Number of residues:%d\n",    atomCount[1]);
      fprintf(stdout, "Number of chains:%d\n",      atomCount[2]);
      fprintf(stdout, "Number of helices:%d\n",     atomCount[3]);
      fprintf(stdout, "Number of sheets:%d\n",      atomCount[4]);
    }
  fclose(pdbdata);
  
  return atomCount; 
}

void  readpdb_v23(pdb_v23 *info, CAcoord *atom, char *inFileName)
{
  //To read both waters and protein atoms in a constant column format
  int i=0;
  int numberofResid=0;
  
  char c_buffer[9]; //Stands for (c)haracter buffer!!!
  memset(c_buffer,'\0', 9);

  char line[100];
  memset(line,'\0', 100);

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  while(1)
    {
      //      line_ptr=;
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+6,   5);
	  info[i].serial=atoi(c_buffer);

	  memset(info[i].name,'\0', 5);
	  strncpy(info[i].name,      line+12,  4);

	  info[i].altLoc=line[16];

	  memset(info[i].resName,'\0', 4);
	  strncpy(info[i].resName,   line+17,  3);

	  info[i].chainID=line[21];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+22,  4);
	  info[i].resSeq=atoi(c_buffer);

	  info[i].iCode=line[26];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+30,  8);
	  info[i].x=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+38,  8);
	  info[i].y=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+46,  8);
	  info[i].z=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+54,  6);
	  info[i].occupancy=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+60,  6);
	  info[i].tempFactor=atof(c_buffer);

	  memset(info[i].element,'\0', 2);
	  strncpy(info[i].element,   line+76,  2);

	  memset(info[i].charge,'\0', 2);
	  strncpy(info[i].charge,    line+78,  2);

	  
	  //If there is an alternative location, I select the first one. This selection may not be healthy depending on the job one is working on!!	  
	  if(((strncmp(" CA ", line+12, 4))==0 ) && ((info[i].altLoc==' ') || (info[i].altLoc=='A')  ) )
	    //	  if((strncmp(" CA", line+12, 3))==0)
	    {
	      strncpy(atom[numberofResid].residname, info[i].resName, 3);
	      atom[numberofResid].residname[3]='\0'; //Add null character to the end of the string!
 	      atom[numberofResid].chain=info[i].chainID;
	      atom[numberofResid].resNo=info[i].resSeq;
	      atom[numberofResid].x=info[i].x;
	      atom[numberofResid].y=info[i].y;
	      atom[numberofResid].z=info[i].z;
	      atom[numberofResid].occ=info[i].occupancy;
	      atom[numberofResid].beta=info[i].tempFactor;
	      if(0)	
		fprintf(stdout,"ATOM   %d %s %c %s %c %d %c %lf %lf %lf %lf %lf %s %s\n", \
			info[i].serial, info[i].name, info[i].altLoc, info[i].resName, info[i].chainID, info[i].resSeq, info[i].iCode, 
			info[i].x, info[i].y, info[i].z, info[i].occupancy, info[i].tempFactor, info[i].element, info[i].charge);
	      
    
	      numberofResid+=1;
	    }
	  
	  i+=1; 
	}
    }
  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }
  fclose(pdbdata);
}

double distanceSquared(double a[3], double b[3])
{
  double dis_components[3]={0.0, 0.0, 0.0};
  dis_components[0]=(a[0]-b[0]);
  dis_components[1]=(a[1]-b[1]);
  dis_components[2]=(a[2]-b[2]);
    
  return ( (dis_components[0])*(dis_components[0]) + (dis_components[1])*(dis_components[1]) + (dis_components[2])*(dis_components[2]) );
}

/* double distanceSquared_v2(double a[3], double b[3], double R_cutoff_squared) */
/* { */
/*   double offset=10.0; */
/*   double dis_squared=0.0; */
/*   double dis_components[3]={0.0, 0.0, 0.0}; */
/*   disSquare_component[0]=(a[0]-b[0])*(a[0]-b[0]); */
/*   if(disSquare_component[0]>=R_cutoff_squared) */
/*     { */
/*     dis_squared+=(R_cutoff_squared+offset); */
/*     return dis_squared; */
/*     } */
/*   else */
/*     { */
/*       disSquare_component[1]=(a[1]-b[1])*(a[1]-b[1]); */
/*       dis_squared=(disSquare_component[0]+disSquare_component[1]); */
/*       if((dis_squared>=R_cutoff_squared) */
/*     } */

    
/*   return ( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) ); */

/* } */

double getMatElement_parabolicPotential(CAcoord *atomSet, int pa, int d1, int ir, int d2, double R_cutoff_squared)
{
	double	fac=1.0;
	double xyz1[3], xyz2[3];

	xyz1[0] = atomSet[pa].x;
	xyz1[1] = atomSet[pa].y;
	xyz1[2] = atomSet[pa].z;
	xyz2[0] = atomSet[ir].x;
	xyz2[1] = atomSet[ir].y;
	xyz2[2] = atomSet[ir].z;
	//	printf("xyz1[0]=%lf, xyz1[1]=%lf, xyz1[2]=%lf\n", xyz1[0], xyz1[1], xyz1[2]);
	//	printf("xyz2[0]=%lf, xyz2[1]=%lf, xyz2[2]=%lf\n", xyz2[0], xyz2[1], xyz2[2]);
	
	double dis_squared = distanceSquared(xyz1, xyz2); 
	// if(abs(pa-ir)==1)
	// {
	// fac=100;
	// }
	// else
	//   fac=1;

	if ( dis_squared >= R_cutoff_squared )
	return 0;
	else
	return (fac * ( xyz1[d1] - xyz2[d1] ) * ( xyz1[d2] - xyz2[d2] ) / (dis_squared));
}

void getHessian_parabolicPotential(int N/*System Size*/, CAcoord *atomSet, double **Hessian, double **forceConstantsMatrix, double R_cutoff, int printDetails)
{

  double R_cutoff_squared=R_cutoff*R_cutoff;
  int i=0, j=0, k=0, l=0;
  //First of all, zero all hessian elements!
  for(i=0; i<(3*N); i++)
    for(j=0;j<(3*N); j++)
      Hessian[i][j]=0.0;

  //Now, you can calculate all offdiagonal hessian blocks. 
  for(i=0; i<N; ++i)         
    for(j=0; j<i; ++j)    
      for(k=0; k<3; ++k)
	for(l=0; l<3; ++l) 
	  {
	    double offdiagonal = -(forceConstantsMatrix[i][j])*getMatElement_parabolicPotential(atomSet, i, k, j, l, R_cutoff_squared); 
	    Hessian[3*i+k][3*j+l] = offdiagonal;
	    //printf("Hessian[%d][%d]=%lf\n", (3*i+k), (3*j+l), Hessian[3*i+k][3*j+l]);
	    Hessian[3*j+l][3*i+k] = offdiagonal;
	    //printf("Hessian[%d][%d]=%lf\n",  (3*j+l), (3*i+k), Hessian[3*j+l][3*i+k]);
	  }	
  
  //Here, we are calculating diagonal blocks.
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      if(i!=j)
	for(k=0; k<3; ++k)
	  for(l=0; l<3; ++l)
	    {
	      Hessian[3*i+k][3*i+l] -= Hessian[3*i+k][3*j+l];
	      //printf("Hessian[%d][%d]=%lf\n",  (3*i+k), (3*i+l), Hessian[3*i+k][3*i+l]);
	    }
  
  fprintf(stdout, "\nCalculated hessian matrix.\n");
  if(printDetails==1)
    {
      FILE *MATRIX_FILE=fopen("matrix.txt", "w");
      if(MATRIX_FILE==NULL)
	{
	  fprintf(stderr, "ERROR: I can't create matrix.txt file!\n"); 
	  exit(EXIT_FAILURE);
	}
      for(i=0; i<3*N; i++)
	{
	  for(j=0; j<3*N; j++)
	    fprintf(MATRIX_FILE, "%lf\t", Hessian[i][j]);
	  fprintf(MATRIX_FILE, "\n");
	}
      fclose(MATRIX_FILE);
  }
}
void writeNMDformat(int N/*Number of CA atoms*/, CAcoord *atomSet, double *W, double *A, char *NMD_file_name, int num_of_modes, double scaleAmplitudes)
{
  //Purpose: NMWiz is a plugin for normal mode analysis implemented in VMD.
  //         NMD is file format for visualization of normal modes. NMWiz and NMD are described in
  //         Reference: "ProDy: Protein Dynamics Inferred from Theory and Experiments 
  //                    Bakan A, Meireles LM, Bahar I 2011 Bioinformatics 27(11):1575-1577"

/*   NMD Format: Taken from "http://csb.pitt.edu/ProDy/tutorials/nmwiz_tutorial/nmwiz.html#nmd-format"                                      */
/*   NMD files (extension .nmd) are plain text files that contain at least normal mode and coordinate data.                                 */
/*   In addition to PCA, EDA, NMA, ANM, or GNM data, arbitrary vectors can be stored in NMD files.                                          */
/*   Following data fields are recognized:                                                                                                  */
/*   coordinates: Coordinates must be provided in one line as a list of decimal numbers.                                                    */
/*                Number of atoms in the system is deduced from size of this data line.                                                     */
/*   mode       : Normal mode array. Each normal mode array must be provided in one line as a list of decimal numbers.                      */
/*                Mode array may be preceded by mode index and mode length (square root of variance or inverse frequency).                  */
/*   title      : A title for the dataset.                                                                                                  */
/*   names      : Atom names. Default is “CA” for all atoms.                                                                                */
/*   resnames   : Residue names. Default value is “GLY”.                                                                                    */
/*   chainids   : Chain identifiers. Default value is “A”.                                                                                  */
/*   resids     : Residue numbers. If this data line if not found, residue numbers are started from 1 and incremented by one for each atom. */
/*   betas      : Beta factors. Default value is 0 (zero). B-factors are used to color the protein representation.                          */

  int i=0, j=0;
  FILE *NMD_FILE=fopen(NMD_file_name, "w");
  if(NMD_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not open %s file!\n", NMD_file_name);
      exit(EXIT_FAILURE);
    }

  fprintf(NMD_FILE, "title %s\n", "test.pdb");
  fprintf(NMD_FILE, "names ");
  for(i=0; i<N; i++) fprintf(NMD_FILE, "CA ");

  fprintf(NMD_FILE, "\nresnames ");
  for(i=0; i<N; i++) fprintf(NMD_FILE, "%s ", atomSet[i].residname);

  fprintf(NMD_FILE, "\nchids ");
  for(i=0; i<N; i++) fprintf(NMD_FILE, "%c ", atomSet[i].chain);

  fprintf(NMD_FILE, "\nresnums ");
  for(i=0; i<N; i++) fprintf(NMD_FILE, "%d ", atomSet[i].resNo);

  fprintf(NMD_FILE, "\nbetas ");
  for(i=0; i<N; i++) fprintf(NMD_FILE, "%.2lf ", atomSet[i].beta);

  fprintf(NMD_FILE, "\ncoordinates ");
  for(i=0; i<N; i++) fprintf(NMD_FILE, "%.3lf %.3lf %.3lf ", atomSet[i].x, atomSet[i].y, atomSet[i].z);

  double coeff[3*N];
  double sum=0.0;

  for(i=0;i<(3*N);i++)
    {
      sum=0.0;
      for(j=0;j<(3*N);j++)
	{
	  sum+=A[3*N*i+j]*A[3*N*i+j];
	}
      coeff[i]=sqrt(N*sum)*scaleAmplitudes;
    }
  if( ((num_of_modes+6)>(3*N-6)) || ((num_of_modes)>(3*N-6)) ) 
    {
      fprintf(stderr, "WARNING: Number of normal modes can not exceed 3*N-6, where N is number of CA atoms!\n");
      fprintf(stderr, "WARNING: Number of normal modes has been reset to (3*N-6)!\n");
      num_of_modes=(3*N-6);
    }
  for(i=6;i<(num_of_modes+6);i++) 
    {
      fprintf(NMD_FILE, "\nmode %d %.5lf ", (i-5), W[i]);
      for(j=0;j<N;j++)
	{
	  fprintf(NMD_FILE, "%.3lf %.3lf %.3lf ", (coeff[i])*A[3*N*i+3*j], (coeff[i])*A[3*N*i+3*j+1], (coeff[i])*A[3*N*i+3*j+2]);
	}
    }
  
  fclose(NMD_FILE);
}


void write_nm_allatom2pdb(int N/*System size*/, pdb_v23 *info, double *W, double *A, char *nma_file, int num_of_modes, \
			  int printDetails, bool writeInitialPdb, double scaleAmplitudes)
{
  int i=0, j=0, k=0;

  double coeff[3*N];
  double sum=0.0;
  
  for(i=0;i<(3*N);i++)
    {
      sum=0.0;
      for(j=0;j<(3*N);j++)
	{
	  sum+=A[3*N*i+j]*A[3*N*i+j];
	}
      coeff[i]=sqrt(N/sum)*scaleAmplitudes;
      //      coeff[i]=sqrt(1.0/(W[i]*sum))*scaleAmplitudes;
      //      printf("coeff[%d]=%lf\n", i, coeff[i]);
    }
  
  FILE *NMA_FILE=fopen(nma_file, "w");
  if(NMA_FILE==NULL)
    {
      fprintf(stderr, "ERROR: I can not create %s file!", nma_file);
      exit(EXIT_FAILURE);
    }
  
  if(writeInitialPdb==true)
    {
      //Write original conformation we started with!
      fprintf(NMA_FILE,"MODEL\n");
      for(j=0;j<N;j++)
	{
	  while(info[k].resSeq==(j+1)) //If there is a missing residue, this will not work!!
	    {
	      fprintf(NMA_FILE, "ATOM  %5d %-5s%3s %1c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf  %10.1s\n", \
		      info[k].serial, info[k].name, info[k].resName, info[k].chainID, info[k].resSeq, \
		      info[k].x, info[k].y, info[k].z, info[k].occupancy, info[k].tempFactor, info[k].element); 
	      k++;
	    }
	}
      fprintf(NMA_FILE,"ENDMDL\n\n");
    }

  //Modes 0, 1, and 2 are translations. Modes 3, 4 and 5 are rotations. 
  //Therefore, we start to write from mode 6. 
  for(i=6;i<(num_of_modes+6);i++) 
    {
      k=0;
      fprintf(NMA_FILE,"MODEL\n");
      for(j=0;j<N;j++)
	{
	  while(info[k].resSeq==(j+1))
	    {
	      fprintf(NMA_FILE, "ATOM  %5d %-5s%3s %1c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf  %10.1s\n", \
		      info[k].serial, info[k].name, info[k].resName, info[k].chainID, info[k].resSeq, \
		      ( info[k].x + (coeff[i]*A[3*N*i+3*j])   ), \
		      ( info[k].y + (coeff[i]*A[3*N*i+3*j+1]) ), \
		      ( info[k].z + (coeff[i]*A[3*N*i+3*j+2]) ), \
		      info[k].occupancy, info[k].tempFactor, info[k].element); 
	      k++;
	    }
	}
      fprintf(NMA_FILE,"ENDMDL\n\n");
    }
  fclose(NMA_FILE);
}

void solveSymEigens(int N/*System size*/, double **hess_ENM, double *W, double *A, int printDetails)
{
  int i=0, j=0;
  for(i=0; i<(3*N); i++)
    for(j=0; j<(3*N); j++)
      {
	A[3*N*i+j]=hess_ENM[i][j];
	if(printDetails==1)
	  fprintf(stdout, "A[%d]=%lf\n", (3*N*i+j), A[3*N*i+j]);
      }
     
  char JOBZ = 'V';
  char UPLO = 'U';
  int  n = 3*N;
  int  INFO;
  int  LDA = n;
   
  double *WORK=(double *)malloc(1*sizeof(double));
  int    LWORK=-1;

  /* Workspace size computation */
  /* Call DSYEV with the same arguments as the "real' computation */
  dsyev_(&JOBZ, &UPLO, &n, A, &LDA, W, WORK, &LWORK, &INFO);  
  if(INFO!=0)
    printf("Error in workspace allocation");
  LWORK=(int)WORK[0];
  free(WORK);

  WORK=(double *)malloc(LWORK*sizeof(double));
  if(WORK==NULL)
    {
      fprintf(stderr, "ERROR: I can not allocate memory for workspace!\n");
      exit(EXIT_FAILURE);
    }
  /* Call DSYEV for the "real' computation with the optimal block size */
  dsyev_(&JOBZ, &UPLO, &n, A, &LDA, W, WORK, &LWORK, &INFO);
  if(INFO!=0)
    fprintf(stderr, "ERROR: Error in workspace allocation!\n");

  if(printDetails==1)
    {
      FILE *EIGENVAL_FILE=fopen("eigenval.txt", "w");
      if(EIGENVAL_FILE==NULL)
	{
	  fprintf(stderr, "ERROR: I can't create eigenval.txt file!\n"); 
	  exit(EXIT_FAILURE);
	}
      FILE *EIGENVEC_FILE=fopen("eigenvec.txt", "w");
      if(EIGENVEC_FILE==NULL)
	{
	  fprintf(stderr, "ERROR: I can't create eigenvec.txt file!\n"); 
	  exit(EXIT_FAILURE);
	}
      
      for (i = 0; i < n; i++)
	fprintf (EIGENVAL_FILE, "%.4lf\n",W[i]);
      for (i = 0; i < (n*n); i++)
	fprintf (EIGENVEC_FILE, "%.4lf\n", A[i]);
      
      fclose(EIGENVAL_FILE);
      fclose(EIGENVEC_FILE);
    }
  //  free(W);
  free(WORK);
}

double diffclock(clock_t clock1, clock_t clock2)
{
  double diffticks=clock1-clock2;
  double diffms=diffticks/CLOCKS_PER_SEC;
  return diffms;
}

double calculateIsotropicDisplacementSquared(int N, double *displacementSquared, double *A, int mode_no)
{
  int i=mode_no;
  int j=0;
  double totalAreaTheoretical=0.0;
  for(j=0;j<N;j++)
    {
      displacementSquared[j]=0.0;
      displacementSquared[j]=( (A[3*N*i+3*j]*A[3*N*i+3*j]) + (A[3*N*i+3*j+1]*A[3*N*i+3*j+1]) + (A[3*N*i+3*j+2]*A[3*N*i+3*j+2]) );
      totalAreaTheoretical+=displacementSquared[j];
    }
  return totalAreaTheoretical;
}
double betaFactorCrossCorrelation(int num_of_C_alphas, double *betaFactorsExperimental, double *betaFactorsTheoretical)
{
  int i=0;
  double normExp=0.0, normTheo=0.0, scalarProduct=0.0;

  for(i=0; i<num_of_C_alphas; i++)
    {
      normExp+=(betaFactorsExperimental[i]*betaFactorsExperimental[i]);
      normTheo+=(betaFactorsTheoretical[i]*betaFactorsTheoretical[i]);
      scalarProduct+=(betaFactorsExperimental[i]*betaFactorsTheoretical[i]);
    }
  fprintf(stdout, "CROSS CORELATION=%.3lf\n", acos(scalarProduct/sqrt(normTheo*normExp)));
  return(acos(scalarProduct/sqrt(normTheo*normExp)));
}

void calculateBetaFactors4CA(int N, CAcoord *atomSet, double *W, double *A, int num_modes, char *betaFactorFile, bool fit2Experimental)
{
  int i=0;
  int j=0;
  double displacementSquared[N];
  double averageDisplacementSquared[N];

  FILE *BETA_FILE=fopen(betaFactorFile, "w");
  if(BETA_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Beta factors file can not be produced!\n");
      exit(EXIT_FAILURE);
    }
  //Calculate normalized experimental beta factor.
  //Normalization is such that the areas enclosed by the curves are equal to unity. This method has been used in reference below:
  //"Vibrational dynamics of proteins: Significance of slow and fast modes in relation to function and stability, 
  //I. Bahar, A. R. Atilgan, M. C. Demirel, & B. Erman, Phys. Rev. Lett. 80, 2733-2736, 1998."
  
  //For area calculation, I will assume that each short side of histogram bar has 1 unit length. Long side has the value of beta factor. 
  //By using this method, I will calculate total area and then divide each value with the total area.

  //Now, we are calculating total area here!
  double totalAreaExperimental=0.0;
  double betaFactorsExperimental[N];
  double eightPISquaredby3=8.0*M_PI*M_PI/3.0*100;
  
  for(j=0;j<N;j++)
    {
      betaFactorsExperimental[j]=0.0;
      betaFactorsExperimental[j]=atomSet[j].beta;
      totalAreaExperimental+=atomSet[j].beta;
    }
  
  if( ((num_modes+6)>(3*N-6)) || ((num_modes)>(3*N-6)) ) 
    {
      fprintf(stderr, "WARNING: Number of normal modes can not exceed 3*N-6, where N is number of CA atoms!\n");
      num_modes=(3*N-6);
    }
  double totalAreaTheoretical[num_modes+1]; //the last plus one is for average coefficient.
  //double allDisplacementSquared[num_of_modes][N];
  //fprintf(stderr, "HERE I AM\n");
  //Here, I am trying to get the sum of squares of displacement for different modes!
  for(j=0;j<N;j++)
    {
      averageDisplacementSquared[j]=0.0;
    }


  for(i=6;i<(num_modes+6);i++) 
    {
      calculateIsotropicDisplacementSquared(N, displacementSquared, A, i);
      totalAreaTheoretical[i-6]=0.0;

      for(j=0;j<N;j++)
	{
	  totalAreaTheoretical[i-6]+=displacementSquared[j];
	    //averageDisplacementSquared[j]+=((1.0/pow(W[i], 0.5))*displacementSquared[j]);
      //averageDisplacementSquared[j]+=((3*KBT/W[i])*displacementSquared[j]);
      averageDisplacementSquared[j]+=((KBT/W[i])*displacementSquared[j]);
	  //allDisplacementSquared[i-6][j]=displacementSquared[j];
	}
    }

  //Here, I am averaging sum of squares of displacements for each CA atom and calculating total area under this curve!
  totalAreaTheoretical[num_modes]=0.0;
  for(j=0;j<N;j++)
    {
      
      //      averageDisplacementSquared[j]=averageDisplacementSquared[j]/((double)num_of_modes);
      totalAreaTheoretical[num_modes]+=averageDisplacementSquared[j];
    }
  
  for(j=0;j<N;j++)
    {
      fprintf(BETA_FILE, "%d\t%.2lf\t", atomSet[j].resNo, atomSet[j].beta );
      // for(i=6;i<(num_of_modes+6);i++) 
	// {
	//   fprintf(BETA_FILE, "%.2lf\t", allDisplacementSquared[i-6][j]*(totalAreaExperimental/totalAreaTheoretical[i-6]) );
	//   //	  fprintf(BETA_FILE, "%.2lf\t", allDisplacementSquared[i-6][j]);
	// }
    if(fit2Experimental)
      fprintf(BETA_FILE, "%.2lf\n", averageDisplacementSquared[j]*(totalAreaExperimental/totalAreaTheoretical[num_modes]));
    else
      fprintf(BETA_FILE, "%.2lf\n", eightPISquaredby3*averageDisplacementSquared[j]);
    }
   
  betaFactorCrossCorrelation(N, betaFactorsExperimental, averageDisplacementSquared);

/*   //Print results to file */
/*   for(j=0;j<N;j++) */
/*     { */
/*       //Not normalized, raw results! */
/*       //fprintf(BETA_FILE, "%d\t%.2lf\t%.2lf\n", atomSet[j].resNo, (atomSet[j].beta), (displacementSquared[j])); */

/*       //Normalized based on total area to be 1 unit. */
/*       //fprintf(BETA_FILE, "%d\t%.2lf\t%.2lf\n", atomSet[j].resNo, (atomSet[j].beta/totalAreaExperimental), (displacementSquared[j]/totalAreaTheoretical)); */

/*       //Scaled based on total area rate. */
/*   fprintf(BETA_FILE, "%d\t%.2lf\t%.2lf\n", atomSet[j].resNo, atomSet[j].beta, averageDisplacementSquared[j]*(totalAreaExperimental/totalAreaTheoretical) ); */
/*     } */
  fclose(BETA_FILE);
}
void calculateCrossCorrelationsCA(int N, double *lambdas, double *A, int num_modes, char *crossCorrelationsFile,\
                    double **forceConstantsMatrix, bool normalized, bool save_matrix_on)
{

  int i=0;
  int j=0;
  int k=0;
  int ind_3i=0;
  int ind_3j=0;
  // double totalAreaTheoretical=0.0;
  //double constant = 3.0*KBT;
  double constant = 100.0;
  //double R_i_dot_R_j = 0.0;
  double inv_eigval = 0.0;

  double **ccMatrix = (double**)malloc(N*sizeof(double*));
  if(ccMatrix==NULL)
    {
      fprintf(stderr, "Malloc cannot allocate memory for ccMatrix array");
      exit(EXIT_FAILURE);
    }
  
  for (i=0; i<N; i++)
    {
      ccMatrix[i] = (double*)calloc((N),sizeof(double));
      if(ccMatrix[i]==NULL)
        {
          fprintf(stderr, "Calloc cannot allocate memory for ccMatrix[i] array");
          exit(EXIT_FAILURE);
        }
    }
  ///////////////////////////////////////////////////////////////
  
  for(k=6; k<num_modes; k++)
  {  //# Start with the first nonzero mode, namely, k=6.
    //fprintf(stderr, "HERE I AM%d\n", k);   
    inv_eigval=(1.0/lambdas[k]);
    //
    for (i=0; i<N; i++)
    {    ind_3i=3*i;
        for (j=0; j<N; j++)
        {
            ind_3j=3*j;
            //#Do the calculation for one half of the matrix
            ccMatrix[i][j]+= (constant*(1.0/(forceConstantsMatrix[i][j]))*inv_eigval*( (A[3*N*k+ind_3i]*A[3*N*k+ind_3j]) + \
                                              (A[3*N*k+ind_3i+1]*A[3*N*k+ind_3j+1]) + \
                                              (A[3*N*k+ind_3i+2]*A[3*N*k+ind_3j+2])));
        }
    }
  
  }
  //fprintf(stderr, "HERE I AM ALSO%d\n", k);
  //fprintf(stderr, "HERE I AM\n"); 
  //# Complete the other half after all calculations
  // for (i=0; i<N; i++)
  //   for (j=(i+1); j< N; j++)
  //     ccMatrix[j][i]=ccMatrix[i][j];



  FILE *CCFILE=NULL;
  
  

  if(save_matrix_on==true)
  {
    if(normalized==true)
    {
      CCFILE=fopen("cross_corr_norm.txt", "w");
         //# Start with the first nonzero mode, namely, k=6.
      for (i=0; i<N; i++)
      {   
          for (j=0; j<N; j++)
          {
              fprintf(CCFILE, "%.6lf\t", ccMatrix[i][j]/(sqrt((ccMatrix[i][i])*(ccMatrix[j][j]) ))); ;
          }
          fprintf(CCFILE, "\n");
      }
   
    }
    else
    {
      CCFILE=fopen("cross_corr_non_norm.txt", "w");
      for (i=0; i<N; i++)
      {    
          for (j=0; j<N; j++)
          {
              fprintf(CCFILE, "%.6lf\t", ccMatrix[i][j]); ;
          }
          fprintf(CCFILE, "\n");
      }
    }
    if(CCFILE==NULL)
    {
      fprintf(stderr, "Could not create cross correlation matrix!\n");
      exit(EXIT_FAILURE);
    }
  
  }

  fclose(CCFILE);
  for (i=0; i<N; i++)
  {
    free(ccMatrix[i]);
  }
  free(ccMatrix);
  //return totalAreaTheoretical;
}


double aa2Mass(char *residueName)
{

  //Put the reference name here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Purpose: This function returns the mass of each residue. 
  //         'char *residueName' is upper case three letter code
  //         and return value is mass in atomic mass units. 
  if (strncmp(residueName,"ALA", 3)==0) return  (71.0779);     /* A */
  if (strncmp(residueName,"ARG", 3)==0) return (156.1857);     /* R */
  if (strncmp(residueName,"ASN", 3)==0) return (114.1026);     /* N */
  if (strncmp(residueName,"ASP", 3)==0) return (115.0874);     /* D */
  if (strncmp(residueName,"CYS", 3)==0) return (103.1429);     /* C */
  if (strncmp(residueName,"GLN", 3)==0) return (128.1292);     /* Q */
  if (strncmp(residueName,"GLU", 3)==0) return (129.1140);     /* E */
  if (strncmp(residueName,"GLY", 3)==0) return  (57.0513);     /* G */
  if (strncmp(residueName,"HIS", 3)==0) return (137.1393);     /* H */
  if (strncmp(residueName,"HSD", 3)==0) return (137.1393);     /* H */
  if (strncmp(residueName,"ILE", 3)==0) return (113.1576);     /* I */
  if (strncmp(residueName,"LEU", 3)==0) return (113.1576);     /* L */
  if (strncmp(residueName,"LYS", 3)==0) return (128.1723);     /* K */
  if (strncmp(residueName,"MET", 3)==0) return (131.1961);     /* M */
  if (strncmp(residueName,"PHE", 3)==0) return (147.1739);     /* F */
  if (strncmp(residueName,"PRO", 3)==0) return  (97.1152);     /* P */
  if (strncmp(residueName,"SER", 3)==0) return  (87.0773);     /* S */
  if (strncmp(residueName,"THR", 3)==0) return (101.1039);     /* T */
  if (strncmp(residueName,"TRP", 3)==0) return (186.2099);     /* W */
  if (strncmp(residueName,"TYR", 3)==0) return (163.1733);     /* Y */
  if (strncmp(residueName,"VAL", 3)==0) return  (99.1311);     /* V */
  //Up to this point, all the data was take from David Whitford's Proteins book.
  //Rest of the data for metals are taken from
  //http://www.iupac.org/publications/pac/2003/pdf/7508x1107.pdf
  if (strncmp(residueName," FE", 3)==0) return  (55.8450);     
  if (strncmp(residueName," MG", 3)==0) return  (24.3050);     
  if (strncmp(residueName," MN", 3)==0) return  (54.9380);     
  if (strncmp(residueName," ZN", 3)==0) return  (65.4090);     
  else
    {
      fprintf(stderr,"ERROR: %s is an unknown residue!\n", residueName);
      return(0.0);
      exit(EXIT_FAILURE);
    } 
}

void alanineMassScanning(int N, int residueNumber, char *mutation, CAcoord *atomSet)
{
  int i=0;
  for(i=0;i<N;i++)
    {
      if(i==residueNumber)
	{
	  atomSet[i].mass=aa2Mass(mutation);
	}
      else
	{
	  atomSet[i].mass=aa2Mass(atomSet[i].residname);
	}
    }
}

void applyMassWeighting(int N, double **hess_ENM, CAcoord *atomSet)
{
  int i=0, j=0, a=0, b=0;
  for(i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
	{
	  for(a=0; a<3; a++)
	    {
	      for(b=0; b<3; b++)
		{
		  //fprintf(stdout, "%.1lf ", sqrt((atomSet[i].mass)*(atomSet[j].mass)));
		  hess_ENM[3*i+a][3*j+b]=(hess_ENM[3*i+a][3*j+b])*(1.0/sqrt((atomSet[i].mass)*(atomSet[j].mass)));
		  //		      hess_ENM[3*i+a][3*j+b]=(hess_ENM[3*i+a][3*j+b])*(1.0/pow((atomSet[i].mass),1.0));
		  // hess_ENM[3*i+a][3*j+b]=(hess_ENM[3*i+a][3*j+b])*(1.0/sqrt((atomSet[i].mass)*(atomSet[j].mass)));
		  //      hess_ENM[3*i+a][3*j+b]=(hess_ENM[3*i+a][3*j+b])*(1.0/((atomSet[i].residname)*(atomSet[j].mass)));
		}	    
	    }	    
	}
      //fprintf(stdout,"\n");
    }
}

void revertMassWeighting(int N, double *A, CAcoord *atomSet, int printDetails)
{
  int i=0, j=0;
  for(i=0; i<(3*N); i++)
    for(j=0; j<(3*N); j++)
      {
	//Now, lets go back to cartesian coordinates from mass-weighted cartesian coordinates!
	A[3*N*i+j]=A[3*N*i+j]*(1.0/sqrt(atomSet[j/3].mass));
	//A[3*N*i+j]=A[3*N*i+j]*(  1.0/(atomSet[j/3].mass)  );
	//A[3*N*i+j]=A[3*N*i+j]*(1.0/(atomSet[j/3].mass)*(atomSet1[j/3].mass) );
	
	if(printDetails==1)
	  fprintf(stdout, "A[%d]=%lf\n", (3*N*i+j), A[3*N*i+j]);
      }
}
void usage()
{
    fputs("\nUsage: ./cpu_enm_nma.exe -i prt1.pdb -o normalModes.nmd -R 15.0 -n 10 -s 1.0 -m 0\n", stdout);

    fputs("-i: Name of input file. (File has to be in pdb format.)                \n\n", stdout);
    fputs("-o: Name of output file.                                                 \n", stdout);
    fputs("    Output file format is deduced implicitly based on extension of       \n", stdout);
    fputs("    output file name. There are two options: nmd or pdb. If output is    \n", stdout);
    fputs("    in nmd format, like normalModes.nmd, you can visualize it with       \n", stdout);
    fputs("    Normal Mode Wizard plugin of VMD program.                          \n\n", stdout);

    fputs("-R: Cutoff radius in Angstrom units (10 Angstrom is default value).      \n", stdout);
    fputs("    Cutoff radius generally take values between 7 and 20 Angstroms.    \n\n", stdout);

    fputs("-n: Number of normal modes to be calculated (10 is default value).       \n", stdout);
    fputs("    If n is set to '-1', 5% of total number of modes will be used.       \n", stdout);
    fputs("    for beta factor calculations. For details, look at reference:        \n", stdout);
    fputs("    http://dx.doi.org/10.1016/j.bpj.2010.03.027                        \n\n", stdout);

    fputs("-s: Scale eigenvectors by a double type value (1.0 is default value).  \n\n", stdout);

    fputs("-m: 0 means no mass weighting hessian.                                   \n", stdout);
    fputs("    1 means linear mass weighting.                                       \n", stdout);
    fputs("    2 means squared mass weighting.                                    \n\n", stdout);

    fputs("-b: Name of file for theoretical beta factors.                            \n", stdout);
    fputs("    The first column of this file is experimental CA beta factors.        \n", stdout);
    fputs("    The last column of this file is average theoretical CA beta factors.\n\n", stdout);

    fputs("-c: Name of file for dynamical cross-coreelations file.                            \n", stdout);
    fputs("    The first column of this file is experimental CA beta factors.        \n", stdout);
    fputs("    The last column of this file is average theoretical CA beta factors.\n\n", stdout);
}
void readForceConstantsMatrixHANM(char *fcFile, double **forceConstantsMatrix)
{
  //Read a custom force constants matrix.
  FILE *fcdata=fopen(fcFile,"r");
  if(fcdata==NULL)
  {
    fprintf(stderr, "No such file:%s\n", fcFile); 
    exit(EXIT_FAILURE); 
  }
  //////////////////////////////////////////////////////////
  //char c_buffer[20]; //Stands for (c)haracter buffer!!!
  char *c_buffer;
  //memset(c_buffer,'\0', 20);

  char line[100];
  memset(line,'\0', 100);

  char *array[3];
  //Just a conversion factor from kJ/mol-nm^2 to kcal/mol-Angstrom^2
  double conversionFactor=23.9E-2;

  while(1)
  {
    //      line_ptr=;
  if(fgets(line, sizeof(line),fcdata)==NULL) break;
  
  c_buffer=strtok(line, " ");
  int i=0;
  while(c_buffer != NULL)
    {

      array[i++]=c_buffer;
      
      c_buffer = strtok (NULL, " ");
    }
  //fprintf(stdout, "%s\t%s\t%s\n", array[0], array[1], array[2]);
    //forceConstantsMatrix[atoi(array[0])-1][atoi(array[1])-1]=conversionFactor*atof(array[2]);
    //forceConstantsMatrix[atoi(array[1])-1][atoi(array[0])-1]=conversionFactor*atof(array[2]);
    forceConstantsMatrix[atoi(array[0])-1][atoi(array[1])-1]=atof(array[2]);
    forceConstantsMatrix[atoi(array[1])-1][atoi(array[0])-1]=atof(array[2]);


  }

  fclose(fcdata);
}
int main (int argc, char **argv)
{
  fputs("================================================================\n", stdout);
  fputs(" ======  =     =   ==     ==       o     o  o       o      o    \n", stdout);
  fputs(" =       = =   =   = =   = =       o o   o  oo     oo     o o   \n", stdout);
  fputs(" ====    =  =  =   =  = =  =  ***  o  o  o  o o   o o    ooooo  \n", stdout);
  fputs(" =       =   = =   =   =   =       o   o o  o  o o  o   o     o \n", stdout);
  fputs(" ======  =     =   =       =       o     o  o   o   o  o       o\n", stdout);
  fputs("================================================================\n", stdout);
  fputs("            Protein Normal Mode Calculation Program             \n", stdout);
  fputs("                              by                                \n", stdout);
  fputs("                     TechPinar Scientific                     \n\n", stdout);

  int option=0;
  int num_modes=10;
  int mass_weight_type=0;
  double scaleAmplitudes=1.0;
  //I know, it is not a good idea to hard code file name length. However, I always use my code on Unix systems and 
  //maximum file name length is 255 bytes. 
  char  *inputfile=(char *)malloc(255*sizeof(char));
  char  *outputfile=(char *)malloc(255*sizeof(char));
  char  *betafile=(char *)malloc(255*sizeof(char));
  char  *crossCorrelationsFile=(char *)malloc(255*sizeof(char));
  char  *forceConstantsFile=(char *)malloc(255*sizeof(char));
  //Now, lets clean inside filenames.
  memset (inputfile,'\0',255);
  memset (outputfile,'\0',255);
  memset (betafile,'\0',255);
  memset (crossCorrelationsFile, '\0', 255);
  if( (inputfile==NULL) || (outputfile==NULL) || (betafile==NULL) || (crossCorrelationsFile==NULL) || (forceConstantsFile==NULL))
    {
      fprintf(stderr, "Warning:  Can not allocate space for file names!\n");
    }

  double R_cutoff=15.0;                    
  
  while ((option = getopt(argc, argv,"hi:o:R:n:s:m:b:c:f:")) != -1) 
    {
      switch (option) 
      {
      case 'i' : inputfile =optarg;
        break;
      case 'o' : outputfile=optarg;
        break;
      case 'R' : R_cutoff  = atof(optarg); 
        break;
      case 'n' : num_modes = atoi(optarg); 
        break;
      case 's' : scaleAmplitudes=atof(optarg); 
        break;
      case 'm' : mass_weight_type=atoi(optarg); 
        break;
      case 'b' : betafile=optarg; 
        break;
      case 'c' : crossCorrelationsFile=optarg; 
        break;
      case 'f' : forceConstantsFile=optarg; 
        break;
      case 'h' : usage(); 
        exit(EXIT_FAILURE);
        break;

      default: usage(); 
        exit(EXIT_FAILURE);
      }
    }
  //betafile="test.dat";

  if((R_cutoff<0.0))
    {
      fprintf(stderr, "ERROR: Cutoff radius can not be negative!\n");
      usage();
      exit(EXIT_FAILURE);
    }
  if(num_modes<-1)
    {
      fprintf(stderr, "ERROR: Wrong selection for number of modes!\n");
      usage();
      exit(EXIT_FAILURE);
    }

  if((mass_weight_type!=0) && (mass_weight_type!=1) && (mass_weight_type!=2))
    {
      fprintf(stderr, "ERROR: Mass weight type has to be 0, 1 or 2!\n");
      usage();
      exit(EXIT_FAILURE);
    }


  fprintf(stdout, "Input file     : %s            \n", inputfile);
  fprintf(stdout, "Output file    : %s            \n", outputfile);
  fprintf(stdout, "Beta file      : %s            \n", betafile);
  fprintf(stdout, "Cutoff radius  : %.1lf Angstrom\n", R_cutoff);
  char extension[4]={'\0', '\0', '\0', '\0'};
  int length=strlen(outputfile);
  extension[0]=outputfile[length-3];
  extension[1]=outputfile[length-2];
  extension[2]=outputfile[length-1];
  fprintf(stdout, "File extension : %s\n", extension);
  if((strncmp (extension, "pdb", 3)!=0) && (strncmp (extension, "nmd", 3)!=0))
    {
      fprintf(stderr, "ERROR: %s is an unknown file format!\n", extension);
      fprintf(stderr, "Select pdb or nmd as output file format!\n");
      usage();
      exit(EXIT_FAILURE);
    }

  fprintf(stdout, "\nStarting calculation of normal modes for  %s.....\n", inputfile);

  int atomCount1[5]={99999/*Atoms*/, 999/*Residues*/, 99/*Chains*/, 100/*Helices*/, 99/*Sheets*/};
  int i=0;//My dear counter. 

  FILE *pdbdata=fopen(inputfile,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "No such file: %s\n", inputfile); 
      exit(EXIT_FAILURE); 
    }
  
  scanpdb(atomCount1, inputfile);
  printf("Number of residues:%d\n",  atomCount1[1]);
  
  //--Get data from pdb file--
  pdb_v23 *info1 = (pdb_v23*)malloc((atomCount1[0])*sizeof(pdb_v23));
  if (info1 == NULL) 
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb\n");
      exit(EXIT_FAILURE);
    }
  
  CAcoord *atomSet1 = (CAcoord*)malloc((atomCount1[1])*sizeof(CAcoord));
  if (atomSet1 == NULL) 
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord\n");
      exit(EXIT_FAILURE);
    }
    
  readpdb_v23(info1, atomSet1, inputfile);
  
  int N=atomCount1[1];//Number residues

  //Since we know number of amino acids in protein now, we can determine 5% of normal modes now. 
  if(num_modes==(-1))
    {
      num_modes=(int)( round((3*N-6)*5.0/100.0) );
      fprintf(stdout, "\n\nNumber of modes: %d            \n\n", num_modes);
    }

  //Assign masses to each amino acid. 
  for(i=0;i<N;i++)
    {
      atomSet1[i].mass=aa2Mass(atomSet1[i].residname);
    }

  clock_t begin=clock();
  //  int k=0;
  //  for(k=0; k<N; k++)
    {

      //      alanineMassScanning(N, k, "ALA", atomSet1);

  double **hess_ENM = (double**)malloc(3*N*sizeof(double*));
  if(hess_ENM==NULL)
    {
      fprintf(stderr, "Malloc cannot allocate memory for hess_ENM array");
      exit(EXIT_FAILURE);
    }
      
  for (i=0; i<3*N; i++)
	  {
	  hess_ENM[i] = (double*)calloc((3*N),sizeof(double));
	  if(hess_ENM[i]==NULL)
	    {
	      fprintf(stderr, "Calloc cannot allocate memory for hess_ENM[i] array");
	      exit(EXIT_FAILURE);
	    }
	  }
      
  double* A = (double*)calloc((3*N*3*N), sizeof(double));
  if(A==NULL)
	{
	  fprintf(stderr, "ERROR: I can not allocate memory for eigenvectors!\n");
	  exit(EXIT_FAILURE);
	}
      double *W=(double *)calloc((3*N), sizeof(double));
      if(W==NULL)
	{
	  fprintf(stderr, "ERROR: I can not allocate memory for eigenvalues!\n");
	  exit(EXIT_FAILURE);
	}
      
  int    printDetails=0;       //Which means obviously no!
  //This is just a theoretical point but I need to deviate coordinates from equillibrium(pdb) coordinates.
  //So, I am deviating all distances 0.001 Angstrom here. For parabolic potential, it didn't not change results.
  //========================================================================================================
  /* double deviateDist=0.001;  */
  /* for(i=0; i<N; i++) */
  /* 	{ */
  /* 	  atomSet1[i].x=(atomSet1[i].x-deviateDist); */
  /* 	  atomSet1[i].y=(atomSet1[i].y-deviateDist); */
  /* 	  atomSet1[i].z=(atomSet1[i].z-deviateDist); */
  /* 	} */
  //========================================================================================================
  
      
      
  //========================================================================================================
  //Get force constants matrix
  double **forceConstantsMatrix = (double**)malloc(N*sizeof(double*));
  if(forceConstantsMatrix==NULL)
    {
      fprintf(stderr, "Malloc cannot allocate memory for forceConstantsMatrix array");
      exit(EXIT_FAILURE);
    }
  
  for (i=0; i<N; i++)
    {
      forceConstantsMatrix[i] = (double*)calloc((N),sizeof(double));
      if(forceConstantsMatrix[i]==NULL)
        {
          fprintf(stderr, "Calloc cannot allocate memory for forceConstantsMatrix[i] array");
          exit(EXIT_FAILURE);
        }
    }
  int j=0;
  //Set all diagonal force constants to 1.0. This is just a trieck to avoid nan values for on diagonal elements!
  for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++)
      {
        forceConstantsMatrix[i][j]=1.0;
        
      }
    }
  //Read force constants matrix
  readForceConstantsMatrixHANM("rc15/fc_all_final_25_1p5_1p0_40_100.xvg", forceConstantsMatrix);
  //========================================================================================================

      //Build hessian!
      //========================================================================================================
  getHessian_parabolicPotential(N, atomSet1, hess_ENM, forceConstantsMatrix, R_cutoff, printDetails);                            
      //========================================================================================================
      
      //Here, we are testing normal mode calculation and beta factor calculation with mass weighted coordinates!
  if(mass_weight_type)
	  {
	  applyMassWeighting(N, hess_ENM, atomSet1);
	  }
      
      //Solve eigenvalues and eigenvectors of hessian!
      //========================================================================================================
  solveSymEigens(N, hess_ENM, W, A, printDetails);
      //========================================================================================================
      
  if(mass_weight_type)
	  {
	  revertMassWeighting(N, A, atomSet1, printDetails);
	  }
      
      fprintf(stdout, "\nCalculated eigenvalues and eigenvectors.\n\n");
      //========================================================================================================
      
      //Unshift CA coordinates back to original values!
      //========================================================================================================
      /* for(i=0; i<N; i++) */
      /* 	{ */
      /* 	  atomSet1[i].x=(atomSet1[i].x+deviateDist); */
      /* 	  atomSet1[i].y=(atomSet1[i].y+deviateDist); */
      /* 	  atomSet1[i].z=(atomSet1[i].z+deviateDist); */
      /* 	} */
      //========================================================================================================

      bool fit2Experimental = false;
      calculateBetaFactors4CA(N, atomSet1, W, A, num_modes, betafile, fit2Experimental);
      calculateCrossCorrelationsCA(N, W, A, num_modes, crossCorrelationsFile, forceConstantsMatrix, true, true);
      if(strncmp (extension, "pdb", 3)==0)
      {
        //Write results to an all atom file!
        //========================================================================================================
        bool writeInitialPdb=true;
        write_nm_allatom2pdb(N, info1, W, A, outputfile, num_modes, printDetails, writeInitialPdb, scaleAmplitudes);
        //========================================================================================================
      }
        if(strncmp (extension, "nmd", 3)==0)
          {
            //Write results in NMD format so that NMWiz can analyze them.
            //========================================================================================================
            writeNMDformat(N, atomSet1, W, A, outputfile, num_modes, scaleAmplitudes);
            //========================================================================================================
          }
      free(W);
      free(A);
      for (i=0; i<3*N; i++)
	      {
          free(hess_ENM[i]);
        }
      free(hess_ENM);
      for (i=0; i<N; i++)
        {
          free(forceConstantsMatrix[i]);
        }
      free(forceConstantsMatrix);
    }

  free(atomSet1);
  free(info1);  
  clock_t end=clock();
  fprintf(stdout, "Calculation finished succesfully within %lf s.\n", diffclock(end, begin));
  return 0;
}
