//Purpose: aa_form_factors.c contains functions related to coarse-grained
//         residue form factors. For some of them, we calculated them from a NMR structure set.
//         Some of the data has been taken from literature. However, all of the source code
//         belongs to the author. 
//Author: Mustafa Tekpinar
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
#include <string.h>
#include <math.h>
#include <assert.h>
//Personal header files.
#include <structures.h>
#include <form_factors.h>
#include <distance_functions.h>
void readYangFF(double F_CG[21][100], char *residDataFile)
{

/*
S=0  //Serine        (SER)
F=1  //Phenylalanine (PHE)
T=2  //Threonine     (THR)
N=3  //Asparagine    (ASN)
K=4  //Lysine        (LYS)
E=5  //Glutamic Acid (GLU)
Y=6  //Tyrosine      (TYR)
V=7  //Valine        (VAL)
Q=8  //Glutamine     (GLN)
M=9  //Methionine    (MET)
C=10 //Cysteine      (CYS)
L=11 //Leucine       (LEU)
A=12 //Alanine       (ALA)
W=13 //Tryptophan    (TRP)
P=14 //Proline       (PRO)
H=15 //Histidine     (HIS)
D=16 //Aspartic Acid (ASP)
R=17 //Arginine      (ARG)
I=18 //Isoleucine    (ILE)
G=19 //Glysine       (GLY)
3=20 //Water         (TIP or HOH)
*/
  FILE *scat_data=fopen(residDataFile, "r");
  if(scat_data==NULL)
    {
      fprintf(stderr, "ERROR: File %s not found!\n", residDataFile);
      exit(EXIT_FAILURE);
    }
  int i=0, j=0;
  int fscanfCheck=0;
  while(!feof(scat_data))
    {
      for(i=0; i<21; i++)
	{
	  fscanfCheck=fscanf(scat_data, "%lf\n", &F_CG[i][j]);
	  if(fscanfCheck==0)
	    {
	      fprintf(stderr, "ERROR: Data in %s file not read properly!\n", residDataFile);
	      exit(EXIT_FAILURE);
	    }
	}
      j++;
    }
  if(0)
    for(j=0; j<100; j++)
      fprintf(stdout, "%lf\t%lf\n", (double)j/1000.0, F_CG[20][j]);
  fclose(scat_data);
}
void readYangFF_v2(int num_data_points, double F_CG[21][num_data_points], char *residDataFile)
{
  //Whats new: Read extended form factors which goes until 1.8(2*pi/d)!
/*
S=0  //Serine        (SER)
F=1  //Phenylalanine (PHE)
T=2  //Threonine     (THR)
N=3  //Asparagine    (ASN)
K=4  //Lysine        (LYS)
E=5  //Glutamic Acid (GLU)
Y=6  //Tyrosine      (TYR)
V=7  //Valine        (VAL)
Q=8  //Glutamine     (GLN)
M=9  //Methionine    (MET)
C=10 //Cysteine      (CYS)
L=11 //Leucine       (LEU)
A=12 //Alanine       (ALA)
W=13 //Tryptophan    (TRP)
P=14 //Proline       (PRO)
H=15 //Histidine     (HIS)
D=16 //Aspartic Acid (ASP)
R=17 //Arginine      (ARG)
I=18 //Isoleucine    (ILE)
G=19 //Glysine       (GLY)
3=20 //Water         (TIP or HOH)
*/
  FILE *scat_data=fopen(residDataFile, "r");
  if(scat_data==NULL)
    {
      fprintf(stderr, "ERROR: File %s not found!\n", residDataFile);
      exit(EXIT_FAILURE);
    }
  int i=0, j=0;
  int fscanfCheck=0;

  for(i=0; i<21; i++)
    {
      for(j=0; j<num_data_points; j++)
	F_CG[i][j]=0.0;
    }
  j=0;

  while(!feof(scat_data))
    {
      for(i=0; i<21; i++)
	{
	  fscanfCheck=fscanf(scat_data, "%lf\n", &F_CG[i][j+9]);
	  F_CG[i][j+8]=F_CG[i][j+9];
	  F_CG[i][j+7]=F_CG[i][j+9];
	  F_CG[i][j+6]=F_CG[i][j+9];

	  if(fscanfCheck==0)
	    {
	      fprintf(stderr, "ERROR: Data in %s file not read properly!\n", residDataFile);
	      exit(EXIT_FAILURE);
	    }
	}
      j+=4;
    }

  for(i=0; i<21; i++)
    {
      F_CG[i][5]=F_CG[i][6];
      F_CG[i][4]=F_CG[i][6];
      F_CG[i][3]=F_CG[i][6];
      F_CG[i][2]=F_CG[i][6];
      F_CG[i][1]=F_CG[i][6];
      F_CG[i][0]=F_CG[i][6];
    }
 

  if(0)
    for(j=0; j<num_data_points; j++)
      {
	//fprintf(stdout, "%.4lf\t%.4lf\n", (double)j*2*M_PI/10000, F_CG[0][j]);
	fprintf(stdout, "%d\t%.4lf\n", (j+1), F_CG[0][j]);
      }
  fclose(scat_data);
}

double resid2YangCoefficient(char *resName, double F_CG[21][100], double q/*q value of form factor*/, double waterWeight)
{
  //Purpose: To obtain residue weights for SAXS specified in: 
  //         A rapid coarse residue-based computational method 
  //         for X-ray solution scattering characterization of 
  //         protein folds and multiple conformational states 
  //         of large protein complexes, Biophysical Journal 96,
  //         4449-4463 (2009). Sichun Yang, Sanghyun Park, 
  //         Lee Makowski, and Benoît Roux. 

  //Keep in mind that you cant go beyond q=0.628 in this table!!!!!
  //For now, I will just take q=0 value.
  //  printf ("%.3s\n", resName);
  //At first, determine the bin: 
  //Divide it by 2*PI
  
  double q_by_2PI=((q*1000)/(2*M_PI)); //1000 just comes to convert it properly to an array index;
  int qIndex=round(q_by_2PI);
  if(qIndex==0) qIndex=1; //Since in the file there is not q=0 value, I am interpolating q=0.00000 -> q=0.00628 

  qIndex=qIndex-1;

  double difference=((double)qIndex - q_by_2PI)/1000.0;
  if(qIndex!=0) assert(abs(difference)<0.0005); //---------Check if you are supposed to use half value??????????
  if(qIndex>99) 
    {
      qIndex=99; 
      fprintf(stderr, "ERROR: You don't have q values greater than 0.628!!!!!!!!!");
      exit(EXIT_FAILURE);
    }

  //For now, I will just take q=0 value.
  if(0) printf ("Residue Name:%.3s\tQ index= %d\tq_by_2PI=%lf\tq=%lf\n", resName, qIndex, q_by_2PI, q);
  if((strncmp(resName, "SER", 3)==0)) return (F_CG[0][qIndex]); // S=0  //Serine        (SER)
  if((strncmp(resName, "PHE", 3)==0)) return (F_CG[1][qIndex]); // F=1  //Phenylalanine (PHE)
  if((strncmp(resName, "THR", 3)==0)) return (F_CG[2][qIndex]); // T=2  //Threonine     (THR)
  if((strncmp(resName, "ASN", 3)==0)) return (F_CG[3][qIndex]); // N=3  //Asparagine    (ASN)
  if((strncmp(resName, "LYS", 3)==0)) return (F_CG[4][qIndex]); // K=4  //Lysine        (LYS)
  if((strncmp(resName, "GLU", 3)==0)) return (F_CG[5][qIndex]); // E=5  //Glutamic Acid (Glu)
  if((strncmp(resName, "TYR", 3)==0)) return (F_CG[6][qIndex]); // Y=6  //Tyrosine      (TYR)
  if((strncmp(resName, "VAL", 3)==0)) return (F_CG[7][qIndex]); // V=7  //Valine        (VAL)
  if((strncmp(resName, "GLN", 3)==0)) return (F_CG[8][qIndex]); // Q=8  //Glutamine     (GLN)
  if((strncmp(resName, "MET", 3)==0)) return (F_CG[9][qIndex]); // M=9  //Methionine    (MET)
  if((strncmp(resName, "CYS", 3)==0)) return (F_CG[10][qIndex]);// C=10 //Cysteine      (CYS)
  if((strncmp(resName, "LEU", 3)==0)) return (F_CG[11][qIndex]);// L=11 //Leucine       (LEU)
  if((strncmp(resName, "ALA", 3)==0)) return (F_CG[12][qIndex]);// A=12 //Alanine       (ALA)
  if((strncmp(resName, "TRP", 3)==0)) return (F_CG[13][qIndex]);// W=13 //Tryptophan    (TRP)
  if((strncmp(resName, "PRO", 3)==0)) return (F_CG[14][qIndex]);// P=14 //Proline       (PRO)

  // H=15 //Histidine(HIS) and all protonation states of it: HSD, HSE, HSP 
  if((strncmp(resName, "HIS", 3)==0) || (strncmp(resName, "HSD", 3)==0) \
  || (strncmp(resName, "HSE", 3)==0) || (strncmp(resName, "HSP", 3)==0)) return (F_CG[15][qIndex]); 

  if((strncmp(resName, "ASP", 3)==0)) return (F_CG[16][qIndex]);// D=16 //Aspartic Acid (ASP)
  if((strncmp(resName, "ARG", 3)==0)) return (F_CG[17][qIndex]);// R=17 //Arginine      (ARG)
  if((strncmp(resName, "ILE", 3)==0)) return (F_CG[18][qIndex]);// I=18 //Isoleucine    (ILE)
  if((strncmp(resName, "GLY", 3)==0)) return (F_CG[19][qIndex]);// G=19 //Glysine       (GLY)
  if((strncmp(resName, "TIP", 3)==0) || (strncmp(resName, "HOH", 3)==0)) return (waterWeight*F_CG[20][qIndex]);// 3=20 //Water         (TIP or HOH)
  else
    {
      fprintf(stderr, "Unknown residue name\n");
      exit(EXIT_FAILURE);
      return (0.0);
    }
}

double resid2YangFF_Long(int num_data_points, char *resName, double F_CG[21][num_data_points],\
			 double q/*q value of form factor*/, double waterWeight)
{
  //Purpose: To obtain residue weights for SAXS specified in: 
  //         A rapid coarse residue-based computational method 
  //         for X-ray solution scattering characterization of 
  //         protein folds and multiple conformational states 
  //         of large protein complexes, Biophysical Journal 96,
  //         4449-4463 (2009). Sichun Yang, Sanghyun Park, 
  //         Lee Makowski, and Benoît Roux. 

  //Keep in mind that you cant go beyond q=0.628 in this table!!!!!
  //For now, I will just take q=0 value.
  //  printf ("%.3s\n", resName);
  //At first, determine the bin: 
  //Divide it by 2*PI

  //Whats new: This version supersedes previous version since it reads amino acid form factors up to higher q(1.8-2pi/d)  

  double q_by_2PI=((q*10000)/(2*M_PI)); //1000 just comes to convert it properly to an array index;
  int qIndex=round(q_by_2PI);
  if(qIndex==0) qIndex=1;
  qIndex=qIndex-1;

  double difference=((double)qIndex - q_by_2PI);
  //  printf("Difference is %lf\n", difference);
  if(qIndex!=0) assert(abs(difference)<=1.0); 
  if(qIndex>2997) 
    {
      qIndex=2997; 
      fprintf(stderr, "ERROR: You don't have q values greater than 1.8\n!!!!!!!!!");
      exit(EXIT_FAILURE);
    }

  //For now, I will just take q=0 value.
  if(0) printf ("Residue Name:%.3s\tQ index= %d\tq_by_2PI=%lf\tq=%lf\n", resName, qIndex, q_by_2PI, q);
  if((strncmp(resName, "SER", 3)==0)) return (F_CG[0][qIndex]); // S=0  //Serine        (SER)
  if((strncmp(resName, "PHE", 3)==0)) return (F_CG[1][qIndex]); // F=1  //Phenylalanine (PHE)
  if((strncmp(resName, "THR", 3)==0)) return (F_CG[2][qIndex]); // T=2  //Threonine     (THR)
  if((strncmp(resName, "ASN", 3)==0)) return (F_CG[3][qIndex]); // N=3  //Asparagine    (ASN)
  if((strncmp(resName, "LYS", 3)==0)) return (F_CG[4][qIndex]); // K=4  //Lysine        (LYS)
  if((strncmp(resName, "GLU", 3)==0)) return (F_CG[5][qIndex]); // E=5  //Glutamic Acid (Glu)
  if((strncmp(resName, "TYR", 3)==0)) return (F_CG[6][qIndex]); // Y=6  //Tyrosine      (TYR)
  if((strncmp(resName, "VAL", 3)==0)) return (F_CG[7][qIndex]); // V=7  //Valine        (VAL)
  if((strncmp(resName, "GLN", 3)==0)) return (F_CG[8][qIndex]); // Q=8  //Glutamine     (GLN)
  if((strncmp(resName, "MET", 3)==0)) return (F_CG[9][qIndex]); // M=9  //Methionine    (MET)
  if((strncmp(resName, "CYS", 3)==0)) return (F_CG[10][qIndex]);// C=10 //Cysteine      (CYS)
  if((strncmp(resName, "LEU", 3)==0)) return (F_CG[11][qIndex]);// L=11 //Leucine       (LEU)
  if((strncmp(resName, "ALA", 3)==0)) return (F_CG[12][qIndex]);// A=12 //Alanine       (ALA)
  if((strncmp(resName, "TRP", 3)==0)) return (F_CG[13][qIndex]);// W=13 //Tryptophan    (TRP)
  if((strncmp(resName, "PRO", 3)==0)) return (F_CG[14][qIndex]);// P=14 //Proline       (PRO)

  // H=15 //Histidine(HIS) and all protonation states of it: HSD, HSE, HSP 
  if((strncmp(resName, "HIS", 3)==0) || (strncmp(resName, "HSD", 3)==0) \
  || (strncmp(resName, "HSE", 3)==0) || (strncmp(resName, "HSP", 3)==0)) return (F_CG[15][qIndex]); 

  if((strncmp(resName, "ASP", 3)==0)) return (F_CG[16][qIndex]);// D=16 //Aspartic Acid (ASP)
  if((strncmp(resName, "ARG", 3)==0)) return (F_CG[17][qIndex]);// R=17 //Arginine      (ARG)
  if((strncmp(resName, "ILE", 3)==0)) return (F_CG[18][qIndex]);// I=18 //Isoleucine    (ILE)
  if((strncmp(resName, "GLY", 3)==0)) return (F_CG[19][qIndex]);// G=19 //Glysine       (GLY)
  if((strncmp(resName, "TIP", 3)==0) || (strncmp(resName, "HOH", 3)==0)) return (waterWeight*F_CG[20][qIndex]);// 3=20 //Water         (TIP or HOH)
  else
    {
      fprintf(stderr, "Unknown residue name\n");
      exit(EXIT_FAILURE);
      return (0.0);
    }
}
void readStovgaardFF(double F_CG[21][100], char *residDataFile)
{
  // Purpose: Read data of residue forma factors from 
  //          "Calculation of Accurate Small Angle X-ray Scattering Curves from Coarse-Grained Protein Models. 
  //           BMC Bioinformatics. 2010; 11: 429. Stovgaard K, Andreetta C, Hamelryck T."

/*
S=0  //Serine        (SER)
F=1  //Phenylalanine (PHE)
T=2  //Threonine     (THR)
N=3  //Asparagine    (ASN)
K=4  //Lysine        (LYS)
E=5  //Glutamic Acid (GLU)
Y=6  //Tyrosine      (TYR)
V=7  //Valine        (VAL)
Q=8  //Glutamine     (GLN)
M=9  //Methionine    (MET)
C=10 //Cysteine      (CYS)
L=11 //Leucine       (LEU)
A=12 //Alanine       (ALA)
W=13 //Tryptophan    (TRP)
P=14 //Proline       (PRO)
H=15 //Histidine     (HIS)
D=16 //Aspartic Acid (ASP)
R=17 //Arginine      (ARG)
I=18 //Isoleucine    (ILE)
G=19 //Glysine       (GLY)
3=20 //Water         (TIP or HOH)*/
  FILE *ff_data=fopen(residDataFile, "r");
  if(ff_data==NULL)
    {
      fprintf(stderr, "ERROR: File %s not found!\n", residDataFile);
      exit(EXIT_FAILURE);
    }
  int j=0;
  int fscanfCheck=0;
  double q=0.000;
  while(!feof(ff_data))
    {
      fscanfCheck=fscanf(ff_data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
			 &q,&F_CG[12][j], &F_CG[17][j], &F_CG[3][j],  &F_CG[16][j], &F_CG[10][j], &F_CG[8][j], &F_CG[5][j], \
	    		   &F_CG[19][j], &F_CG[15][j], &F_CG[18][j], &F_CG[11][j], &F_CG[4][j],  &F_CG[9][j], &F_CG[1][j], \
			   &F_CG[14][j], &F_CG[0][j],  &F_CG[2][j],  &F_CG[13][j], &F_CG[6][j],  &F_CG[7][j]);
      if(fscanfCheck<=0)
	{
	  fprintf(stderr, "ERROR: Data in %s file not read properly!\n", residDataFile);
	  exit(EXIT_FAILURE);
	}
      j++;
    }
  fprintf(stdout, "Number of Data points in file %s is %d\n", residDataFile, j);
  if(0)
    for(j=0; j<51; j++)
      fprintf(stdout, "%d\t%lf\n", j, F_CG[19][j]);
  fclose(ff_data);
}

double resid2StovgaardCoefficient(char *resName, double F_CG[21][100], double q/*q value of form factor*/, double waterWeight)
{
  // Purpose: Read data of residue forma factors from 
  //          "Calculation of Accurate Small Angle X-ray Scattering Curves from Coarse-Grained Protein Models. 
  //           BMC Bioinformatics. 2010; 11: 429. Stovgaard K, Andreetta C, Hamelryck T."
  
  int qIndex=round(q/(0.015)); //Since intervals are 0.015, I divide it to that number and round it to obtain the corresponding integer.
  double difference=(((double)qIndex)*0.015 - q);
  assert(abs(difference)<0.015); //Check if your difference is smaller than 0.015. 
  
  //For now, I will just take q=0 value.
  if(0) printf ("Residue Name:%.3s\tQ index= %d\tq=%lf\n", resName, qIndex, q);
  if((strncmp(resName, "SER", 3)==0)) return (F_CG[0][qIndex]); // S=0  //Serine        (SER)
  if((strncmp(resName, "PHE", 3)==0)) return (F_CG[1][qIndex]); // F=1  //Phenylalanine (PHE)
  if((strncmp(resName, "THR", 3)==0)) return (F_CG[2][qIndex]); // T=2  //Threonine     (THR)
  if((strncmp(resName, "ASN", 3)==0)) return (F_CG[3][qIndex]); // N=3  //Asparagine    (ASN)
  if((strncmp(resName, "LYS", 3)==0)) return (F_CG[4][qIndex]); // K=4  //Lysine        (LYS)
  if((strncmp(resName, "GLU", 3)==0)) return (F_CG[5][qIndex]); // E=5  //Glutamic Acid (Glu)
  if((strncmp(resName, "TYR", 3)==0)) return (F_CG[6][qIndex]); // Y=6  //Tyrosine      (TYR)
  if((strncmp(resName, "VAL", 3)==0)) return (F_CG[7][qIndex]); // V=7  //Valine        (VAL)
  if((strncmp(resName, "GLN", 3)==0)) return (F_CG[8][qIndex]); // Q=8  //Glutamine     (GLN)
  if((strncmp(resName, "MET", 3)==0)) return (F_CG[9][qIndex]); // M=9  //Methionine    (MET)
  if((strncmp(resName, "CYS", 3)==0)) return (F_CG[10][qIndex]);// C=10 //Cysteine      (CYS)
  if((strncmp(resName, "LEU", 3)==0)) return (F_CG[11][qIndex]);// L=11 //Leucine       (LEU)
  if((strncmp(resName, "ALA", 3)==0)) return (F_CG[12][qIndex]);// A=12 //Alanine       (ALA)
  if((strncmp(resName, "TRP", 3)==0)) return (F_CG[13][qIndex]);// W=13 //Tryptophan    (TRP)
  if((strncmp(resName, "PRO", 3)==0)) return (F_CG[14][qIndex]);// P=14 //Proline       (PRO)

  // H=15 //Histidine(HIS) and all protonation states of it: HSD, HSE, HSP 
  if((strncmp(resName, "HIS", 3)==0) || (strncmp(resName, "HSD", 3)==0) \
  || (strncmp(resName, "HSE", 3)==0) || (strncmp(resName, "HSP", 3)==0)) return (F_CG[15][qIndex]); 

  if((strncmp(resName, "ASP", 3)==0)) return (F_CG[16][qIndex]);// D=16 //Aspartic Acid (ASP)
  if((strncmp(resName, "ARG", 3)==0)) return (F_CG[17][qIndex]);// R=17 //Arginine      (ARG)
  if((strncmp(resName, "ILE", 3)==0)) return (F_CG[18][qIndex]);// I=18 //Isoleucine    (ILE)
  if((strncmp(resName, "GLY", 3)==0)) return (F_CG[19][qIndex]);// G=19 //Glysine       (GLY)
  if((strncmp(resName, "TIP", 3)==0) || (strncmp(resName, "HOH", 3)==0)) return (waterWeight*F_CG[20][qIndex]);// 3=20 //Water         (TIP or HOH)
  else
    {
      fprintf(stderr, "Unknown residue name\n");
      exit(EXIT_FAILURE);
      return (0.0);
    }
}
void determineResLengths(int *atomCount, pdb_v23 *info, int *residLengths)
{
  //Purpose: Given an all atom structure, one sometimes needs to know number of atoms in each residue. 
  //         This function will return an integer array which contains that information! 
  //         If there are alternative locations for atoms in pdb, it may cause counting problems. 
  
  int i=0, j=0, k=0; //i is atom counter, j is internal counter for residue length, k is internal counter for residues.

  //Find the beginning and end of each residue and save them in an array that has elements as much as number of residues.
  for(i=0; i<(atomCount[0]-1); i++)
    {      //If their residue sequence number is same, all atoms are in the same residue.
      if( (info[i].resSeq==info[i+1].resSeq)&&(i!=(atomCount[0]-1))&&	\
	  ((info[i].altLoc=='A') || (info[i].altLoc==' ')) && ((info[i+1].altLoc=='A') || (info[i+1].altLoc==' ')))
	{
	  j++;
	  //	  printf("End of residue:%d at atom %d\n", info[i].resSeq, i);
	}
      else if((info[i].resSeq!=info[i+1].resSeq)&&(i!=(atomCount[0]-1))&& \
	      ((info[i].altLoc=='A') || (info[i].altLoc==' ')) && ((info[i+1].altLoc=='A') || (info[i+1].altLoc==' ')))
	{
	  residLengths[k]=(j+1); //Except the last residue, we have lengths of all residues.
	  k++;
	  //	  printf("End of residue:%d at atom %d. Residue length=%d\n", info[i].resSeq, i, (j+1));
	  j=0;
	}
    }

  //Lets calculate length of the last residue also:
  //Doing it in this way prevent wrong counts due to missing residues and alternate location indicators!!!
  j=0; //Zero your counter!
  for(i=(atomCount[0]-1); i>0; i--)
    {      
      if( (info[i].resSeq==info[i-1].resSeq) &&
	  ((info[i].altLoc=='A') || (info[i].altLoc==' ')) && ((info[i-1].altLoc=='A') || (info[i-1].altLoc==' ')))
	{
	  j++;
	  //	  printf("End of residue:%d at atom %d\n", info[i].resSeq, i);
	}
      else 
	{
	  j=j+1;
	  break;
	}
    }

  residLengths[atomCount[1]-1]=j;

  if(0)
    {
      for(k=0; k<atomCount[1]; k++)
      	printf("Residue %d Length %d\n", k , residLengths[k]);
    }
}

double calculateResidueFF_v3( pdb_v23 *info, int *atomCount, double s/*Scattering vector*/,\
			     double rho_s/*Solvent density=0.334(for water)*/, int isSolventOn/*0:NO - 1:YES*/, double *F_CG_per_q, int *badRes)
{
  //What's new: In CHARMM pdb, there is not 'element' data section. Instead, there is segment name. Therefore, I modified where I get
  //element information of a residue in this v2.
  //Keep in mind that s=sin(theta)/lambda
  //What's new:05/17/2011
  //This time the purpose is to calculate form factors for one pdb file including end of chain residues. 

  int residLengths[atomCount[1]];
  int i=0, j=0, k=0; //i is atom counter, j is internal counter for residue length, k is internal counter for residues.

  //=====================================================================================================================
  determineResLengths(atomCount, info, residLengths);
  //=====================================================================================================================

  //Count through all residues.
  int res_end_atom  =(-1);
  int res_begin_atom=0;
  int howManyRes=0;
  double F_CG_temp=0.0;


  for(k=0; k<atomCount[1]; k++)
    {
      res_end_atom+=(residLengths[k]);
      res_begin_atom=(res_end_atom - residLengths[k] + 1);

      double F_CG_res=0.0;
      //      printf("Residue %d begining=%d and ending=%d\n", k, res_begin_atom, res_end_atom);

      //      if(residLengths[k]==(aa_atomNumbers(residueName)-3) )//Check if all atoms of residue exists.
      if(badRes[k]) //No need to check integrity of residue here. badRes tagging will be used for control purposes only. 
	  {           
	    for(i=res_begin_atom; i<=res_end_atom; i++)
	      {
		for(j=res_begin_atom; j<=i; j++)
		  {
		    char element_i=' ', element_j=' ';
		    if(info[i].name[0]==' ') element_i=info[i].name[1];
		    else element_i=info[i].name[0];
		    
		    //		    element_i=info[i].name[0];
		    //		    element_j=info[j].name[0];
		    
		    if(info[j].name[0]==' ') element_j=info[j].name[1];
		    else element_j=info[j].name[0];
		    
		    if(0)
		      printf("Name 1:%s\tElement 1:%c\tLength 1:%d\tName 2:%s\tElement 2:%c\tLength 2:%d\n", \
			     info[i].name, element_i, (int)strlen(info[i].name), info[j].name, element_j, (int)strlen(info[j].name));
		    
		    if((info[i].altLoc=='A') || (info[i].altLoc==' '))   //Select one of alternative locations i'th atom.
		      if((info[j].altLoc=='A') || (info[j].altLoc==' ')) //Select one of alternative locations i'th atom.
			{
			  //			  double f_prime_i=atom_form_factors(element_i, s);
			  //			  double f_prime_j=atom_form_factors(element_j, s);
			  double f_prime_i=calculateAtomicFF_v2(element_i, s, rho_s, isSolventOn);
			  double f_prime_j=calculateAtomicFF_v2(element_j, s, rho_s, isSolventOn);
			  double f_iXf_j=f_prime_i*f_prime_j;
			  if((s==0.000)&&(i==j))
			    F_CG_res+=(f_iXf_j);
			  else if((s==0.000)&&(i!=j))
			    F_CG_res+=(2*f_iXf_j);
			  else if((s!=0.000) && (i==j))
			    F_CG_res+=(f_iXf_j);
			  else if((s!=0.000) && (i!=j))
			    {
			      double r_ij=0.0;
			      r_ij=disAll_pdb_v23(info, i, j);
			      double q=s;
			      double qXr_ij=q*r_ij;
			      F_CG_res+=(2*f_iXf_j*sin(qXr_ij)/qXr_ij);
			      //			      printf("(%d, %d)\n", i, j);
			    }
			}
		  }
		if (0) printf("Atom name: %s\n", info[i].name);
	      }
	    F_CG_per_q[k]=sqrt(F_CG_res);
	    //	    F_CG_temp+=F_CG_res;
	    //	    printf("Name=%s\tResidue No=%d F_CG=%lf\n", info[res_begin_atom].resName, k, F_CG_per_q[k]);
	    //	    printf("%s(%d) %lf =%lf\n", info[res_begin_atom].resName,  k, s, F_CG_res);
	    howManyRes++;
	  }
    }
  //  printf("Number of residues in file=%d \n", howManyRes);
  //  printf("F_CG_temp=%lf\n", F_CG_temp);
      //  resCount[0]=howManyRes;
  //Do the averaging all over pdb set, not just one residue set.
  return (F_CG_temp);
}
