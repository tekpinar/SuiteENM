//Purpose: enm_functions.c contains core functions of my PhD research. 
//         All functions related to elastic network model are in this file. 

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
#include <assert.h>
#include <math.h>
//Personal header files
//#include <defines.h>
#include <structures.h>
#include <distance_functions.h>
#define ENM_DEBUG 0
#define CL_CON 10.0                        /* If close contact(neighbor residues) force constant is 10.                */
#define FR_CON 1.0                         /* If far contact(non-neighbor residues) force constant is 1                */
#define C_COL 1000.0                         /* If far contact(non-neighbor residues) approach each other too much!!!!   */
#define SS_FORCE_CONSTANT       1.0        //SS: (S)econdary (S)tructure-This is used in SAXS flexible fitting!!

#define INVALID_CRD 999.999
#define CA_DIST 3.9


double R_CUTOFF=10.0;                      /* Don't forget to change R_CUTOFF_SQUARED if you change this value!!       */
double R_CUTOFF_SQUARED=100.0;             /* Don't forget to change R_CUTOFF if you change this value!!               */
#define VAN_DER_WAALS 0.0                  /* Force constant for van der Walls interactions                            */
int isCoordValid(CAcoord *atom, int i)
{
  if((atom[i].x==INVALID_CRD) || (atom[i].y==INVALID_CRD) || (atom[i].z==INVALID_CRD))
    return 0;
  else 
    return 1;
}

int isContact(CAcoord *atom, int i, int j, double *dstnc_ij)
{
  /*Purpose: To check if two CA atoms are in contact distance which is defined by R_CUTOFF*/
  //  printf("R_CUTOFF=%lf\tR_CUTOFF_SQUARED=%lf\n", R_CUTOFF, R_CUTOFF_SQUARED);
  double diff_x=fabs(atom[i].x-atom[j].x);
  if(diff_x>=R_CUTOFF)
    return 0;
  
      double diff_y=fabs(atom[i].y-atom[j].y);
      if(diff_y>=R_CUTOFF)
	return 0;
      
      double diff_z=fabs(atom[i].z-atom[j].z);
      if(diff_z>=R_CUTOFF)
	return 0;
      
      double rSquared=((diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z) );
      
      if(rSquared >=  R_CUTOFF_SQUARED)
	return 0;
      
      else
	{
	  *dstnc_ij=sqrt(rSquared);
	  //      printf("Distance[%d][%d]=%lf\n",i, j,  *dstnc_ij);
	  return 1;
	}

}

int isContact_v2(CAcoord *atom, int i, int j, double *dstnc_ij, double R_cutoff)
{

  /*Purpose: To check if two CA atoms are in contact distance which is defined by R_CUTOFF*/
  //  printf("R_CUTOFF=%lf\tR_CUTOFF_SQUARED=%lf\n", R_CUTOFF, R_CUTOFF_SQUARED);
  //In this version, R_CUTOFF is not a 'defined' value. Instead, it is an argument of the function. 
  double R_cutoff_sqrd=(R_cutoff*R_cutoff);
  double diff_x=fabs(atom[i].x-atom[j].x);
  if(diff_x>=R_cutoff)
    return 0;
  
  double diff_y=fabs(atom[i].y-atom[j].y);
  if(diff_y>=R_cutoff)
    return 0;
  
  double diff_z=fabs(atom[i].z-atom[j].z);
  if(diff_z>=R_cutoff)
    return 0;
  
  double rSquared=((diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z) );
  
  if(rSquared >=  R_cutoff_sqrd)
    return 0;
      
  else
    {
      *dstnc_ij=sqrt(rSquared);
      //      printf("Distance[%d][%d]=%lf\n",i, j,  *dstnc_ij);
      return 1;
    }
  
}
void getPairInfo_new(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, int *map1)
{
  /*Purpose: To get the information of pairs in contact and put them into 'pairInfoNew' structure*/
  int i=0, j=0, l=0; /*Pair counter*/
  
  /*   for(k=0; k< N; k++) */
  /*     fprintf(stdout, "Atom1: %dCorresponding Atom2: %d\n", k, map1[k]); */
  
  double dstnc_ij[1]; 
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)         
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 1) )
	      {
		if((abs(i-j))==1) pair[l].k=CL_CON;//Since j<i, no need to put abs here!!
		
		else pair[l].k=FR_CON;
		
		pair[l].dis= dstnc_ij[0];
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)	  
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }
	    else
	      {
		if((isContact(atomFixed, i, j, dstnc_ij)) == 1)
		  if(abs(i-j)==1)
		    {
		      pair[l].k=CL_CON;//Since j<i, no need to put abs here!!
		      pair[l].dis= CA_DIST;;
		      pair[l].posI=i;
		      pair[l].posJ=j;
		      FILE* f=fopen("UnmatchedAtomPairs.txt", "a");
		      fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		      fclose(f);
		      l++;
		    }
	      }
	  }
      }
  
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}

//========
void getPairInfo_v2(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, int *map1)
{
  /*Purpose: To get the information of pairs in contact and put them into 'pairInfoNew' structure*/
  int i=0, j=0, l=0; /*Pair counter*/
  
  /*   for(k=0; k< N; k++) */
  /*     fprintf(stdout, "Atom1: %dCorresponding Atom2: %d\n", k, map1[k]); */
  
  double dstnc_ij[1]; 
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)         
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 1) )
	      {
		if((abs(i-j))==1) pair[l].k=CL_CON;//Since j<i, no need to put abs here!!
		
		else pair[l].k=FR_CON;
		
		pair[l].dis= dstnc_ij[0];
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)	  
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }
	  }
	else
	  {
	    //Sometimes, after Dali alignment there are not some atoms corresponding to each other in two
	    //conformations. 
	    //This part of code ensures that even if some coordinates are invalid after alignment using Dali
	    //we consider the bonding energies for those missing atoms!!!
	    if((isContact(atomFixed, i, j, dstnc_ij)) == 1)
	      if(abs(i-j)==1)
		{
		  pair[l].k=CL_CON;//Since j<i, no need to put abs here!!
		  pair[l].dis= CA_DIST;;
		  pair[l].posI=i;
		  pair[l].posJ=j;
		  FILE* f=fopen("UnmatchedAtomPairs.txt", "a");
		  fprintf(f, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  fclose(f);
		  l++;
		}
	  }
      }

  
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}
//========
void getPairInfo_back(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, int *map1)
{
  /*Purpose: To get the information of pairs in contact and put them into 'pairInfoNew' structure*/
  int i=0, j=0, l=0; /*Pair counter*/
  double dstnc_ij[1]; 
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)         
    for(j=0; j<i; ++j)
      {    
	{
	  if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 1) )
	    {
	      if((abs(i-j))==1) pair[l].k=CL_CON;//Since j<i, no need to put abs here!!
	      
	      else pair[l].k=FR_CON;
	      
	      pair[l].dis= dstnc_ij[0];
	      
	      pair[l].posI=i;
	      pair[l].posJ=j;
	      
	      if(ENM_DEBUG)	  
		{
		  fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		  fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		  fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		}
	      l++;
	    }
	}
	
      }
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}



void getAllPairInfo(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo)
{
  /*Purpose: To get the information of pairs so that one can construct pair distribution function and 
    put that information into 'pairInfoNew' structure*/
  //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations.
  //In order just to get rid of ENM hessian calculation for pairs not in contact, I will just make their k==0
  int i=0, j=0, l=0; /*Pair counter*/
    double dstnc_ij[1]; 
  printf("Amino acid number=%d\n", N);

  for(i=0; i<N; ++i)         
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 1) )
	      {
		if((abs(i-j))==1) pair[l].k=CL_CON;//Since j<i, no need to put abs here!!
		
		else pair[l].k=FR_CON;
		
		pair[l].dis= dstnc_ij[0];
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)	  
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }

	    else if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 0) )
	      {
		pair[l].k=0.0;
		pair[l].dis= dis(atomFixed, i, j);
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)	  
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }

	  }
      }
  fprintf(stdout, "Sen almak istiyor duj, sen odeyecek 50 dolar daha\n");  
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}

void getAllPairInfo_v2(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo)
{
  /*Purpose: To get the information of pairs so that one can construct pair distribution function and
    put that information into 'pairInfoNew' structure*/
  //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations.
  //In order just to get rid of ENM hessian calculation for pairs not in contact, I will just make their k==0
  //Check residue numbers instead of your internal numbering!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int i=0, j=0, l=0; /*Pair counter*/
  double dstnc_ij[1];
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 1) )
	      {
		
		//If their residue numbers are consequtive and they are in the same chain, they are covalently bonded!!
		if(( abs(atomFixed[i].resNo - atomFixed[j].resNo)==1) && (atomFixed[i].chain==atomFixed[j].chain ))
		  {
		    if(    (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		      )
		    pair[l].k=VAN_DER_WAALS;
		    else
		    pair[l].k=CL_CON;
		  }
		else
		  { 
		    if(    (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		      )
		      pair[l].k=VAN_DER_WAALS;
		    else
		      pair[l].k=FR_CON;
		  }
		pair[l].dis= dstnc_ij[0];
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(0)
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i, j, l, pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }
	    //You dont need to check it again!!!!!!!!!!!!!!!!!!!!!!!
	    else if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 0) )
	      {
		pair[l].k=0.0;
		pair[l].dis= dis(atomFixed, i, j);
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }

	  }
      }
  
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}

void getAllPairInfo_v3(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, double R_cutoff)
{
  /*Purpose: To get the information of all pairs so that one can construct pair distribution function and
    put that information into 'pairInfoNew' structure*/
  //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations.
  //In order just to get rid of ENM hessian calculation for pairs not in contact, I will just make their k==0
  //Check residue numbers instead of your internal numbering!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Whats new: In version 2, I was using 'isContact'. In this version isContact_v2 is used. I added a new argument called R_cutoff, 
  //and it makes the program more flexible!
  int i=0, j=0, l=0; /*Pair counter*/
  double dstnc_ij[1];
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    if( ( (isContact_v2(atomFixed, i, j, dstnc_ij, R_cutoff)) == 1) )
	      {
		
		//If their residue numbers are consequtive and they are in the same chain, they are covalently bonded!!
		if(( abs(atomFixed[i].resNo - atomFixed[j].resNo)==1) && (atomFixed[i].chain==atomFixed[j].chain ))
		  {
		    if(    (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		      )
		    pair[l].k=VAN_DER_WAALS;
		    else
		    pair[l].k=CL_CON;
		  }
		else
		  { 
		    if(    (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		      )
		      pair[l].k=VAN_DER_WAALS;
		    else
		      pair[l].k=FR_CON;
		  }
		pair[l].dis= dstnc_ij[0];
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(0)
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i, j, l, pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }
	    else
	      {
		pair[l].k=0.0;
		pair[l].dis= dis(atomFixed, i, j);
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;

/* 		//Just to see CA water distances.  */
/* 		if((strncmp(atomFixed[i].residname, "TIP", 3)==0) &&(strncmp(atomFixed[j].residname, "TIP", 3)!=0)) */
/* 		  fprintf(stdout, "Water - CA pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
/* 		else if(((strncmp(atomFixed[i].residname, "TIP", 3)!=0) &&(strncmp(atomFixed[j].residname, "TIP", 3)==0))) */
/* 		  fprintf(stdout, "CA-water pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
		
	      }

	  }
      }
  
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}
void getAllNonBondedPairsInfo(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, double R_cutoff)
{
  /*Purpose: To get the information of all nonbonded pairs to use in modiefies ENM for flexible saxs fitting algorithm.*/
  //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations.

  int i=0, j=0, l=0; /*Pair counter*/
  double dstnc_ij[1];
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    if( ( (isContact_v2(atomFixed, i, j, dstnc_ij, R_cutoff)) == 1) )
	      {
		//If their residue numbers are consequtive and they are in the same chain, they are covalently bonded!!
		if(( abs(atomFixed[i].resNo - atomFixed[j].resNo)==1) && (atomFixed[i].chain==atomFixed[j].chain ))
		  {
		    if(    (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		      )
		      pair[l].k=VAN_DER_WAALS;
		    else
		      pair[l].k=0.0; //If they are waters, they are loosely connected. Else, they are not a part of nonbonded residue set!
 		  }
		else
		  { 
		    if(    (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		       ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		      )
		      pair[l].k=VAN_DER_WAALS;
		    else
		      {
			//Following our paper: Zheng, Tekpinar, Accurate Flexible Fitting ..., 2011, Biophysical Journal. 
			pair[l].k=(16.0/(dstnc_ij[0]*dstnc_ij[0]));
		      }
		  }
		pair[l].dis= dstnc_ij[0];
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(0)
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i, j, l, pair[l].k, l, pair[l].dis);
		  }
		l++;
	      }
	    else
	      {
		pair[l].k=0.0;
		pair[l].dis= dis(atomFixed, i, j);
		pair[l].posI=i;
		pair[l].posJ=j;
		
		if(ENM_DEBUG)
		  {
		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j));
		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
		  }
		l++;

/* 		//Just to see CA water distances.  */
/* 		if((strncmp(atomFixed[i].residname, "TIP", 3)==0) &&(strncmp(atomFixed[j].residname, "TIP", 3)!=0)) */
/* 		  fprintf(stdout, "Water - CA pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
/* 		else if(((strncmp(atomFixed[i].residname, "TIP", 3)!=0) &&(strncmp(atomFixed[j].residname, "TIP", 3)==0))) */
/* 		  fprintf(stdout, "CA-water pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
		
	      }

	  }
      }
  
  pairNo[0]=l;
  printf("Total number of residue pairs are %d\n", l);
}
void getSSandBondedPairsInfo(char *ssArray, pairInfoNew *pairSS, int N, CAcoord *atomFixed, int *pairSSNo, double *weight, int is_ss_on)
{
  /*Purpose: To get the information of pairs that are bonded or in a secondary structure!*/
  //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations.
  //In order just to get rid of ENM hessian calculation for pairs not in contact, I will just make their k==0
  int i=0, j=0, l=0; /*Pair counter*/

  printf("Amino acid number=%d\n", N);
  //  for(i=0; i<N; ++i)         
  //    fprintf(stdout, "SSarray[%d]=%c\n", i, ssArray[i]);
  //Count bonded pairs at first
  for(i=0; i<N; ++i)         
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1))
	  {
	    //Count just bonded pairs
	    //==========================
	    int posDiff=0;
	    // Bonded:      posDiff=1
	    // Strand:      posDiff=2
	    // 3_10 helix:  posDiff=3
	    // Alpha helix: posDiff=4
	    // Pi heliz:    posDiff=5
	    
	    posDiff= abs(atomFixed[i].resNo - atomFixed[j].resNo);
	    int isInSameChain=0;
	    if(atomFixed[i].chain==atomFixed[j].chain ) isInSameChain=1;
	    else                                        isInSameChain=0;
	    
	    if((posDiff==1) && (isInSameChain))
	      {
		if(   (strncmp(atomFixed[i].residname, "TIP", 3)==0) \
		  ||  (strncmp(atomFixed[i].residname, "HOH", 3)==0) \
		  ||  (strncmp(atomFixed[j].residname, "TIP", 3)==0) \
		  ||  (strncmp(atomFixed[j].residname, "HOH", 3)==0) \
		  )
		  pairSS[l].k=0.0;
		//pairSS[l].k=VAN_DER_WAALS;
		else
		  pairSS[l].k=CL_CON;

		pairSS[l].dis= dis(atomFixed, i, j);
		pairSS[l].posI=i;
		pairSS[l].posJ=j;
		pairSS[l].Pweight=weight[i]*weight[j];

		if(0)	  
		  {
		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\tpair[%d].Pweight=%lf\n", \
			    i, j, l, pairSS[l].k, l, pairSS[l].dis, l, pairSS[l].Pweight);
		  }
	      }

	    else 
	      {
		pairSS[l].k=0.0;
		pairSS[l].dis= dis(atomFixed, i, j);
		pairSS[l].posI=i;
		pairSS[l].posJ=j;
		pairSS[l].Pweight=weight[i]*weight[j];
		
		if((posDiff>=2) && (posDiff<=5) && (isInSameChain ))
		  {
		    int t=0;
		    int isContinuous=0; //Check if there is any missing residue.
		    //Since j is always smaller that i, this will work at j_min=0 and j_max=i-1.
		    for(t=j; t<(j+posDiff); t++)
		      {
			if(abs(atomFixed[t+1].resNo - atomFixed[t].resNo)==1) isContinuous=1;
			else
			  {  
			    isContinuous=0;
			    break;
			  }
		      }
		    
		    int isHelixORStrand =0;
		    //Dogrusu bu strand olayi cok sacma geliyor!
		    for(t=j; t<=(j+posDiff); t++)
		      {
			if((posDiff==2)&&(ssArray[t]=='E')) isHelixORStrand=1;
			else if((posDiff==3)&&(ssArray[t]=='G')) isHelixORStrand=1;
			else if((posDiff==4)&&(ssArray[t]=='H')) isHelixORStrand=1;
			else if((posDiff==5)&&(ssArray[t]=='I')) isHelixORStrand=1;
			
			else
			  {  
			    isHelixORStrand=0;
			    break;
			  }
		      }
		    if((isHelixORStrand==1)&&(isContinuous))
		      {
			if(is_ss_on)
			  pairSS[l].k=SS_FORCE_CONSTANT; //To keep secondary structure intact, I need a number bigger than 1.0.
			if(0)
			  fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\tpair[%d].Pweight=%lf\n", \
				  i, j, l, pairSS[l].k, l, pairSS[l].dis, l, pairSS[l].Pweight);
		      }
		  }
	      }
	  }
	l++;
      }
  //  fprintf(stdout, "Number of bonded pairs is %d\n", l);

  pairSSNo[0]=l;
  printf("Total SS number of residue pairs are %d\n", l);
}

void getAllPairInfo_v4(pairInfoNew *pair, int N/*Number of CA*/, CAcoord *protCA, int *pairNo, double exponent, double coeff, int printDetails/*1:yes, 0:No*/)
{
  /*Purpose: To get the information of all pairs so that one can construct pair distribution function and
    put that information into 'pairInfoNew' structure*/
  //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations.
  //In order just to get rid of ENM hessian calculation for pairs not in contact, I will just make their k==0

  //Whats new: In version 2, I was using 'isContact'. In this version isContact_v2 is used. I added a new argument called R_cutoff, 
  //and it makes the program more flexible!

  //Whats new: In version 4,  I am trying to a distance based force constant system. (19 December 2011)

  int i=0, j=0, l=0; /*Pair counter*/
  double dstnc_ij[1];
  printf("Amino acid number=%d\n", N);
  for(i=0; i<N; ++i)
    for(j=0; j<i; ++j)
      {
	if((isCoordValid(protCA, i)==1) && (isCoordValid(protCA, j)==1) )
	  {

	    pair[l].dis= dis(protCA, i, j);
	    //If residues are valid and they are not waters, I will include them.
	    if(    (strncmp(protCA[i].residname, "TIP", 3)==0)		\
		   ||  (strncmp(protCA[i].residname, "HOH", 3)==0)	\
		   ||  (strncmp(protCA[j].residname, "TIP", 3)==0)	\
		   ||  (strncmp(protCA[j].residname, "HOH", 3)==0)	)
	      pair[l].k=0.0; //Zero weight for waters
	    else
	      pair[l].k=(coeff/(pow ( pair[l].dis,  exponent ))); 
	    
	    pair[l].posI=i;
	    pair[l].posJ=j;
	    
	    if(printDetails)
	      {
		fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(protCA, i, j));
		fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i, j, l, pair[l].k, l, pair[l].dis);
	      }
	    l++;
	  }
	else
	  {
	    pair[l].k=0.0;
	    pair[l].dis= 3.9;  //Correct the case where one of the coordinates is invalid!
	    //pair[l].dis= dis(protCA, i, j);
	    pair[l].posI=i;
	    pair[l].posJ=j;
	    
	    if(printDetails)
	      {
		fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(protCA, i, j));
		fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis);
		fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis);
	      }
	    //	    l++;
	    
	    /* 		//Just to see CA water distances.  */
	    /* 		if((strncmp(protCA[i].residname, "TIP", 3)==0) &&(strncmp(protCA[j].residname, "TIP", 3)!=0)) */
	    /* 		  fprintf(stdout, "Water - CA pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
	    /* 		else if(((strncmp(protCA[i].residname, "TIP", 3)!=0) &&(strncmp(protCA[j].residname, "TIP", 3)==0))) */
	    /* 		  fprintf(stdout, "CA-water pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
	    
	  }
 
      }
  
  pairNo[0]=l;
  if(pairNo[0]!=(N*(N-1)/2))
    {
      fprintf(stderr, "ERROR: There are invalid coordinates!\n");
      exit(EXIT_FAILURE);
    }
  printf("Total number of residue pairs are %d\n", l);
}


double  **constructHessian_new(int conformationNo, int N/*Number of residues*/, double **hessian, CAcoord *atomMoving, int *pairNo, pairInfoNew *pair)
{
  //Purpose: A more efficient way of constructing ENM hessian using pair information
  //Next   : Use a sparse matrix format instead of usual array
  int i=0, j=0, l=0; /*Counters for atom pairs                */
  int a=0, b=0;      /*Counters for atom cartesian coordinates*/

  for (i=0; i<3*N; i++) /*To zero input hessian so that previous calculation have nothing to do w/ this one*/
    for (j=0; j<3*N; j++)
      hessian[i][j]=0.0;
  
  for(l=0; l<pairNo[0]; ++l)         
    {
      double disMobile_ij=0.0;
      disMobile_ij= dis(atomMoving, pair[l].posI, pair[l].posJ);
      
      double x_i[3]={999.0, 999.0, 999.0};
      double x_j[3]={999.0, 999.0, 999.0};
      
      x_i[0]=atomMoving[pair[l].posI].x;
      x_i[1]=atomMoving[pair[l].posI].y;
      x_i[2]=atomMoving[pair[l].posI].z;
      
      x_j[0]=atomMoving[pair[l].posJ].x;
      x_j[1]=atomMoving[pair[l].posJ].y;
      x_j[2]=atomMoving[pair[l].posJ].z;

      double disSqMob_ij=0.0;
      disSqMob_ij=(disMobile_ij*disMobile_ij);
      
      double coeff1=0.0, coeff2=0.0;
      coeff1= pair[l].dis /disMobile_ij;
      coeff2= (pair[l].k)*coeff1;
      
      for(a=0; a<3; ++a)
	for(b=0; b<3; ++b) 
	  {
	    if(a!=b)
	      {
	       	double offdiagonal = -( coeff2 )   *(   ( x_i[a] - x_j[a] ) * ( x_i[b] - x_j[b])/   disSqMob_ij    ); 
		//double offdiagonal = -pair[l].k*   (    ( coeff1)    *(   ( x_i[a] - x_j[a] ) * ( x_i[b] - x_j[b])/   (  disSqMob_ij   )  )); 
		//double offdiagonal = -pair[l].k*(coeff1 * ( x_i[a] - x_j[a] ) * ( x_i[b] - x_j[b] ) / disSqMob_ij)  ; 
		//	  printf("offdiag: %lf\n", offdiagonal);    
		
		hessian[3*pair[l].posI+a][3*pair[l].posJ+b]=offdiagonal;
		hessian[3*pair[l].posJ+b][3*pair[l].posI+a]=offdiagonal;
		//	  printf("pair[%d].posI: %d\n", l, pair[l].posI);
		//	  printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		//	  printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), hessian[3*pair[l].posI+a][3*pair[l].posJ+b]);
	      }
	    else if (a==b)
	      {
		double ondiagonal = pair[l].k*(        (  coeff1   )          *(1- (          ( x_i[a] - x_j[a] ) * ( x_i[b] - x_j[b])/ disSqMob_ij    )   )-1);
		//	 printf("ondiag: %lf\n", ondiagonal);    
		hessian[3*pair[l].posI+a][3*pair[l].posJ+b]=ondiagonal;
		hessian[3*pair[l].posJ+b][3*pair[l].posI+a]=ondiagonal;
		//	  printf("pair[%d].posI: %d\n", l, pair[l].posI);
		//	  printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		//	  printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), hessian[3*pair[l].posI+a][3*pair[l].posJ+b]);
		
	      }
	  }	
    }
  
  
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      hessian[3*i+a][3*i+b] -= hessian[3*i+a][3*j+b]; 
 	      //printf("Hessian[%d][%d]=%lf\n",  (3*i+a), (3*i+b), Hessian[3*i+a][3*i+b]); 
 	    } 
  // printf("\nCalculated Hessian Matrix\n");

  if(ENM_DEBUG)
    {
      char temp0[15];
      char outFileName[14]="mtrlConf";
      
      sprintf(temp0, "%d", conformationNo);
      strncat(outFileName, temp0, 1);
      strncat(outFileName, ".txt",4);
      // printf("Output:%s\n", outFileName);

      FILE *matrix=fopen(outFileName, "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "%s file could not be produced\n", outFileName);
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
 	{ 
 	  for(j=0; j<i; j++)
	    for(a=0; a<3; a++) 
	      for(b=0; b<3; b++) 
		if(hessian[3*i+a][3*j+b]!=0.0)
		  fprintf(matrix, "%d\t%d\t%lf\n", 3*i+a, 3*j+b, hessian[3*i+a][3*j+b]); 
	  
	  //fprintf(matrix, "\n"); 
 	} 

      for(i=0; i<N; i++) 
 	{ 
 	  for(j=0; j<N; j++)
	    if(i==j)
	      {
		for(a=0; a<3; a++) 
		  for(b=0; b<3; b++) 
		    if(a<=b)
		      if(hessian[3*i+a][3*j+b]!=0.0)
			fprintf(matrix, "%d\t%d\t%lf\n",  3*i+a, 3*j+b,  hessian[3*i+a][3*j+b]); 
		
		//fprintf(matrix, "\n"); 
	      }
 	} 
      
      fclose(matrix);
    }
  
  return hessian;
}

//======Gradient and hessian ENM together
double gradientANDhessian_ENM(int conformationNo, int N/*Number of residues*/, 
			      double *grad_ENM, double **H_ENM, 
			      CAcoord *atomMoving, 
			      int *pairNo, pairInfoNew *pair,
			      int *pairInitialNo, pairInfoNew *pairInitial)
{
  int i=0, j=0, l=0;        /*Counters for atom pairs                */
  int a=0, b=0;             /*Counters for atom cartesian coordinates*/
  
  double fx=0.0, fy=0.0, fz=0.0;
  double coeff1=0.0, coeff2=0.0, coeff3=0.0;
  double disMobile_ij=0.0, disSqMob_ij=0.0;
  double coordDif[3]={0.0, 0.0, 0.0};

  for (i=0; i< (3*N); i++) /*To zero input hessian and gradient so that previous calculation have nothing to do w/ this one*/
    {
      for (j=0; j<=i; j++)
	{
	  H_ENM[i][j]=0.0;
	  H_ENM[j][i]=0.0;
	}
      grad_ENM[i]=0.0;
    }
  double E_ENM=0.0;
  double dis_diff=0.0;
  for (l=0;l<pairNo[0]; ++l) 
    {
      if(pair[l].k!=0.0)
	{      
	  dis_diff=(pair[l].dis-pairInitial[l].dis);
	  E_ENM+=pair[l].k*dis_diff*dis_diff;
	  //disMobile_ij=dis(atomMoving, (pair[l].posI), (pair[l].posJ));
	  disMobile_ij=pairInitial[l].dis;
	  disSqMob_ij=(disMobile_ij*disMobile_ij);
	  
	  coeff1= pair[l].dis /disMobile_ij;
	  coeff2= (pair[l].k)*coeff1;
	  coeff3= (pair[l].k)*(1- coeff1 );
	  coordDif[0]=(atomMoving[pair[l].posI].x-atomMoving[pair[l].posJ].x);
  	  coordDif[1]=(atomMoving[pair[l].posI].y-atomMoving[pair[l].posJ].y);
	  coordDif[2]=(atomMoving[pair[l].posI].z-atomMoving[pair[l].posJ].z);

	  fx=  coeff3*(coordDif[0]);
	  fy=  coeff3*(coordDif[1]);
	  fz=  coeff3*(coordDif[2]);
	  
	  int Ind_I=3*(pair[l].posI);
	  int Ind_J=3*(pair[l].posJ);

	  //	  ptr=3*(pair[l].posI);
	  grad_ENM[Ind_I]+=fx;
	  grad_ENM[Ind_I+1]+=fy;
	  grad_ENM[Ind_I+2]+=fz;
	  
	  //	  ptr=3*(pair[l].posJ);
	  grad_ENM[Ind_J]-=fx;
	  grad_ENM[Ind_J+1]-=fy;
	  grad_ENM[Ind_J+2]-=fz;
	  
	    //========Copy hessian here
	    //fprintf(stdout, "Pair Initial Distance=%lf\tDistance calculated=%lf\n", pairInitial[l].dis, disMobile_ij);
	  for(a=0; a<3; ++a)
	    for(b=0; b<3; ++b) 
	      {
		if(a!=b)
		  {
		    double offdiagonal = -( coeff2 )   *(   ( coordDif[a] ) * ( coordDif[b])/   disSqMob_ij    ); 
		    //double offdiagonal = -pair[l].k*   (    ( coeff1)    *(   ( x_i[a] - x_j[a] ) * ( x_i[b] - x_j[b])/   (  disSqMob_ij   )  )); 
		    //double offdiagonal = -pair[l].k*(coeff1 * ( x_i[a] - x_j[a] ) * ( x_i[b] - x_j[b] ) / disSqMob_ij)  ; 
		    //printf("offdiag: %lf\n", offdiagonal);    
		    
		    H_ENM[Ind_I+a][Ind_J+b]=offdiagonal;
		    H_ENM[Ind_J+b][Ind_I+a]=offdiagonal;
		    //printf("pair[%d].posI: %d\n", l, pair[l].posI);
		    //printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		    //printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), H_ENM[3*pair[l].posI+a][3*pair[l].posJ+b]);
		  }
		else// if (a==b)
		  {
		    double ondiagonal = pair[l].k*(        (  coeff1   )          *(1- (  (coordDif[a])*(coordDif[a])/disSqMob_ij    )   )-1);
		    //printf("ondiag: %lf\n", ondiagonal);    
		    H_ENM[Ind_I+a][Ind_J+b]=ondiagonal;
		    H_ENM[Ind_J+b][Ind_I+a]=ondiagonal;
		    //printf("pair[%d].posI: %d\n", l, pair[l].posI);
		    //printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		    //printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), hessian[3*pair[l].posI+a][3*pair[l].posJ+b]);
		  }
	      }	
	}
    }

  //Build diagonal blocks of hessian
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      H_ENM[3*i+a][3*i+b] -= H_ENM[3*i+a][3*j+b]; 
 	      //printf("Hessian[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_ENM[3*i+a][3*i+b]); 
 	    } 
  // printf("\nCalculated Hessian Matrix\n");
  
  if(ENM_DEBUG)
    {
      //Print gradient
      FILE *f=fopen("gradient_ENM.txt", "w");
      if(f==NULL)
	{
	  fprintf(stderr, "Can not produce gradient.txt file");
	  exit(EXIT_FAILURE);
	}
      for(l=0; l<3*N; l++)
	fprintf(f, "Gradient[%d]=%lf\n", l, grad_ENM[l]);
      fclose(f);
      //Print lower triangle of hessian
      char outFileName[14];
      memset(outFileName, '\0', 14);
      sprintf(outFileName, "%.8s%d%.4s", "mtrlConf", conformationNo, ".txt");
      // printf("Output:%s\n", outFileName);
      
      FILE *matrix=fopen(outFileName, "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "%s file could not be produced\n", outFileName);
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
 	{ 
 	  for(j=0; j<=i; j++)
	    for(a=0; a<3; a++) 
	      for(b=0; b<3; b++) 
		if(H_ENM[3*i+a][3*j+b]!=0.0)
		  fprintf(matrix, "%d\t%d\t%lf\n", 3*i+a, 3*j+b, H_ENM[3*i+a][3*j+b]); 
	  
	  //fprintf(matrix, "\n"); 
 	} 
      
      fclose(matrix);
    }

  //  fprintf(stdout, "ENM energy=%lf\n", 0.5*E_ENM);
  return (0.5*E_ENM);
}

double gradientANDhessian_mENM(int conformationNo, int N/*Number of residues*/, 
			       double *grad_ENM, double **H_ENM, 
			       CAcoord *atomMoving, 
			       int *pairNo, pairInfoNew *pair,
			       int *pairInitialNo, pairInfoNew *pairInitial, double exponent)
{

  //Purpose: This function is based on modified Elastic Network Model explained in 
  //         Zheng W. Accurate Flexible Fitting of High-Resolution Protein Structures
  //         into Cryo-Electron Microscopy Maps Using Coarse-Grained Pseudo-Energy Minimization.
  //         Biophys. J. 100, 478-488 (2011).

  int i=0, j=0, l=0;        /*Counters for atom pairs                */
  int a=0, b=0;             /*Counters for atom cartesian coordinates*/
  
  double fx=0.0, fy=0.0, fz=0.0;
  double coeff1=0.0, coeff2=0.0, coeff3=0.0;
  double d_ij_sqrd=0.0, inv_d_ij_sqrd=0.0, d_ij_pow=0.0;
  double d_ij0_sqrd=0.0, d_ij0_pow=0.0;
  double coordDif[3]={0.0, 0.0, 0.0};

  for (i=0; i< (3*N); i++) /*To zero input hessian and gradient so that previous calculation have nothing to do w/ this one*/
    {
      for (j=0; j<=i; j++)
	{
	  H_ENM[i][j]=0.0;
	  H_ENM[j][i]=0.0;
	}
      grad_ENM[i]=0.0;
    }
  double E_ENM=0.0;
  double oneMinusDisRate=0.0;
  double disRate=0.0;
  double invExp    =(1.0/exponent);
  double invExpSqrd=(1.0/(exponent*exponent));
  for (l=0;l<pairNo[0]; ++l) 
    {
      if(pair[l].k!=0.0)
	{      
	  d_ij0_sqrd=(pair[l].dis*pair[l].dis);                 //
	                                          
	  d_ij0_pow=pow(pair[l].dis, exponent);
	  d_ij_pow=pow(pairInitial[l].dis, exponent);

	  disRate=(d_ij0_pow/d_ij_pow);
	  oneMinusDisRate=(1-disRate);

	  //Just test with this value: This is distance based force constant and it is obsolete.
	  //pair[l].k=(16.0/d_ij0_sqrd);

	  E_ENM+=(pair[l].k*d_ij0_sqrd*invExpSqrd*oneMinusDisRate*oneMinusDisRate);
	  
	  //disMobile_ij=dis(atomMoving, (pair[l].posI), (pair[l].posJ));
	  d_ij_sqrd=(   (pairInitial[l].dis)*(pairInitial[l].dis)   );
	  inv_d_ij_sqrd=(    1.0/(  (pairInitial[l].dis)*(pairInitial[l].dis)  )    );
	  
	  coeff1= (pair[l].k)*invExp*disRate*d_ij0_sqrd*inv_d_ij_sqrd;
	  coeff2= oneMinusDisRate*coeff1;
	  coeff3= (pair[l].k)*(1.0 - coeff1 );

	  coordDif[0]=(atomMoving[pair[l].posI].x-atomMoving[pair[l].posJ].x);
  	  coordDif[1]=(atomMoving[pair[l].posI].y-atomMoving[pair[l].posJ].y);
	  coordDif[2]=(atomMoving[pair[l].posI].z-atomMoving[pair[l].posJ].z);

	  fx=  coeff2*(coordDif[0]);
	  fy=  coeff2*(coordDif[1]);
	  fz=  coeff2*(coordDif[2]);
	  
	  int Ind_I=3*(pair[l].posI);
	  int Ind_J=3*(pair[l].posJ);

	  grad_ENM[Ind_I]+=fx;
	  grad_ENM[Ind_I+1]+=fy;
	  grad_ENM[Ind_I+2]+=fz;
	  
	  grad_ENM[Ind_J]-=fx;
	  grad_ENM[Ind_J+1]-=fy;
	  grad_ENM[Ind_J+2]-=fz;
	  
	    //fprintf(stdout, "Pair Initial Distance=%lf\tDistance calculated=%lf\n", pairInitial[l].dis, disMobile_ij);
	  for(a=0; a<3; ++a)
	    for(b=0; b<3; ++b) 
	      {
		double offdiagonal = coeff1 * inv_d_ij_sqrd * ( (exponent + 2.0) - (2.0*exponent + 2.0)*disRate ) * coordDif[a] * coordDif[b]; 
		if(a!=b)
		  {
		    H_ENM[Ind_I+a][Ind_J+b]=offdiagonal;
		    H_ENM[Ind_J+b][Ind_I+a]=offdiagonal;
		    //printf("pair[%d].posI: %d\n", l, pair[l].posI);
		    //printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		    //printf("offdiag: %lf\n", offdiagonal);    
		    //printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), H_ENM[3*pair[l].posI+a][3*pair[l].posJ+b]);
		  }
		else// if (a==b)
		  {
		    double ondiagonal = -coeff2 + offdiagonal;
		    //printf("ondiag: %lf\n", ondiagonal);
		    H_ENM[Ind_I+a][Ind_J+b]=ondiagonal;
		    H_ENM[Ind_J+b][Ind_I+a]=ondiagonal;
		    //printf("pair[%d].posI: %d\n", l, pair[l].posI);
		    //printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		    //printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), hessian[3*pair[l].posI+a][3*pair[l].posJ+b]);
		  }
	      }	
	}
    }

  //Build diagonal blocks of hessian
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      H_ENM[3*i+a][3*i+b] -= H_ENM[3*i+a][3*j+b]; 
 	      //printf("Hessian[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_ENM[3*i+a][3*i+b]); 
 	    } 
  // printf("\nCalculated Hessian Matrix\n");
  
  if(ENM_DEBUG)
    {
      //Print gradient
      FILE *f=fopen("gradient_ENM.txt", "w");
      if(f==NULL)
	{
	  fprintf(stderr, "Can not produce gradient.txt file");
	  exit(EXIT_FAILURE);
	}
      for(l=0; l<3*N; l++)
	fprintf(f, "Gradient[%d]=%lf\n", l, grad_ENM[l]);
      fclose(f);
      //Print lower triangle of hessian
      char outFileName[14];
      memset(outFileName, '\0', 14);
      sprintf(outFileName, "%.8s%d%.4s", "mtrlConf", conformationNo, ".txt");
      // printf("Output:%s\n", outFileName);
      
      FILE *matrix=fopen(outFileName, "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "%s file could not be produced\n", outFileName);
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
 	{ 
 	  for(j=0; j<=i; j++)
	    for(a=0; a<3; a++) 
	      for(b=0; b<3; b++) 
		if(H_ENM[3*i+a][3*j+b]!=0.0)
		  fprintf(matrix, "%d\t%d\t%lf\n", 3*i+a, 3*j+b, H_ENM[3*i+a][3*j+b]); 
	  
	  //fprintf(matrix, "\n"); 
 	} 
      
      fclose(matrix);
    }

  //  fprintf(stdout, "ENM energy=%lf\n", 0.5*E_ENM);
  return (0.5*E_ENM);
}

double gradientANDhessian_Collision(int conformationNo, int N/*Number of residues*/, 
				    double *grad_col, double **H_col, 
				    CAcoord *atomMoving, 
				    int *pairNo, pairInfoNew *pair,
				    int *pairInitialNo, pairInfoNew *pairInitial, 
				    double R_collision)
{
  int i=0, j=0, l=0; /*Counters for atom pairs                */
  int a=0, b=0;             /*Counters for atom cartesian coordinates*/
  int ind_I=0, ind_J=0;     /*For temporary internal indices*/
  double fx=0.0, fy=0.0, fz=0.0;
  double coeff1=0.0, coeff2=0.0, coeff3=0.0;
  double d_ij=0.0, d_ij_sqrd=0.0, inv_d_ij_sqrd=0.0;
  double coordDif[3]={0.0, 0.0, 0.0};

  for (i=0; i< (3*N); i++) /*To zero input hessian and gradient so that previous calculation have nothing to do w/ this one*/
    {
      for (j=0; j<=i; j++)
	{
	  H_col[i][j]=0.0;
	  H_col[j][i]=0.0;
	}
      grad_col[i]=0.0;
    }
  double E_col=0.0;

  double dis_diff=0.0;
  for (l=0;l<pairNo[0]; l++) 
    {
      if(pair[l].k==FR_CON) //This condition is set to just check non-bonded pairs
	{
	  //	  if(abs(pair[l].posI-pair[l].posJ)>1)
	  if(pairInitial[l].dis<R_collision)
	    {
	      dis_diff=(R_collision - pairInitial[l].dis);
	      E_col+=C_COL*dis_diff*dis_diff;
	      
	      d_ij=pairInitial[l].dis;
	      d_ij_sqrd=(d_ij*d_ij);
	      inv_d_ij_sqrd=(1.0/d_ij_sqrd);
	      
	      coeff1= R_collision/d_ij;
	      coeff2= (C_COL)*coeff1;
	      coeff3= (C_COL)*(1- coeff1 );
	      
	      coordDif[0]=(atomMoving[pair[l].posI].x-atomMoving[pair[l].posJ].x);
	      coordDif[1]=(atomMoving[pair[l].posI].y-atomMoving[pair[l].posJ].y);
	      coordDif[2]=(atomMoving[pair[l].posI].z-atomMoving[pair[l].posJ].z);
	      
	      fx=  coeff3*(coordDif[0]);
	      fy=  coeff3*(coordDif[1]);
	      fz=  coeff3*(coordDif[2]);
	      
	      ind_I=3*(pair[l].posI);
	      grad_col[ind_I]+=fx;
	      grad_col[ind_I+1]+=fy;
	      grad_col[ind_I+2]+=fz;
	      
	      ind_J=3*(pair[l].posJ);
	      grad_col[ind_J]-=fx;
	      grad_col[ind_J+1]-=fy;
	      grad_col[ind_J+2]-=fz;
	      
	      //========Copy hessian here
	      //fprintf(stdout, "Pair Initial Distance=%lf\tDistance calculated=%lf\n", pairInitial[l].dis, disMobile_ij);
	      for(a=0; a<3; ++a)
		for(b=0; b<3; ++b) 
		  {
		    if(a!=b)
		      {
			double offdiagonal = -( coeff2 )   *(   ( coordDif[a] ) * ( coordDif[b]) *inv_d_ij_sqrd); 
			H_col[ind_I+a][ind_J+b]=offdiagonal;
			H_col[ind_J+b][ind_I+a]=offdiagonal;
			//printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), H_col[3*pair[l].posI+a][3*pair[l].posJ+b]);
		      }
		    else//if(a==b)
		      {
			double ondiagonal = C_COL*(    coeff1*(1.0 - (  (coordDif[a])*(coordDif[a])*inv_d_ij_sqrd)  ) - 1.0);
			H_col[ind_I+a][ind_J+b]=ondiagonal;
			H_col[ind_J+b][ind_I+a]=ondiagonal;
			//printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), hessian[3*pair[l].posI+a][3*pair[l].posJ+b]);
		      }
		  }
	    }
	}
    }
  
  //Build diagonal blocks of hessian
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      H_col[3*i+a][3*i+b] -= H_col[3*i+a][3*j+b]; 
 	      //printf("Hessian[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_col[3*i+a][3*i+b]); 
 	    } 
  // printf("\nCalculated Hessian Matrix\n");
  
  if(ENM_DEBUG)
    {
      //Print gradient
      FILE *f=fopen("gradient_E_col.txt", "w");
      if(f==NULL)
	{
	  fprintf(stderr, "Can not produce gradient.txt file");
	  exit(EXIT_FAILURE);
	}
      for(l=0; l<3*N; l++)
	fprintf(f, "Gradient[%d]=%lf\n", l, grad_col[l]);
      fclose(f);
      //Print lower triangle of hessian
      char outFileName[14];
      memset(outFileName, '\0', 14);
      sprintf(outFileName, "%.8s%d%.4s", "mtrlConf", conformationNo, ".txt");
      printf("Output:%s\n", outFileName);
      
      FILE *matrix=fopen(outFileName, "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "%s file could not be produced\n", outFileName);
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++)
 	{
 	  for(j=0; j<=i; j++)
	    for(a=0; a<3; a++)
	      for(b=0; b<3; b++)
		if(H_col[3*i+a][3*j+b]!=0.0)
		  fprintf(matrix, "%d\t%d\t%lf\n", 3*i+a, 3*j+b, H_col[3*i+a][3*j+b]);
	  
	  //fprintf(matrix, "\n");
 	}
      
      fclose(matrix);
      fprintf(stdout, "Collision energy=%lf\n", 0.5*E_col);
    }
    
  return (0.5*E_col);
}
