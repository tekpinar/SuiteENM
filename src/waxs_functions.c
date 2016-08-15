//Purpose: waxs_functions.c contains functions I developed to calculate WAXS 
//         intensities of proteins. This calculation method is explained in 
//         Park, S., J. P. Bardhan, B. Roux, and L. Makowski (2009) Simulated X-Ray 
//         Scattering of Protein Solutions Using Explicit-Solvent Molecular Dynamics. 
//         J. Chem. Phys. 130, 134114. PMID: 19355724 
//Author : Mustafa Tekpinar
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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <string.h>
//Personal header files
#include <structures.h>
#include <pdb_io.h>
#include <form_factors.h>

void q_x(int J_max, double q, double q_x_array[J_max])
{
  //Purpose: To calculate x component of q vector for J_max data points. 
  //Instead of calculating those values, lets keep them in an array and use whenever needed!
  int    j=0;
  double t=0;
  double coeff=sqrt(M_PI*J_max);
  double J_maxPlus1=(1.0+(double)J_max);
  for(j=1; j<=J_max; j++)
    {
      //      t=(  ((double)(2*j-1-J_max)) / ((double)(J_max))  );
      t=(  ((double)(2*j - J_maxPlus1)) / ((double)(J_max))  );
      
      q_x_array[j-1]=(  q*sin( acos(t)  )*cos( coeff*asin(t) )  );
      
      if(0)
	fprintf(stdout, "q_x[%d]=%lf\n", (j-1), q_x_array[j-1]);
    }

}

void q_y(int J_max, double q, double q_y_array[J_max])
{
  //Purpose: To calculate y component of q vector for J_max data points. 
  //Instead of calculating those values, lets keep them in an array and use whenever needed!
  int    j=0;
  double t=0.0;
  double coeff=sqrt(M_PI*J_max);
  double J_maxPlus1=(1.0+(double)J_max);
  for(j=1; j<=J_max; j++)
    {
      //      t=(  ((double)(2*j-1-J_max)) / ((double)(J_max))  );
      t=(  ((double)(2*j-J_maxPlus1)) / ((double)(J_max))  );
      q_y_array[j-1]=(  q*sin( acos(t)  )*sin( coeff*asin(t) )  );

      if(0)
	fprintf(stdout, "q_y[%d]=%lf\n", (j-1), q_y_array[j-1]);
    }

}

void q_z(int J_max, double q, double q_z_array[J_max])
{
  //Purpose: To calculate z component of q vector for J_max data points. 
  //Instead of calculating those values, lets keep them in an array and use whenever needed!
  int    j=0;
  double t=0.0;
  double J_maxPlus1=(1.0+(double)J_max);
  for(j=1; j<=J_max; j++)
    {
      //      t=(  ((double)(2*j-1-J_max)) / ((double)(J_max))  );
      t=(  ((double)(2*j-J_maxPlus1)) / ((double)(J_max))  );
      q_z_array[j-1]=(q*t) ;

      if(0)
	fprintf(stdout, "q_z[%d]=%lf\n", (j-1), q_z_array[j-1]);
    }

}

gsl_complex X1_q_tilde_v1(int N/*Number of atoms*/, pdb_v23 *info, int k/*A replacement for j.*/, int J_max, double q, 
			  double q_x_array[J_max], double q_y_array[J_max], double q_z_array[J_max], double *f_l)
{
  //Normally j starts from 1 to J_max. The k here I use starts from 0 to J_max-1. 
  int i=0; //My dear counter. 
  gsl_complex allAtoms; //This will be used for sum.
  
  allAtoms=gsl_complex_rect(0.0, 0.0); //Zero whatever inside!
  
  gsl_complex oneAtom; //And this will be used to calculate per atom. 
  gsl_complex exponent;//This is the complex number that will be used as exponent of exponential function in eq. 30 of reference cited at the top.
  double tempImag;     //Just a temporary double variable.
  for(i=0; i<N; i++)
    {
      tempImag=(-1.0*( q_x_array[k]*info[i].x  + q_y_array[k]*info[i].y + q_z_array[k]*info[i].z)  );
      
      exponent=gsl_complex_rect(0.0, tempImag);
      
      oneAtom=gsl_complex_rect(0.0, 0.0); //Zero whatever inside!
      oneAtom=gsl_complex_mul_real(gsl_complex_exp(exponent), f_l[i]);
      allAtoms=gsl_complex_add(allAtoms, oneAtom);
    }
  return allAtoms;  
}

void X1_q_tilde_allS_j_v1(int totalFrameNum, char framesList[][255], int J_max, \
			  double q, double waterWeight, double rho_s, int isSolventOn, \
			  double q_x_array[J_max], double q_y_array[J_max], double q_z_array[J_max], gsl_complex **X1_S)
{
  //Lets try to read all snapshots in a for loop!!!
  int i=0, k=0;
  int atomCount[5]={99999/*Atoms*/, 999/*Residues*/, 99/*Chains*/, 100/*Helices*/, 99/*Sheets*/};
  for(i=0; i<totalFrameNum; i++)
    {
      //--Scan first conformation
      FILE *pdbdata=fopen(framesList[i],"r");
      if(pdbdata==NULL)
	{
	  fprintf(stderr, "No such file: %s\n", framesList[i]);
	  exit(EXIT_FAILURE);
	}
      fprintf(stdout, "%s\n", framesList[i]);
      scanpdb_CAandOH2(atomCount, framesList[i]);
      printf("Number of residues:%d\n",  atomCount[1]);
      
      //--Get data from first pdb file--
      pdb_v23 *info = (pdb_v23*)malloc((atomCount[0])*sizeof(pdb_v23));
      if (info == NULL)
	{
	  fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb\n");
	  exit(EXIT_FAILURE);
	}
      
      CAcoord *protAndWater = (CAcoord*)malloc((atomCount[1])*sizeof(CAcoord));
      if (protAndWater == NULL)
	{
	  fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord\n");
	  exit(EXIT_FAILURE);
	}
      
      readpdb_CAandOH2(info, protAndWater, framesList[i]);
      fclose(pdbdata);
      //Now, lets calculate form factors for a certain q on both of pdb files.
      double *f_l= (double*) calloc ((atomCount[0]), sizeof(double));
      if (f_l==NULL)
	{
	  fprintf(stderr, "No memory callocation for f_l calculation\n");
	  exit (EXIT_FAILURE);
	}
      
      formFactorAllAtoms(atomCount[0], info, waterWeight, f_l, q, rho_s, isSolventOn);
      for(k=0; k<J_max; k++)
	{ 
	  //Here, i is the value for snapshot and k is the value at a spesific j according to paper. 
	  //However, I have to note again that k starts from 0 and j from 1.        
	  X1_S[i][k]=X1_q_tilde_v1(atomCount[0], info, k/*A replacement for j.*/, J_max, q, q_x_array, q_y_array, q_z_array, f_l);
	  if(0) 
	    printf("Structure: %d - J value=%d\t\tReal part= %lf and complex part= %lf\n", i, k, GSL_REAL(X1_S[i][k]), GSL_IMAG(X1_S[i][k]));
	}
      
      free(f_l);
      free(protAndWater);
      free(info);
    }

}

gsl_complex X1_q_tilde_v2(int frameNo, int frm_beg, int frm_end, float *crdX, float *crdY, float *crdZ,\
			  int k/*A replacement for j.*/, int J_max, double q, \
			  double q_x_array[J_max], double q_y_array[J_max], double q_z_array[J_max], double *f_l)
{
  //Normal j starts from 1 to J_max. The k here I use starts from 0 to J_max-1. 
  int i=0; //My dear counter. 
  gsl_complex allAtoms; //This will be used for sum.
  
  allAtoms=gsl_complex_rect(0.0, 0.0); //Zero whatever inside!
  
  double temp_arg=0.0;     //Just a temporary double variable.
  //  double temp_cos=0.0;     //Just a temporary double variable.
  //  double temp_sin=0.0;     //Just a temporary double variable.
  double sumReal=0.0;
  double sumImag=0.0;
  for(i=frm_beg; i<frm_end; i++)
    {
      temp_arg=( ( q_x_array[k]*crdX[i]  + q_y_array[k]*crdY[i] + q_z_array[k]*crdZ[i])  );
      //      temp_cos=cos(temp_arg);      
      //      temp_sin=sqrt(1.0 -  (temp_cos*temp_cos));      
      
      //      sumReal+=(f_l[i]*temp_cos);
      //      sumImag+=(f_l[i]*temp_sin);
      
      sumReal+=(f_l[i]*cos(temp_arg));
      sumImag+=(f_l[i]*sin(temp_arg));
    }

  allAtoms=gsl_complex_rect(sumReal, (-sumImag)); 
  return allAtoms;  
}

/* void parseSWAXSexperimental(int n_exp, char *experimentalFile, double *exp_q, double *exp_I_q, double *sigma_q) */
/* { */
/*   int i=0; */

/*   char line[255]; */
/*   FILE *EXP_FILE=fopen(experimentalFile, "r"); */
/*   if(EXP_FILE==NULL) */
/*     { */
/*       fprintf(stderr, "ERROR: Experimental data file %s not found!\n", experimentalFile); */
/*       exit(EXIT_FAILURE); */
/*     } */
/*   memset(line, '\0', 255); */
/*   fprintf(stdout, "Number of data points in experimental file is %d.\n", i); */
/*   rewind(EXP_FILE); */
/*   memset(line, '\0', 255); */
/*   i=0; */
/*   //Read experimental data to three arrays */
/*   while(1) */
/*     { */
/*       if(fgets(line, sizeof(line), EXP_FILE)==NULL) */
/* 	{ */
/* 	  break; */
/* 	} */
/*       else if (line[0]=='#') //Do not read comment line in experimental data file! */
/* 	{ */
/* 	  fprintf(stderr, "ERROR: This file contanins non data lines:\n%s\n", line); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
/*       else */
/* 	{ */
/* 	  //	  printf("%s", line); */
/* 	  sscanf(line, "%lf %lf %lf\n", &exp_q[i], &exp_I_q[i], &sigma_q[i]); */
/* 	  printf("%lf\t%lf\t%lf\n", exp_q[i], exp_I_q[i], sigma_q[i]); */
/* 	  i++; */
/* 	} */
/*     } */
/*   fclose(EXP_FILE); */
/* } */
