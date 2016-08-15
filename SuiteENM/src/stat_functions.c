//Purpose: stat_functions.c contains statistical analysis functions for various purposes.
//Author : Mustafa Tekpinar
//Date   : 04/22/2011
//Revisions: 
//     v.00: average(), standardDeviation(), percentError() functions added. 
//     v.01: chiSquareTest(), scalingCoefficient(), and linterpolate() added at 10 May 2011. 
//     v.02: chiSquareLogTest() added at 17 August 2011. 
//     v.03: chi_of_data_files_v2() function was added to calculate chi between experimental
//           and theoretical data. Theoretical data points has been produced at exactly experimental
//           q data points. 11 October 2011.
//     v.04: double chi_of_data_files(char *experimentalFile, char *theoreticalFile,  int scaleTheoretical/*0: NO, 1:YES*/)
//           This function calculates chi value no by producing theoretical values at experimental data points!
//           It has been added to library at 2 January 2012!
//     v.05: double chi_of_data_files_v3(char *experimentalFile, char *theoreticalFile,  double min_range, double max_range, 
//			                 int scaleTheoretical/*0: NO, 1:YES*/, int printDetails/*0: NO, 1:YES*/, FILE *FIT_FILE)
//           function is added to make fitting procedure much more flexible. 12 January 2012

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
double average(int n/*Number of elements in array*/, double *values)
{
  int i=0;//My dear regular counter.
  double average=0.0;
  for(i=0; i<n; i++)
    {
      average+=(values[i]);
    }
  return (average/((double)n));
}

double standardDeviation(int n/*Number of elements in array*/, double *values)
{
  double ave=average(n, values);
  double stdDev=0.0;
  int i=0;
  for(i=0; i<n; i++)
    {
      stdDev+=( (values[i]-ave) * (values[i]-ave) );
    }
  return (sqrt ( stdDev/(double)n));
}

double percentError(double measured, double actual)
{
  double error=0.0;
  error=fabs((measured-actual)/actual);
  return error;
}

double chiSquareTest(int n_exp, double *experimental, double *theoretical, double *sigma, double c)
{
  //Purpose: To calculate chiSquare value between experimental and theoretical data using Equation 15 of CRYSOL paper. 
  //         int n_exp           : Number of experimantal data points. 
  //         double *experimental: The experimental data array. 
  //         double *theoretical : The theoretical data array. 
  //         double *sigma       : The experimental error array. 
  //         double c            : A scaling coefficient. 

  int    i=0;
  double chiSquare=0.0; 
  for(i=0; i<n_exp; i++)
    {
      double diff=( ( experimental[i]-c*theoretical[i] ) / sigma[i] );
      chiSquare+= ( diff*diff );
    }
  chiSquare=(chiSquare/(double)n_exp);
  return chiSquare;
}

double chiSquareLogTest(int n_exp, double *experimental, double *theoretical, double *sigma)
{
  //Purpose: To calculate chiSquare value between experimental and theoretical data using Equation 15 of CRYSOL paper.
  //         int n_exp           : Number of experimantal data points.
  //         double *experimental: The experimental data array.
  //         double *theoretical : The theoretical data array.
  //         double *sigma       : The experimental error array.
  //         double delta        : Normalization shift amount.

  int    i=0;
  double chiSquare=0.0;
  double delta=fabs( experimental[0] - log10(theoretical[0]));
  printf("Delta=%lf\n", delta);
  for(i=0; i<n_exp; i++)
    {
      printf("%d\t%lf\t%lf\n", i, experimental[i], (log10(theoretical[i])  - delta));
      double diff=( ( experimental[i] - log10(theoretical[i])  - delta ) / sigma[i] );
      chiSquare+= ( diff*diff );
    }
  //It can be divided by n_exp or one can skip it!
  chiSquare=(chiSquare/(double)n_exp);
  return chiSquare;
}

double scalingCoefficient(int n_exp, double *experimental, double *theoretical, double *sigma)
{
  //Purpose: To calculate scaling coefficient described in CRYSOL paper, equation 16.
  //         int n_exp           : Number of experimantal data points. 
  //         double *experimantal: The experimental data array. 
  //         double *theoretical : The theoretical data array. 
  //         double *sigma       : The experimental error array. 
  //         double c            : A scaling coefficient. 
  int    i=0;
  double c=0.0;
  double tempSum1=0.0; 
  double tempSum2=0.0; 
  for(i=0; i<n_exp; i++)
    {
      double sigmaSqrd=(sigma[i]*sigma[i]);
      tempSum1+=( ( experimental[i]*theoretical[i] ) / sigmaSqrd ) ;
      tempSum2+=( ( theoretical[i]*theoretical[i]  ) / sigmaSqrd ) ;
    }

  c=(tempSum1/tempSum2);
  return c;
}

double linterpolate(double y_a, double y_b, double x_a, double x, double x_b)
{
  //Purpose: Given a function y(x) at points a and b,
  //         find the linear interpolation value specied at x which is between a and b.
  double y=0;
  y= (   y_a + ( (y_b - y_a)*(x - x_a)/(x_b-x_a) )   );
  return y;
}
double chi_of_data_files(char *experimentalFile, char *theoreticalFile,  int scaleTheoretical/*0: NO, 1:YES*/)
{
  int i=0, j=0;
  //========================================================================================================================
  char line[255];
  FILE *EXP_FILE=fopen(experimentalFile, "r");
  if(EXP_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Experimental data file %s not found!\n", experimentalFile);
      exit(EXIT_FAILURE);
    }
  memset(line, '\0', 255);
  
  int n_exp=0; //Number of experimental data points. 
  //Scan file to determine number experimental dat points!
  while(1)
    {
      if(fgets(line, sizeof(line), EXP_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in experimental data file!
	    {
	      // printf("%s", line);
	      i++;
	    }
	}
    }
  n_exp=i;
  fprintf(stdout, "Number of data points in experimental file is %d.\n", i);
  double expArray[n_exp][3]; //First Column: q, Second Column: I(q), Third Column: Sigma(q)
  rewind(EXP_FILE);
  memset(line, '\0', 255);
  i=0;
  //Read experimental data to an array!
  while(1)
    {
      if(fgets(line, sizeof(line), EXP_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in experimental data file!
	    {
	      //	  printf("%s", line);
	      sscanf(line, "%lf %lf %lf\n", &expArray[i][0], &expArray[i][1], &expArray[i][2]);
	      printf("%lf\t%lf\t%lf\n", expArray[i][0], expArray[i][1], expArray[i][2]);
	      i++;
	    }
	}
    }

  fclose(EXP_FILE);
  //========================================================================================================================
  //Copy data to one dimensional arrays for later use!
  double exp_q[n_exp];
  double exp_I_q[n_exp];
  double sigma[n_exp];
  for(i=0; i<n_exp; i++)
    {
      exp_q[i]=expArray[i][0];
      exp_I_q[i]=expArray[i][1];
      sigma[i]=expArray[i][2];
    }

  int n_theo=0;
  //Open theoretical data file!
  FILE *THEO_FILE=fopen(theoreticalFile, "r");
  if(THEO_FILE==NULL)
    {
      fprintf(stderr, "ERROR: My dear, I couldn't read %s file!\n", theoreticalFile);
      exit(EXIT_FAILURE);
    }
  memset(line, '\0', 255);
  i=0;
  //Scan file to determine number of theoretical data points!
  while(1)
    {
      if(fgets(line, sizeof(line), THEO_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in theoretical data file!
	    {
	      //	  printf("%s", line);
	      //	  sscanf(line, "%lf %lf\n", &theoArray[i][0], &theoArray[i][1]);
	      //	  printf("%lf\t%lf\n", theoArray[i][0], theoArray[i][1]);
	      i++;
	    }
	}
    }
  n_theo=i;
  fprintf(stdout, "Number of data points in theoretical file is %d.\n", i);
  double theoArray[n_theo][2];   //Theoretical data: 1st  Column: q, 2nd Column: I(q) in solution.
  double theo_I_q[n_exp];
  //  At first, initialize it to zero.
  for(i=0; i<n_theo; i++)
    {
      theoArray[i][0]=0.0;
      theoArray[i][1]=0.0;
    }

  rewind(THEO_FILE);
  memset(line, '\0', 255);
  i=0;
  //Read theoretical data to an array!
  while(1)
    {
      if(fgets(line, sizeof(line), THEO_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in theoretical data file!
	    {
	      //	  printf("%s", line);
	      sscanf(line, "%lf %lf\n", &theoArray[i][0], &theoArray[i][1]);
	      //	  printf("%lf\t%lf\n", theoArray[i][0], theoArray[i][1]);
	      i++;
	    }
	}
    }
  fclose(THEO_FILE);

  /*   //Now, lets interpolate linearly my results to get values at experimental q data points. */
  
  for(i=0; i<n_exp; i++)
    {
      //Our 'x's are q values of experimental data.
      double x=exp_q[i];//This is the point where we want to find theoretical value!
      //Find x_a and x_b(or q_a and q_b). Lets find in which interval it is.
      for(j=0; j<n_theo; j++)
	if((x>=theoArray[j][0])&&(x<=theoArray[j+1][0]))
	  {
	    //	    printf("Exp data point %lf is between %lf  and  %lf\n",exp_q[i], theoArray[j][0], theoArray[j+1][0]);
	    break;
	  }
      theo_I_q[i]=linterpolate(theoArray[j][1], theoArray[j+1][1], theoArray[j][0], x, theoArray[j+1][0]);
    }

  //Calculate scaling coefficient and chi square value for Theoretical vs Experiment.
  double chiSqrdTestValue=0.0;
  double coeff=0.0;

  coeff=scalingCoefficient(n_exp, exp_I_q, theo_I_q, sigma);
  chiSqrdTestValue=chiSquareTest(n_exp, exp_I_q, theo_I_q, sigma, coeff/*c=1.0=no scaling*/);
  
  //Just in case if you want to see the fit visually!!
  if(scaleTheoretical)
    for(i=0; i<n_exp; i++)
      {
	theo_I_q[i]=coeff*(theo_I_q[i]);
	
	printf("%lf\t%lf\n", exp_q[i], theo_I_q[i]);
	// crysArray[i][1]=log10(crysArray[i][1]);
	// printf("%lf\t%lf\n", crysArray[i][0], crysArray[i][1]);
      }

  printf("Chi Value=%lf\tScaling coefficient:%lf\n", sqrt(chiSqrdTestValue), coeff);
  return (sqrt(chiSqrdTestValue));
}

double chi_of_data_files_v2(char *experimentalFile, char *theoreticalFile,  int scaleTheoretical/*0: NO, 1:YES*/)
{

  //Whats new: In this version, instead of interpolating theoretical data, I calculate intensities exactly at experimental data points. 
  //It is hoped that it will give me a better fit to experimental!!!

  int i=0;
  //========================================================================================================================
  char line[255];
  FILE *EXP_FILE=fopen(experimentalFile, "r");
  if(EXP_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Experimental data file %s not found!\n", experimentalFile);
      exit(EXIT_FAILURE);
    }
  memset(line, '\0', 255);
  
  int n_exp=0; //Number of experimental data points. 
  //Scan file to determine number experimental dat points!
  while(1)
    {
      if(fgets(line, sizeof(line), EXP_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment lines in experimental data file!
	    {
	      // printf("%s", line);
	      i++;
	    }
	}
    }
  n_exp=i;
  fprintf(stdout, "Number of data points in experimental file is %d.\n", i);
  
  rewind(EXP_FILE);
  memset(line, '\0', 255);
  i=0;
  //Read experimental data to three arrays!
  //First Column: q, Second Column: I(q), Third Column: Sigma(q)
  double exp_q[n_exp];
  double exp_I_q[n_exp];
  double sigma[n_exp];

  while(1)
    {
      if(fgets(line, sizeof(line), EXP_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in experimental data file!
	    {
	      //	  printf("%s", line);
	      //	      sscanf(line, "%lf %lf %lf\n", &expArray[i][0], &expArray[i][1], &expArray[i][2]);
	      sscanf(line, "%lf %lf %lf\n", &exp_q[i], &exp_I_q[i], &sigma[i]);
	      //	      printf("%lf\t%lf\t%lf\n", expArray[i][0], expArray[i][1], expArray[i][2]);
	      i++;
	    }
	}
    }

  fclose(EXP_FILE);
  //========================================================================================================================

  int n_theo=0;
  //Open theoretical data file!
  FILE *THEO_FILE=fopen(theoreticalFile, "r");
  if(THEO_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Dear user, I couldn't read %s file!\n", theoreticalFile);
      exit(EXIT_FAILURE);
    }
  memset(line, '\0', 255);
  i=0;
  //Scan file to determine number of theoretical data points!
  while(1)
    {
      if(fgets(line, sizeof(line), THEO_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in theoretical data file!
	    {
	      //	  printf("%s", line);
	      //	  sscanf(line, "%lf %lf\n", &theoArray[i][0], &theoArray[i][1]);
	      //	  printf("%lf\t%lf\n", theoArray[i][0], theoArray[i][1]);
	      i++;
	    }
	}
    }
  n_theo=i;

  fprintf(stdout, "Number of data points in theoretical file is %d.\n", i);

  if(n_theo!=n_exp)
    {
      fprintf(stderr, "ERROR: Number of data points in experimental and theoretical files is supposed to be same!");
      exit(EXIT_FAILURE);
    }

  //Theoretical data: 1st  Column: q, 2nd Column: I(q) in solution.
  double theo_I_q[n_exp];
  double theo_q[n_exp];
  //  At first, initialize it to zero.
  for(i=0; i<n_theo; i++)
    {
      theo_q[i]=0.0;
      theo_I_q[i]=0.0;
    }

  rewind(THEO_FILE);
  memset(line, '\0', 255);
  i=0;
  //Read theoretical data to arrays!
  while(1)
    {
      if(fgets(line, sizeof(line), THEO_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in theoretical data file!
	    {
	      //	  printf("%s", line);
	      sscanf(line, "%lf %lf\n", &theo_q[i], &theo_I_q[i]);
	      //	  printf("%lf\t%lf\n", theo_q[i], theo_I_q[i]);
	      i++;
	    }
	}
    }
  fclose(THEO_FILE);
  
  //================================================================================================
  //Calculate scaling coefficient and chi square value for Theoretical vs Experiment.
  double chiSqrdTestValue=0.0;
  double coeff=0.0;

  coeff=scalingCoefficient(n_exp, exp_I_q, theo_I_q, sigma);
  chiSqrdTestValue=chiSquareTest(n_exp, exp_I_q, theo_I_q, sigma, coeff/*c=1.0=no scaling*/);
  
  //Just in case if you want to see the fit visually!!
  if(scaleTheoretical)
    {
      FILE *THEO_FIT_FILE=fopen("I_q_makowski_waxs_fit.dat", "w");
      if(THEO_FIT_FILE==NULL)
	{
	  fprintf(stderr, "I_q_makowski_waxs_fit.dat file can not be produced!\n");
	  exit(EXIT_FAILURE);
	}

      for(i=0; i<n_exp; i++)
	{
	  theo_I_q[i]=coeff*(theo_I_q[i]);
	  fprintf(THEO_FIT_FILE, "%lf\t%lf\n", exp_q[i], theo_I_q[i]);
	}
      fclose(THEO_FIT_FILE);
    }
  printf("Chi Value=%lf\tScaling coefficient:%lf\n", sqrt(chiSqrdTestValue), coeff);
  return (sqrt(chiSqrdTestValue));
}
double chi_of_data_files_v3(char *experimentalFile, char *theoreticalFile,  double min_range, double max_range, 
			    int scaleTheoretical/*0: NO, 1:YES*/, int printDetails/*0: NO, 1:YES*/, FILE *FIT_FILE)
{
  int i=0, j=0;
  //========================================================================================================================
  char line[255];
  FILE *EXP_FILE=fopen(experimentalFile, "r");
  if(EXP_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Experimental data file %s not found!\n", experimentalFile);
      exit(EXIT_FAILURE);
    }
  memset(line, '\0', 255);
  
  int n_exp=0; //Number of experimental data points. 
  //Scan file to determine number experimental dat points!
  while(1)
    {
      if(fgets(line, sizeof(line), EXP_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in experimental data file!
	    {
	      // printf("%s", line);
	      i++;
	    }
	}
    }
  n_exp=i;
  fprintf(stdout, "Number of data points in experimental file is %d.\n", i);
  double expArray[n_exp][3]; //First Column: q, Second Column: I(q), Third Column: Sigma(q)
  rewind(EXP_FILE);
  memset(line, '\0', 255);
  i=0;
  //Read experimental data to an array!
  while(1)
    {
      if(fgets(line, sizeof(line), EXP_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in experimental data file!
	    {
	      //	  printf("%s", line);
	      sscanf(line, "%lf %lf %lf\n", &expArray[i][0], &expArray[i][1], &expArray[i][2]);
	      if(printDetails)
		printf("Original experimental data: q: %lf\tI(q): %lf\tSigma(q): %lf\n", expArray[i][0], expArray[i][1], expArray[i][2]);
	      i++;
	    }
	}
    }

  fclose(EXP_FILE);
  //========================================================================================================================
  //Copy data to one dimensional arrays for later use!
  double exp_q[n_exp];
  double exp_I_q[n_exp];
  double sigma[n_exp];
  int reCounter=0;

  fprintf(stdout, "WARNING: Data between %lf and %lf will be fit!\n", min_range, max_range);

  for(i=0; i<n_exp; i++)
    {
      //To fit just a certain range of data! 
      if( (expArray[i][0]>min_range) && (expArray[i][0]<max_range) )
	{
	  if(printDetails)
	    printf("Data point %lf has been selected for fitting!\n", expArray[i][0]);
	      
	  exp_q[reCounter]=expArray[i][0];
	  exp_I_q[reCounter]=expArray[i][1];
	  sigma[reCounter]=expArray[i][2];
	  reCounter++;
	}
    }


  n_exp=reCounter; //Reassign number of data points!

  if(printDetails)
    printf("Number of data points to be fitted is %d\n", reCounter);
  int n_theo=0;
  //Open theoretical data file!
  FILE *THEO_FILE=fopen(theoreticalFile, "r");
  if(THEO_FILE==NULL)
    {
      fprintf(stderr, "ERROR: My dear, I couldn't read %s file!\n", theoreticalFile);
      exit(EXIT_FAILURE);
    }
  memset(line, '\0', 255);
  i=0;
  //Scan file to determine number of theoretical data points!
  while(1)
    {
      if(fgets(line, sizeof(line), THEO_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in theoretical data file!
	    {
	      //	  printf("%s", line);
	      //	  sscanf(line, "%lf %lf\n", &theoArray[i][0], &theoArray[i][1]);
	      //	  printf("%lf\t%lf\n", theoArray[i][0], theoArray[i][1]);
	      i++;
	    }
	}
    }
  n_theo=i;
  fprintf(stdout, "Number of data points in theoretical file is %d.\n", i);
  double theoArray[n_theo][2];   //Theoretical data: 1st  Column: q, 2nd Column: I(q) in solution.
  double theo_I_q[n_exp];
  //  At first, initialize it to zero.
  for(i=0; i<n_theo; i++)
    {
      theoArray[i][0]=0.0;
      theoArray[i][1]=0.0;
    }

  rewind(THEO_FILE);
  memset(line, '\0', 255);
  i=0;
  //Read theoretical data to an array!
  while(1)
    {
      if(fgets(line, sizeof(line), THEO_FILE)==NULL)
	break;
      else
	{
	  if(line[0]!='#') //Do not read comment line in theoretical data file!
	    {
	      //	  printf("%s", line);
	      sscanf(line, "%lf %lf\n", &theoArray[i][0], &theoArray[i][1]);
	      if(printDetails)	      
		printf("Theoretical data: q: %lf\tI(q): %lf\n", theoArray[i][0], theoArray[i][1]);
	      i++;
	    }
	}
    }
  fclose(THEO_FILE);

  /*   //Now, lets interpolate linearly my results to get values at experimental q data points. */  
  for(i=0; i<n_exp; i++)
    {
      //Our 'x's are q values of experimental data.
      double x=exp_q[i];//This is the point where we want to find theoretical value!
      //Find x_a and x_b(or q_a and q_b). Lets find in which interval it is.
      for(j=0; j<n_theo; j++)
	if((x>=theoArray[j][0])&&(x<=theoArray[j+1][0]))
	  {
	    //	    printf("Exp data point %lf is between %lf  and  %lf\n",exp_q[i], theoArray[j][0], theoArray[j+1][0]);
	    break;
	  }
      theo_I_q[i]=linterpolate(theoArray[j][1], theoArray[j+1][1], theoArray[j][0], x, theoArray[j+1][0]);
    }

  //Calculate scaling coefficient and chi square value for Theoretical vs Experiment.
  double chiSqrdTestValue=0.0;
  double coeff=0.0;

  coeff=scalingCoefficient(n_exp, exp_I_q, theo_I_q, sigma);
  chiSqrdTestValue=chiSquareTest(n_exp, exp_I_q, theo_I_q, sigma, coeff/*c=1.0=no scaling*/);
  
  //Just in case if you want to see the fit visually!!
  if(scaleTheoretical)
    for(i=0; i<n_exp; i++)
      {
	theo_I_q[i]=coeff*(theo_I_q[i]);
	
	fprintf(FIT_FILE, "%lf\t%lf\n", exp_q[i], theo_I_q[i]);
	// crysArray[i][1]=log10(crysArray[i][1]);
	// printf("%lf\t%lf\n", crysArray[i][0], crysArray[i][1]);
      }

  printf("Chi Value=%lf\tScaling coefficient:%lf\n", sqrt(chiSqrdTestValue), coeff);
  return (sqrt(chiSqrdTestValue));
}
