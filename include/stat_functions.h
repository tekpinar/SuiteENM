//Purpose: To write statistical analysis functions and to put them all here for later use.
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

#ifndef  STAT_FUNCTIONS_H
#define  STAT_FUNCTIONS_H

double average(int n/*Number of elements in array*/, double *values);
double standardDeviation(int n/*Number of elements in array*/, double *values);
double percentError(double measured, double actual);
double chiSquareTest(int n_exp, double *experimental, double *theoretical, double *sigma, double c);
double chiSquareLogTest(int n_exp, double *experimental, double *theoretical, double *sigma);
double scalingCoefficient(int n_exp, double *experimental, double *theoretical, double *sigma);
double linterpolate(double y_a, double y_b, double x_a, double x, double x_b);
double chi_of_data_files(char *experimentalFile, char *theoreticalFile,  int scaleTheoretical/*0: NO, 1:YES*/);
double chi_of_data_files_v2(char *experimentalFile, char *theoreticalFile,  int scaleTheoretical/*0: NO, 1:YES*/);
double chi_of_data_files_v3(char *experimentalFile, char *theoreticalFile,  double min_range, double max_range, 
			    int scaleTheoretical/*0: NO, 1:YES*/, int printDetails/*0: NO, 1:YES*/, FILE *FIT_FILE);
#endif
