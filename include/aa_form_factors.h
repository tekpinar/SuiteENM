#ifndef AA_FORM_FACTORS_H
#define AA_FORM_FACTORS_H
void   readYangFF(double F_CG[21][100], char *residDataFile);
void readYangFF_v2(int num_data_points, double F_CG[21][num_data_points], char *residDataFile);
double resid2YangCoefficient(char *resName, double F_CG[21][100], double q/*q value of form factor*/, double waterWeight);
double resid2YangFF_Long(int num_data_points, char *resName, double F_CG[21][num_data_points],\
			 double q/*q value of form factor*/, double waterWeight);
void   readStovgaardFF(double F_CG[21][100], char *residDataFile);
double resid2StovgaardCoefficient(char *resName, double F_CG[21][100], double q/*q value of form factor*/, double waterWeight);
void determineResLengths(int *atomCount, pdb_v23 *info, int *residLengths);
double calculateResidueFF_v3(pdb_v23 *info, int *atomCount, double s/*Scattering vector*/, double rho_s/*Solvent density=0.334(for water)*/,\
			     int isSolventOn/*0:NO - 1:YES*/,  double *F_CG_per_q, int *badRes);
#endif
