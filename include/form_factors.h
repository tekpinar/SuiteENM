//Purpose: To keep all atomic and residue small angle x-ray or wide angle x-ray form factor calculation functions in a file. 
//Author : Mustafa Tekpinar
//Date   : 08/17/2011
//Revisions:
//     v.00: calculateAtomicFF_v2() and calculateAtomicFF_v3() functions added. 
//     v.01: calculateAtomicFF_v4() and formFactorAllAtoms() functions added on 26 October 2011.  
#ifndef  FORM_FACTORS_H
#define  FORM_FACTORS_H
double calculateAtomicFF_v2(char element, double q/*Scattering vector*/, double rho_s/*Solvent density=0.334(water)*/, int isSolventOn/*0:NO - 1:YES*/);
double calculateAtomicFF_v3(char element, double q/*Scattering vector*/);
double calculateAtomicFF_v4(char element, double q/*Scattering vector*/, double rho_s/*Solvent density=0.334(water)*/, int isSolventOn/*0:NO - 1:YES*/);
void formFactorAllAtoms(int N/*Total number of atoms*/, pdb_v23 *info, double waterWeight, double *f_l, double q, double rho_s, int isSolventOn);
void formFactorAllAtoms_v3(int frameNo, int frm_beg, int frm_end, char *element, double waterWeight, double *f_l, double q, double rho_s, int isSolventOn);
void formFactorAllAtoms_v4(int frameNo, int frm_beg, int frm_end, char *element, double waterWeight, double *f_l, double q, double rho_s, int isSolventOn);
#endif


