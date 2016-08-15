#ifndef DEFINES_H
#define DEFINES_H
//====================GLOBAL VARIABLES=================================================
#define PI 3.1415926535
//#define R_CUTOFF 10          /* Don't forget to change R_CUTOFF_SQUARED if you change this value!!       */
//#define R_CUTOFF_SQUARED 100 /* Don't forget to change R_CUTOFF if you change this value!!               */
#define R_COLLISION 4.0        /* Collision radius to prevent steric collisions in ENM                     */
#define CL_CON 10.0            /* If close contact(neighbor residues) force constant is 10.                */
#define FR_CON 1.0             /* If far contact(non-neighbor residues) force constant is 1                */
//#define VAN_DER_WAALS 0.0      /* Force constant for van der Walls interactions                            */
#define C_COL 10.0             /* If far contact(non-neighbor residues) approach each other too much!!!!   */
#define totalSteps 500         /*This number limits total number of steps in case of slow convergence      */
#define DEBUGMODE 0
#define INVALID_CRD 999.999
#define DENSE_FORMAT 1         /*Use dense arrays instead of sparse. It is faster but requires more memory */
#define CA_DIST 3.9
#define SIGMA 2.0              /*The full width at half maximum (FWHM) defined based on Gorba, C., Tama, F.*/
#define SIGMA_SQRD 4.0         /*The full width at half maximum (FWHM) defined based on Gorba, C., Tama, F.*/
#define SIGMA_CBD  8.0         /*The full width at half maximum (FWHM) defined based on Gorba, C., Tama, F.*/
float erfTable[500000];

//double R_CUTOFF=10.0;                      /* Don't forget to change R_CUTOFF_SQUARED if you change this value!!       */
//double R_CUTOFF_SQUARED=100.0;             /* Don't forget to change R_CUTOFF if you change this value!!               */

//double N_pairSum=0.0;
#endif
