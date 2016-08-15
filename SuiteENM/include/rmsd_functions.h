#ifndef RMSD_FUNCTIONS_H
#define RMSD_FUNCTIONS_H
extern void u3b_(double *W, double *X, double *Y, int *N,int *MODE, double *RMS, double *U, double *T,int  *IER);
void   superpose_v0(int N/*Number of mapped residues*/, CAcoord *atom1, CAcoord *atom2, int *map1, FILE *FID_debug);
double *superpose_v1(int N/*Number of mapped residues*/, CAcoord *atomUpdated, double *resultVec, FILE *FID_debug);
double superpose(int N/*Number of mapped residues*/, CAcoord *atom1, CAcoord *atom2, int *map1, double *w/*Weight for RMSD*/);
double *superpose_v2(int N/*Number of mapped residues*/, CAcoord *atomUpdated, double *resultVec, double *w);
double calculateRMS_v1(int N/*System size*/, double *CAcoordSet1, CAcoord *atomUpdated, FILE *FID_debug);
double calculateRMS_v2(int N/*System size*/, CAcoord *atomSet1, CAcoord *atomSet2, double *w/*Weight for each atom*/);
double superimpose(int N/*System size*/, CAcoord *atomSet1, CAcoord *atomSet2, double *w/*Weight for each atom*/, int mode);
double superimpose_allatoms(int N/*# of all atoms*/, pdb_v23 *info1, pdb_v23 *info2, double *w/*Weight for each atom*/, int mode);
#endif
