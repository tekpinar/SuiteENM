#ifndef PDB_IO_H
#define PDB_IO_H
//Version: 0.01- I added xyz2pdb() and pdb2seqANDxyz() functions. February, 12, 2012. 
//Version: 0.02- I added scanpdb_trajectory() functions. November 16, 2014. 
int  *scanpdb(int *atomCount, char *inFileName);
void  readpdb(pdb *info, CAcoord *atom, char *inFileName);
void  readpdb_v23(pdb_v23 *info, CAcoord *atom, char *inFileName);
int  *scanpdb_CAandOH2(int *atomCount, char *inFileName);
int  *scanpdb_CAandOH2_v1(int *atomCount, char *inFileName);
void  readpdb_CAandOH2(pdb_v23 *info, CAcoord *atom, char *inFileName);
void  readHelices(char *inFileName, helixInfo *helix);
void  readSheets(char *inFileName, sheetInfo *sheet);
int   scan4Helices(char *inFileName);
int   scan4Sheets(char *inFileName);
//int write2pdb_v2(int N, char *inFile, CAcoord *atomUpdated,  CAcoord *atomsFixed1, double RMSD, char **argv, double scaler, double *w/*Weight for RMSD*/);
void write2pdb_v23_CA(int N/*System size*/, CAcoord *atoms, char *outFile, int modelNo);
int  *scanpdb_mm(int *atomCount, char *inFileName);
void  readpdb_mm(pdb_v23 *info, CAcoord *atom, char *inFileName, int *atomCount, int modelNo);
void writeCA2pdb(int N/*Number of CA*/, FILE *OUTFILE, CAcoord *atoms);
int xyz2pdb(char *seqFileName, char *xyzFileName, char *outFileName, int writeWater, int writeProtein, int printDetails);
void pdb2seqANDxyz(char *pdbFileName, char *seqFileName,  char *xyzFileName, int printDetails);
void combineALLtoCG(char *allAtomPdb, char *coarsePdb, char *allProtPlusCGwatPdb);
CAcoord *readPdbReturnCA(char *fileName, int *protCounts);
CAcoord *readPdbReturnCAandOH2(char *fileName, int *protCounts);
int  *scanpdb_trajectory(int *atomCount, int *resCount, int *frameCount, char *inFileName, int printDetails);
void  readpdb_trajectory(pdb_v23 *info, CAcoord *atom, FILE *PDBDATA, int *atomCount, int printDetails);
#endif
