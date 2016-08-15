//Author of all of these functions except 'blastAll' is Wenjun Zheng. 
//The functions have been copied during our study with his permission.  
#ifndef  BLAST_FUNCTIONS_H
#define  BLAST_FUNCTIONS_H
//#include <stdio.h>
//Personal header files
//#include "/user/tekpinar/mylib/defines.h"
//#include "/user/tekpinar/mylib/structures.h"
//FILE*   FID_debug=NULL;
//int     LIGAND_ON=0;  /* read ligand coords from pdb files */
//char    aaMap[28]="ACDEFGHIKLMNPQRSTVWYXACGTUX";
int  read_line(int lenth_max, char* line, FILE* fid) ;
int  parse_alignment(char* fname_in, char* fname_out);
void showMap(int* seqAA1, int* seqAA2, int len1, int len2,int* map);
int* BLAST2seq_char (char* name, char*    seq1, char*    seq2, int len1, int len2, int* length_mapped);
int* BLAST2seq(char *name , int*    seqAA1, int*    seqAA2, int len1, int len2, int* length_mapped);
int  is_ligand(char* resname);
int  aa2int(char *s);
int  blastAll(char *name, CAcoord *atom1, CAcoord *atom2, int lengthProtein1, int lengthProtein2);
int  blastAll_v2(char *name, CAcoord *atom1, CAcoord *atom2, int lengthProtein1, int lengthProtein2);
int* blastAll_v3(CAcoord *atom1, CAcoord *atom2, int lengthProtein1, int lengthProtein2, \
		char *file1, char *file2, int printDetails, int printMap2File);
#endif
