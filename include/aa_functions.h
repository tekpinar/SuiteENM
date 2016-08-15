//Purpose: I noticed that I need some functions to convert 
//         3-letter amino codes 1-letter codes and vice versa. 
//Author: Mustafa Tekpinar
//Date: 31/01/2011: 
//Revision history:
//v0.01: Initial Version
//v0.02: I decided to collect all functions that contain basic information about 
//       amino acids in this module. Therefore, I added aa2Mass() and aa2Area()
//       functions here. 03/18/2011
#ifndef  AA_FUNCTIONS_H
#define  AA_FUNCTIONS_H

char   aa_3letter_to_1letter(char *resName);
char  *aa_1letter_to_3letter(char oneLetterResName);
double aa2Mass(char *residueName);
double aa2Area(char *residueName);
int    aa_atomNumbers(char *residueName);
#endif
