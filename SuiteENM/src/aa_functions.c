//Purpose: I noticed that I need some functions to convert 
//         3-letter amino acid codes 1-letter codes and vice versa. 
//Author: Mustafa Tekpinar
//Email : tekpinar@buffalo.edu
//Date: 31/01/2011: 
//Revision history:
//v0.01: Initial Version
//v0.02: I decided to collect all functions that contain basic information about 
//       amino acids in this module. Therefore, I added aa2Mass() and aa2Area()
//       functions here. 03/18/2011
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
  From: http://imgt.cines.fr/textes/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html

  Amino acid     Abbreviations Molecular mass (Da) Number of atoms Volume (A3) Hydropathy index
  Alanine        Ala A         89                  13              88.6        1.8
  Arginine       Arg R         174                 26              173.4       -4.5
  Asparagine     Asn N         132                 17              114.1       -3.5
  Aspartic acid  Asp D         133                 16              111.1       -3.5
  Asparagine or 
  Aspartic acid Asx B 
  Cysteine      Cys C          121                 14              108.5        2.5
  Glutamine     Gln Q          146                 20              143.8       -3.5
  Glutamic Acid Glu E          147                 19              138.4       -3.5
  Glutamine or 
  Glutamic acid Glx Z 
  Glycine       Gly G          75                  10              60.1        -0.4
  Histidine     His H         155                  20              153.2       -3.2
  Isoleucine    Ile I         131                  22              166.7        4.5
  Leucine       Leu L         131                  22              166.7        3.8
  Lysine        Lys K         146                  24              168.6       -3.9
  Methionine    Met M         149                  20              162.9        1.9
  Phenylalanine Phe F         165                  23              189.9        2.8
  Proline       Pro P         115                  17              112.7       -1.6
  Serine        Ser S         105                  14              89.0        -0.8
  Threonine     Thr T         119                  17             116.1        -0.7
  Tryptophan    Trp W         204                  27             227.8        -0.9
  Tyrosine      Tyr Y         181                  24             193.6        -1.3
  Valine        Val V         117                  19             140.0         4.2
*/

char aa_3letter_to_1letter(char *resName)
{
  //Purpose: Obviously, it converts 3-letter amino acid codes to 1-letter codes.
  //         It works only for upper case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if((strncmp(resName, "SER", 3)==0)) return 'S'; // S=0  //Serine        (SER)
  if((strncmp(resName, "PHE", 3)==0)) return 'F'; // F=1  //Phenylalanine (PHE)
  if((strncmp(resName, "THR", 3)==0)) return 'T'; // T=2  //Threonine     (THR)
  if((strncmp(resName, "ASN", 3)==0)) return 'N'; // N=3  //Asparagine    (ASN)
  if((strncmp(resName, "LYS", 3)==0)) return 'K'; // K=4  //Lysine        (LYS)
  if((strncmp(resName, "GLU", 3)==0)) return 'E'; // E=5  //Glutamic Acid (Glu)
  if((strncmp(resName, "TYR", 3)==0)) return 'Y'; // Y=6  //Tyrosine      (TYR)
  if((strncmp(resName, "VAL", 3)==0)) return 'V'; // V=7  //Valine        (VAL)
  if((strncmp(resName, "GLN", 3)==0)) return 'Q'; // Q=8  //Glutamine     (GLN)
  if((strncmp(resName, "MET", 3)==0)) return 'M'; // M=9  //Methionine    (MET)
  if((strncmp(resName, "CYS", 3)==0)) return 'C'; // C=10 //Cysteine      (CYS)
  if((strncmp(resName, "LEU", 3)==0)) return 'L'; // L=11 //Leucine       (LEU)
  if((strncmp(resName, "ALA", 3)==0)) return 'A'; // A=12 //Alanine       (ALA)
  if((strncmp(resName, "TRP", 3)==0)) return 'W'; // W=13 //Tryptophan    (TRP)
  if((strncmp(resName, "PRO", 3)==0)) return 'P'; // P=14 //Proline       (PRO)
  if((strncmp(resName, "HIS", 3)==0)) return 'H'; // H=15 //Histidine     (HIS)
  if((strncmp(resName, "HSD", 3)==0)) return 'H'; // H=15 //Histidine     (HSD)
  if((strncmp(resName, "HSE", 3)==0)) return 'H'; // H=15 //Histidine     (HSE)
  if((strncmp(resName, "HSP", 3)==0)) return 'H'; // H=15 //Histidine     (HSP)
  if((strncmp(resName, "ASP", 3)==0)) return 'D'; // D=16 //Aspartic Acid (ASP)
  if((strncmp(resName, "ARG", 3)==0)) return 'R'; // R=17 //Arginine      (ARG)
  if((strncmp(resName, "ILE", 3)==0)) return 'I'; // I=18 //Isoleucine    (ILE)
  if((strncmp(resName, "GLY", 3)==0)) return 'G'; // G=19 //Glysine       (GLY)
  else
    {
      fprintf(stderr, "ERROR: Unknown residue:%s!\n", resName);
      exit(EXIT_FAILURE);
      return 'X';
    }
}

char *aa_1letter_to_3letter(char oneLetterResName)
{
  //Purpose: This function converts 1-letter amino acid names to 1-letter codes.
  //         It works only for upper case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //         Fuction returns an upper case string WITHOUT a null character. 

  if(oneLetterResName=='S') return "SER"; // S=0  //Serine        (SER)
  if(oneLetterResName=='F') return "PHE"; // F=1  //Phenylalanine (PHE)
  if(oneLetterResName=='T') return "THR"; // T=2  //Threonine     (THR)
  if(oneLetterResName=='N') return "ASN"; // N=3  //Asparagine    (ASN)
  if(oneLetterResName=='K') return "LYS"; // K=4  //Lysine        (LYS)
  if(oneLetterResName=='E') return "GLU"; // E=5  //Glutamic Acid (Glu)
  if(oneLetterResName=='Y') return "TYR"; // Y=6  //Tyrosine      (TYR)
  if(oneLetterResName=='V') return "VAL"; // V=7  //Valine        (VAL)
  if(oneLetterResName=='Q') return "GLN"; // Q=8  //Glutamine     (GLN)
  if(oneLetterResName=='M') return "MET"; // M=9  //Methionine    (MET)
  if(oneLetterResName=='C') return "CYS"; // C=10 //Cysteine      (CYS)
  if(oneLetterResName=='L') return "LEU"; // L=11 //Leucine       (LEU)
  if(oneLetterResName=='A') return "ALA"; // A=12 //Alanine       (ALA)
  if(oneLetterResName=='W') return "TRP"; // W=13 //Tryptophan    (TRP)
  if(oneLetterResName=='P') return "PRO"; // P=14 //Proline       (PRO)
  if(oneLetterResName=='H') return "HIS"; // H=15 //Histidine     (HIS)
  if(oneLetterResName=='D') return "ASP"; // D=16 //Aspartic Acid (ASP)
  if(oneLetterResName=='R') return "ARG"; // R=17 //Arginine      (ARG)
  if(oneLetterResName=='I') return "ILE"; // I=18 //Isoleucine    (ILE)
  if(oneLetterResName=='G') return "GLY"; // G=19 //Glysine       (GLY)
  else
    {
      fprintf(stderr, "ERROR: Unknown residue: %c!\n", oneLetterResName);
      return "XXX";
      exit(EXIT_FAILURE);
    }
}

double aa2Mass(char *residueName)
{

  //Put the reference name here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Purpose: This function returns the mass of each residue. 
  //         'char *residueName' is upper case three letter code
  //         and return value is mass in atomic mass units. 
  if (strncmp(residueName,"ALA", 3)==0) return  (71.0779);     /* A */
  if (strncmp(residueName,"ARG", 3)==0) return (156.1857);     /* R */
  if (strncmp(residueName,"ASN", 3)==0) return (114.1026);     /* N */
  if (strncmp(residueName,"ASP", 3)==0) return (115.0874);     /* D */
  if (strncmp(residueName,"CYS", 3)==0) return (103.1429);     /* C */
  if (strncmp(residueName,"GLN", 3)==0) return (128.1292);     /* Q */
  if (strncmp(residueName,"GLU", 3)==0) return (129.1140);     /* E */
  if (strncmp(residueName,"GLY", 3)==0) return  (57.0513);     /* G */
  if (strncmp(residueName,"HIS", 3)==0) return (137.1393);     /* H */
  if (strncmp(residueName,"ILE", 3)==0) return (113.1576);     /* I */
  if (strncmp(residueName,"LEU", 3)==0) return (113.1576);     /* L */
  if (strncmp(residueName,"LYS", 3)==0) return (128.1723);     /* K */
  if (strncmp(residueName,"MET", 3)==0) return (131.1961);     /* M */
  if (strncmp(residueName,"PHE", 3)==0) return (147.1739);     /* F */
  if (strncmp(residueName,"PRO", 3)==0) return  (97.1152);     /* P */
  if (strncmp(residueName,"SER", 3)==0) return  (87.0773);     /* S */
  if (strncmp(residueName,"THR", 3)==0) return (101.1039);     /* T */
  if (strncmp(residueName,"TRP", 3)==0) return (186.2099);     /* W */
  if (strncmp(residueName,"TYR", 3)==0) return (163.1733);     /* Y */
  if (strncmp(residueName,"VAL", 3)==0) return  (99.1311);     /* V */
  else
    {
      fprintf(stderr,"ERROR: %s is unknown residue encountered!\n", residueName);
      return(0.0);
      exit(EXIT_FAILURE);
    } 
}

double aa2Area(char *residueName)
{
  //Surface area data in (Angstrom^2) taken from 'Proteins, Structure and Functions' by
  //David Whitford. 
  if (strncmp(residueName,"ALA", 3)==0) return (115.0);     /* A */
  if (strncmp(residueName,"ARG", 3)==0) return (225.0);     /* R */
  if (strncmp(residueName,"ASN", 3)==0) return (160.0);     /* N */
  if (strncmp(residueName,"ASP", 3)==0) return (150.0);     /* D */
  if (strncmp(residueName,"CYS", 3)==0) return (135.0);     /* C */
  if (strncmp(residueName,"GLN", 3)==0) return (180.0);     /* Q */
  if (strncmp(residueName,"GLU", 3)==0) return (190.0);     /* E */
  if (strncmp(residueName,"GLY", 3)==0) return  (75.0);     /* G */
  if (strncmp(residueName,"HIS", 3)==0) return (195.0);     /* H */
  if (strncmp(residueName,"ILE", 3)==0) return (175.0);     /* I */
  if (strncmp(residueName,"LEU", 3)==0) return (170.0);     /* L */
  if (strncmp(residueName,"LYS", 3)==0) return (200.0);     /* K */
  if (strncmp(residueName,"MET", 3)==0) return (185.0);     /* M */
  if (strncmp(residueName,"PHE", 3)==0) return (210.0);     /* F */
  if (strncmp(residueName,"PRO", 3)==0) return (145.0);     /* P */
  if (strncmp(residueName,"SER", 3)==0) return (115.0);     /* S */
  if (strncmp(residueName,"THR", 3)==0) return (140.0);     /* T */
  if (strncmp(residueName,"TRP", 3)==0) return (255.0);     /* W */
  if (strncmp(residueName,"TYR", 3)==0) return (230.0);     /* Y */
  if (strncmp(residueName,"VAL", 3)==0) return (155.0);     /* V */

  else
    {
      fprintf(stderr,"ERROR: %s is unknown residue\n",residueName);
      return(0.0);                          
      exit(EXIT_FAILURE);
    }
}

int aa_atomNumbers(char *residueName)
{
  //This function returns number of atoms for each residue:
  if (strncmp(residueName,"ALA", 3)==0) return (13);     /* A */
  if (strncmp(residueName,"ARG", 3)==0) return (27);     /* R One must put a caution remark here!-I added one more atom to chemical formula*/
  if (strncmp(residueName,"ASN", 3)==0) return (17);     /* N */
  if (strncmp(residueName,"ASP", 3)==0) return (15);     /* D: In all of the pdb files total number is 15 except terminal residues. One is subtracted!*/
  if (strncmp(residueName,"CYS", 3)==0) return (14);     /* C */
  if (strncmp(residueName,"GLN", 3)==0) return (20);     /* Q */
  if (strncmp(residueName,"GLU", 3)==0) return (18);     /* E In all of the pdb files total number is 15 except terminal residues. One is subtracted!*/
  if (strncmp(residueName,"GLY", 3)==0) return (10);     /* G */
  if (strncmp(residueName,"HIS", 3)==0) return (20);     /* H */
  if (strncmp(residueName,"HSD", 3)==0) return (20);     /* H */
  if (strncmp(residueName,"HSE", 3)==0) return (20);     /* H */
  if (strncmp(residueName,"HSP", 3)==0) return (20);     /* H */
  if (strncmp(residueName,"ILE", 3)==0) return (22);     /* I */
  if (strncmp(residueName,"LEU", 3)==0) return (22);     /* L */
  if (strncmp(residueName,"LYS", 3)==0) return (25);     /* K: One must put a caution remark here!-I added one more atom to chemical formula*/
  if (strncmp(residueName,"MET", 3)==0) return (20);     /* M */
  if (strncmp(residueName,"PHE", 3)==0) return (23);     /* F */
  if (strncmp(residueName,"PRO", 3)==0) return (17);     /* P */
  if (strncmp(residueName,"SER", 3)==0) return (14);     /* S */
  if (strncmp(residueName,"THR", 3)==0) return (17);     /* T */
  if (strncmp(residueName,"TRP", 3)==0) return (27);     /* W */
  if (strncmp(residueName,"TYR", 3)==0) return (24);     /* Y */
  if (strncmp(residueName,"VAL", 3)==0) return (19);     /* V */
  else
    {
      fprintf(stderr,"ERROR: %s is unknown residue\n",residueName);
      return(0);                          
      exit(EXIT_FAILURE);
    }
}
//Main program has been written just for test purposes!
/* int main() */
/* { */
/* /\*   char resName[4]={'\0', '\0', '\0', '\0'}; *\/ */
/* /\*   printf("Enter three letter code of an amino acid:\n"); *\/ */
/* /\*   scanf("%s", resName); *\/ */
/* /\*   printf("Amino acid=%s\n", resName); *\/ */
/* /\*   printf("One letter code for %s is %c\n", resName, aa_3letter_to_1letter(resName)); *\/ */

/*   char oneLetterResName='\0'; */
/*   printf("Enter one letter code of an amino acid:\n"); */
/*   scanf("%c", &oneLetterResName); */
/*   printf("Amino acid=%c\n", oneLetterResName); */
/*   printf("Mass for  %c is %lf \n", oneLetterResName, aa2Mass(aa_1letter_to_3letter(oneLetterResName))); */

/*   return 0; */
/* } */
