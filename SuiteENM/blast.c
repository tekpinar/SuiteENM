//Author of all of these functions except 'blastAll' is Wenjun Zheng. 
//The functions have been copied during our study with his permission.  
#include <stdio.h>
//Personal header files
#include "/user/tekpinar/mylib/defines.h"
#include "/user/tekpinar/mylib/structures.h"

FILE*   FID_debug=NULL;
int     LIGAND_ON=0;  /* read ligand coords from pdb files */
char    aaMap[28]="ACDEFGHIKLMNPQRSTVWYXACGTUX";
int read_line(int lenth_max, char* line, FILE* fid) 
{
  /*Purpose: To read lines in a text file(Blast output) for parsing*/
  int     j=0;
  char    ch;
  
  while ((!feof(fid))&&('\n'!=(ch=fgetc(fid))))
    if (j<lenth_max)
      {
	line[j]=ch;
	++j;
      }
  if (feof(fid)) j-=1;
  
  line[j]=0;
  if (0) printf("%d:%s\n",j,line);
  
  return (j);
}

int parse_alignment(char* fname_in, char* fname_out)
{
  /*Purpose: Parse blast sequence alignment data*/
  FILE    *fid1=fopen(fname_in,"r");
  FILE    *fid2=fopen(fname_out,"w");
  char    line[201];
  if (fid1==NULL) return (0);
  while (!feof(fid1))
    {
      if (0<read_line(200, line, fid1))
	{
	  if  (strstr(line,"Query:")==line)
	    {
	      fprintf(fid2,"%s\n",line);
	      read_line(200, line, fid1);
	      fprintf(fid2,"%s\n",line);
	      read_line(200, line, fid1);
	      assert(strstr(line,"Sbjct")==line);
	      fprintf(fid2,"%s\n",line);
	    }
	}
    }
  
  fclose(fid1);
  fclose(fid2);
  
  return (1);
}

void showMap(int* seqAA1, int* seqAA2, int len1, int len2,int* map)
{
  /*Purpose: Show the alignment between two sequences with gap inserted*/
  int     len_tot=(len1+len2)*2, ptr1=0,ptr2=0,ptr=0;
  char    Seq1[len_tot], Seq2[len_tot];
  
  while (((ptr1<len1)||(ptr2<len2))&&(ptr<len_tot-2))
    {
      if (ptr%100==0)
	{ /* 100 chars per line */
	  Seq1[ptr]=Seq2[ptr]='\n';
	  ++ptr;
	  continue;
	}
      if (ptr1<len1)  assert((seqAA1[ptr1]>=0)&&(seqAA1[ptr1]<25));
      if (ptr2<len2)  assert((seqAA2[ptr2]>=0)&&(seqAA2[ptr2]<25));
      
      if (ptr1<len1)
	{
	  if (map[ptr1]<0)
	    { /* a gap in seq2 */
	      Seq1[ptr]=aaMap[seqAA1[ptr1]];
	      Seq2[ptr]='.';
	      ++ptr1;
	      ++ptr;
	    }
	  else if (ptr2<map[ptr1])
	    { /* a gap in seq1 */
	      Seq1[ptr]='.';
	      Seq2[ptr]=aaMap[seqAA2[ptr2]];
	      ++ptr2;
	      ++ptr;
	    }
	  else if (ptr2==map[ptr1])
	    { /* a match */
	      Seq1[ptr]=aaMap[seqAA1[ptr1]];
	      Seq2[ptr]=aaMap[seqAA2[ptr2]];
	      ++ptr1;
	      ++ptr2;
	      ++ptr;
	    }
	  else
	    {
	      printf("Error in showMap at %d,%d\n",ptr1,ptr2);
	      assert(0);
	    }
	}
      else if (ptr2<len2)
	{ /* done with seq1 but not seq2 */
	  Seq1[ptr]='.';
	  Seq2[ptr]=aaMap[seqAA2[ptr2]];
	  ++ptr2;
	  ++ptr;
	}
    }
  Seq1[ptr]='\n';
  Seq2[ptr]='\n';
  Seq1[ptr+1]=Seq2[ptr+1]=0;
  
  if (DEBUGMODE)
    {
      fprintf(FID_debug,"Sequence Alignment as shown below:\n");
      
      ptr1=ptr2=1;
      while ((ptr1<ptr)&&(ptr2<ptr))
	{
	  fprintf(FID_debug,"1:");
	  while (Seq1[ptr1]!='\n')
	    {
	      fprintf(FID_debug,"%c",Seq1[ptr1]);
	      ++ptr1;
	    }
	  fprintf(FID_debug,"%c",Seq1[ptr1]);
	  ++ptr1;
	  fprintf(FID_debug,"2:");
	  while (Seq2[ptr2]!='\n')
	    {
	      fprintf(FID_debug,"%c",Seq2[ptr2]);
	      ++ptr2;
	    }
	  fprintf(FID_debug,"%c",Seq2[ptr2]);
	  ++ptr2;
	  fprintf(FID_debug,"\n");
	}
    }
}

int*    BLAST2seq_char (char* name, char*    seq1, char*    seq2, int len1, int len2, int* length_mapped)
{
  /*Purpose: Call bl2seq to align two sequences */

  int*    map = (int*)malloc(sizeof(int)*len1);
  FILE    *fid1, *fid2;
  int     i,j,ptr1[2],ptr2[2],ptr1_old=-1,ptr2_old=-1;
  char    line[201], s1[200], s2[200], sim[200];
  length_mapped[0]=0;
  assert(map!=NULL);
  if (DEBUGMODE)   fprintf(FID_debug,"\nBLAST2seq: Run bl2seq\n");
  for (i=0;i<len1;++i)    map[i] = -1;
  /* write 2 seq files */
  sprintf(s1,"seq1_%s",name);
  sprintf(s2,"seq2_%s",name);
  assert((fid1=fopen(s1,"w"))!=NULL);
  assert((fid2=fopen(s2,"w"))!=NULL);
  for (i=0;i<len1;++i)    fprintf(fid1,"%c",seq1[i]);
  for (i=0;i<len2;++i)    fprintf(fid2,"%c",seq2[i]);
  fclose(fid1);
  fclose(fid2);
  
  /* call bl2seq */
  sprintf(line,"/user/tekpinar/cholmodsolver/blast-2.2.21/bin/bl2seq -i seq1_%s -j seq2_%s -p blastp -o out_blast_%s -F F -e 0.001 > junk",name,name,name);  /* default E=0.001 */
  system(line);
  sprintf(s1,"out_blast_%s",name);
  sprintf(s2,"out_align_%s",name);
  if (0==parse_alignment(s1, s2))
    {
      free(map);
      return (NULL);
    }

  fid1=fopen(s2,"r");
  assert(fid1!=NULL);
  
  /* parse output from bl2seq */
  while (!feof(fid1)) {
    j=read_line(200, line, fid1);   /* read Query line */
    if (feof(fid1)) break;
    assert(strstr(line,"Query")==line);
    sscanf(line+7,"%d %s %d",&(ptr1[0]),s1,&(ptr1[1]));
    if (ptr1_old>=ptr1[0])  break;
    else    ptr1_old=ptr1[1];
    
    j=read_line(200, line, fid1);   /* read Similarity line */
    if (feof(fid1)) break;
    assert(memcmp(line,"           ",6)==0);
    strcpy(sim,line+11);
    
    j=read_line(200, line, fid1);   /* read Sbjct line */
    assert(strstr(line,"Sbjct")==line);
    sscanf(line+7,"%d %s %d",&(ptr2[0]),s2,&(ptr2[1]));
    if (ptr2_old>=ptr2[0])  break;
    else    ptr2_old=ptr2[1];
    
    if (DEBUGMODE)   fprintf(FID_debug,"Q: %s %d-%d\n",s1,ptr1[0],ptr1[1]);
    if (DEBUGMODE)   fprintf(FID_debug,"   %s\n",sim);
    if (DEBUGMODE)   fprintf(FID_debug,"S: %s %d-%d\n\n",s2,ptr2[0],ptr2[1]);
    
    assert(strlen(s1)==strlen(s2));
    /* assert(strlen(s1)==strlen(sim)); */
    
    for (j=0;j<strlen(s1);++j) {
      if ((s1[j]!='-')&&(s2[j]!='-')) {
	/* if ((sim[j]=='+')||(s1[j]==s2[j]))   map[ptr1[0]-1]=ptr2[0]-1; */
	map[ptr1[0]-1]=ptr2[0]-1;
	ptr1[0]+=1;
	ptr2[0]+=1;
	if (s1[j]==s2[j])       length_mapped[0]+=1;
	
      } else if ((s1[j]!='-')&&(s2[j]=='-')) {
	ptr1[0]+=1;
      } else if ((s1[j]=='-')&&(s2[j]!='-')) {
	ptr2[0]+=1;
      } else assert(0);
    }
    assert((ptr1[0]==ptr1[1]+1)&&(ptr2[0]==ptr2[1]+1));
  }
  fclose(fid1);
  
  if (1) printf("Sequence alignment: Map sequence 1(%d) to sequence 2(%d) with %d residues matched (identity=%3.2f)\n", len1, len2, length_mapped[0], \
		1.0*length_mapped[0]/len1);
  /* assert(length_mapped[0]>0); */
  
  /* clean up */
  sprintf(line,"rm -f seq1_%s seq2_%s out_blast_%s out_align_%s  junk", name, name, name, name);
  if (1) system(line);
  return (map);
}

int *BLAST2seq(char *name , int*    seqAA1, int*    seqAA2, int len1, int len2, int* length_mapped)
{
  int     i, *map;
  char    seq1[len1], seq2[len2];
  
  for (i=0;i<len1;++i)
    seq1[i]=aaMap[seqAA1[i]];
  
  for (i=0;i<len2;++i)
    seq2[i]=aaMap[seqAA2[i]];
  
  map=BLAST2seq_char (name, seq1, seq2, len1, len2, length_mapped);
  
  if (DEBUGMODE)
    showMap(seqAA1, seqAA2, len1, len2, map);
  return (map);
}

int is_ligand(char* resname)
{
  if (LIGAND_ON==0) return (0);
  
  if ((memcmp(resname,"MG ",3)==0)||(memcmp(resname," MG",3)==0)) return (1);    /* ligand Mg */
  if (memcmp(resname,"ADP", 3)==0) return (1);    /* ligand ADP */
  if (memcmp(resname,"ANP", 3)==0) return (1);    /* ligand ANP */
  if (memcmp(resname,"ATP", 3)==0) return (1);    /* ligand ATP */
  if (memcmp(resname,"ACP", 3)==0) return (1);    /* ligand ACP */
  if (memcmp(resname,"VO4", 3)==0) return (1);    /* ligand VO4 */
  if (memcmp(resname,"PIH", 3)==0) return (1);    /* ligand Pi  */
  if (memcmp(resname,"TIP3",3)==0) return (1);    /* ligand water */
  if (memcmp(resname,"HOH", 3)==0) return (1);    /* ligand water */
  return (0);
}

int aa2int(char *s)
{                         /* aa #heavy-atom #all-atom */
  if (strncmp(s,"ALA", 3)==0) return (0);     /* ala A: 6  13 */
  if (strncmp(s,"CYS", 3)==0) return (1);     /* cys C: 7  14 */
  if (strncmp(s,"ASP", 3)==0) return (2);     /* asp D: 9  15 */
  if ((strncmp(s,"GLU", 3)==0)||(strcmp(s,"PCA")==0))  return (3);    /* glu E: 10 18 */
  if (strncmp(s,"PHE", 3)==0) return (4);     /* phe F: 12 23 */
  if (strncmp(s,"GLY", 3)==0) return (5);     /* gly G: 5  10 */
  if (strncmp(s,"HIS", 3)==0) return (6);     /* HIS H: 11 20 */
  if (strncmp(s,"HSD", 3)==0) return (6);
  if (strncmp(s,"ILE", 3)==0) return (7);     /* ile I: 9  22 */
  if (strncmp(s,"LYS", 3)==0) return (8);     /* lys K: 10 25 */
  if (strncmp(s,"LEU", 3)==0) return (9);     /* leu L: 9  22 */
  if ((strncmp(s,"MET", 3)==0)||(strcmp(s,"MSE")==0)||(strcmp(s,"FME")==0)) return (10);      /* MET M: 9  20 */
  if (strncmp(s,"ASN", 3)==0) return (11);    /* asn N: 9  17 */
  if (strncmp(s,"PRO", 3)==0) return (12);    /* pro P: 8  17 */
  if (strncmp(s,"GLN", 3)==0) return (13);    /* gln Q: 10 20 */
  if (strncmp(s,"ARG", 3)==0) return (14);    /* arg R: 12 27 */
  if ((strncmp(s,"SER", 3)==0)||(strcmp(s,"OSE")==0)) return (15);    /* ser S: 7  14 */
  if (strncmp(s,"THR", 3)==0) return (16);    /* thr T: 8  17 */
  if (strncmp(s,"VAL", 3)==0) return (17);    /* val V: 8  19 */
  if (strncmp(s,"TRP", 3)==0) return (18);    /* trp W: 15 27 */
  if (strncmp(s,"TYR", 3)==0) return (19);    /* tyr Y: 13 24 */
  if (strncmp(s,"  A", 3)==0) return (21);    /* A */
  if (strncmp(s," DA", 3)==0) return (21);    /* d-na, see 2ktq */
  if (strncmp(s,"  C", 3)==0) return (22);    /* C */
  if (strncmp(s," DC", 3)==0) return (22);
  if (strncmp(s,"  G", 3)==0) return (23);    /* G */
  if (strncmp(s," DG", 3)==0) return (23);
  if (strncmp(s,"  T", 3)==0) return (24);    /* T */
  if (strncmp(s," DT", 3)==0) return (24);
  if (strncmp(s,"  U", 3)==0) return (25);    /* U */
  if (strncmp(s," DU", 3)==0) return (25);
  if (is_ligand(s))       return (26);    /* ligand */

  if (0&&DEBUGMODE) fprintf(FID_debug,"Warning:%s is unknown AA\n",s);
  return(20);                             /* unknown and invalid aa name */
}

int blastAll(char *name, CAcoord *atom1, CAcoord *atom2, int lengthProtein1, int lengthProtein2)
{
  /*Purpose: Align two sequences and write aligned coordinates to CAcoordConf1.txt*/
  /*         and CAcoordConf2.txt files.                                          */
  int i=0, j=0;
  
  int *seqAA1=(int *)malloc(lengthProtein1*sizeof(int));
  if(seqAA1==NULL)
    {
      fprintf(stderr, "Memory not allocated for seqAA1\n");
      exit(EXIT_FAILURE);
    }

  int *seqAA2=(int *)malloc(lengthProtein2*sizeof(int));
  if(seqAA2==NULL)
    {
      fprintf(stderr, "Memory not allocated for seqAA2\n");
      exit(EXIT_FAILURE);
    }
  
  for(i=0;i<lengthProtein1;i++)
    seqAA1[i]=aa2int(atom1[i].residname);

  for(i=0;i<lengthProtein2;i++)
    seqAA2[i]=aa2int(atom2[i].residname);
  
  int *map;
  int length_mapped=0;
  // BLAST2seq();
  // FID_debug=fopen("test.txt", "w");//If(DEBUG==)
  map=BLAST2seq (name, seqAA1, seqAA2, lengthProtein1, lengthProtein2, &length_mapped);
  
  for(i=0;i<lengthProtein1;i++)
    {
    map[i]=i;
    printf("map[%d]=%d\n",i, map[i]);
    }
  //printf("Length mapped=%d\n", length_mapped);
  FILE *CAcoord1=fopen("CAcoordConf1.txt", "w");
  assert(CAcoord1!=NULL);
  FILE *CAcoord2=fopen("CAcoordConf2.txt", "w");
  assert(CAcoord2!=NULL);
  FILE *residueInfo1=fopen("map1.pdb", "w");
  assert(residueInfo1!=NULL);
  FILE *residueInfo2=fopen("map2.pdb", "w");
  assert(residueInfo2!=NULL);
  
  for(i=0;i<lengthProtein1;i++)
    {
      if(map[i]!=(-1))
	{
	  fprintf(CAcoord1, "%.3lf\t%.3lf\t%.3lf\n", atom1[i].x,  atom1[i].y, atom1[i].z);
	  fprintf(CAcoord2, "%.3lf\t%.3lf\t%.3lf\n", atom2[map[i]].x,  atom2[map[i]].y, atom2[map[i]].z);
	  
	  //"ATOM  %5d %4s %3.3s %c %4d   %8.3lf %8.3lf %8.3lf %6.2lf %6.2lf\n" old!!!

	  fprintf(residueInfo1,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n" , (i+1), " CA ", \
		  atom1[i].residname,  atom1[i].chain, atom1[i].resNo, atom1[i].x,  atom1[i].y, atom1[i].z, 1.00, 0.00);
	  fprintf(residueInfo2, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (map[i]+1), " CA ", \
		  atom2[map[i]].residname,  atom2[map[i]].chain, atom2[map[i]].resNo, atom2[map[i]].x,  atom2[map[i]].y, atom2[map[i]].z, 1.00, 0.00);
	  j+=1;
	}
    }
  
  printf("Number of residues mapped=%d\n", j);
  fclose(CAcoord1);
  fclose(CAcoord2);
  fclose(residueInfo1);
  fclose(residueInfo2);
  return (j);
  
}
int blastAll_v2(CAcoord *atom1, CAcoord *atom2, int lengthProtein1, int lengthProtein2, 
		char *file1, char *file2, int printDetails)
{
  /*Purpose: Align two sequences and write aligned coordinates to two pdb files whose name are file1.pdb and file2.pdb.*/
  
  int i=0, j=0;
  int *map;
  int length_mapped=0;

  int *seqAA1=(int *)malloc(lengthProtein1*sizeof(int));
  if(seqAA1==NULL)
    {
      fprintf(stderr, "Memory not allocated for seqAA1\n");
      exit(EXIT_FAILURE);
    }
  
  int *seqAA2=(int *)malloc(lengthProtein2*sizeof(int));
  if(seqAA2==NULL)
    {
      fprintf(stderr, "Memory not allocated for seqAA2\n");
      exit(EXIT_FAILURE);
    }
  
  for(i=0;i<lengthProtein1;i++)
    seqAA1[i]=aa2int(atom1[i].residname);

  for(i=0;i<lengthProtein2;i++)
    seqAA2[i]=aa2int(atom2[i].residname);
  
  char *name="test";
  map=BLAST2seq (name, seqAA1, seqAA2, lengthProtein1, lengthProtein2, &length_mapped);
  
  if(printDetails)
    for(i=0;i<lengthProtein1;i++)
      {
	map[i]=i;
	printf("map[%d]=%d\n",i, map[i]);
      }
  FILE *residueInfo1=fopen(file1, "w");
  assert(residueInfo1!=NULL);
  FILE *residueInfo2=fopen(file2, "w");
  assert(residueInfo2!=NULL);
  
  for(i=0;i<lengthProtein1;i++)
    {
      if(map[i]!=(-1))
	{
	  fprintf(CAcoord1, "%.3lf\t%.3lf\t%.3lf\n", atom1[i].x,  atom1[i].y, atom1[i].z);
	  fprintf(CAcoord2, "%.3lf\t%.3lf\t%.3lf\n", atom2[map[i]].x,  atom2[map[i]].y, atom2[map[i]].z);
	  
	  fprintf(residueInfo1,"ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n" , (i+1), " CA ", \
		  atom1[i].residname,  atom1[i].chain, atom1[i].resNo, atom1[i].x,  atom1[i].y, atom1[i].z, 1.00, 0.00);
	  fprintf(residueInfo2, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (map[i]+1), " CA ", \
		  atom2[map[i]].residname,  atom2[map[i]].chain, atom2[map[i]].resNo, atom2[map[i]].x,  atom2[map[i]].y, \
		  atom2[map[i]].z, 1.00, 0.00);
	  j+=1;
	}
    }
  
  printf("Number of residues mapped=%d\n", j);
  fclose(residueInfo1);
  fclose(residueInfo2);
  return (j);
  
}
