//Purpose: pdb_io.c contains all routines to read and write pdb files. 
/* Copyright (C) 2012  Mustafa Tekpinar

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; 
version 2.1 of the License.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
02110-1301  USA

Additionally, you can obtain a copy of the licence from 
http://www.gnu.org/licenses/lgpl-2.1.html

*/

//Email: tekpinar@buffalo.edu

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//Personal header file
#include <structures.h>
#include <rmsd_functions.h>
#include <aa_functions.h>
//Version: 0.01- I added xyz2pdb() and pdb2seqANDxyz() functions. February, 12, 2012. 

int  *scanpdb(int *atomCount, char *inFileName)
{
/*Purpose: To determine the total number of atoms and residues for memory allocation*/
/*In future I may use it to keep protein chain and segment infos*/

  int i=0;
  int numberofResids=0;
  int numberofChains=1; /*Remember that if there is just one chain, you will not see any */
                        /*TER signal to count chain numbers.                             */
  int numberofHelices=0;
  int numberofSheets=0;
  char line[100];
  
  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }

  while(1)
    {
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  if(((strncmp(" CA", line+12, 3))==0))

	    numberofResids+=1;
	  
	  i+=1; 
	}
      if((strncmp("TER   ", line, 6))==0)
	{
	  numberofChains+=1;
	}

      if((strncmp("HELIX ", line, 6))==0)
	{
	  numberofHelices+=1;
	}

      if((strncmp("SHEET ", line, 6))==0)
	{
	  numberofSheets+=1;
	}
    }
  atomCount[0]=i;               //atomCount[0] is the total number of atoms
  atomCount[1]=numberofResids;  //atomCount[1] is the total number of residues
  atomCount[2]=numberofChains;  //atomCount[2] is the total number of chains
  atomCount[3]=numberofHelices; //atomCount[3] is the total number of helices
  atomCount[4]=numberofSheets;  //atomCount[4] is the total number of sheets

  if((atomCount[0]==0) || (atomCount[1]==0))
    {
      fprintf(stderr, "Error: Can not read atoms and residues\n");
      exit(EXIT_FAILURE);
    }
  else
    {
      fprintf(stdout, "Total number of atoms:%d\n", atomCount[0]);
      fprintf(stdout, "Number of residues:%d\n",    atomCount[1]);
      fprintf(stdout, "Number of chains:%d\n",      atomCount[2]);
      fprintf(stdout, "Number of helices:%d\n",     atomCount[3]);
      fprintf(stdout, "Number of sheets:%d\n",      atomCount[4]);
    }
  fclose(pdbdata);
  
  return atomCount; 
}

void  readpdb(pdb *info, CAcoord *atom, char *inFileName)
{
  /*Purpose: To read a pdb file and keep the data into the 'pdb' and 'CAcoord' structures*/
  int i=0, k=0;
  int numberofResid=0;

  char line[100];
  char *line_ptr;
  char CAtest[5]={'\0', '\0', '\0', '\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  while(1)
    {
      line_ptr=fgets(line, sizeof(line),pdbdata);
      if(line_ptr==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  
	  sscanf(line,"%s %d %s %s %c %d %lf %lf %lf %lf %lf           %c",	   \
		 info[i].text, &info[i].atnum, info[i].attyp, info[i].resname,     \
		 &info[i].chain, &info[i].resno, &info[i].coordX, &info[i].coordY, \
		 &info[i].coordZ, &info[i].occ, &info[i].temp, &info[i].element );
	  

	  if(0)	
	    fprintf(stdout,"%s %d %s %s %c %d %lf %lf %lf %lf %lf           %c\n", \
		    info[i].text, info[i].atnum, info[i].attyp, info[i].resname,   \
		    info[i].chain, info[i].resno, info[i].coordX, info[i].coordY,  \
		    info[i].coordZ, info[i].occ, info[i].temp, info[i].element );
	  
	  for(k=12;k<16;k++)
	    CAtest[k-12]=line[k];

	  //	  if((strncmp(" CA", line+12, 3))==0)
	  if((strncmp(" CA", CAtest, 3))==0)
	    {
	      strncpy(atom[numberofResid].residname, info[i].resname, 3);
	      atom[numberofResid].chain=info[i].chain;
	      atom[numberofResid].resNo=info[i].resno;
	      atom[numberofResid].x=info[i].coordX;
	      atom[numberofResid].y=info[i].coordY;
	      atom[numberofResid].z=info[i].coordZ;
	      
	      numberofResid+=1;
	    }
	  
	  i+=1; 
	}
    }
  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }
  fclose(pdbdata);
}

/*
PDB FORMAT Version 2.3 ATOM Record

COLUMNS      DATA TYPE        FIELD      DEFINITION
------------------------------------------------------
 1 -  6      Record name      "ATOM    "
 7 - 11      Integer          serial     Atom serial number.
13 - 16      Atom             name       Atom name.
17           Character        altLoc     Alternate location indicator.
18 - 20      Residue name     resName    Residue name.
22           Character        chainID    Chain identifier.
23 - 26      Integer          resSeq     Residue sequence number.
27           AChar            iCode      Code for insertion of residues.
31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
                                         Angstroms
39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
                                         Angstroms
47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
                                         Angstroms
55 - 60      Real(6.2)        occupancy  Occupancy.
61 - 66      Real(6.2)        tempFactor Temperature factor.
77 - 78      LString(2)       element    Element symbol, right-justified.
79 - 80      LString(2)       charge     Charge on the atom.

*/
void  readpdb_v23(pdb_v23 *info, CAcoord *atom, char *inFileName)
{
  //To read both waters and protein atoms in a constant column format
  int i=0;
  int numberofResid=0;
  
  char c_buffer[9]; //Stands for (c)haracter buffer!!!
  memset(c_buffer,'\0', 9);

  char line[100];
  memset(line,'\0', 100);

  //  char *line_ptr;
  //  char CAtest[5]={'\0', '\0', '\0', '\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  while(1)
    {
      //      line_ptr=;
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+6,   5);
	  info[i].serial=atoi(c_buffer);

	  memset(info[i].name,'\0', 5);
	  strncpy(info[i].name,      line+12,  4);

	  info[i].altLoc=line[16];

	  memset(info[i].resName,'\0', 4);
	  strncpy(info[i].resName,   line+17,  3);

	  info[i].chainID=line[21];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+22,  4);
	  info[i].resSeq=atoi(c_buffer);

	  info[i].iCode=line[26];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+30,  8);
	  info[i].x=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+38,  8);
	  info[i].y=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+46,  8);
	  info[i].z=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+54,  6);
	  info[i].occupancy=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+60,  6);
	  info[i].tempFactor=atof(c_buffer);

	  memset(info[i].element,'\0', 2);
	  strncpy(info[i].element,   line+76,  2);

	  memset(info[i].charge,'\0', 2);
	  strncpy(info[i].charge,    line+78,  2);

	  if(0)	
	    fprintf(stdout,"ATOM   %d %s %c %s %c %d %c %lf %lf %lf %lf %lf %s %s\n", \
		    info[i].serial, info[i].name, info[i].altLoc, info[i].resName, info[i].chainID, info[i].resSeq, info[i].iCode, 
		    info[i].x, info[i].y, info[i].z, info[i].occupancy, info[i].tempFactor, info[i].element, info[i].charge);
	  
	  if((strncmp(" CA", line+12, 3))==0)
	    {
	      strncpy(atom[numberofResid].residname, info[i].resName, 3);
	      atom[numberofResid].residname[3]='\0'; //Add null character to the end of the string!
 	      atom[numberofResid].chain=info[i].chainID;
	      atom[numberofResid].resNo=info[i].resSeq;
	      atom[numberofResid].x=info[i].x;
	      atom[numberofResid].y=info[i].y;
	      atom[numberofResid].z=info[i].z;
	      atom[numberofResid].occ=info[i].occupancy;
	      atom[numberofResid].beta=info[i].tempFactor;
	      
	      numberofResid+=1;
	    }
	  
	  i+=1; 
	}
    }
  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }
  fclose(pdbdata);
}

int  *scanpdb_CAandOH2(int *atomCount, char *inFileName)
{
/*Purpose: To determine the total number of atoms and residues for memory allocation*/
/*In future I may use it to keep protein chain and segment infos*/

  int i=0;
  int numberofResids=0;
  int numberofChains=0;                            
  int numberofHelices=0;
  int numberofSheets=0;
  char line[100];
  
  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }

  while(1)
    {
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  if(((strncmp(" CA ", line+12, 4))==0 ) && ((line[16]==' ') || (line[16]=='A')  ) )
	    numberofResids+=1;
	  
	  else if((strncmp(" OH2", line+12, 4))==0)
	    numberofResids+=1;
	  	  
	  i+=1; 
	}
      if((strncmp("TER   ", line, 6))==0)
	{
	  numberofChains+=1;
	}

      if((strncmp("HELIX ", line, 6))==0)
	{
	  numberofHelices+=1;
	}

      if((strncmp("SHEET ", line, 6))==0)
	{
	  numberofSheets+=1;
	}
    }

  atomCount[0]=i;               //atomCount[0] is the total number of atoms
  atomCount[1]=numberofResids;  //atomCount[1] is the total number of residues
  atomCount[2]=numberofChains;  //atomCount[2] is the total number of chains
  atomCount[3]=numberofHelices; //atomCount[3] is the total number of helices
  atomCount[4]=numberofSheets;  //atomCount[4] is the total number of sheets
  if(atomCount[2]==0)   atomCount[2]=1;  /*Remember that if there is just one chain, you will not see any */
                                         /*TER signal to count chain numbers.*/
  if((atomCount[0]==0) || (atomCount[1]==0))
    {
      fprintf(stderr, "Error: Can not read atoms and residues\n");
      exit(EXIT_FAILURE);
    } 
  else
    {
      fprintf(stdout, "Total number of atoms:%d\n", atomCount[0]);
      fprintf(stdout, "Number of residues:%d\n",    atomCount[1]);
      fprintf(stdout, "Number of chains:%d\n",      atomCount[2]);
      fprintf(stdout, "Number of helices:%d\n",     atomCount[3]);
      fprintf(stdout, "Number of sheets:%d\n",      atomCount[4]);
    }
  fclose(pdbdata);
  
  return atomCount; 
}
int  *scanpdb_CAandOH2_v1(int *atomCount, char *inFileName)
{
/*Purpose: To determine the total number of atoms and residues for memory allocation*/
/*In future I may use it to keep protein chain and segment infos*/
/*This function returns number of waters and CA also. atomCount array has to have 7 elements.*/

  int i=0;
  int numberofResids=0;
  int numberofCA=0;
  int numberofOH2=0;
  int numberofChains=0;                            
  int numberofHelices=0;
  int numberofSheets=0;
  char line[100];
  
  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }

  while(1)
    {
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  if(((strncmp(" CA", line+12, 3))==0 ) && ((line[16]==' ') || (line[16]=='A')  ) )
	    numberofCA+=1;
	  
	  else if((strncmp(" OH2", line+12, 4))==0)
	    numberofOH2+=1;
	  	  
	  i+=1; 
	}
      if((strncmp("TER   ", line, 6))==0)
	{
	  numberofChains+=1;
	}

      if((strncmp("HELIX ", line, 6))==0)
	{
	  numberofHelices+=1;
	}

      if((strncmp("SHEET ", line, 6))==0)
	{
	  numberofSheets+=1;
	}
    }

  numberofResids=numberofCA+numberofOH2;
  atomCount[0]=i;               //atomCount[0] is the total number of atoms
  atomCount[1]=numberofResids;  //atomCount[1] is the total number of residues
  atomCount[2]=numberofChains;  //atomCount[2] is the total number of chains
  atomCount[3]=numberofHelices; //atomCount[3] is the total number of helices
  atomCount[4]=numberofSheets;  //atomCount[4] is the total number of sheets
  atomCount[5]=numberofCA;      //atomCount[5] is the total number of CA
  atomCount[6]=numberofOH2;     //atomCount[6] is the total number of waters
  if(atomCount[2]==0)   atomCount[2]=1;  /*Remember that if there is just one chain, you will not see any */
                                         /*TER signal to count chain numbers.*/
  if((atomCount[0]==0) || (atomCount[1]==0))
    {
      fprintf(stderr, "Error: Can not read atoms and residues\n");
      exit(EXIT_FAILURE);
    } 
  else
    {
      fprintf(stdout, "Total number of atoms:%d\n", atomCount[0]);
      fprintf(stdout, "Number of residues:%d\n",    atomCount[1]);
      fprintf(stdout, "Number of chains:%d\n",      atomCount[2]);
      fprintf(stdout, "Number of helices:%d\n",     atomCount[3]);
      fprintf(stdout, "Number of sheets:%d\n",      atomCount[4]);
      fprintf(stdout, "Number of CA:%d\n",          atomCount[5]);
      fprintf(stdout, "Number of waters:%d\n",      atomCount[6]);
    }
  fclose(pdbdata);
  
  return atomCount; 
}

void  readpdb_CAandOH2(pdb_v23 *info, CAcoord *atom, char *inFileName)
{
  //To read both waters and protein atoms in a constant column format
  int i=0;
  int numberofResid=0;
  
  char c_buffer[9]; //Stands for (c)haracter buffer!!!
  memset(c_buffer,'\0', 9);

  char line[100];
  memset(line,'\0', 100);

  //  char *line_ptr;
  //  char CAtest[5]={'\0', '\0', '\0', '\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  int numCA=0;
  int numOH2=0;
  while(1)
    {
      //      line_ptr=;
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+6,   5);
	  info[i].serial=atoi(c_buffer);

	  memset(info[i].name,'\0', 5);
	  strncpy(info[i].name,      line+12,  4);

	  info[i].altLoc=line[16];

	  memset(info[i].resName,'\0', 4);
	  strncpy(info[i].resName,   line+17,  3);

	  info[i].chainID=line[21];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+22,  4);
	  info[i].resSeq=atoi(c_buffer);

	  info[i].iCode=line[26];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+30,  8);
	  info[i].x=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+38,  8);
	  info[i].y=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+46,  8);
	  info[i].z=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+54,  6);
	  info[i].occupancy=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+60,  6);
	  info[i].tempFactor=atof(c_buffer);

	  memset(info[i].element,'\0', 2);
	  strncpy(info[i].element,   line+76,  2);

	  memset(info[i].charge,'\0', 2);
	  strncpy(info[i].charge,    line+78,  2);

	  if(0)	
	    fprintf(stdout,"ATOM   %d %s %c %s %c %d %c %lf %lf %lf %lf %lf %s %s\n", \
		    info[i].serial, info[i].name, info[i].altLoc, info[i].resName, info[i].chainID, info[i].resSeq, info[i].iCode, 
		    info[i].x, info[i].y, info[i].z, info[i].occupancy, info[i].tempFactor, info[i].element, info[i].charge);

	  //If there is an alternative location, I select the first one. This selection may not be healthy depending on the job one is working on!!	  
	  if(((strncmp(" CA ", line+12, 4))==0 ) && ((info[i].altLoc==' ') || (info[i].altLoc=='A')  ) )
	    {
	      strncpy(atom[numberofResid].residname, info[i].resName, 3);
	      atom[numberofResid].residname[3]='\0'; //Add null character to the end of the string!
	      atom[numberofResid].chain=info[i].chainID;
	      atom[numberofResid].resNo=info[i].resSeq;
	      atom[numberofResid].x=info[i].x;
	      atom[numberofResid].y=info[i].y;
	      atom[numberofResid].z=info[i].z;
	      if(0)	
		fprintf(stdout, "%s\t%c\t%d\t%.3lf\t%.3lf\t%.3lf\n",atom[numberofResid].residname, \
			atom[numberofResid].chain, atom[numberofResid].resNo, atom[numberofResid].x,\
			atom[numberofResid].y, atom[numberofResid].z );
	      
	      numCA+=1;
	      numberofResid+=1;
	    }

	  if(((strncmp(" OH2", line+12, 4))==0))
	    {
	      strncpy(atom[numberofResid].residname, info[i].resName, 3);
	      atom[numberofResid].chain=info[i].chainID;
	      atom[numberofResid].resNo=info[i].resSeq;
	      atom[numberofResid].x=info[i].x;
	      atom[numberofResid].y=info[i].y;
	      atom[numberofResid].z=info[i].z;
	      numOH2+=1;
	      numberofResid+=1;
	    }

	  i+=1; 
	}
    }
  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }
  printf("Number of CA=%d\tNumber of waters=%d\n", numCA, numOH2);
  fclose(pdbdata);
}
void readpdbV23_v1(char *inFileName, pdb_v23 *infoAll, pdb_v23 *subSet, CAcoord *atom, \
		   int readCA/*1:YES, 0:NO*/, int readWater/*1:YES, 0:NO*/, int readHetAtom/*1:YES, 0:NO*/)
{
  //To read both waters and protein atoms in a constant column format
  int i=0;
  int numberofResid=0;
  
  char c_buffer[9]; //Stands for (c)haracter buffer!!!
  memset(c_buffer,'\0', 9);

  char line[100];
  memset(line,'\0', 100);

  //  char *line_ptr;
  //  char CAtest[5]={'\0', '\0', '\0', '\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  int numCA=0;
  int numOH2=0;
  while(1)
    {
      //      line_ptr=;
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if(  ((strncmp("ATOM  ", line, 6))==0) || ( readHetAtom && ((strncmp("HETATM", line, 6))==0) )  )
	{
	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+6,   5);
	  infoAll[i].serial=atoi(c_buffer);

	  memset(infoAll[i].name,'\0', 5);
	  strncpy(infoAll[i].name,      line+12,  4);

	  infoAll[i].altLoc=line[16];

	  memset(infoAll[i].resName,'\0', 4);
	  strncpy(infoAll[i].resName,   line+17,  3);

	  infoAll[i].chainID=line[21];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+22,  4);
	  infoAll[i].resSeq=atoi(c_buffer);

	  infoAll[i].iCode=line[26];

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+30,  8);
	  infoAll[i].x=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+38,  8);
	  infoAll[i].y=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+46,  8);
	  infoAll[i].z=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+54,  6);
	  infoAll[i].occupancy=atof(c_buffer);

	  memset(c_buffer,'\0', 9);
	  strncpy(c_buffer,         line+60,  6);
	  infoAll[i].tempFactor=atof(c_buffer);

	  memset(infoAll[i].element,'\0', 2);
	  strncpy(infoAll[i].element,   line+76,  2);

	  memset(infoAll[i].charge,'\0', 2);
	  strncpy(infoAll[i].charge,    line+78,  2);

	  if(0)	
	    fprintf(stdout,"ATOM   %d %s %c %s %c %d %c %lf %lf %lf %lf %lf %s %s\n", \
		    infoAll[i].serial, infoAll[i].name, infoAll[i].altLoc, infoAll[i].resName, infoAll[i].chainID, infoAll[i].resSeq, infoAll[i].iCode, 
		    infoAll[i].x, infoAll[i].y, infoAll[i].z, infoAll[i].occupancy, infoAll[i].tempFactor, infoAll[i].element, infoAll[i].charge);

	  //If there is an alternative location, I select the first one. This selection may not be healthy depending on the job one is working on!!	  
	  if(  (readCA) && ((strncmp(" CA", line+12, 3))==0 ) && ( (infoAll[i].altLoc==' ') || (infoAll[i].altLoc=='A') )  )
	    {
	      strncpy(atom[numberofResid].residname, infoAll[i].resName, 3);
	      atom[numberofResid].chain=infoAll[i].chainID;
	      atom[numberofResid].resNo=infoAll[i].resSeq;
	      atom[numberofResid].x=infoAll[i].x;
	      atom[numberofResid].y=infoAll[i].y;
	      atom[numberofResid].z=infoAll[i].z;
	      if(0)	
		fprintf(stdout, "%s\t%c\t%d\t%.3lf\t%.3lf\t%.3lf\n",atom[numberofResid].residname, \
			atom[numberofResid].chain, atom[numberofResid].resNo, atom[numberofResid].x,\
			atom[numberofResid].y, atom[numberofResid].z );
	      
	      numCA+=1;
	      numberofResid+=1;
	    }

	  if(  (readWater) && ((strncmp(" OH2", line+12, 4))==0)  )
	    {
	      strncpy(atom[numberofResid].residname, infoAll[i].resName, 3);
	      atom[numberofResid].chain=infoAll[i].chainID;
	      atom[numberofResid].resNo=infoAll[i].resSeq;
	      atom[numberofResid].x=infoAll[i].x;
	      atom[numberofResid].y=infoAll[i].y;
	      atom[numberofResid].z=infoAll[i].z;
	      numOH2+=1;
	      numberofResid+=1;
	    }

	  i+=1; 
	}
    }
  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }
  printf("Number of CA=%d\tNumber of waters=%d\n", numCA, numOH2);
  fclose(pdbdata);

}

/* void selectSubSet(int *atomCount, pdb_v23 *infoAll, pdb_v23 *subSet) */
/* { */
/*   int i=0; */

/*   for(i=0; i<atomCount[0]; i++) */
/*     { */


/*     } */
/* } */

void readHelices(char *inFileName, helixInfo *helix)
{
  int i=0, j=0;//My infamous regular counter
  //All C_ prefixed variables are integers and they will be converted
  //after reading them properly!
  char  f_C_serNum[4]=    {'\0', '\0', '\0', '\0'};
  char  f_helixID[4]=     {'\0', '\0', '\0', '\0'};
  char  f_initResName[4]= {'\0', '\0', '\0', '\0'};
  char  f_initChainID[2]= {'\0', '\0'};
  char  f_C_initSeqNum[5]={'\0','\0', '\0', '\0', '\0'};
  char  f_initICode[2]=   {'\0', '\0'};
  char  f_endResName[4]=  {'\0', '\0', '\0', '\0'};
  char  f_endChainID[2]=  {'\0', '\0'};
  char  f_C_endSeqNum[5]= {'\0','\0', '\0', '\0', '\0'};
  char  f_endICode[2]=    {'\0', '\0'};
  char  f_C_helixClass[3]={'\0', '\0', '\0'};
  char  f_comment[31];
  for(j=0; j<31; j++)   f_comment[j]='\0';
  char  f_C_length[6]=    {'\0','\0', '\0', '\0', '\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  char *line_ptr=NULL;
  char  line[81];
  //Why do I use strncpy?: The format of pdb file is a constant column format!
  //Why do I use strcpy after strncpy?: Just not to mess up w/ null terminating character in the string!
  while(1)
    {
      line_ptr=fgets(line, sizeof(line),pdbdata);
      if(line_ptr==NULL) break;
      
      if((strncmp("HELIX ", line, 6))==0)
	{

	  strncpy(f_C_serNum,     line+7,   3);
	  helix[i].serNum=atoi(f_C_serNum);

	  strncpy(f_helixID,      line+11,  3);
	  strcpy(helix[i].helixID, f_helixID);

	  strncpy(f_initResName,  line+15,  3);
	  strcpy(helix[i].initResName, f_initResName);

	  strncpy(f_initChainID,  line+19,  1);
	  strcpy(helix[i].initChainID, f_initChainID);

	  strncpy(f_C_initSeqNum, line+21,  4);
	  helix[i].initSeqNum=atoi(f_C_initSeqNum);

	  strncpy(f_initICode,    line+25,  1);
	  strcpy(helix[i].initICode, f_initICode);

	  strncpy(f_endResName,   line+27,  3);
	  strcpy(helix[i].endResName, f_endResName);

	  strncpy(f_endChainID,   line+31,  1);
	  strcpy(helix[i].endChainID, f_endChainID);

	  strncpy(f_C_endSeqNum,  line+33,  4);
	  helix[i].endSeqNum=atoi(f_C_endSeqNum);

	  strncpy(f_endICode,     line+37,  1);
	  strcpy(helix[i].endICode, f_endICode);

	  strncpy(f_C_helixClass, line+38,  2);
	  helix[i].helixClass=atoi(f_C_helixClass);

	  strncpy(f_comment,      line+40, 30);
	  strcpy(helix[i].comment, f_comment);

	  strncpy(f_C_length,     line+71,  5);
	  helix[i].length=atoi(f_C_length);
	  if(0)
	    {
	      fprintf(stdout, "%s\n", line);
	      fprintf(stdout, "HELIX  %s %s %s %s %s%s %s %s %s%s%s%s %s\n", \
		      f_C_serNum, f_helixID, f_initResName, f_initChainID, \
		      f_C_initSeqNum, f_initICode, f_endResName, f_endChainID, \
		      f_C_endSeqNum, f_endICode, f_C_helixClass, f_comment, f_C_length);
	      
	      fprintf(stdout, "HELIX  %3d %s %s %s %4d%s %s %s %4d%s%2d%s %5d\n", \
		      helix[i].serNum, helix[i].helixID, helix[i].initResName, \
		      helix[i].initChainID, helix[i].initSeqNum, helix[i].initICode, \
		      helix[i].endResName, helix[i].endChainID, helix[i].endSeqNum, \
		      helix[i].endICode, helix[i].helixClass, helix[i].comment, helix[i].length);
	    }
	  
	  i+=1; 
	}
    }

  fclose(pdbdata);

}

void readSheets(char *inFileName, sheetInfo *sheet)
{
  int i=0;//My infamous regular counter
  //All C_ prefixed variables are integers and they will be converted
  //after reading them properly!
  char f_C_strand[4]=    {'\0', '\0', '\0', '\0'};
  char f_sheetID[4]=     {'\0', '\0', '\0', '\0'};
  char f_C_numStrands[3]={'\0', '\0', '\0'};

  char f_initResName[4]= {'\0', '\0', '\0', '\0'};
  char f_initChainID[2]= {'\0', '\0'};
  char f_C_initSeqNum[5]={'\0','\0', '\0', '\0', '\0'};
  char f_initICode[2]=   {'\0', '\0'};

  char f_endResName[4]=  {'\0', '\0', '\0', '\0'};
  char f_endChainID[2]=  {'\0', '\0'};
  char f_C_endSeqNum[5]= {'\0','\0', '\0', '\0', '\0'};
  char f_endICode[2]=    {'\0', '\0'};
  
  char f_C_sense[3]=     {'\0', '\0', '\0'};

  char f_curAtom[5]=     {'\0','\0', '\0', '\0', '\0'};
  char f_curResName[4]=  {'\0', '\0', '\0', '\0'};
  char f_curChainID[2]=  {'\0', '\0'};
  char f_C_curResSeq[5]= {'\0','\0', '\0', '\0', '\0'};
  char f_curICode[2]=    {'\0', '\0'};

  char f_prevAtom[5]=     {'\0','\0', '\0', '\0', '\0'};
  char f_prevResName[4]=  {'\0', '\0', '\0', '\0'};
  char f_prevChainID[2]=  {'\0', '\0'};
  char f_C_prevResSeq[5]= {'\0','\0', '\0', '\0', '\0'};
  char f_prevICode[2]=    {'\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  char *line_ptr=NULL;
  char  line[81];
  //Why do I use strncpy?: The format of pdb file is a constant column format!
  //Why do I use strcpy after strncpy?: Just not to mess up w/ null terminating character in the string!
  while(1)
    {
      line_ptr=fgets(line, sizeof(line),pdbdata);
      if(line_ptr==NULL) break;
      
      if((strncmp("SHEET ", line, 6))==0)
	{

	  strncpy(f_C_strand,     line+7,   3);
	  sheet[i].strand=atoi(f_C_strand);

	  strncpy(f_sheetID,      line+11,  3);
	  strcpy(sheet[i].sheetID, f_sheetID);

	  strncpy(f_C_numStrands,  line+14,   2);
	  sheet[i].numStrands=atoi(f_C_numStrands);

	  strncpy(f_initResName,  line+17,  3);
	  strcpy(sheet[i].initResName, f_initResName);

	  strncpy(f_initChainID,  line+21,  1);
	  strcpy(sheet[i].initChainID, f_initChainID);

	  //Dont forget that if the string is NULL, the integer number is 0
	  strncpy(f_C_initSeqNum, line+22,  4);
	  sheet[i].initSeqNum=atoi(f_C_initSeqNum);

	  strncpy(f_initICode,    line+26,  1);
	  strcpy(sheet[i].initICode, f_initICode);

	  strncpy(f_endResName,   line+28,  3);
	  strcpy(sheet[i].endResName, f_endResName);

	  strncpy(f_endChainID,   line+32,  1);
	  strcpy(sheet[i].endChainID, f_endChainID);

	  //Dont forget that if the string is NULL, the integer number is 0
	  strncpy(f_C_endSeqNum,  line+33,  4);
	  sheet[i].endSeqNum=atoi(f_C_endSeqNum);

	  strncpy(f_endICode,     line+37,  1);
	  strcpy(sheet[i].endICode, f_endICode);


	  strncpy(f_C_sense, line+38,  2);
	  sheet[i].sense=atoi(f_C_sense);


	  strncpy(f_curAtom,   line+41,  4);
	  strcpy(sheet[i].curAtom, f_curAtom);

	  strncpy(f_curResName,   line+45,  3);
	  strcpy(sheet[i].curResName, f_curResName);

	  strncpy(f_curChainID,   line+49,  1);
	  strcpy(sheet[i].curChainID, f_curChainID);

	  //Dont forget that if the string is NULL, the integer number is 0
	  strncpy(f_C_curResSeq,  line+50,  4);
	  sheet[i].curResSeq=atoi(f_C_curResSeq);

	  strncpy(f_curICode,     line+54,  1);
	  strcpy(sheet[i].curICode, f_curICode);


	  strncpy(f_prevAtom,   line+56,  4);
	  strcpy(sheet[i].prevAtom, f_prevAtom);

	  strncpy(f_prevResName,   line+60,  3);
	  strcpy(sheet[i].prevResName, f_prevResName);

	  strncpy(f_prevChainID,   line+64,  1);
	  strcpy(sheet[i].prevChainID, f_prevChainID);

	  //Dont forget that if the string is NULL, the integer number is 0
	  strncpy(f_C_prevResSeq,  line+65,  4);
	  sheet[i].prevResSeq=atoi(f_C_prevResSeq);

	  strncpy(f_prevICode,     line+69,  1);
	  strcpy(sheet[i].prevICode, f_prevICode);

	  if(0)
	    {
	      fprintf(stdout, "%s\n", line);
	      fprintf(stdout, "SHEET  %s %s%s %s %s%s%s %s %s%s%s%s %s%s %s%s%s %s%s %s%s%s\n",	\
		      f_C_strand, f_sheetID, f_C_numStrands, f_initResName, f_initChainID,     \
		      f_C_initSeqNum, f_initICode, f_endResName, f_endChainID, f_C_endSeqNum,  \
		      f_endICode, f_C_sense, f_curAtom, f_curResName, f_curChainID,            \
		      f_C_curResSeq, f_curICode, f_prevAtom, f_prevResName, f_prevChainID,     \
		      f_C_prevResSeq, f_prevICode);

	      fprintf(stdout, "SHEET  %3d %s%2d %s %s%4d%s %s %s%4d%s%2d %s%s %s%4d%s %s%s %s%4d%s\n",   \
		      sheet[i].strand, sheet[i].sheetID, sheet[i].numStrands, sheet[i].initResName,      \
		      sheet[i].initChainID, sheet[i].initSeqNum, sheet[i].initICode, sheet[i].endResName,\
		      sheet[i].endChainID, sheet[i].endSeqNum, sheet[i].endICode, sheet[i].sense,        \
		      sheet[i].curAtom, sheet[i].curResName, sheet[i].curChainID, sheet[i].curResSeq,    \
		      sheet[i].curICode, sheet[i].prevAtom, sheet[i].prevResName, sheet[i].prevChainID,	 \
		      sheet[i].prevResSeq, sheet[i].prevICode);
	    }
	  
	  i+=1; 
	}
    }

  fclose(pdbdata);
}
int scan4Helices(char *inFileName)
{
  int numberofHelices=0;

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file to read\n", inFileName);
      exit(EXIT_FAILURE);
    }
  char *line_ptr=NULL;
  char  line[81];

  while(1)
    {
      line_ptr=fgets(line, sizeof(line),pdbdata);
      if(line_ptr==NULL) break;
      
      if((strncmp("HELIX ", line, 6))==0)
	{
	  numberofHelices+=1; 
	}
    }

  fclose(pdbdata);
  return numberofHelices;
}

int scan4Sheets(char *inFileName)
{
  int numberofSheets=0;

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  char *line_ptr=NULL;
  char  line[81];

  while(1)
    {
      line_ptr=fgets(line, sizeof(line),pdbdata);
      if(line_ptr==NULL) break;
      
      if((strncmp("SHEET ", line, 6))==0)
	{
	  numberofSheets+=1; 
	}
    }

  fclose(pdbdata);
  return numberofSheets;
}
/* int write2pdb_v2(int N, char *inFile,  */
/* 		 CAcoord *atomUpdated,    */
/* 		 CAcoord *atomsFixed1, */
/* 		 double RMSD, */
/* 		 char **argv, */
/* 		 double scaler, double *w/\*Weight for RMSD*\/) */
/* { */
/*   //Purpose: This function writes pdb trajectory files without removal of translations and rotations */

/*   //  if(RMSD<0.01) */
/*   //    { */
/*   //    return 0; */
/*   //    } */
/*   //  else */
/*   //    { */
/*   char outFile[255]="path_"; */
/*   strncat(outFile, argv[1], 4); */
/*   strncat(outFile, "_", 1); */
/*   strncat(outFile, argv[2], 4); */
/*   strncat(outFile, "_Rc", 3); */
/*   strcat(outFile, argv[3]); */
/*   strncat(outFile, "_Scal", 5); */
/*   strcat(outFile, argv[4]); */
/*   strncat(outFile, "_expli", 6); */
/*   strcat(outFile, argv[5]); */
/*   strncat(outFile, ".pdb", 4); */
/*   //  sprintf(); */
/*   int  i=0, k=0; */
/*   int numberofResid=0; */
  
/*   //Instead of keeping pdb data on memory, it reads pdb file again to get pdb data */
/*   char line[256]; */
/*       char* line_ptr; */
/*       char comp2[4]="CA"; */
      
/*       /\*Allocate memory for structures*\/ */
/*       pdb *info = (pdb*)malloc(N*sizeof(pdb));/\*map1.pdb has just CA residues. That's the reason of N*\/ */
/*       if (info == NULL)  */
/* 	{ */
/* 	  fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb info\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       CAcoord *atom = (CAcoord*)malloc(N*sizeof(CAcoord)); */
/*       if (atom == NULL)  */
/* 	{ */
/* 	  fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord *atom\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
      
/*       FILE *pdbdata=fopen(inFile,"r"); */
/*       if(pdbdata==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "No such file:%s\n", inFile);  */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       while(1) */
/* 	{ */
/* 	  line_ptr=fgets(line, sizeof(line),pdbdata); */
/* 	  if(line_ptr==NULL) break; */
	  
/* 	  if((strncmp("ATOM", line, 4))==0) */
/* 	    { */
/* 	      sscanf(line,"%s %d %s%s %c %d %lf %lf %lf %lf %lf",	\ */
/* 		     info[i].text, &info[i].atnum, info[i].attyp, info[i].resname,  */
/* 		     &info[i].chain, &info[i].resno, &info[i].coordX, &info[i].coordY,  */
/* 		     &info[i].coordZ, &info[i].occ, &info[i].temp); */
	      
/* 	      if((strncmp(comp2, info[i].attyp, 3))==0) */
/* 		{ */
		  
/* 		  (atom[numberofResid].chain)=(info[i].chain); */
/* 		  atom[numberofResid].resNo=info[i].resno; */
/* 		  strncpy(atom[numberofResid].residname, info[i].resname, 3); */
/* 		  numberofResid+=1; */
/* 		} */
/* 	      i+=1; */
/* 	    } */
/* 	} */
/*       fclose(pdbdata); */
/*       free(info); */
/*       //-------------------------------------------------------------------------- */
/*       //This part is used to superimpose each frame to the final structure so that */
/*       //a smooth transition path can be observed.  */
      
/*       double *resultVec = (double*)calloc((3*N), sizeof(double)); */
/*       if (resultVec==NULL) */
/* 	{     */
/* 	  fprintf(stderr, "Program could not allocate memory for \"resultVec\"\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       resultVec=superpose_v2(N, atomUpdated, resultVec, w); */
/*       //---------------------------------------------------------------------------   */
/*       FILE *modespdb=fopen(outFile, "a"); */
/*       if(modespdb==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "%s has not been produced\n", outFile); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       k=0; */
/*       while(k<N) */
/* 	{ */
/* 	  fprintf(modespdb, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (k+1), " CA ",  */
/* 		  atom[k].residname,  atom[k].chain, atom[k].resNo, resultVec[3*k], resultVec[3*k+1], resultVec[3*k+2], 1.00, 0.00);  */
/* /\* 	  if((k!=0)&&(atom[k].chain!=atom[k-1].chain)) *\/ */
/* /\* 	    fprintf(modespdb, "TER\n"); *\/ */
/* 	  k++; */
/* 	} */
/*       fprintf(modespdb,"ENDMDL\n\n"); */
/*       free(atom); */
/*       free(resultVec); */
/*       fclose(modespdb); */
/*       return 1; */
/*       // } */
/* } */

void write2pdb_v23_CA(int N/*System size*/, CAcoord *atom, char *outFile, int modelNo)
{
  //Purpose: To write CA coordinates of "N" atoms specified in "CAcoord *atoms"
  //         structure in pdb v2.3 format to outFile. 
  //Note   : modelNo starts from 1. 
  int i=0;
  FILE *pdbdata=fopen(outFile, "w");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "%s has not been produced\n", outFile);
      exit(EXIT_FAILURE);
    }
  
  fprintf(pdbdata, "MODEL     %4d\n", modelNo);
  for(i=0; i<N; i++)
    {
      fprintf(pdbdata, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " CA ", 
	      atom[i].residname,  atom[i].chain, atom[i].resNo, atom[i].x, atom[i].y, atom[i].z, atom[i].occ, atom[i].beta); 
      /* 	  if((i!=0)&&(atom[i].chain!=atom[i-1].chain)) */
      /* 	    fprintf(modespdb, "TER\n"); */
      //      i++;
    }
  fprintf(pdbdata, "ENDMDL\n\n");
  fclose(pdbdata);
}


int  *scanpdb_mm(int *atomCount, char *inFileName)
{
/*Purpose: To determine the total number of atoms and residues for memory allocation*/
/*In future I may use it to keep protein chain and segment infos*/
/*What's new: I tried to extend it so that it can check multiple models. 04/04/2011*/
/*04/26/2011: I added a new argument(int modelNo) so that it will get just one models data.*/
  //This function read only one model from a  multimodel file

  int i=0;
  int numberofResids=0;
  int numberofChains=0;                            
  int numberofHelices=0;
  int numberofSheets=0;
  int numberofModels=0;
  char line[100];
  
  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  while(1)
    {
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  //Only CA atoms with null or 'A' alternate location indicator are selected. 
	  if(((strncmp(" CA", line+12, 3))==0 ) && ((line[16]==' ') || (line[16]=='A')  ) )
	    numberofResids+=1;
	  
	  else if((strncmp(" OH2", line+12, 4))==0)
	    numberofResids+=1;
	  
	  i+=1; 
	}
      if((strncmp("TER   ", line, 6))==0)
	{
	  numberofChains+=1;
	}
      
      if((strncmp("HELIX ", line, 6))==0)
	{
	  numberofHelices+=1;
	}
      
      if((strncmp("SHEET ", line, 6))==0)
	{
	  numberofSheets+=1;
	}
      if((strncmp("MODEL ", line, 6))==0)
	{
	  numberofModels+=1;
	}
    }
  atomCount[0]=i;               //atomCount[0] is the total number of atoms
  atomCount[1]=numberofResids;  //atomCount[1] is the total number of residues
  atomCount[2]=numberofChains;  //atomCount[2] is the total number of chains
  atomCount[3]=numberofHelices; //atomCount[3] is the total number of helices
  atomCount[4]=numberofSheets;  //atomCount[4] is the total number of sheets
  atomCount[5]=numberofModels;  //atomCount[5] is the total number of models
  if(atomCount[2]==0)   atomCount[2]=1;  /*Remember that if there is just one chain, you will not see any TER signal to count chain numbers.*/
  if(atomCount[5]==0)   atomCount[5]=1;  /*Remember that if there is just one model, you will not see any MODEL word to count model numbers.*/

  if((atomCount[0]==0) || (atomCount[1]==0))
    {
      fprintf(stderr, "Error: Can not read atoms and residues\n");
      exit(EXIT_FAILURE);
    } 
  else
    {
      fprintf(stdout, "Total number of atoms:%d\n", atomCount[0]);
      fprintf(stdout, "Total number of residues:%d\n",    atomCount[1]);
      fprintf(stdout, "Total number of chains:%d\n",      atomCount[2]);
      fprintf(stdout, "Total number of helices:%d\n",     atomCount[3]);
      fprintf(stdout, "Total number of sheets:%d\n",      atomCount[4]);
      fprintf(stdout, "Total number of models:%d\n",      atomCount[5]);
    }
  fclose(pdbdata);
  
  return atomCount; 
}

void  readpdb_mm(pdb_v23 *info, CAcoord *atom, char *inFileName, int *atomCount, int modelNo)
{
  //To read both waters and protein atoms in a constant column format
  /*What's new: I tried to extend it so that it can read on of multiple models. 04/04/2011*/
  //This function read only one model from a  multimodel file

  int i=0;
  int numberofResid=0;
  
  if(atomCount[5]==1)
    printf("WARNING: This function is NOT appropriate to read one model pdb files!\n");
  char c_buffer[9]; //Stands for (c)haracter buffer!!!
  memset(c_buffer,'\0', 9);

  char line[100];
  memset(line,'\0', 100);

  //  char *line_ptr;
  //  char CAtest[5]={'\0', '\0', '\0', '\0', '\0'};

  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  int numCA=0;
  int numOH2=0;
  while(1)
    {
      //      line_ptr=;
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if(((strncmp("MODEL ", line, 6))==0) && (atoi(strncpy(c_buffer, line+10, 4))==modelNo))
	{
	  while( fgets(line, sizeof(line),pdbdata)!=NULL)
	    if((strncmp("ATOM  ", line, 6))==0)
	      {
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+6,   5);
		info[i].serial=atoi(c_buffer);
		
		memset(info[i].name,'\0', 5);
		strncpy(info[i].name,      line+12,  4);
		
		info[i].altLoc=line[16];
		
		memset(info[i].resName,'\0', 4);
		strncpy(info[i].resName,   line+17,  3);
		
		info[i].chainID=line[21];
		
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+22,  4);
		info[i].resSeq=atoi(c_buffer);
		
		info[i].iCode=line[26];
		
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+30,  8);
		info[i].x=atof(c_buffer);
		
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+38,  8);
		info[i].y=atof(c_buffer);
		
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+46,  8);
		info[i].z=atof(c_buffer);
		
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+54,  6);
		info[i].occupancy=atof(c_buffer);
		
		memset(c_buffer,'\0', 9);
		strncpy(c_buffer,         line+60,  6);
		info[i].tempFactor=atof(c_buffer);
		
		memset(info[i].element,'\0', 2);
		strncpy(info[i].element,   line+76,  2);
		
		memset(info[i].charge,'\0', 2);
		strncpy(info[i].charge,    line+78,  2);
		
		if(0)	
		  fprintf(stdout,"ATOM   %d %s %c %s %c %d %c %lf %lf %lf %lf %lf %s %s\n", \
			  info[i].serial, info[i].name, info[i].altLoc, info[i].resName, info[i].chainID, info[i].resSeq, info[i].iCode, 
			  info[i].x, info[i].y, info[i].z, info[i].occupancy, info[i].tempFactor, info[i].element, info[i].charge);
		
		//If there is an alternative location, I select the first one. This selection may not be healthy depending on the job one is working on!!	  
		if(((strncmp(" CA", line+12, 3))==0 ) && ((info[i].altLoc==' ') || (info[i].altLoc=='A')  ) )
		  {
		    strncpy(atom[numberofResid].residname, info[i].resName, 3);
		    atom[numberofResid].chain=info[i].chainID;
		    atom[numberofResid].resNo=info[i].resSeq;
		    atom[numberofResid].x=info[i].x;
		    atom[numberofResid].y=info[i].y;
		    atom[numberofResid].z=info[i].z;
		    numCA+=1;
		    numberofResid+=1;
		  }
		
		if(((strncmp(" OH2", line+12, 4))==0))
		  {
		    strncpy(atom[numberofResid].residname, info[i].resName, 3);
		    atom[numberofResid].chain=info[i].chainID;
		    atom[numberofResid].resNo=info[i].resSeq;
		    atom[numberofResid].x=info[i].x;
		    atom[numberofResid].y=info[i].y;
		    atom[numberofResid].z=info[i].z;
		    numOH2+=1;
		    numberofResid+=1;
		  }
		i+=1; 
	      }
	    else if((strncmp("ENDMDL", line, 6))==0)
	      break;
	}
    }
  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }

  atomCount[0]=i;
  atomCount[1]=numberofResid;
  printf("Number of CA=%d\tNumber of waters=%d\n", numCA, numOH2);
  fclose(pdbdata);
}
void writeCA2pdb(int N/*Number of CA*/, FILE *OUTFILE, CAcoord *atoms)
{
  //Purpose: This function writes CA atoms in pdb format given in the structure 'CAcoord *atoms'!
  //         I use this function to wrote a iENM trajectory frame by frame. 
  //         Keep in mind it ends pdb file with ENDMDL line!
  int k=0;
  for(k=0; k<N; k++)
    {
      fprintf(OUTFILE, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (k+1), " CA ", \
	      atoms[k].residname,  atoms[k].chain, atoms[k].resNo, atoms[k].x, atoms[k].y, atoms[k].z, atoms[k].occ, atoms[k].beta); 

      //If the structure contains multiple chains, we had better put a TER signal between chains!
      if((k!=0)&&(atoms[k].chain!=atoms[k-1].chain)) 
	fprintf(OUTFILE, "TER\n");
    }
  fprintf(OUTFILE,"ENDMDL\n\n");
}

int xyz2pdb(char *seqFileName, char *xyzFileName, char *outFileName, int writeWater, int writeProtein, int printDetails)
{

  //Purpose: If we are given a sequence file (1 letter amino acid codes) and a xyz coordinate file,
  //         this function will convert them into pdb format!
  //  int printDetails : 0 or 1
  //  int writeWater   : 0 or 1
  //  int writeProtein : 0 or 1


  if(printDetails) printf("%s\n", seqFileName);
  if(printDetails) printf("%s\n", xyzFileName);
  if(printDetails) printf("%s\n", outFileName);

  //Read sequence file. 
  FILE *SEQ_FILE=fopen(seqFileName, "r");
  if(SEQ_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not open sequence file: %s\n", seqFileName);
      exit(EXIT_FAILURE);
    }

  //Read xyz file. 
  FILE *XYZ_FILE=fopen(xyzFileName, "r");
  if(XYZ_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not open xyz file: %s\n", xyzFileName);
      exit(EXIT_FAILURE);
    }

  //Write results to outFile
  FILE *OUT_FILE=fopen(outFileName, "w");
  if(OUT_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not produce output file: %s\n", outFileName);
      exit(EXIT_FAILURE);
    }

  int  N=0;  /*System size: Number of CA plus number of waters(or OH2)*/
  int  i=0;  /*My regular counter.*/

  while( (!feof(XYZ_FILE)) || (!feof(SEQ_FILE)))
    {
      int    retValueXYZ=0, retValueSEQ=0;
      char   resName=' ';
      double x=0.000, y=0.000, z=0.000;

      retValueXYZ=fscanf(XYZ_FILE, "%lf %lf %lf\n", &x, &y, &z);
      retValueSEQ=fscanf(SEQ_FILE, "%c\n", &resName);      
      if(0) fprintf(stdout, "%c\t%d\n", resName, retValueSEQ);
      if((retValueXYZ<=0) || (retValueSEQ<=0))
	{
	  fprintf(stderr, "Error happened while reading:\nReturn value for xyz: %c!\nReturn value for seq: %c", retValueXYZ, retValueSEQ);
	  exit(EXIT_FAILURE);
	}
      else
	{


	  if((resName!='3') && writeProtein)
	    {
	      fprintf(OUT_FILE, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " CA ", \
		      aa_1letter_to_3letter(resName),  'X', (i+1), x, y, z, 1.00, 0.00);
	      i++;
	    }
	  else if(writeWater)
	    {
	      fprintf(OUT_FILE, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (i+1), " OH2", \
		      "TIP3",  'Y', (i+1), x, y, z, 1.00, 0.00); 
	      i++;
	    }
	  else
	    {
	      fprintf(stderr, "ERROR: Unknown parameter!\n");
	      exit(EXIT_FAILURE);
	    }
	}
    }

  N=i;
  fprintf(stdout, "Number of residues read=%d\n", N);
  fclose(SEQ_FILE);
  fclose(XYZ_FILE);
  fclose(OUT_FILE);
  return N;
}

void pdb2seqANDxyz(char *pdbFileName, char *seqFileName,  char *xyzFileName, int printDetails)
{

  //Purpose: This function obtains a sequence file and xyz coordinate file from a pdb file

  if(printDetails)
    printf("Sequence file is %s.\n", seqFileName);

  if(printDetails)
    printf("xyz file is %s.\n", xyzFileName);

  //Read pdb file
  FILE *PDB_FILE=fopen(pdbFileName, "r");
  if(PDB_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not open pdb file: %s\n", pdbFileName);
      exit(EXIT_FAILURE);
    }

  //Write CA sequence file. 
  FILE *SEQ_FILE=fopen(seqFileName, "w");
  if(SEQ_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not produce sequence file: %s\n", seqFileName);
      exit(EXIT_FAILURE);
    }

  //Write CA xyz file. 
  FILE *XYZ_FILE=fopen(xyzFileName, "w");
  if(XYZ_FILE==NULL)
    {
      fprintf(stderr, "ERROR: Can not produce  xyz file: %s\n", xyzFileName);
      exit(EXIT_FAILURE);
    }

  char   line[80];
  char   c_buffer[9]; //Stands for (c)haracter buffer!!!
  double x=0.0, y=0.0, z=0.0;
  memset(c_buffer,'\0', 9);

  memset(line, '\0', 80);
  while(!feof(PDB_FILE))
    {
      if(fgets(line, sizeof(line), PDB_FILE)==NULL) break;
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  //If there is an alternative location, I select the first one. This selection may not be healthy depending on the job one is working on!!	  
	  if(  ((strncmp(" CA ", line+12, 4))==0 ) && ( (line[16]==' ') || (line[16]=='A')  ) )
	    {
	      fprintf(SEQ_FILE, "%c\n", aa_3letter_to_1letter( (line+17) ));

	      memset(c_buffer,'\0', 9);
	      strncpy(c_buffer, line+30,  8);
	      x=atof(c_buffer);
	      
	      memset(c_buffer,'\0', 9);
	      strncpy(c_buffer, line+38,  8);
	      y=atof(c_buffer);
	      
	      memset(c_buffer,'\0', 9);
	      strncpy(c_buffer, line+46,  8);
	      z=atof(c_buffer);

	      fprintf(XYZ_FILE, "%8.3lf%8.3lf%8.3lf\n", x, y, z);
	      //	      strncpy(atom[numberofResid].residname, info[i].resName, 3);
	      //	      numberofResid+=1;
	    }
	  //	  i+=1; 
	}
      memset(line, '\0', 80);
    }
  fclose(PDB_FILE);
  fclose(SEQ_FILE);
  fclose(XYZ_FILE);

  if(1) printf("%s has been produced successfully!\n", seqFileName);
  if(1) printf("%s has been produced successfully!\n", xyzFileName);
}


void combineALLtoCG(char *allAtomPdb, char *coarsePdb, char *allProtPlusCGwatPdb)
{
  //Purpose: First file contains just all atom pdb and the second file contains a coarse grained
  //         water and CA model. The purpose here is to replace CA model of protein with all atom
  //         model so that it can be used for other purposes. Here, we assume that there is not any rotation.
  //         This assumption is based on the information of how equillibrated water is added around protein in
  //         addwat function of Yang et al. 
  //=============================================================================================
  //Open all atom pdb!
  int atomCount1[5]={99999/*Atoms*/, 999/*Residues*/, 99/*Chains*/, 100/*Helices*/, 99/*Sheets*/}; 
  //--Scan first conformation
  FILE *pdbdata1=fopen(allAtomPdb,"r");
  if(pdbdata1==NULL)
    {
      fprintf(stderr, "No such file:%s\n", allAtomPdb);
      exit(EXIT_FAILURE);
    }
  fprintf(stdout, "All atom pdb file is %s\n", allAtomPdb);
  scanpdb(atomCount1, allAtomPdb);

  printf("Number of residues:%d\n",  atomCount1[1]);
  
  //--Get data from all atom pdb file to the structures below--
  pdb_v23 *info1 = (pdb_v23*)malloc((atomCount1[0])*sizeof(pdb_v23));
  if (info1 == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb\n");
      exit(EXIT_FAILURE);
    }
  
  CAcoord *prot1atoms = (CAcoord*)malloc((atomCount1[1])*sizeof(CAcoord));
  if (prot1atoms == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord\n");
      exit(EXIT_FAILURE);
    }
  
  readpdb_v23(info1, prot1atoms, allAtomPdb);
  fclose(pdbdata1);
  //=============================================================================================

  //=============================================================================================
  //Open coarse pdb!
  int atomCount2[7]={99999/*Atoms*/, 999/*Residues*/, 99/*Chains*/, 100/*Helices*/, 99/*Sheets*/, 0/*CAs*/, 0/*Waters*/}; 
  //--Scan first conformation
  FILE *pdbdata2=fopen(coarsePdb,"r");
  if(pdbdata2==NULL)
    {
      fprintf(stderr, "No such file:%s\n", coarsePdb);
      exit(EXIT_FAILURE);
    }
  fprintf(stdout, "Coarse pdb file is %s\n", coarsePdb);
  scanpdb_CAandOH2_v1(atomCount2, coarsePdb);

  printf("Number of residues:%d\n",  atomCount2[1]);
  
  //--Get data from all atom pdb file to the structures below--
  pdb_v23 *info2 = (pdb_v23*)malloc((atomCount2[0])*sizeof(pdb_v23));
  if (info2 == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb\n");
      exit(EXIT_FAILURE);
    }
  
  CAcoord *prot2atoms = (CAcoord*)malloc((atomCount2[1])*sizeof(CAcoord));
  if (prot2atoms == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord\n");
      exit(EXIT_FAILURE);
    }
  
  readpdb_CAandOH2(info2, prot2atoms, coarsePdb);
  fclose(pdbdata2);
  //=============================================================================================
  //While addwat immerses protein into an equillibrated water box, it moves the geometrical center to zero.
  //Lets find the translation vector so that we can translate waters to original pdb coordinates! 
  double translationVector[3]={0.00, 0.00, 0.00};

  translationVector[0]= (prot2atoms[0].x-prot1atoms[0].x);
  translationVector[1]= (prot2atoms[0].y-prot1atoms[0].y);
  translationVector[2]= (prot2atoms[0].z-prot1atoms[0].z);

  fprintf(stdout, "Translation vector: x=%lf\ty=%lf\tz=%lf\n", translationVector[0], translationVector[1], translationVector[2]);
  
  FILE *OUT_FILE=fopen(allProtPlusCGwatPdb, "w");
  if(OUT_FILE==NULL)
    {
      fprintf(stderr, "ERROR: File %s can not be produced!!\n", allProtPlusCGwatPdb);
      exit(EXIT_FAILURE);
    }

  int i=0;

  //Now, rewrite translated coordinates and combine them with coarse grained water!
  for(i=0; i<atomCount1[0]; i++)
    {
      fprintf(OUT_FILE, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", info1[i].serial, info1[i].name, \
	      info1[i].resName, info1[i].chainID, info1[i].resSeq, info1[i].x, info1[i].y, info1[i].z, info1[i].occupancy, info1[i].tempFactor);
    }
  fprintf(OUT_FILE, "TER\n");
  //Use the number of waters in coarse grained model!

  for(i=atomCount2[5]; i<atomCount2[0]; i++)
    {
      prot2atoms[i].x=(prot2atoms[i].x-translationVector[0]);
      prot2atoms[i].y=(prot2atoms[i].y-translationVector[1]);
      prot2atoms[i].z=(prot2atoms[i].z-translationVector[2]);

      //All number tricks are for the sake of counting according to pdb format. 
      fprintf(OUT_FILE, "ATOM  %5d %s %3.3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", (atomCount1[0]+i-atomCount2[5]+1), " OH2", \
	      "TIP3",  'W', (i+1), prot2atoms[i].x, prot2atoms[i].y, prot2atoms[i].z, 1.00, 0.00); 
    }
  fprintf(OUT_FILE, "END\n");

  fclose(OUT_FILE);
  free(info1);
  free(prot1atoms);
  free(info2);
  free(prot2atoms);
}

CAcoord *readPdbReturnCA(char *fileName, int *protCounts)
{
  FILE *pdbData=fopen(fileName,"r");
  if(pdbData==NULL)
    {
      fprintf(stderr, "No such file:%s\n", fileName);
      exit(EXIT_FAILURE);
    }
  fprintf(stdout, "%s\n", fileName);
  scanpdb(protCounts, fileName);
  printf("Number of residues:%d\n",  protCounts[1]);
  //--Get data from first pdb file--
  pdb_v23 *info = (pdb_v23*)malloc((protCounts[0])*sizeof(pdb_v23));
  if (info == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb\n");
      exit(EXIT_FAILURE);
    }
  CAcoord *protCA = (CAcoord*)malloc((protCounts[1])*sizeof(CAcoord));
  if (protCA == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord\n");
      exit(EXIT_FAILURE);
    }
  readpdb_v23(info, protCA, fileName);
  fclose(pdbData);
  free(info);
  return protCA;
}

CAcoord *readPdbReturnCAandOH2(char *fileName, int *protCounts)
{
  FILE *pdbData=fopen(fileName,"r");
  if(pdbData==NULL)
    {
      fprintf(stderr, "No such file:%s\n", fileName);
      exit(EXIT_FAILURE);
    }
  fprintf(stdout, "%s\n", fileName);
  scanpdb_CAandOH2(protCounts, fileName);
  printf("Number of residues:%d\n",  protCounts[1]);
  //--Get data from first pdb file--
  pdb_v23 *info = (pdb_v23*)malloc((protCounts[0])*sizeof(pdb_v23));
  if (info == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct pdb\n");
      exit(EXIT_FAILURE);
    }
  CAcoord *protCA = (CAcoord*)malloc((protCounts[1])*sizeof(CAcoord));
  if (protCA == NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for struct CAcoord\n");
      exit(EXIT_FAILURE);
    }
  
  readpdb_CAandOH2(info, protCA, fileName);
  fclose(pdbData);
  free(info);
  return protCA;
}
int  *scanpdb_trajectory(int *atomCount, int *resCount, int *frameCount, char *inFileName, int printDetails)
{
/*Purpose: To determine the total number of atoms and residues and frames in a pdb trajectory file for memory allocation*/
//Date added: November 16, 2014
  int i=0;
  int numberofAtoms=0;
  int numberofResids=0;
  int numberofFrames=0;                            
  char line[100];
  
  FILE *pdbdata=fopen(inFileName,"r");
  if(pdbdata==NULL)
    {
      fprintf(stderr, "Can not open %s file\n", inFileName);
      exit(EXIT_FAILURE);
    }
  while(1)
    {
      if(fgets(line, sizeof(line),pdbdata)==NULL) break;
      
      if((strncmp("ATOM  ", line, 6))==0)
	{
	  //Only CA atoms with null or 'A' alternate location indicator are selected. 
	  if(((strncmp(" CA", line+12, 3))==0 ) && ((line[16]==' ') || (line[16]=='A')  ) )
	    numberofResids+=1;
	  
	  numberofAtoms+=1;
	  i+=1; 
	}
      
      if(((strncmp("ENDMDL", line, 6))==0) || ((strncmp("END", line, 3))==0))
	{
	  numberofFrames+=1;
	}
      
    }
  atomCount[0]=numberofAtoms;   //atomCount[0] is the total number of atoms
  resCount[0]=numberofResids;   //resCount[0] is the total number of residues
  frameCount[0]=numberofFrames; //frameCount[0] is the total number of helices

  if((atomCount[0]==0) || (resCount[0]==0))
    {
      fprintf(stderr, "Error: Can not read atoms and residues\n");
      exit(EXIT_FAILURE);
    } 
  if(printDetails)
    {
      fprintf(stdout, "Total number of atoms:%d\n", atomCount[0]);
      fprintf(stdout, "Total number of residues:%d\n",   resCount[0]);
      fprintf(stdout, "Total number of frames:%d\n",      frameCount[0]);
    }
  fclose(pdbdata);
  
  return atomCount; 
}

void  readpdb_trajectory(pdb_v23 *info, CAcoord *atom, FILE *PDBDATA, int *atomCount, int printDetails)
{
  //This function read a frame from a multiframe pdb file and returns coordinates into info and atom data structures. 
  //Date: November 17, 2014

  int i=0;
  int numberofResid=0;
  
  char c_buffer[9]; //Stands for (c)haracter buffer!!!
  memset(c_buffer,'\0', 9);

  char line[100];
  memset(line,'\0', 100);

  int numCA=0;
  while(1)
    {
      if(fgets(line, sizeof(line),PDBDATA)==NULL) 
	{
	  fprintf(stderr, "ERROR: This file does not contain any data!\n");
	  break;
	}
      else if((strncmp("ATOM  ", line, 6))==0) 
	{
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+6, 5); info[i].serial=atoi(c_buffer);
	memset(info[i].name,'\0', 5)   ; strncpy(info[i].name,      line+12,  4);
	info[i].altLoc=line[16];
	memset(info[i].resName,'\0', 4); strncpy(info[i].resName,   line+17,  3);
	info[i].chainID=line[21];
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+22,  4); info[i].resSeq=atoi(c_buffer);
	info[i].iCode=line[26];
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+30,  8); info[i].x=atof(c_buffer);
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+38,  8); info[i].y=atof(c_buffer);
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+46,  8); info[i].z=atof(c_buffer);
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+54,  6); info[i].occupancy=atof(c_buffer);
	memset(c_buffer,'\0', 9)       ; strncpy(c_buffer,         line+60,  6); info[i].tempFactor=atof(c_buffer);
	memset(info[i].element,'\0', 2); strncpy(info[i].element,   line+76,  2);
	memset(info[i].charge,'\0', 2) ; strncpy(info[i].charge,    line+78,  2);
	
	if(printDetails)	
	  fprintf(stdout,"ATOM   %d %s %c %s %c %d %c %lf %lf %lf %lf %lf %s %s\n", \
		  info[i].serial, info[i].name, info[i].altLoc, info[i].resName, info[i].chainID, info[i].resSeq, info[i].iCode, 
		  info[i].x, info[i].y, info[i].z, info[i].occupancy, info[i].tempFactor, info[i].element, info[i].charge);
	
	//If there is an alternative location, I select the first one. This selection may not be healthy depending on the job one is working on!!	  
	if(((strncmp(" CA", line+12, 3))==0 ) && ((info[i].altLoc==' ') || (info[i].altLoc=='A')  ) )
	  {
	    strncpy(atom[numberofResid].residname, info[i].resName, 3);
	    atom[numberofResid].chain=info[i].chainID;
	    atom[numberofResid].resNo=info[i].resSeq;
	    atom[numberofResid].x=info[i].x;
	    atom[numberofResid].y=info[i].y;
	    atom[numberofResid].z=info[i].z;
	    numCA+=1;
	    numberofResid+=1;
	  }
	i+=1; 
	}
      else if(  ((strncmp("ENDMDL", line, 6))==0) || ((strncmp("END", line, 3))==0)  )
	break;
    }

  if(numberofResid==0)
    {
      fprintf(stderr, "Error: Can not read residues\n");
      exit(EXIT_FAILURE);
    }
  
  atomCount[0]=i;
  atomCount[1]=numberofResid;

  if(printDetails)
    printf("Number of CA=%d\n", numCA);
}
