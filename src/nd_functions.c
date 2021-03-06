//Purpose: nd_functions.c contains a few functions to calculate normalized distance (ND).
//         ND is used to analyze conformational transitions between
//         two protein conformations. It gives percentage of completing action per residue. 
//         ND for 'one residue' can be formulated as below:
//         ND_{i}=(d_{b,i}) / (d_{b,i}+d_{i,e})
//         where b stands for beginning structure, e stands for end structure and i stands for intermediate structure.
//         d denotes distance. ND is calculated for all residues and one can make a color plot to details of conformational 
//         transition for all amino acids. 
//         All functions related to elastic network model are in this file. 

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
#include <math.h>
//Personal header files!
#include <structures.h>

void printNormalizedDistanceEachAtom(int N/*System Size*/, CAcoord *atom_beg, CAcoord *atom_cur, CAcoord *atom_end,
				     int frameNo, int printDetails, char *outFileName)
{
  double distanceForAtom[N];
  int i=0;

  FILE *RESULTS_FILE=fopen(outFileName, "a");
  if(RESULTS_FILE==NULL)
    {
      fprintf(stderr, "ERROR: %s file can not be opened!\n", outFileName);
      exit(EXIT_FAILURE);
    }


  //Zero all elements!
  for(i=0; i<N; i++)
    {
      distanceForAtom[i]=0.0;
    }

  double disMax=0.0;
  for(i=0; i<N; i++)
    {
      disMax=0.00;
      distanceForAtom[i]=sqrt(  ((atom_beg[i].x - atom_cur[i].x)*(atom_beg[i].x - atom_cur[i].x))    + \
			        ((atom_beg[i].y - atom_cur[i].y)*(atom_beg[i].y - atom_cur[i].y))    + \
				((atom_beg[i].z - atom_cur[i].z)*(atom_beg[i].z - atom_cur[i].z)) );
      
      disMax            =sqrt(  ((atom_beg[i].x - atom_end[i].x)*(atom_beg[i].x - atom_end[i].x))    + \
				((atom_beg[i].y - atom_end[i].y)*(atom_beg[i].y - atom_end[i].y))    + \
				((atom_beg[i].z - atom_end[i].z)*(atom_beg[i].z - atom_end[i].z)));

      distanceForAtom[i]=(distanceForAtom[i]/disMax);
    }


  for(i=0; i<N; i++)
    {
      fprintf(RESULTS_FILE, "%d\t%d\t%lf\n", frameNo, i, (distanceForAtom[i]));
    }
  fprintf(RESULTS_FILE, "\n");
  
  fclose(RESULTS_FILE);
}

void printNormalizedDistanceEachAtom_v2(int N/*System Size*/, CAcoord *atom_beg, CAcoord *atom_cur, CAcoord *atom_end,
					int frameNo, int printDetails, char *outFileName)
{
  //Whats new in v2: It assigns the percentages to beta factors in trajectory file!
  int i=0;

  FILE *RESULTS_FILE=fopen(outFileName, "a");
  if(RESULTS_FILE==NULL)
    {
      fprintf(stderr, "ERROR: %s file can not be opened!\n", outFileName);
      exit(EXIT_FAILURE);
    }

  double diff_beg_curX=0.0, diff_beg_curY=0.0, diff_beg_curZ=0.0;
  double diff_cur_endX=0.0, diff_cur_endY=0.0, diff_cur_endZ=0.0;
  double temp1=0.0, temp2=0.0;
  double result=0.0;
  for(i=0; i<N; i++)
    {
      diff_beg_curX=0.0, diff_beg_curY=0.0, diff_beg_curZ=0.0;
      diff_cur_endX=0.0, diff_cur_endY=0.0, diff_cur_endZ=0.0;

      result=0.0;
      diff_beg_curX  =(atom_beg[i].x - atom_cur[i].x);
      diff_beg_curY  =(atom_beg[i].y - atom_cur[i].y);
      diff_beg_curZ  =(atom_beg[i].z - atom_cur[i].z);

      diff_cur_endX  =(atom_cur[i].x - atom_end[i].x);
      diff_cur_endY  =(atom_cur[i].y - atom_end[i].y);
      diff_cur_endZ  =(atom_cur[i].z - atom_end[i].z);

      temp1=0.0, temp2=0.0;

      temp1=((diff_beg_curX*diff_beg_curX)+(diff_beg_curY*diff_beg_curY)+(diff_beg_curZ*diff_beg_curZ));
      
      temp2=((diff_cur_endX*diff_cur_endX)+(diff_cur_endY*diff_cur_endY)+(diff_cur_endZ*diff_cur_endZ)); 
      result=( temp1/(temp1+temp2));
      atom_cur[i].beta=(result);
      fprintf(RESULTS_FILE, "%d\t%d\t%lf\n", frameNo, atom_beg[i].resNo, fabs(sqrt(result)));

    }
  fprintf(RESULTS_FILE, "\n");
  fclose(RESULTS_FILE);
}

void printNormalizedDistanceEachAtom_v3(int N/*System Size*/, CAcoord *atom_beg, CAcoord *atom_cur, CAcoord *atom_end, int frameNo, char *outFileName, int range_beg, int range_end, double coeff, int printDetails)
{
  //Whats new in v3: Two features have been added to this version:
  //                 i) I added a range of residues to ND plot. This feature is useful if you want to compare movements of two domains!
  //                ii) Sometimes, it is useful to use time instead of frame number. Therefore, I added a coefficient to multiply with frame number!
  int i=0;
  if(printDetails)
    {
      fprintf(stdout, "Scaling coefficient for frames: %lf\n", coeff);
    }

  FILE *RESULTS_FILE=fopen(outFileName, "a");
  if(RESULTS_FILE==NULL)
    {
      fprintf(stderr, "ERROR: %s file can not be opened!\n", outFileName);
      exit(EXIT_FAILURE);
    }

  double diff_beg_curX=0.0, diff_beg_curY=0.0, diff_beg_curZ=0.0;
  double diff_cur_endX=0.0, diff_cur_endY=0.0, diff_cur_endZ=0.0;
  double temp1=0.0, temp2=0.0;
  double result=0.0;
  for(i=0; i<N; i++)
    {
      if((atom_beg[i].resNo>=range_beg) && (atom_beg[i].resNo<=range_end))
	{

	  diff_beg_curX=0.0, diff_beg_curY=0.0, diff_beg_curZ=0.0;
	  diff_cur_endX=0.0, diff_cur_endY=0.0, diff_cur_endZ=0.0;
	  
	  result=0.0;
	  diff_beg_curX  =(atom_beg[i].x - atom_cur[i].x);
	  diff_beg_curY  =(atom_beg[i].y - atom_cur[i].y);
	  diff_beg_curZ  =(atom_beg[i].z - atom_cur[i].z);
	  
	  diff_cur_endX  =(atom_cur[i].x - atom_end[i].x);
	  diff_cur_endY  =(atom_cur[i].y - atom_end[i].y);
	  diff_cur_endZ  =(atom_cur[i].z - atom_end[i].z);
	  
	  temp1=0.0, temp2=0.0;
	  
	  temp1=((diff_beg_curX*diff_beg_curX)+(diff_beg_curY*diff_beg_curY)+(diff_beg_curZ*diff_beg_curZ));
	  
	  temp2=((diff_cur_endX*diff_cur_endX)+(diff_cur_endY*diff_cur_endY)+(diff_cur_endZ*diff_cur_endZ));
	  result=( temp1/(temp1+temp2));
	  atom_cur[i].beta=(result);
	  
	  fprintf(RESULTS_FILE, "%.2lf\t%d\t%lf\n", (coeff*frameNo), atom_beg[i].resNo, fabs(sqrt(result)));
	}
    }
  fprintf(RESULTS_FILE, "\n");
  fclose(RESULTS_FILE);
}
