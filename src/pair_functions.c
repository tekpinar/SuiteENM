//Purpose: pair_functions.c contains functions related to pair distribution functions
//         in our SAXS-based flexible fitting project.

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
#include <assert.h>
#include <math.h>
//Personal header files
#include <defines.h>
#include <structures.h>
#include <distance_functions.h>

#define _GNU_SOURCE


void pairDistanceDistro(int N, int *pairNo)
{
  //Purpose: This 
  int N_pair=pairNo[0];
  int i=0, j=0;          // Regular counters
  double d_I=1.0;        // Data interval parameter
  double d_ij=2.0;       // Test case, you will get them from pdb file
  double P_dI=0.0;       // P(d_I) pair distance distrubion function

  FILE *pairdata=fopen("pairdata.txt", "w");
  double constant=sqrt(2*PI)*SIGMA*N_pair;
  for(i=0; i<100; i++)
    {
      for(j=0; j<N_pair; j++)
	{
	  P_dI+=(1/constant)*exp(-((d_ij-d_I)*(d_ij-d_I))/(2*SIGMA*SIGMA));
	}
      fprintf(pairdata, "%lf\t%lf\n", P_dI, d_I);
      d_I+=1.0;
    }

  fclose(pairdata);
}
void pairCorrPlot1(int N, int *pairNo, pairInfoNew* pair)
{
  //Purpose: This fucntion produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional for is delta function form.
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  int i=0; //Regular counters
  int pairInBin=0;

  double binIndex=0.25;
  double maxDistance=0.0;

  FILE *pairdata=fopen("pairCorrPlot1.txt", "w");
  for(i=0; i<pairNo[0]; i++)
    {

      if(pair[i].dis>=maxDistance)
	maxDistance=pair[i].dis;
    }

  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //I will put my data into 1/4 Angstrom bins
  //Count pairs whose distances are between than 1 A
  while(binIndex<(maxDistance + 0.25))
    {
      for(i=0; i<pairNo[0]; i++)
	{
	  if( (pair[i].dis>=(binIndex-0.25)) && (pair[i].dis<binIndex) )
	    {
	      pairInBin+=1;
	    }
	}
      double normPairInBin=(double) pairInBin/pairNo[0];
      fprintf(pairdata, "%lf\t%lf\n", binIndex, normPairInBin);    
      pairInBin=0;
      binIndex+=0.25;
    }

  fclose(pairdata);
}

void pairCorrPlot2(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair)
{
  //Purpose: This function produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional form is gaussian so that we can take derivatives analytically. 
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  int i=0; //Regular counters

  double disMobile_ij=0.0;
  double binIndex=1.0;
  double maxDistance=0.0;
  FILE *pairdata=fopen("pairCorrPlot2.txt", "w");
  for(i=0; i<pairNo[0]; i++)
    {
      if(pair[i].dis>=maxDistance)
	maxDistance=pair[i].dis;
    }
  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //I will put my data into 1 Angstrom bins
  //Count pairs whose distances are between that 1 A
  while(binIndex<(maxDistance + 1.0))
    {
      double P_dI=0.0;
      for(i=0; i<pairNo[0]; i++)
	{
	  disMobile_ij=dis(atom, (pair[i].posI), (pair[i].posJ));
	  double subt=(disMobile_ij-binIndex);
	  if(subt<3*SIGMA)
	    P_dI+=(1/(pairNo[0]*2.506628274631*SIGMA))*exp(-subt*subt/(2*SIGMA_SQRD));
	}
      
      fprintf(pairdata, "%lf\t%lf\n", binIndex, P_dI);    
      binIndex+=1.0;
    }
  
  fclose(pairdata);
}
void pairCorrPlot3(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair)
{
  //Purpose: This function produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional form is error function so that we can take derivatives analytically. 
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  int i=0,  l=0; //Regular counters

  double disMobile_ij=0.0;
  double binSize=1.0;
  double binIndex=1.0;
  double maxDistance=0.0;
  FILE *pairdata=fopen("pairCorrPlot_Erf.txt", "w");
  for(i=0; i<pairNo[0]; i++)
    {
      if(pair[i].dis>=maxDistance)
	maxDistance=pair[i].dis;
    }
  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //I will put my data into 1 Angstrom bins
  //Count pairs whose distances are between that 1 A
  double coeff1=(1/(pairNo[0]*2*binSize));
  //  double coeff2=1.41421356*SIGMA;
  double coeff2=2.8284271212474619;
  double halfBin=binSize/2;
  while(binIndex<=(ceil(maxDistance)))
    {
      double P_dI=0.0;
      for(l=0; l<pairNo[0]; l++)
	{
	  disMobile_ij=pair[l].dis;
	  double subt=(binIndex-disMobile_ij);
	  //Now, lets put the distance restriction off for now.
	  //      if(subt<3*SIGMA)
	  P_dI+=(   ( erf (  (subt+halfBin)/(coeff2 )  )   )   -  ( erf (  (subt-halfBin)/(coeff2)  )   )   );
	}
      
      fprintf(pairdata, "%lf\t%.15lf\n", binIndex, P_dI*coeff1);    
      binIndex+=1.0;
    }
  
  fclose(pairdata);
}
void pairCorrPlot4(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair)
{
  //Purpose: This function produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional form is error function. In this case, I am making an approximation to avoid some calculations!!
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  int i=0,  l=0; //Regular counters

  double disMobile_ij=0.0;
  double binSize=1.0;
  double binIndex=1.0;
  double maxDistance=0.0;
  FILE *pairdata=fopen("pairCorrPlot_Erf4.txt", "w");
  for(i=0; i<pairNo[0]; i++)
    {
      if(pair[i].dis>=maxDistance)
	maxDistance=pair[i].dis;
    }
  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //I will put my data into 1 Angstrom bins
  //Count pairs whose distances are between that 1 A
  double coeff1=(1/(pairNo[0]*2*binSize));
  //  double coeff2=1.41421356*SIGMA;
  double coeff2=2.8284271212474619;
  double coeff3=0.0;
  double coeff4=0.0;
  double halfBin=binSize/2;
  while(binIndex<=(ceil(maxDistance)))
    {
      double P_dI=0.0;
      for(l=0; l<pairNo[0]; l++)
	{
	  disMobile_ij=pair[l].dis;
	  double subt=(binIndex-disMobile_ij);
	  //Now, lets put the distance restriction off for now.
	  //==============
      coeff3=(subt+halfBin)/coeff2;
      coeff4=(subt-halfBin)/coeff2;
      //Now, lets put the distance restriction off for now.
      //      if(subt<3*SIGMA)
      //Error function is almost 1 (or -1) after 3.5(-3.5)!!!!!!!
      if((fabs(coeff3)>(6.0))&&(fabs(coeff4)>(6.0)))
	{
        if((coeff3>0)&&(coeff4<0))
	  P_dI+= 2.0;
        if((coeff3<0)&&(coeff4>0))
	  P_dI+=-2.0;
	/* 	else */
	/* 	  { */
	/* 	    //	if((coeff3>0)&&(coeff4>0)) */
	/* 	  P_dI+= 0.0; */
	/*         if((coeff3<0)&&(coeff4<0)) */
	/* 	  P_dI+= 0.0; */
	//	  }
	}
      else
	P_dI+=(   ( erf ( coeff3  )   )   -  ( erf ( coeff4  )   )   );
      //=================
	}
      
      fprintf(pairdata, "%lf\t%.15lf\n", binIndex, P_dI*coeff1);    
      binIndex+=1.0;
    }
  
  fclose(pairdata);
}

double pairDistroGaussian1(int N, CAcoord *atomMoving, int *pairNo,  pairInfoNew *pair,  double binIndex/*Initially use 1.0*/)
{
  //Purpose: This function produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional form is gaussian so that we can take derivatives analytically. 
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  int i=0; //Regular counters
  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //  double sigmaSqrd=sigma*sigma;
  double disMobile_ij=0.0;
  double P_dI=0.0;
  double coeff=(1/(pairNo[0]*2.506628274631*SIGMA));
  for(i=0; i<pairNo[0]; i++)
    {
      disMobile_ij=dis(atomMoving, (pair[i].posI), (pair[i].posJ));
      double subt=(disMobile_ij-binIndex);
      //Now, lets put the distance restriction on.
      if(subt<3*SIGMA)
	P_dI+=coeff*exp(-subt*subt/(2*SIGMA_SQRD));
    }
  //  fclose(pairdata);
  return P_dI;
}
double pairDistroGaussian2(int N, CAcoord *atomMoving, int *pairNo,  pairInfoNew *pair, double binIndex/*Initially use 1.0*/, double increment,  int atomNo, int component)
{
  //Purpose: This function produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional form is gaussian so that we can take derivatives analytically. 
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  //This is the most inefficient way of doing it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int i=0, j=0; //Regular counters

  double disMobile_ij=0.0;
  double dis_x=0.0;
  double dis_y=0.0;
  double dis_z=0.0;

  //  fprintf(stdout, "Once: AtomMoving[%d].y=%lf\n", atomNo, atomMoving[atomNo].y);
  if(component==0)
    atomMoving[atomNo].x=atomMoving[atomNo].x+increment;
  else  if(component==1)
    atomMoving[atomNo].y=atomMoving[atomNo].y+increment;
  else  if(component==2)
    atomMoving[atomNo].z=atomMoving[atomNo].z+increment;

  //  fprintf(stdout, "Sonra: AtomMoving[%d].y=%lf\n", atomNo, atomMoving[atomNo].y);
  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  double P_dI=0.0;
  for(i=0; i<N; i++)
    for(j=0; j<i; j++)
    {
      dis_x=atomMoving[i].x-atomMoving[j].x;
      dis_y=atomMoving[i].y-atomMoving[j].y;
      dis_z=atomMoving[i].z-atomMoving[j].z;
      //      disMobile_ij=dis(atomMoving, i, j);
      disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
      double subt=(disMobile_ij-binIndex);
      //Now, lets put the distance restriction on.
      //      if(subt<3*sigma)
      P_dI+=(2/(N*(N-1)*2.506628274631*SIGMA))*exp(-subt*subt/(2*SIGMA_SQRD));
    }
  //  fclose(pairdata);

  if(component==0)
    atomMoving[atomNo].x=atomMoving[atomNo].x-increment;
  else  if(component==1)
    atomMoving[atomNo].y=atomMoving[atomNo].y-increment;
  else  if(component==2)
    atomMoving[atomNo].z=atomMoving[atomNo].z-increment;

  return P_dI;
}
double pairDistroGaussian3(int N, CAcoord *atomMoving, 
			   int *pairNo,  pairInfoNew *pair, 
			   int *pairInitialNo,  pairInfoNew *pairInitial, 
			   double binIndex/*Initially use 1.0*/)
{
  //Purpose: This function produces data of pair distribution function for a coarse grained model using all pairs.
  //         The functional form is gaussian so that we can take derivatives analytically. 
  //How many pairs do I have?
  //If there is not any restriction: N*(N-1)/2 
  //where N is my CA atom number for m simple model.
  int i=0; //Regular counters
  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //  double sigmaSqrd=sigma*sigma;
  double disMobile_ij=0.0;
  double P_dI=0.0;
  double coeff=(1/(pairNo[0]*2.506628274631*SIGMA));
  for(i=0; i<pairNo[0]; i++)
    {
      disMobile_ij=pairInitial[i].dis;
      //      disMobile_ij=dis(atomMoving, (pair[i].posI), (pair[i].posJ));
      double subt=(disMobile_ij-binIndex);
      //Now, lets put the distance restriction on.
      if(subt<3*SIGMA)
	P_dI+=coeff*exp(-subt*subt/(2*SIGMA_SQRD));
    }
  //  fclose(pairdata);
  return P_dI;
}
double pairDistroGaussian3_I(int N, CAcoord *atomMoving, 
			     int *pairNo,  pairInfoNew *pair, 
			     int *pairInitialNo,  pairInfoNew *pairInitial, 
			     double binIndex, double binSize)
{
  //Purpose: This function produces data of pair distribution function  for a coarse grained model using all pairs.
  //         The functional form is error function.
  //In equation sheet: binIndex=d_I, binSize= u

  //If there is not any restriction: N*(N-1)/2 pairs exist
  //where N is my CA atom number for m simple model.

  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData_I.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //  double sigmaSqrd=sigma*sigma;
  int l=0; //Regular counters
  double disMobile_ij=0.0;
  double P_dI=0.0;
  double coeff1=(1/(pairNo[0]*2*binSize));
  //  double coeff2=1.41421356*SIGMA;
  double coeff2=2.8284271212474619;
  double halfBin=binSize/2;
  for(l=0; l<pairNo[0]; l++)
    {
      disMobile_ij=pairInitial[l].dis;
      double subt=(binIndex-disMobile_ij);
      //Now, lets put the distance restriction off for now.
      //      if(subt<3*SIGMA)
	P_dI+=(   ( erf (  (subt+halfBin)/(coeff2 )  )   )   -  ( erf (  (subt-halfBin)/(coeff2)  )   )   );
    }
  //  fclose(pairdata);
  //  return P_dI;
  return P_dI*coeff1;
}

double pairDistroGaussian3_I_water(int N, CAcoord *atomMoving, 
				   int *pairNo,  pairInfoNew *pair, 
				   int *pairInitialNo,  pairInfoNew *pairInitial, 
				   double binIndex, double binSize, double N_pairSum)
{
  //Purpose: This function produces data of pair distribution function  for a coarse grained model using all pairs.
  //         The functional form is error function.
  //In equation sheet: binIndex=d_I, binSize= u

  //If there is not any restriction: N*(N-1)/2 pairs exist
  //where N is my CA atom number for m simple model.

  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData_I.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //  double sigmaSqrd=sigma*sigma;
  int i=0, j=0;
  //  int l=0; //Regular counters
  double disMobile_ij=0.0;
  double P_dI=0.0;
  //double coeff1=(1/(pairNo[0]*2*binSize));
  double coeff1=(1/(N_pairSum*2*binSize));
  //  printf("N_pairSum=%lf\n", N_pairSum);
 //  double coeff2=1.41421356*SIGMA;
  double coeff2=2.8284271212474619;
  double halfBin=binSize/2;
  double b1=1.0;  
  double subt=0.0;
  for(i=0; i<N; i++)
    {
      for(j=0; j<i; j++)
	{
	  /* 	  if( ( (strncmp(atomMoving[i].residname, "TIP", 3)!=0) || (strncmp(atomMoving[i].residname, "HOH", 3)!=0)  ) \ */
	  /* 	  &&  ( (strncmp(atomMoving[j].residname, "TIP", 3)!=0) || (strncmp(atomMoving[j].residname, "HOH", 3)!=0)  ) ) */
	  if( ( (strncmp(atomMoving[i].residname, "TIP", 3)!=0)  )  &&  ( (strncmp(atomMoving[j].residname, "TIP", 3)!=0)   ) )

	    {
	      b1=1.0;
	      if(0)	
		fprintf(stdout, "I found %s and %s\n", atomMoving[i].residname, atomMoving[j].residname);		  
	      
	      disMobile_ij=dis(atomMoving, i, j);
	      subt=(disMobile_ij-binIndex);
	      
	      P_dI+=b1*(  ( erf (  (subt+halfBin)/(coeff2 )  )  )  -  ( erf (  (subt-halfBin)/(coeff2)  )  )  );
	      
	    }
	  else if( (  strncmp(atomMoving[i].residname, "TIP", 3)==0 )    &&  (strncmp(atomMoving[j].residname, "TIP", 3)==0)       )
	    {
	      b1=0.49;
	      if(0)	
		fprintf(stdout, "I found %s and %s\n", atomMoving[i].residname, atomMoving[j].residname);		  
	      
	      disMobile_ij=dis(atomMoving, i, j);
	      subt=(disMobile_ij-binIndex);
	      
	      P_dI+=b1*(  ( erf (  (subt+halfBin)/(coeff2 )  )  )  -  ( erf (  (subt-halfBin)/(coeff2)  )  )  );
	      P_dI+=0.0;
	      
	    }
	  else
	    {
	      b1=0.7;
	      if(0)	
		fprintf(stdout, "I found %s and %s\n", atomMoving[i].residname, atomMoving[j].residname);		  
	      
	      disMobile_ij=dis(atomMoving, i, j);
	      subt=(disMobile_ij-binIndex);
	      
	      P_dI+=b1*(  ( erf (  (subt+halfBin)/(coeff2 )  )  )  -  ( erf (  (subt-halfBin)/(coeff2)  )  )  );
	    }
	  
	}
      
    }

  /*   for(l=0; l<pairNo[0]; l++) */
  /*     { */
  /*       disMobile_ij=pairInitial[l].dis; */
  /*       double subt=(binIndex-disMobile_ij); */
  /*       //Now, lets put the distance restriction off for now. */
  /*       //      if(subt<3*SIGMA) */
  /* 	P_dI+=(   ( erf (  (subt+halfBin)/(coeff2 )  )   )   -  ( erf (  (subt-halfBin)/(coeff2)  )   )   ); */
  /*     } */
  //  fclose(pairdata);
  //  return P_dI;
  return P_dI*coeff1;
}


//--New Pair Distribution function
double pairDistroGaussian4_I(int N, CAcoord *atomMoving, 
			     int *pairNo,  pairInfoNew *pair, 
			     int *pairInitialNo,  pairInfoNew *pairInitial, 
			     double binIndex, double binSize)
{
  //Purpose: This function produces data of pair distribution function  for a coarse grained model using all pairs.
  //         The functional form is error function.
  //In equation sheet: binIndex=d_I, binSize= u

  //If there is not any restriction: N*(N-1)/2 pairs exist
  //where N is my CA atom number for m simple model.

  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData_I.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);    
  //  double sigmaSqrd=sigma*sigma;
  int l=0; //Regular counters
  double disMobile_ij=0.0;
  double P_dI=0.0;
  double coeff1=(1/(pairNo[0]*2*binSize));
  //double coeff2=sqrt(2.0)*SIGMA;
  double coeff2=2.8284271212474619;
  double coeff3=0.0;
  double coeff4=0.0;
  double halfBin=binSize/2;
  //  double halfBinByCoeff2=halfBin/coeff2;
  for(l=0; l<pairNo[0]; l++)
    {
      disMobile_ij=pairInitial[l].dis;
      double subt=(binIndex-disMobile_ij);
      coeff3=(subt+halfBin)/coeff2;
      coeff4=(subt-halfBin)/coeff2;
      //Now, lets put the distance restriction off for now.
      //      if(subt<3*SIGMA)
      //Error function is almost 1 (or -1) after 3.5(-3.5)!!!!!!!
      if((fabs(coeff3)>(6.0))&&(fabs(coeff4)>(6.0)))
	{
        if((coeff3>0)&&(coeff4<0))
	  P_dI+= 2.0;
        if((coeff3<0)&&(coeff4>0))
	  P_dI+=-2.0;
	/* 	else */
	/* 	  { */
	/* 	    //	if((coeff3>0)&&(coeff4>0)) */
	/* 	  P_dI+= 0.0; */
	/*         if((coeff3<0)&&(coeff4<0)) */
	/* 	  P_dI+= 0.0; */
	//	  }

	}
      else
	P_dI+=(   ( erf ( coeff3  )   )   -  ( erf ( coeff4  )   )   );
    }
  //  fclose(pairdata);
  //  return P_dI;
  return P_dI*coeff1;
}

//---------------------------------
double pairDistroGaussian5_I(int N, CAcoord *atomMoving,
			     int *pairNo,  pairInfoNew *pair,
			     int *pairInitialNo,  pairInfoNew *pairInitial,
			     double binIndex, double binSize)
{
  //Purpose: This function produces data of pair distribution function  for a coarse grained model using all pairs.
  //         The functional form is error function lookup table.
  //In equation sheet: binIndex=d_I, binSize= u

  //If there is not any restriction: N*(N-1)/2 pairs exist
  //where N is my CA atom number for m simple model.

  //  double maxDistance=0.0;
  //  FILE *pairdata=fopen("pairFuncData_I.txt", "w");
  /*   for(i=0; i<pairNo[0]; i++) */
  /*     { */
  /*       if(pair[i].dis>=maxDistance) */
  /* 	maxDistance=pair[i].dis; */
  /*     } */
  //  fprintf(stdout, "Max distance is %lf\n", maxDistance);
  //  double sigmaSqrd=sigma*sigma;
  int l=0; //Regular counters
  double disMobile_ij=0.0;
  double P_dI=0.0;
  double coeff1=(1/(pairNo[0]*2*binSize));
  //double coeff2=sqrt(2.0)*SIGMA;
  double coeff2=2.8284271212474619;
  double coeff3=0.0;
  double coeff4=0.0;
  int    coeff3Int=0;
  int    coeff4Int=0;

  double halfBin=binSize/2;
  //  double halfBinByCoeff2=halfBin/coeff2;
  for(l=0; l<pairNo[0]; l++)
    {
      disMobile_ij=pairInitial[l].dis;
      double subt=(binIndex-disMobile_ij);
      coeff3=(subt+halfBin)/coeff2;
      coeff4=(subt-halfBin)/coeff2;
      //Now, lets put the distance restriction off for now.
      //      if(subt<3*SIGMA)
      //Error function is 1 (or -1) after 3.5(-3.5)!!!!!!!
      if(  (fabs(coeff3)>(5.0))  &&  (fabs(coeff4)>(5.0))  )
	{
        if((coeff3>0)&&(coeff4<0))
	  P_dI+= 2.0;
        if((coeff3<0)&&(coeff4>0))
	  P_dI+=-2.0;
	}
      else if(  (fabs(coeff3)<(5.0))  &&  (fabs(coeff4)<(5.0))  )
	{
	  coeff3Int=(int)(coeff3*100000);
	  coeff4Int=(int)(coeff4*100000);
	  if((coeff3>0)&&(coeff4>0))
	    {
	      //	   printf("Nerde patladi?l=%d\tcoeff3=%lf\t\tcoeff4=%lf\n", l, coeff3, coeff4);
	      //     fprintf(stdout, "erf(%lf)=%lf\terfTable[(%d)]=%lf\n", coeff3, erf(coeff3), (int)(coeff3*100000), erfTable[(int)(coeff3*100000)]);
	      //     P_dI+=(   ( erfTable [ (int)(round(coeff3*10000))]   )   -  ( erfTable [ (int)(round(coeff4*10000)) ]   )   );
	      P_dI+=(   ( erfTable [coeff3Int]   )   -  ( erfTable [coeff4Int ]   )   );
	    }

	  if((coeff3<0)&&(coeff4<0))
	    {
	      //	printf("Nerde patladi?l=%d\tcoeff3=%lf\tcoeff4=%lf\n", l, coeff3, coeff4);
	      //       P_dI+=(   ( -erfTable [ (int)(round((-1.0)*coeff3*10000))]   )   -  ( -erfTable [ (int)(round((-1.0)*coeff4*10000)) ]   )   );
	      P_dI+=(   ( -erfTable [ -coeff3Int ]   )   -  ( -erfTable [ -coeff4Int ]   )   );
	    }

	  if((coeff3>0)&&(coeff4<0))
	    {
	      //	      printf("Nerde patladi?l=%d\tcoeff3=%lf\tcoeff4=%lf\n", l, coeff3, coeff4);

	      //	      P_dI+=(   ( erfTable [ (int)(round(coeff3*10000))]   )   -  ( -erfTable [ (int)(round((-1.0)*coeff4*10000)) ]   )   );
	      P_dI+=(   ( erfTable [ coeff3Int ]   )   -  ( -erfTable [ -coeff4Int ]   )   );
	    }
	  
	  if((coeff3<0)&&(coeff4>0))
	    {
	      //	      printf("Nerde patladi?l=%d\tcoeff3=%lf\tcoeff4=%lf\n", l, coeff3, coeff4);
	      //	      P_dI+=(   ( -erfTable [ (int)(round((-1.0)*coeff3*10000))]   )   -  ( erfTable [ (int)(round(coeff4*10000)) ]   )   );
	      P_dI+=(   ( -erfTable [ -coeff3Int]   )   -  ( erfTable [ (int) (coeff4*100000) ]   )   );
	    }
	}
      else 
	{
	  //	  printf("Nerde patladi?l=%d\tcoeff3=%lf\tcoeff4=%lf\n", l, coeff3, coeff4);
	  P_dI+=(   ( erf ( coeff3  )   )   -  ( erf ( coeff4  )   )   );
	}
    }
  //  fclose(pairdata);
  //  return P_dI;
  return P_dI*coeff1;
}


void plotModelandTarget(int N, CAcoord *atomsFixed1, pairInfoNew *pair1, int *pair1No, 
			       CAcoord *atomsFixed2, pairInfoNew *pair2, int *pair2No, 
			       double increment, int atomNo, int component)
{
  int i=0;
  double maxDist1=0.0, maxDist2=0.0;
  double maxBigger=0.0;
  double binIndex=1.0;
  
  for(i=0; i<pair1No[0]; i++)
    {
      if(pair1[i].dis>=maxDist1)
	maxDist1=pair1[i].dis;
    }
  
  for(i=0; i<pair2No[0]; i++)
    {
      if(pair2[i].dis>=maxDist2)
	maxDist2=pair2[i].dis;
    }

  FILE *pairdata=fopen("pairFuncData2.txt", "w");
  if(maxDist1>=maxDist2)
    maxBigger=maxDist1;
  else
    maxBigger=maxDist2;
  fprintf(stdout, "Max distance is %lf\n", maxBigger);      
  while(binIndex<(maxBigger + 1.0))
    {
      double P_model=0.0;
      double P_target=0.0;
      
      P_model  = pairDistroGaussian2(N, atomsFixed1, pair1No, pair1, binIndex, increment, atomNo, component);
      P_target = pairDistroGaussian2(N, atomsFixed2, pair2No, pair2, binIndex, increment, atomNo, component);
      fprintf(pairdata, "%lf\t%lf\t%lf\n", binIndex, P_model, P_target);
      binIndex+=1.0;
    }
  fclose(pairdata); 

}


void gradient_P_model(int N, double * grad, CAcoord *atomMoving, pairInfoNew *pair, int *pairNo, double binIndex) 

{
  //Purpose: To obtain gradient of pair distribution function 
  //using pair information kept in  pairInfoNew structure
  int	  n=0, ptr=0;
  double  fx=0.0, fy=0.0, fz=0.0;

  //  double sigmaSqrd=sigma*sigma;
  //  double sigmaCbd=sigma*sigma*sigma;
  double disMobile_ij=0.0;
  double coeff1=(1/(pairNo[0]*2.506628274631*SIGMA_CBD));
  double coeff2=0.0;
  for (n=0;n<3*N;++n)
    grad[n]=0.0;
  
  for (n=0;n<pairNo[0]; ++n) 
    {
      disMobile_ij=dis(atomMoving, (pair[n].posI), (pair[n].posJ));
      double subt=(disMobile_ij-binIndex);
      if(subt<3*SIGMA)
	{
	  //Note that 2.506628274631 is sqrt(2*PI);
	  coeff2= coeff1*exp(-(subt*subt)/(2*SIGMA_SQRD))*((binIndex/disMobile_ij ) - 1);

	  fx=  coeff2*(atomMoving[pair[n].posI].x-atomMoving[pair[n].posJ].x);
	  fy=  coeff2*(atomMoving[pair[n].posI].y-atomMoving[pair[n].posJ].y);
	  fz=  coeff2*(atomMoving[pair[n].posI].z-atomMoving[pair[n].posJ].z);
	  
	  ptr=3*(pair[n].posI);
	  grad[ptr]+=fx;
	  grad[ptr+1]+=fy;
	  grad[ptr+2]+=fz;
	  
	  ptr=3*(pair[n].posJ);
	  grad[ptr]-=fx;
	  grad[ptr+1]-=fy;
	  grad[ptr+2]-=fz;
	}
    }
  
  if(1)  
    {
      FILE *f=fopen("gradP_model.txt", "w");
      if(f==NULL)
	fprintf(stderr, "Can not produce gradPmodel.txt file");
      
      for(n=0; n<3*N; n++)
	fprintf(f, "\t%.10lf\n", grad[n]);
      fclose(f);
    }
  
}


void gradient_P_numerical(int N, double *gradP, CAcoord *atomMoving, pairInfoNew *pair1, int *pair1No, double binIndex)
{
  int i=0, j=0;
  double increment=0.001;//Just for test purposes
  double con=(1/(pair1No[0]*2.506628274631*SIGMA));
  double dis_x=0.0;
  double dis_y=0.0;
  double dis_z=0.0;
  double disMobile_ij=0.0;
  double subt=0.0;

  for(i=0; i<N; i++)
    {
      double P_dI1=0.0;
      double P_dI2=0.0;
      for(j=0; j<N; j++)
	if(i!=j)
	  {
	    dis_x=0.0;
	    dis_y=0.0;
	    dis_z=0.0;
	    disMobile_ij=0.0;
	    dis_x= atomMoving[i].x-atomMoving[j].x;
	    dis_y= atomMoving[i].y-atomMoving[j].y;
	    dis_z= atomMoving[i].z-atomMoving[j].z;
	    disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
	    subt=(disMobile_ij-binIndex);
	    P_dI1+=con*exp(-subt*subt/(2*SIGMA_SQRD));
	    
	    dis_x=0.0;
	    dis_y=0.0;
	    dis_z=0.0;
	    disMobile_ij=0.0;
	    dis_x= (atomMoving[i].x+increment)-atomMoving[j].x;
	    dis_y= atomMoving[i].y-atomMoving[j].y;
	    dis_z= atomMoving[i].z-atomMoving[j].z;
	    disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
	    subt=(disMobile_ij-binIndex);
	    P_dI2+=con*exp(-subt*subt/(2*SIGMA_SQRD ));
	  } 
      gradP[3*i]=(P_dI2-P_dI1)/increment;
    } 
  for(i=0; i<N; i++)
    {
      double P_dI1=0.0;
      double P_dI2=0.0;
      for(j=0; j<N; j++)
	if(i!=j)
	  {
	    dis_x=0.0;
	    dis_y=0.0;
	    dis_z=0.0;
	    disMobile_ij=0.0;
	    dis_x= atomMoving[i].x-atomMoving[j].x;
	    dis_y= atomMoving[i].y-atomMoving[j].y;
	    dis_z= atomMoving[i].z-atomMoving[j].z;
	    disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
	    subt=(disMobile_ij-binIndex);
	    P_dI1+=con*exp(-subt*subt/(2*SIGMA_SQRD));
	    
	    dis_x=0.0;
	    dis_y=0.0;
	    dis_z=0.0;
	    disMobile_ij=0.0;
	    dis_x= atomMoving[i].x-atomMoving[j].x;
	    dis_y= atomMoving[i].y+increment-atomMoving[j].y;
	    dis_z= atomMoving[i].z-atomMoving[j].z;
	    disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
	    subt=(disMobile_ij-binIndex);
	    P_dI2+=con*exp(-subt*subt/(2*SIGMA_SQRD));
	  } 
      gradP[3*i+1]=(P_dI2-P_dI1)/increment;
    } 
  for(i=0; i<N; i++)
    {
      double P_dI1=0.0;
      double P_dI2=0.0;
      for(j=0; j<N; j++)
	if(i!=j)
	  {
	    dis_x=0.0;
	    dis_y=0.0;
	    dis_z=0.0;
	    disMobile_ij=0.0;
	    dis_x= atomMoving[i].x-atomMoving[j].x;
	    dis_y= atomMoving[i].y-atomMoving[j].y;
	    dis_z= atomMoving[i].z-atomMoving[j].z;
	    disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
	    subt=(disMobile_ij-binIndex);
	    P_dI1+=con*exp(-subt*subt/(2*SIGMA_SQRD));
	    
	    dis_x=0.0;
	    dis_y=0.0;
	    dis_z=0.0;
	    disMobile_ij=0.0;
	    dis_x= atomMoving[i].x-atomMoving[j].x;
	    dis_y= atomMoving[i].y-atomMoving[j].y;
	    dis_z= atomMoving[i].z+increment-atomMoving[j].z;
	    disMobile_ij=sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z);
	    subt=(disMobile_ij-binIndex);
	    P_dI2+=con*exp(-subt*subt/(2*SIGMA_SQRD));
	  } 
      gradP[3*i+2]=(P_dI2-P_dI1)/increment;
    } 

  FILE *pairdata=fopen("gradient_P_numerical.txt", "w");  
  for(i=0; i<3*N; i++)
    fprintf(pairdata, "%lf\n", gradP[i]);
  fclose(pairdata);
}


double E_p(CAcoord *atomMoving, int atomCount1, int *pair1No,  pairInfoNew *pair1,
	   CAcoord *atomsFixed2, int atomCount2, int *pair2No,  pairInfoNew *pair2, double maxBigger)
{
  double binIndex=1.0;
  fprintf(stdout, "MaxBigger:%lf\n", maxBigger);
  FILE *pairdata=fopen("E_p.txt", "a");  
  double E_p=0.0;
  while(binIndex<(maxBigger + 1.0))
    {
      double P_model=0.0;
      double P_target=0.0;
      //Just for test purposes, lets scale E_p with 10000,
      //Dont forget to remove those 100s in two equations below!!
      P_model  = pairDistroGaussian1(atomCount1, atomMoving, pair1No, pair1, binIndex);
      P_target = pairDistroGaussian1(atomCount2, atomsFixed2, pair2No, pair2, binIndex);
      double subt= (P_model-P_target);
      //      double coeff= ( 1/ ( 0.5*P_target*P_target )   );
      
      E_p+= ((subt)*(subt)/(P_target*P_target));
      //      fprintf(pairdata, "%lf\t%lf\t%lf\n", binIndex, P_model, P_target);
      binIndex+=1.0;
    }
  fprintf(pairdata, "%lf\n", 10*E_p);
  fclose(pairdata);
  return E_p;
}

void gradient_E_p(int N, double *grad_E_p,  CAcoord *atomMoving,  int *pair1No,  pairInfoNew *pair1,
                        		    CAcoord *atomsFixed2, int *pair2No,  pairInfoNew *pair2, double maxBigger) 
{
  int i=0;
  double binIndex=1.0;
  double P_model=0.0;
  double P_target=0.0;

  for (i=0;i<3*N;i++)
    grad_E_p[i]=0.0;
  double *grad_P_model = (double*) calloc ((3*N), sizeof(double));
  if (grad_P_model==NULL)
    {
      fprintf(stderr, "No memory callocation for grad_P_model calculation\n");
      exit (EXIT_FAILURE);
    }
  
  while(binIndex<(maxBigger + 1.0))
    {
      P_model=0.0;
      P_target=0.0;

      P_model  = pairDistroGaussian1(N, atomMoving,  pair1No, pair1, binIndex);
      P_target = pairDistroGaussian1(N, atomsFixed2, pair2No, pair2, binIndex);
      gradient_P_model(N, grad_P_model, atomMoving, pair1, pair1No,  binIndex);
      double subt= (P_model-P_target);
      double coeff= ( 20/ ( 0.5*P_target*P_target )   );
      for(i=0; i<(3*N); i++)
	{    
	  grad_E_p[i]+= coeff*subt*grad_P_model[i];
	} 
      //      fprintf(stdout, "%lf\t%lf\t%lf\n", binIndex, P_model, P_target);
      
      binIndex+=1.0;
    }

  if(1)
    {
      FILE *pairdata=fopen("gradient_E_p.txt", "w");  
      if(pairdata==NULL)
	{
	  fprintf(stderr, "gradient_E_p.txt file has not been produced\n");
	  exit(EXIT_FAILURE);
	}
      for (i=0; i<(3*N); i++)
	fprintf(pairdata, "%d\t%.10lf\n", i, grad_E_p[i]);
      fclose(pairdata);
    }  
  
  //  fprintf(stdout, "Gradient of scoring function has been calculated\n");
  free(grad_P_model);

}

void gradient_E_p_log(int N, double *grad_E_p,  CAcoord *atomMoving,  int *pair1No,  pairInfoNew *pair1,
		                CAcoord *atomsFixed2, int *pair2No,  pairInfoNew *pair2, double maxBigger)
{
  int i=0, n=0;
  double binIndex=1.0;

  for (n=0;n<3*N;++n)
    grad_E_p[n]=0.0;

  double *gradP = (double*) calloc ((3*N), sizeof(double));
  if (gradP==NULL)
    {
      fprintf(stderr, "No memory callocation for gradP calculation\n");
      exit (EXIT_FAILURE);
    }
  
  while(binIndex<(maxBigger + 1.0))
    {
      double P_model=0.0;
      double P_target=0.0;
      P_model  = pairDistroGaussian1(N, atomMoving,  pair1No, pair1, binIndex);
      P_target = pairDistroGaussian1(N, atomsFixed2, pair2No, pair2, binIndex);
      gradient_P_model(N, gradP, atomMoving, pair1, pair1No,  binIndex);
      for(i=0; i<3*N; i++)
	{    
	  grad_E_p[i]+=0.5*((log(P_model)-log(P_target))/P_model)*gradP[i];
	} 
      //      fprintf(stdout, "%lf\t%lf\t%lf\n", binIndex, P_model, P_target);
      
      binIndex+=1.0;
    }

  if(DEBUGMODE)
    {
      FILE *pairdata=fopen("gradient_E_p.txt", "w");  
      for (n=0;n<3*N;++n)
	fprintf(pairdata, "%d\t%.10lf\n", n, grad_E_p[n]);
      fclose(pairdata);
    }  
  //  fprintf(stdout, "Gradient of scoring function has been calculated\n");
  free(gradP);
}

void gradient_E_p_numerical2(int N, double *grad_E_p,  CAcoord *atomMoving, int *pair1No,  pairInfoNew *pair1,
			     CAcoord *atom2Fixed, int *pair2No,  pairInfoNew *pair2, double increment)
{
  int i=0, j=0;
  double maxDist1=0.0, maxDist2=0.0;
  double maxBigger=0.0;
  double binIndex=1.0;

  for (i=0;i<3*N;++i)
    grad_E_p[i]=0;


  for(i=0; i<pair1No[0]; i++)
    {
      if(pair1[i].dis>=maxDist1)
	maxDist1=pair1[i].dis;
    }
  
  for(i=0; i<pair2No[0]; i++)
    {
      if(pair2[i].dis>=maxDist2)
	maxDist2=pair2[i].dis;
    }
  
  if(maxDist1>=maxDist2)
    maxBigger=maxDist1;
  else
    maxBigger=maxDist2;
  FILE *pairdata=fopen("gradient_E_p_numerical.txt", "w");
  //  double increment=0.1;
  //Keep in mind that their pair numbers may not be same if you dont consider all pairs!!
  //Therefore, the formula below may have to change!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  double con=(1/(pair1No[0]*2.506628274631*sigma));
  //  double dis_x=0.0;
  // double dis_y=0.0;
  //  double dis_z=0.0;
  //  double disMobile_ij=0.0;
  //  double subt=0.0;
  double P_model1=0.0;
  double P_model2=0.0;
  double P_target=0.0;

  double E_p1=0.0;
  double E_p2=0.0;	  
  for(i=0; i<N; i++)
    for(j=0; j<3; j++)
      {
	E_p1=0.0;
	E_p2=0.0;
	while(binIndex<(maxBigger + 1.0))
	  {
	    P_model1=0.0;
	    P_model2=0.0;
	    P_target=0.0;
	    
	    P_model1 = pairDistroGaussian2(N, atomMoving, pair1No, pair1, binIndex, 0.0000000, i, j);
	    P_model2 = pairDistroGaussian2(N, atomMoving, pair1No, pair1, binIndex, increment, i, j);
	    
	    P_target = pairDistroGaussian2(N, atom2Fixed, pair2No, pair2, binIndex, 0.0000000, i, j);
	    
	    //    fprintf(stdout, "P_model1=%.8lf\tP_model2=%.8lf\tP_target=%.8lf\n", P_model1, P_model2, P_target);
	    E_p1+=0.5*(P_model1-P_target)*(P_model1-P_target);
	    E_p2+=0.5*(P_model2-P_target)*(P_model2-P_target);
	    binIndex+=1.0;
	  }
	fprintf(stdout, "%.10lf\n", E_p1);
	fprintf(stdout, "%.10lf\n", E_p2);
	binIndex=1.0;
	grad_E_p[3*i+j]=(E_p2-E_p1)/increment;
      }
  
  if(1)
    for (i=0;i<3*N;++i)
      fprintf(pairdata, "%.10lf\n", grad_E_p[i]);
  
  
  fprintf(stdout, "Gradient of scoring function has been calculated numerically.\n");
  fclose(pairdata);
}

void hessian_P_model(int N/*Number of residues in model*/, 
		     double **H_p, CAcoord *atomMoving, 
		     int *pairNo, pairInfoNew *pair, double binSize)
{
  //Purpose: Constructing hessian of pair distribution function(H_p)
  //Next   : Use a sparse matrix format instead of usual array
  int i=0, j=0, l=0;  /*Counters for atom pairs                           */
  int a=0, b=0;       /*Counters for atom cartesian coordinates           */

  for (i=0; i<3*N; i++) /*To zero input hessian so that previous calculation have nothing to do w/ this one*/
    for (j=0; j<=i; j++)
      {
	H_p[i][j]=0.0;
	H_p[j][i]=0.0;
      }
  for(l=0; l<pairNo[0]; ++l)         
    {
      double disMobile_ij=0.0;
      disMobile_ij= dis(atomMoving, pair[l].posI, pair[l].posJ);
     
      double disSqMob_ij=0.0;
      disSqMob_ij=(disMobile_ij*disMobile_ij);
      
      double coeff1=0.0, coeff2=0.0;
      
      double subt1=(disMobile_ij-binSize);

      if(subt1<3*SIGMA)
      	{
	  //Note that 2.506628274631 is sqrt(2*PI)
	  //Dont forget N_pair!!!!!!!!!!!!!!!!!!!!
	  coeff1= (1/ (pairNo[0]*2.506628274631*SIGMA_CBD))*exp(-(subt1*subt1)/(2*SIGMA_SQRD));
	  coeff2= (binSize/disMobile_ij);
	  //      fprintf(stdout, "Subtraction1:%.10lf\tCoefficient1:%.10lf\tCoefficient2:%.10lf\n", subt1, coeff1, coeff2);
	  //      fprintf(stdout, "Exponential Part: %.10lf\n", exp(-(subt1*subt1)/(2*sigmaSqrd)));
	  double x_i[3]={999.0, 999.0, 999.0};
	  double x_j[3]={999.0, 999.0, 999.0};
	  
	  x_i[0]=atomMoving[pair[l].posI].x;
	  x_i[1]=atomMoving[pair[l].posI].y;
	  x_i[2]=atomMoving[pair[l].posI].z;
	  
	  x_j[0]=atomMoving[pair[l].posJ].x;
	  x_j[1]=atomMoving[pair[l].posJ].y;
	  x_j[2]=atomMoving[pair[l].posJ].z;
	  
	  for(a=0; a<3; ++a)
	    for(b=0; b<3; ++b) 
	      {
		double subt2=( x_i[a] - x_j[a]);
		double subt3=( x_i[b] - x_j[b]);
		if(a!=b)
		  {
		    double offdiagonal = /*Dogru mu isaret?*/ +coeff1 * ( (subt2*subt3/disSqMob_ij)*(coeff2- subt1*subt1/SIGMA_SQRD )); 
		    //	  printf("offdiag: %lf\n", offdiagonal);    
		    H_p[3*pair[l].posI+a][3*pair[l].posJ+b]=offdiagonal;
		    H_p[3*pair[l].posJ+b][3*pair[l].posI+a]=offdiagonal;
		    //	  printf("pair[%d].posI: %d\n", l, pair[l].posI);
		    //	  printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		    //	  printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), H_p[3*pair[l].posI+a][3*pair[l].posJ+b]);
		  }
		else if (a==b)
		  {
		    double ondiagonal =coeff1*( 1 - (coeff2) + (subt2*subt3/disSqMob_ij)*(  coeff2 - subt1*subt1/SIGMA_SQRD ) );
		    //	 printf("ondiag: %lf\n", ondiagonal);    
		    H_p[3*pair[l].posI+a][3*pair[l].posJ+b]=ondiagonal;
		    H_p[3*pair[l].posJ+b][3*pair[l].posI+a]=ondiagonal;
		    //	  printf("pair[%d].posI: %d\n", l, pair[l].posI);
		    //	  printf("pair[%d].posJ: %d\n", l, pair[l].posJ);
		    //	  printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[l].posI+a) , (3*pair[l].posJ+b), H_p[3*pair[l].posI+a][3*pair[l].posJ+b]);
		  }
	      }
	}	
    }
  
  
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      H_p[3*i+a][3*i+b] -= H_p[3*i+a][3*j+b]; 
 	      //printf("H_p[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_p[3*i+a][3*i+b]); 
 	    } 
  // printf("\nCalculated Hessian Matrix\n");
  
  if(0)
    {
      FILE *matrix=fopen("hessian_P.txt", "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "hessian_P.txt file could not be produced\n");
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
	for(j=0; j<=i; j++)
	  for(a=0; a<3; a++) 
	    for(b=0; b<3; b++) 
	      if(H_p[3*i+a][3*j+b]!=0.0)
		fprintf(matrix, "%.10lf\n", H_p[3*i+a][3*j+b]); 
      //	      fprintf(matrix, "%d\t%d\t%.10lf\n", 3*i+a, 3*j+b, H_p[3*i+a][3*j+b]); 
      
      //fprintf(matrix, "\n"); 
      //    } 
  //?///////////////////////////////??????????????????????????????????????????????????
/*       for(i=0; i<N; i++) */
/*  	{ */
/*  	  for(j=0; j<N; j++) */
/* 	    if(i==j) */
/* 	      { */
/* 		for(a=0; a<3; a++) */
/* 		  for(b=0; b<3; b++) */
/* 		    if(a<=b) */
/* 		      if(H_p[3*i+a][3*j+b]!=0.0) */
/* 			fprintf(matrix, "%d\t%d\t%lf\n",  3*i+a, 3*j+b,  H_p[3*i+a][3*j+b]); */
		
/* 		//fprintf(matrix, "\n"); */
/* 	      } */
/*  	}  */
      
      fclose(matrix);
    }
  
}
void hessian_P_numerical(int N, double **H_p, CAcoord *atomMoving, int *pairNo, pairInfoNew *pair, double binSize, double increment)
{
  int i=0, a=0, j=0, b=0;
  fprintf(stdout, "Sen istiyor duj, sen verecek 50 dolar daha\n");
  double *gradP1 = (double*) calloc (3*N, sizeof(double));
  if (gradP1==NULL)
    {
      fprintf(stderr, "No memory allocation fo gradP1 calculation\n");
      exit (EXIT_FAILURE);
    }

  double *gradP2 = (double*) calloc (3*N, sizeof(double));
  if (gradP2==NULL)
    {
      fprintf(stderr, "No memory allocation fo gradP2 calculation\n");
      exit (EXIT_FAILURE);
    }
  for(i=0; i<3*N; i++)
    for(j=0; j<=i; j++)
      {
	H_p[i][j]=0.0;
	H_p[j][i]=0.0;
      }
  for(i=0; i<N; i++)
    for(j=0; j<i; j++)
      {
      for(a=0; a<3; a++)
	{

	  gradient_P_model(N, gradP1, atomMoving, pair, pairNo, binSize) ;
	  if(a==0) atomMoving[i].x=atomMoving[i].x+increment;
	  if(a==1) atomMoving[i].y=atomMoving[i].y+increment;
	  if(a==2) atomMoving[i].z=atomMoving[i].z+increment;
	  gradient_P_model(N, gradP2, atomMoving, pair, pairNo, binSize) ;
	  
	  for(b=0; b<3; b++)
	    H_p[3*i+a][3*j+b]= (gradP2[3*j+b]-gradP1[3*j+b])/increment;
	  
	  if(a==0) atomMoving[i].x=atomMoving[i].x-increment;
	  if(a==1) atomMoving[i].y=atomMoving[i].y-increment;
	  if(a==2) atomMoving[i].z=atomMoving[i].z-increment;
	  
	}
      fprintf(stdout, "Sen istiyor duj, sen verecek 50 dolar daha: %d\n", j);
      }
  fprintf(stdout, "Sen istiyor duj, sen verecek 50 dolar daha: SON\n");
  //Upper symmetric offdiagonal part      
  for(i=0; i<N; i++)
    for(j=0; j<i;j++)
      for(a=0; a<3; a++)
	for(b=0; b<3; b++)
	 H_p[3*j+b][3*i+a]=H_p[3*i+a][3*j+b];
  //Diagonal blocks: Bunu sorgula!!!!!!!!!!
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      H_p[3*i+a][3*i+b] -= H_p[3*i+a][3*j+b]; 
 	      //printf("H_p[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_p[3*i+a][3*i+b]); 
 	    } 

  if(DEBUGMODE)
    {
      FILE *matrix=fopen("hessian_P_numerical.txt", "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "hessian_P_numerical.txt file could not be produced\n");
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
	for(j=0; j<=i; j++)
	  for(a=0; a<3; a++) 
	    for(b=0; b<3; b++) 
	      //	      if(H_p[3*i+a][3*j+b]!=0.0)
	      fprintf(matrix, "%.10lf\n", H_p[3*i+a][3*j+b]); 
      //	      fprintf(matrix, "%d\t%d\t%.10lf\n", 3*i+a, 3*j+b, H_p[3*i+a][3*j+b]); 
      
      fclose(matrix);
    }


  free(gradP1);
  free(gradP2);
}
void hessian_E_p_numerical(int N/*System size*/, double **H_E_p,
			   CAcoord *atomMoving,  CAcoord *atom2Fixed,
			   int *pair1No, pairInfoNew *pair1, 
			   int *pair2No, pairInfoNew *pair2, 
			   double increment, double maxBigger)
{
  int i=0, a=0, j=0, b=0;
 
  double *gradE_p1 = (double*) calloc (3*N, sizeof(double));
  if (gradE_p1==NULL)
    {
      fprintf(stderr, "No memory allocation fo gradE_p1 calculation\n");
      exit (EXIT_FAILURE);
    }

  double *gradE_p2 = (double*) calloc (3*N, sizeof(double));
  if (gradE_p2==NULL)
    {
      fprintf(stderr, "No memory allocation fo gradE_p2 calculation\n");
      exit (EXIT_FAILURE);
    }
  for(i=0; i<3*N; i++)
    for(j=0; j<3*N; j++)
      H_E_p[i][j]=0.0;

  FILE *matrix=fopen("hessian_E_p_numerical.txt", "w");
  if (matrix==NULL)
    {
      fprintf(stderr, "hessian_E_p_numerical.txt file could not be produced\n");
      exit(EXIT_FAILURE);
    }
  
  gradient_E_p(N, gradE_p1,  atomMoving, pair1No,  pair1, atom2Fixed, pair2No,  pair2, maxBigger);  
  for(i=0; i<N; i++)
    for(j=0; j<=i; j++)
      {
      for(a=0; a<3; a++)
	{

	  if(a==0) atomMoving[i].x=atomMoving[i].x+increment;
	  if(a==1) atomMoving[i].y=atomMoving[i].y+increment;
	  if(a==2) atomMoving[i].z=atomMoving[i].z+increment;
	  gradient_E_p(N, gradE_p2,  atomMoving, pair1No,  pair1, atom2Fixed, pair2No,  pair2, maxBigger);

	  
	  for(b=0; b<3; b++)
	    {
	    H_E_p[3*i+a][3*j+b]= (gradE_p2[3*j+b]-gradE_p1[3*j+b])/increment;
	    fprintf(matrix, "%.10lf\n", H_E_p[3*i+a][3*j+b]);  
	    }
	  
	  if(a==0) atomMoving[i].x=atomMoving[i].x-increment;
	  if(a==1) atomMoving[i].y=atomMoving[i].y-increment;
	  if(a==2) atomMoving[i].z=atomMoving[i].z-increment;
	  
	}
      fprintf(stdout, "Sen istiyor duj, sen verecek 50 dolar daha: %d\n", j);
      }

/*   //Upper symmetric offdiagonal part       */
/*   for(i=0; i<N; i++) */
/*     for(j=0; j<i;j++) */
/*       for(a=0; a<3; a++) */
/* 	for(b=0; b<3; b++) */
/* 	 H_E_p[3*j+b][3*i+a]=H_E_p[3*i+a][3*j+b]; */
  //Diagonal blocks: Bunu sorgula!!!!!!!!!!
/*   for(i=0;i<N;++i) */
/*     for(j=0;j<N;++j) */
/*       if(i!=j)  */
/*  	for(a=0; a<3; ++a)  */
/*  	  for(b=0; b<3; ++b)  */
/*  	    {  */
/* 	      H_E_p[3*i+a][3*i+b] -= H_E_p[3*i+a][3*j+b];  */
/*  	      //printf("H_E_p[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_E_p[3*i+a][3*i+b]);  */
/*  	    }  */

/*   if(1) */
/*     { */
/*       FILE *matrix=fopen("hessian_E_p_numerical.txt", "w"); */
/*       if (matrix==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "hessian_E_p_numerical.txt file could not be produced\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       for(i=0; i<N; i++)  */
/* 	for(j=0; j<=i; j++) */
/* 	  for(a=0; a<3; a++)  */
/* 	    for(b=0; b<3; b++)  */
/* 	      //	      if(H_E_p[3*i+a][3*j+b]!=0.0) */
/* 	      fprintf(matrix, "%.10lf\n", H_E_p[3*i+a][3*j+b]);  */
/*       //	      fprintf(matrix, "%d\t%d\t%.10lf\n", 3*i+a, 3*j+b, H_E_p[3*i+a][3*j+b]);  */
      
      fclose(matrix);
      //    }

  free(gradE_p1);
  free(gradE_p2);
}

void hessian_E_p(int N/*System size*/, double **H_E_p,
		 CAcoord *atomMoving, CAcoord *atomsFixed2,
		 int *pair1No, pairInfoNew *pair1, 
		 int *pair2No, pairInfoNew *pair2, double maxBigger)
{
 
  int i=0, j=0, a=0, b=0;
  //Zero all values from previous calculation
  for(i=0; i<(3*N); i++)
    for(j=0; j<(3*N); j++)
      H_E_p[i][j]=0.0;

  double binSize=1.0;
  double *grad_P_model = (double*) malloc (3*N*sizeof(double));
  if (grad_P_model==NULL)
    {
      fprintf(stderr, "No memory allocation fo grad_P_model calculation\n");
      exit (EXIT_FAILURE);
    }
  
  double **H_part1 = (double**)malloc(3*N*sizeof(double*));
  if(H_part1==NULL)
    {
      fprintf(stderr, "Calloc cannot allocate memory for H_part1 array");
      exit(EXIT_FAILURE);
    }
  for (i=0; i<3*N; i++)
    {
      H_part1[i] = (double*)calloc(3*N,sizeof(double));
      if(H_part1[i]==NULL)
	{
	  fprintf(stderr, "Calloc cannot allocate memory for H_part1[i] array");
	  exit(EXIT_FAILURE);
	}
    }

  double value=0.0;
  while(binSize<(maxBigger + 1.0))
    {
      double P_model=0.0;
      double P_target=0.0;
      
      P_model  = pairDistroGaussian1(N, atomMoving, pair1No, pair1, binSize);
      P_target = pairDistroGaussian1(N, atomsFixed2, pair2No, pair2, binSize);
      
      hessian_P_model(N, H_part1, atomMoving, pair1No, pair1, binSize);
      gradient_P_model(N, grad_P_model, atomMoving, pair1, pair1No, binSize);
      double subt =(P_model-P_target);
      double coeff= ( 20/ ( 0.5*P_target*P_target )   );
      for(i=0; i<3*N; i++)
	for(j=0; j<3*N; j++)
	  {
	    value=0.0;
	    value=grad_P_model[i]*grad_P_model[j];
	    H_E_p[i][j]+= coeff*(  (subt)*H_part1[i][j] +  value );
	    
	    //  H_E_p[i][j]+= (  (P_model-P_target)*H_part1[i][j] +  H_part2[i][j] );
	    // fprintf(stdout, "%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n", binSize, P_model, P_target, i, j, H_part1[i][j], H_part2[i][j], H_E_p[i][j]);
	  }
      binSize+=1.0;
    }
  if(1)
    {
      FILE *f=fopen("hessian_E_p.txt", "w");
      if(f==NULL)
	{
	fprintf(stderr, "hessian_E_p.txt file has not been produced\n");
	exit(EXIT_FAILURE);
	}
      for(i=0; i<N; i++) 
	for(j=0; j<=i; j++)
	  for(a=0; a<3; a++) 
	    for(b=0; b<3; b++) 
	      fprintf(f, "%.10lf\n",   H_E_p[3*i+a][3*j+b]); 
      fclose(f);
    }
  
  fprintf(stdout, "Hessian of scoring function has been calculated\n");
  free(H_part1);
  free(grad_P_model);
}
void gradientANDhessian_P(int N, double * grad, double **H_p, CAcoord *atomMoving, 
			  pairInfoNew *pair, int *pairNo, 
			  pairInfoNew *pairInitial, int *pairInitialNo, 
			  double binIndex) 
{
  //Purpose: To obtain gradient of pair distribution function 
  //using pair information kept in  pairInfoNew structure
  int	  i=0, j=0, a=0, b=0, n=0, ptr=0;
  double  fx=0.0, fy=0.0, fz=0.0;

  double disMobile_ij=0.0;
  double coeff0=(1/(pairNo[0]*2.506628274631*SIGMA_CBD));
  double coeff1= 0.0, coeff2= 0.0, coeff3= 0.0;

  for (i=0; i<(3*N); i++)
    {
      for (j=0;j<=i; j++)
	{
	  H_p[i][j]=0.0;
	  H_p[j][i]=0.0;
	}
      grad[i]=0.0;
    }

  for (n=0;n<pairNo[0]; ++n) 
    {
      //      disMobile_ij=dis(atomMoving, (pair[n].posI), (pair[n].posJ));
      disMobile_ij=pairInitial[n].dis;
      double disSqMob_ij=0.0;
      disSqMob_ij=(disMobile_ij*disMobile_ij);

      double subt=(disMobile_ij-binIndex);
      double subtSqrd=subt*subt;
      if(subt<3*SIGMA)
	{
	  coeff1=exp(-(subtSqrd)/(2*SIGMA_SQRD));
	  coeff2=(binIndex/disMobile_ij);

	  //Note that 2.506628274631 is sqrt(2*PI);
	  coeff3= coeff0*coeff1*( coeff2 - 1);
	  //---Gradient P computation-----------
	  fx=  coeff3*(atomMoving[pair[n].posI].x-atomMoving[pair[n].posJ].x);
	  fy=  coeff3*(atomMoving[pair[n].posI].y-atomMoving[pair[n].posJ].y);
	  fz=  coeff3*(atomMoving[pair[n].posI].z-atomMoving[pair[n].posJ].z);
	  
	  ptr=3*(pair[n].posI);
	  grad[ptr]+=fx;
	  grad[ptr+1]+=fy;
	  grad[ptr+2]+=fz;
	  
	  ptr=3*(pair[n].posJ);
	  grad[ptr]-=fx;
	  grad[ptr+1]-=fy;
	  grad[ptr+2]-=fz;
	  //---End of gradient P computation---
	  //---Hessian P computation-----------

	  double x_i[3]={999.0, 999.0, 999.0};
	  double x_j[3]={999.0, 999.0, 999.0};

	  x_i[0]=atomMoving[pair[n].posI].x;
	  x_i[1]=atomMoving[pair[n].posI].y;
	  x_i[2]=atomMoving[pair[n].posI].z;
	  
	  x_j[0]=atomMoving[pair[n].posJ].x;
	  x_j[1]=atomMoving[pair[n].posJ].y;
	  x_j[2]=atomMoving[pair[n].posJ].z;
	  for(a=0; a<3; ++a)
	    for(b=0; b<3; ++b) 
	      {
		double subt2=( x_i[a] - x_j[a]);
		double subt3=( x_i[b] - x_j[b]);
		double subt2Xsubt3=subt2*subt3;
		if(a!=b)
		  {
		    double offdiagonal = /*Dogru mu isaret?*/+coeff0*coeff1*((subt2Xsubt3/disSqMob_ij)*(coeff2- subtSqrd/SIGMA_SQRD )); 
		    H_p[3*pair[n].posI+a][3*pair[n].posJ+b]=offdiagonal;
		    H_p[3*pair[n].posJ+b][3*pair[n].posI+a]=offdiagonal;
		    //	  fprintf(stdout, "pair[%d].posI: %d\tOffdiagonal: %lf\n", n, pair[n].posI, offdiagonal);
		    //	  fprintf(stdout, "pair[%d].posJ: %d\n", n, pair[n].posJ);
		    //	  printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[n].posI+a) , (3*pair[n].posJ+b), H_p[3*pair[n].posI+a][3*pair[n].posJ+b]);
		  }
		else if (a==b)
		  {
		    double ondiagonal =coeff0*coeff1*( 1 - coeff2 + (subt2Xsubt3/disSqMob_ij)*(  coeff2 - subtSqrd/SIGMA_SQRD ) );
		    H_p[3*pair[n].posI+a][3*pair[n].posJ+b]=ondiagonal;
		    H_p[3*pair[n].posJ+b][3*pair[n].posI+a]=ondiagonal;
		    //	 fprintf(stdout, "pair[%d].posI: %d\tOndiagonal:%lf\n", n, pair[n].posI, ondiagonal);
		    //	 fprintf(stdout, "pair[%d].posJ: %d\n", n, pair[n].posJ);
		    //	 printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[n].posI+a) , (3*pair[n].posJ+b), H_p[3*pair[n].posI+a][3*pair[n].posJ+b]);
		  }
	      }
	  //---End hessian P computation-------
	}
    }

  
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i!=j) 
 	for(a=0; a<3; ++a) 
 	  for(b=0; b<3; ++b) 
 	    { 
	      H_p[3*i+a][3*i+b] -= H_p[3*i+a][3*j+b]; 
 	      //printf("H_p[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_p[3*i+a][3*i+b]); 
 	    } 
  // printf("\nCalculated Hessian Matrix\n");
  if(DEBUGMODE)  
    {
      FILE *f=fopen("gradP_modelComb.txt", "w");
      if(f==NULL)
	{
	fprintf(stderr, "Can not produce gradPmodelComb.txt file");
	exit(EXIT_FAILURE);
	}
      
      for(n=0; n<3*N; n++)
	fprintf(f, "\t%.10lf\n", grad[n]);
      fclose(f);
    }

  if(DEBUGMODE)
    {
      FILE *matrix=fopen("hessian_P.txt", "w");
      if (matrix==NULL)
	{
	  fprintf(stderr, "hessian_P.txt file could not be produced\n");
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
	for(j=0; j<=i; j++)
	  for(a=0; a<3; a++) 
	    for(b=0; b<3; b++) 
	      if(H_p[3*i+a][3*j+b]!=0.0)
		fprintf(matrix, "%.10lf\n", H_p[3*i+a][3*j+b]); 
      //	      fprintf(matrix, "%d\t%d\t%.10lf\n", 3*i+a, 3*j+b, H_p[3*i+a][3*j+b]); 
      fclose(matrix);
    }

  //hessianpart
}
//======================================================
//=========Gradient and hessian of E_p in the same loops


void gradientANDhessian_E_p(int N/*System size*/, double **H_E_p, double *grad_E_p,
			    CAcoord *atomMoving,  CAcoord *atomsFixed2,
			    int *pair1No,         pairInfoNew *pair1, 
			    int *pair2No,         pairInfoNew *pair2, 
			    int *pairInitialNo,   pairInfoNew *pairInitial, 
			    pdf *target,
			    double maxBigger,     double scaler)
{
 
  int i=0, j=0, a=0, b=0;
  //Zero all values from previous calculation
  for(i=0; i<(3*N); i++)
    {
      for(j=0; j<=i; j++)
	{
	  H_E_p[i][j]=0.0;
	  H_E_p[j][i]=0.0;

	}
      grad_E_p[i]=0.0;
    }


/*   int l=0; */
/*   for(i=0; i<N; i++) */
/*     for(j=0; j<i; j++) */
/*     { */
/*       pairInitial[l].dis=dis(atomMoving, i, j); */
/*       fprintf(stdout, "Dis[%d]=%lf\n", l, pairInitial[l].dis); */
/*       l++; */
/*     } */
  
  double binIndex=1.0;

  double *grad_P_model = (double*) calloc ((3*N), sizeof(double));
  if (grad_P_model==NULL)
    {
      fprintf(stderr, "No memory callocation for grad_P_model calculation\n");
      exit (EXIT_FAILURE);
    }
  
  double **H_part1 = (double**)malloc(3*N*sizeof(double*));
  if(H_part1==NULL)
    {
      fprintf(stderr, "Malloc cannot allocate memory for H_part1 array");
      exit(EXIT_FAILURE);
    }
  for (i=0; i<3*N; i++)
    {
      H_part1[i] = (double*)calloc(3*N,sizeof(double));
      if(H_part1[i]==NULL)
	{
	  fprintf(stderr, "Calloc cannot allocate memory for H_part1[i] array");
	  exit(EXIT_FAILURE);
	}
    }
  double P_model=0.0;
  double P_target=0.0;
  double value=0.0;

/*   FILE *f=fopen("P_t_funct.txt", "a"); */
/*   if(f==NULL) */
/*     { */
/*       fprintf(stderr, "P_t_funct.txt file not produced\n"); */
/*       exit(EXIT_FAILURE); */
/*     } */
     
  while(binIndex<(maxBigger + 1.0))
    {
      P_model=0.0;
      P_target=0.0;
      P_model  = pairDistroGaussian3(N, atomMoving,  pair1No, pair1, pairInitialNo, pairInitial, binIndex);
      //      P_model  = pairDistroGaussian1(N, atomMoving,  pair1No, pair1, binIndex);

      //      P_target = pairDistroGaussian1(N, atomsFixed2, pair2No, pair2, binIndex);
      int test=binIndex-1;
      P_target=target[test].binValue;
      //  fprintf(f,"%lf\n", P_target);      
      //      hessian_P_model (N, H_part1,      atomMoving, pair1No, pair1, binIndex);
      //      gradient_P_model(N, grad_P_model, atomMoving, pair1, pair1No, binIndex);
      gradientANDhessian_P(N, grad_P_model, H_part1, atomMoving, pair1, pair1No, pairInitial, pairInitialNo, binIndex);
      double subt =(P_model-P_target);
      //Logarithmic: double subt =(log(P_model)-log(P_target));

      double coeff= (scaler/ ( 0.5*P_target*P_target )   );
      //Change it to P_target instead of P_target^2
      //      double coeff= (scaler/ ( 0.5*P_target)   );

      for(i=0; i<3*N; i++)
	{
	  for(j=0; j<3*N; j++)
	    {
	      value=0.0;
	      value=grad_P_model[i]*grad_P_model[j];
	      H_E_p[i][j]+= coeff*( subt*H_part1[i][j] +  value );
	      //	      H_E_p[i][j]+= coeff*( value );
	      //Logarithmic: H_E_p[i][j]+=  (  subt/P_model) *H_part1[i][j] +   ( ( 1- subt  )/( P_model * P_model))*value;

	      // fprintf(stdout, "%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n", binIndex, P_model, P_target, i, j, H_part1[i][j], H_part2[i][j], H_E_p[i][j]);
	    }
	  grad_E_p[i]+= coeff*subt*grad_P_model[i];	  
	  //Logarithmic: grad_E_p[i]+=0.5*(( subt )/P_model)*grad_P_model[i];
	}
      binIndex+=1.0;
    }
  //  fclose(f);
  if(0)
    {
      FILE *pairdata=fopen("gradient_E_p.txt", "w");  
      if(pairdata==NULL)
	{
	  fprintf(stderr, "gradient_E_p.txt file has not been produced\n");
	  exit(EXIT_FAILURE);
	}
      
      for (i=0;i<3*N;i++)
	fprintf(pairdata, "%d\t%.10lf\n", i, grad_E_p[i]);
      fclose(pairdata);
    }  
  
  
  if(0)
    {
      FILE *f=fopen("hessian_E_p_old.txt", "w");
      if(f==NULL)
	{
	  fprintf(stderr, "hessian_E_p.txt file has not been produced\n");
	  exit(EXIT_FAILURE);
	}
      
      for(i=0; i<N; i++) 
	for(j=0; j<=i; j++)
	  for(a=0; a<3; a++) 
	    for(b=0; b<3; b++) 
	      fprintf(f, "%.10lf\n",   H_E_p[3*i+a][3*j+b]); 
      fclose(f);
    }
  
  fprintf(stdout, "Hessian and gradient of scoring function has been calculated\n");
  free(H_part1);
  free(grad_P_model);
}
/*From now on, I will use integral form of pair distribution function described in Gorba, Tama.*/
/*All functions in integral form will have a '_I' suffix.                                      */


/* double pairDistroGaussian3_I(int N, CAcoord *atomMoving,  */
/* 			     int *pairNo,  pairInfoNew *pair,  */
/* 			     int *pairInitialNo,  pairInfoNew *pairInitial,  */
/* 			     double binIndex, double binSize) */
/* { */
/*   //Purpose: This function produces data of pair distribution function  for a coarse grained model using all pairs. */
/*   //         The functional form is error function. */
/*   //In equation sheet: binIndex=d_I, binSize= u */


/*   //If there is not any restriction: N*(N-1)/2 pairs exist */
/*   //where N is my CA atom number for m simple model. */
/*   int i=0; //Regular counters */
/*   //  double maxDistance=0.0; */
/*   //  FILE *pairdata=fopen("pairFuncData_I.txt", "w"); */
/*   /\*   for(i=0; i<pairNo[0]; i++) *\/ */
/*   /\*     { *\/ */
/*   /\*       if(pair[i].dis>=maxDistance) *\/ */
/*   /\* 	maxDistance=pair[i].dis; *\/ */
/*   /\*     } *\/ */
/*   //  fprintf(stdout, "Max distance is %lf\n", maxDistance);     */
/*   //  double sigmaSqrd=sigma*sigma; */
/*   double disMobile_ij=0.0; */
/*   double P_dI=0.0; */
/*   double coeff1=(1/(pairNo[0]*2*binSize)); */
/*   double coeff2=1.41421356*SIGMA; */
/*   double halfBin=binSize/2; */
/*   for(i=0; i<pairNo[0]; i++) */
/*     { */
/*       disMobile_ij=pairInitial[i].dis; */
/*       //      disMobile_ij=dis(atomMoving, (pair[i].posI), (pair[i].posJ)); */
/*       double subt=(binIndex-disMobile_ij); */
/*       //Now, lets put the distance restriction off for now. */
/*       //      if(subt<3*SIGMA) */
/*       P_dI+=coeff1* (   ( erf (  (subt+halfBin)/(coeff2 )  )   )   -  ( erf (  (subt-halfBin)/(coeff2)  )   )   ); */
/*     } */
/*   //  fclose(pairdata); */
/*   return P_dI; */
/* } */

/* void gradientANDhessian_P_I(int N, double * grad, double **H_p, CAcoord *atomMoving, */
/* 			    pairInfoNew *pair, int *pairNo,  */
/* 			    pairInfoNew *pairInitial, int *pairInitialNo,  */
/* 			    double binIndex, double binSize)  */
/* { */
/*   //Purpose: To obtain gradient of pair distribution function  */
/*   //using pair information kept in  pairInfoNew structure */
/*   int	  i=0, j=0, a=0, b=0, n=0, ptr=0; */
/*   double  fx=0.0, fy=0.0, fz=0.0; */

/*   double disMobile_ij=0.0; */
/*   double coeff0=(1/(pairNo[0]*2.506628274631*SIGMA*binSize)); */
/*   double coeff1= 0.0, coeff2= 0.0; */
/*   double coeffA=0.0, coeffB=0.0; */

/*   for (i=0; i<(3*N); i++) */
/*     { */
/*       for (j=0;j<=i; j++) */
/* 	{ */
/* 	  H_p[i][j]=0.0; */
/* 	  H_p[j][i]=0.0; */
/* 	} */
/*       grad[i]=0.0; */
/*     } */

/*   for (n=0;n<pairNo[0]; ++n)  */
/*     { */
/*       //      disMobile_ij=dis(atomMoving, (pair[n].posI), (pair[n].posJ)); */
/*       disMobile_ij=pairInitial[n].dis; */
/*       double disSqMob_ij=0.0; */
/*       disSqMob_ij=(disMobile_ij*disMobile_ij); */

/*       double subt1=(disMobile_ij - binIndex + binSize/2); */
/*       double subt2=(disMobile_ij - binIndex - binSize/2); */
/*       double subtSqrd1=subt1*subt1; */
/*       double subtSqrd2=subt2*subt2; */
/*       //Turn it on later!!!!!!  */
/*       //      if(subt<3*SIGMA) */
/* 	{ */
/* 	  coeffA=exp(-(subtSqrd1)/(2*SIGMA_SQRD)); */
/* 	  coeffB=exp(-(subtSqrd2)/(2*SIGMA_SQRD)); */
/* 	  coeff1=(coeffA-coeffB); */

/* 	  //Note that 2.506628274631 is sqrt(2*PI); */
/* 	  coeff2= coeff0*coeff1/disMobile_ij; */
/* 	  //---Gradient P computation----------- */
/* 	  fx=  coeff2*(atomMoving[pair[n].posI].x-atomMoving[pair[n].posJ].x); */
/* 	  fy=  coeff2*(atomMoving[pair[n].posI].y-atomMoving[pair[n].posJ].y); */
/* 	  fz=  coeff2*(atomMoving[pair[n].posI].z-atomMoving[pair[n].posJ].z); */
	  
/* 	  ptr=3*(pair[n].posI); */
/* 	  grad[ptr]+=fx; */
/* 	  grad[ptr+1]+=fy; */
/* 	  grad[ptr+2]+=fz; */
	  
/* 	  ptr=3*(pair[n].posJ); */
/* 	  grad[ptr]-=fx; */
/* 	  grad[ptr+1]-=fy; */
/* 	  grad[ptr+2]-=fz; */
/* 	  //---End of gradient P computation--- */
/* 	  //---Hessian P computation----------- */

/* 	  double x_i[3]={999.0, 999.0, 999.0}; */
/* 	  double x_j[3]={999.0, 999.0, 999.0}; */

/* 	  x_i[0]=atomMoving[pair[n].posI].x; */
/* 	  x_i[1]=atomMoving[pair[n].posI].y; */
/* 	  x_i[2]=atomMoving[pair[n].posI].z; */
	  
/* 	  x_j[0]=atomMoving[pair[n].posJ].x; */
/* 	  x_j[1]=atomMoving[pair[n].posJ].y; */
/* 	  x_j[2]=atomMoving[pair[n].posJ].z; */
/* 	  for(a=0; a<3; ++a) */
/* 	    for(b=0; b<3; ++b)  */
/* 	      { */
/* 		double subt3=( x_i[a] - x_j[a]); */
/* 		double subt4=( x_i[b] - x_j[b]); */
/* 		double subt3Xsubt4=subt3*subt4; */
/* 		double subANDdiv=subt3Xsubt4/disSqMob_ij; */
/* 		if(a!=b) */
/* 		  { */
/* 		    double offdiagonal = coeff0*( coeff1 + (-coeffA*subt1 + coeffB*subt2)/SIGMA_SQRD  )*subANDdiv;  */
/* 		    H_p[3*pair[n].posI+a][3*pair[n].posJ+b]=offdiagonal; */
/* 		    H_p[3*pair[n].posJ+b][3*pair[n].posI+a]=offdiagonal; */
/* 		    //	  fprintf(stdout, "pair[%d].posI: %d\tOffdiagonal: %lf\n", n, pair[n].posI, offdiagonal); */
/* 		    //	  fprintf(stdout, "pair[%d].posJ: %d\n", n, pair[n].posJ); */
/* 		    //	  printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[n].posI+a) , (3*pair[n].posJ+b), H_p[3*pair[n].posI+a][3*pair[n].posJ+b]); */
/* 		  } */
/* 		else if (a==b) */
/* 		  { */
/* 		    //	  coeff3= coeff0*( (coeffA - coeffB)/dis_Mobile_ij ); */
/* 		    double ondiagonal = coeff0*(coeff1*(subANDdiv - 1) + ((-coeffA*subt1 + coeffB*subt2)/SIGMA_SQRD  )*subANDdiv); */
/* 		    H_p[3*pair[n].posI+a][3*pair[n].posJ+b]=ondiagonal; */
/* 		    H_p[3*pair[n].posJ+b][3*pair[n].posI+a]=ondiagonal; */
/* 		    //	 fprintf(stdout, "pair[%d].posI: %d\tOndiagonal:%lf\n", n, pair[n].posI, ondiagonal); */
/* 		    //	 fprintf(stdout, "pair[%d].posJ: %d\n", n, pair[n].posJ); */
/* 		    //	 printf("Row:%d\tColumn:%d\tValue=%lf\n",(3*pair[n].posI+a) , (3*pair[n].posJ+b), H_p[3*pair[n].posI+a][3*pair[n].posJ+b]); */
/* 		  } */
/* 	      } */
/* 	  //---End hessian P computation------- */
/* 	} */
/*     } */

  
/*   for(i=0;i<N;++i) */
/*     for(j=0;j<N;++j) */
/*       if(i!=j)  */
/*  	for(a=0; a<3; ++a)  */
/*  	  for(b=0; b<3; ++b)  */
/*  	    {  */
/* 	      H_p[3*i+a][3*i+b] -= H_p[3*i+a][3*j+b];  */
/*  	      //printf("H_p[%d][%d]=%lf\n",  (3*i+a), (3*i+b), H_p[3*i+a][3*i+b]);  */
/*  	    }  */
/*   // printf("\nCalculated Hessian Matrix\n"); */
/*   if(DEBUGMODE)   */
/*     { */
/*       FILE *f=fopen("gradP_modelComb.txt", "w"); */
/*       if(f==NULL) */
/* 	{ */
/* 	fprintf(stderr, "Can not produce gradPmodelComb.txt file"); */
/* 	exit(EXIT_FAILURE); */
/* 	} */
      
/*       for(n=0; n<3*N; n++) */
/* 	fprintf(f, "\t%.10lf\n", grad[n]); */
/*       fclose(f); */

/*       FILE *matrix=fopen("hessian_P.txt", "w"); */
/*       if (matrix==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "hessian_P.txt file could not be produced\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
/*       //Print only lower triangular part! */
/*       for(i=0; i<N; i++)  */
/* 	for(j=0; j<=i; j++) */
/* 	  for(a=0; a<3; a++)  */
/* 	    for(b=0; b<3; b++)  */
/* 	      if(H_p[3*i+a][3*j+b]!=0.0) */
/* 		fprintf(matrix, "%.10lf\n", H_p[3*i+a][3*j+b]);  */
/*       //	      fprintf(matrix, "%d\t%d\t%.10lf\n", 3*i+a, 3*j+b, H_p[3*i+a][3*j+b]);  */
/*       fclose(matrix); */
/*     } */

/*   //hessianpart */
/* } */
/* void gradientANDhessian_E_p_I(int N/\*System size*\/, double **H_E_p, double *grad_E_p, */
/* 			      CAcoord *atomMoving,  CAcoord *atomsFixed2, */
/* 			      int *pair1No,         pairInfoNew *pair1,  */
/* 			      int *pair2No,         pairInfoNew *pair2,  */
/* 			      int *pairInitialNo,   pairInfoNew *pairInitial,  */
/* 			      pdf *target, */
/* 			      double maxBigger,     double scaler) */
/* { */
 
/*   int i=0, j=0, a=0, b=0; */
/*   double binSize=1.0; */
/*   //Zero all values from previous calculation */
/*   for(i=0; i<(3*N); i++) */
/*     { */
/*       for(j=0; j<=i; j++) */
/* 	{ */
/* 	  H_E_p[i][j]=0.0; */
/* 	  H_E_p[j][i]=0.0; */
	  
/* 	} */
/*       grad_E_p[i]=0.0; */
/*     } */
  
/*   /\*   int l=0; *\/ */
/*   /\*   for(i=0; i<N; i++) *\/ */
/*   /\*     for(j=0; j<i; j++) *\/ */
/*   /\*     { *\/ */
/*   /\*       pairInitial[l].dis=dis(atomMoving, i, j); *\/ */
/*   /\*       fprintf(stdout, "Dis[%d]=%lf\n", l, pairInitial[l].dis); *\/ */
/*   /\*       l++; *\/ */
/*   /\*     } *\/ */
  
/*   double binIndex=1.0; */

/*   double *grad_P_model = (double*) calloc ((3*N), sizeof(double)); */
/*   if (grad_P_model==NULL) */
/*     { */
/*       fprintf(stderr, "No memory callocation for grad_P_model calculation\n"); */
/*       exit (EXIT_FAILURE); */
/*     } */
  
/*   double **H_part1 = (double**)malloc(3*N*sizeof(double*)); */
/*   if(H_part1==NULL) */
/*     { */
/*       fprintf(stderr, "Malloc cannot allocate memory for H_part1 array"); */
/*       exit(EXIT_FAILURE); */
/*     } */
/*   for (i=0; i<3*N; i++) */
/*     { */
/*       H_part1[i] = (double*)calloc(3*N,sizeof(double)); */
/*       if(H_part1[i]==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "Calloc cannot allocate memory for H_part1[i] array"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
/*     } */
/*   double P_model=0.0; */
/*   double P_target=0.0; */
/*   double value=0.0; */

/* /\*   FILE *f=fopen("P_t_funct.txt", "a"); *\/ */
/* /\*   if(f==NULL) *\/ */
/* /\*     { *\/ */
/* /\*       fprintf(stderr, "P_t_funct.txt file not produced\n"); *\/ */
/* /\*       exit(EXIT_FAILURE); *\/ */
/* /\*     } *\/ */
     
/*   while(binIndex<(maxBigger + 1.0)) */
/*     { */
/*       P_model=0.0; */
/*       P_target=0.0; */
/*       P_model  = pairDistroGaussian3_I(N, atomMoving,  pair1No, pair1, pairInitialNo, pairInitial, binIndex, binSize); */
/*       //      P_model  = pairDistroGaussian1(N, atomMoving,  pair1No, pair1, binIndex); */

/*       //      P_target = pairDistroGaussian1(N, atomsFixed2, pair2No, pair2, binIndex); */
/*       int test=binIndex-1; */
/*       P_target=target[test].binValue; */
/*       //  fprintf(f,"%lf\n", P_target);       */
/*       //      hessian_P_model (N, H_part1,      atomMoving, pair1No, pair1, binIndex); */
/*       //      gradient_P_model(N, grad_P_model, atomMoving, pair1, pair1No, binIndex); */
/*       gradientANDhessian_P_I(N, grad_P_model, H_part1, atomMoving, pair1, pair1No, pairInitial, pairInitialNo, binIndex, binSize); */
/*       double subt =(P_model-P_target); */
/*       //Logarithmic: double subt =(log(P_model)-log(P_target)); */

/*       double coeff= (scaler/ ( 0.5*P_target*P_target )   ); */
/*       //Change it to P_target instead of P_target^2 */
/*       //      double coeff= (scaler/ ( 0.5*P_target)   ); */

/*       for(i=0; i<3*N; i++) */
/* 	{ */
/* 	  for(j=0; j<3*N; j++) */
/* 	    { */
/* 	      value=0.0; */
/* 	      value=grad_P_model[i]*grad_P_model[j]; */
/* 	      H_E_p[i][j]+= coeff*( subt*H_part1[i][j] +  value ); */
/* 	      //	      H_E_p[i][j]+= coeff*( value ); */
/* 	      //Logarithmic: H_E_p[i][j]+=  (  subt/P_model) *H_part1[i][j] +   ( ( 1- subt  )/( P_model * P_model))*value; */

/* 	      // fprintf(stdout, "%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n", binIndex, P_model, P_target, i, j, H_part1[i][j], H_part2[i][j], H_E_p[i][j]); */
/* 	    } */
/* 	  grad_E_p[i]+= coeff*subt*grad_P_model[i];	   */
/* 	  //Logarithmic: grad_E_p[i]+=0.5*(( subt )/P_model)*grad_P_model[i]; */
/* 	} */
/*       binIndex+=1.0; */
/*     } */
/*   //  fclose(f); */
/*   if(0) */
/*     { */
/*       FILE *pairdata=fopen("gradient_E_p.txt", "w");   */
/*       if(pairdata==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "gradient_E_p.txt file has not been produced\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       for (i=0;i<3*N;i++) */
/* 	fprintf(pairdata, "%d\t%.10lf\n", i, grad_E_p[i]); */
/*       fclose(pairdata); */

/*       FILE *f=fopen("hessian_E_p_old.txt", "w"); */
/*       if(f==NULL) */
/* 	{ */
/* 	  fprintf(stderr, "hessian_E_p.txt file has not been produced\n"); */
/* 	  exit(EXIT_FAILURE); */
/* 	} */
      
/*       for(i=0; i<N; i++)  */
/* 	for(j=0; j<=i; j++) */
/* 	  for(a=0; a<3; a++)  */
/* 	    for(b=0; b<3; b++)  */
/* 	      fprintf(f, "%.10lf\n",   H_E_p[3*i+a][3*j+b]);  */
/*       fclose(f); */
/*     } */
  
/*   fprintf(stdout, "Hessian and gradient of scoring function has been calculated\n"); */

/*   for (i=0; i<3*N; i++) */
/*     free(H_part1[i]); */

/*   free(H_part1); */
/*   free(grad_P_model); */
/* } */

/* void getAllPairInfo(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo) */
/* { */
/*   /\*Purpose: To get the information of pairs so that one can construct pair distribution function and  */
/*     put that information into 'pairInfoNew' structure*\/ */
/*   //Keep in mind that in this function you dont consider missing residues and unmatching residues in two conformations. */
/*   //In order just to get rid of ENM hessian calculation for pairs not in contact, I will just make their k==0 */
/*   int i=0, j=0, l=0; /\*Pair counter*\/ */
/*     double dstnc_ij[1];  */
/*   printf("Amino acid number=%d\n", N); */
/*   for(i=0; i<N; ++i)          */
/*     for(j=0; j<i; ++j) */
/*       { */
/* 	if((isCoordValid(atomFixed, i)==1) && (isCoordValid(atomFixed, j)==1)) */
/* 	  { */
/* 	    if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 1) ) */
/* 	      { */
/* 		if((abs(i-j))==1) pair[l].k=CL_CON;//Since j<i, no need to put abs here!! */
		
/* 		else pair[l].k=FR_CON; */
		
/* 		pair[l].dis= dstnc_ij[0]; */
/* 		pair[l].posI=i; */
/* 		pair[l].posJ=j; */
		
/* 		if(DEBUGMODE)	   */
/* 		  { */
/* 		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j)); */
/* 		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
/* 		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i, j, l, pair[l].k, l, pair[l].dis); */
/* 		  } */
/* 		l++; */
/* 	      } */

/* 	    else if( ( (isContact(atomFixed, i, j, dstnc_ij)) == 0) ) */
/* 	      { */
/* 		pair[l].k=0.0; */
/* 		pair[l].dis= dis(atomFixed, i, j); */
/* 		pair[l].posI=i; */
/* 		pair[l].posJ=j; */
		
/* 		if(DEBUGMODE)	   */
/* 		  { */
/* 		    fprintf(stdout, "dis_%d%d=%lf\n", i, j, dis(atomFixed, i, j)); */
/* 		    fprintf(stdout, "pair[%d][%d]->dis=%lf\n", i, j,  pair[l].dis); */
/* 		    fprintf(stdout, "posI: %d, posJ:%d, pair[%d].k=%lf, pair[%d].dis=%lf\n", i,j, l,  pair[l].k, l, pair[l].dis); */
/* 		  } */
/* 		l++; */
/* 	      } */

/* 	  } */
/*       } */
  
/*   pairNo[0]=l; */
/*   printf("Total number of residue pairs are %d\n", l); */
/* } */


