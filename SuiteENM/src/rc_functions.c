#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//Personal header files
#include <structures.h>
#include <pdb_io.h>
#include <rmsd_functions.h>

double reactionCoordinate(int N/*System size*/, CAcoord *protInitialCA, CAcoord *protInterCA, CAcoord *protFinalCA)
{
  double *deltaX_S=(double*)malloc(3*N*sizeof(double));
  if(deltaX_S==NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for coordinates of Domain S\n");
      exit(EXIT_FAILURE);
    }

  double *deltaX_S_obs=(double*)malloc(3*N*sizeof(double));
  if(deltaX_S_obs==NULL)
    {
      fprintf (stderr, "Error: malloc failed in memory allocation for coordinates of Domain S\n");
      exit(EXIT_FAILURE);
    }
  //Lets produce RC_global in a few steps!
  //Step 2 - We can produce displacement vectors now. 
  //deltaX_S_obs is displacement between initial and final conformations. 
  int    i=0, ix3=0;
  double normSqrd=0.0;
  double RC=0.0;

  //Step 1 - We need to superimpose structures at first. 
  double w[N];
  for(i=0; i<N; i++)
    {
      if(protInitialCA[i].selTag==1)
	w[i]=1.0;
      else
	w[i]=0.0;
    }
  double rmsd=superimpose(N,  protInterCA, protInitialCA, w, 1); 
  printf("RMSD of Intermediate: %.3lf\n", rmsd);

  rmsd=superimpose(N,  protFinalCA, protInitialCA, w, 1); 
  printf("RMSD of Final: %.3lf\n", rmsd);

  for(i=0; i<N; i++)
    {
      ix3=3*i;
      if(0)
	fprintf(stdout, "%s\t%.3lf\t%.3lf\t%.3lf\n", protInitialCA[i].residname,  \
		protInitialCA[i].x, protInitialCA[i].y, protInitialCA[i].z);
      //Produce deltaX_S_obs vector. 
      if(protInitialCA[i].selTag==1)
	{
	  deltaX_S_obs[ix3]  = (protInitialCA[i].x - protFinalCA[i].x);
	  deltaX_S_obs[ix3+1]= (protInitialCA[i].y - protFinalCA[i].y);
	  deltaX_S_obs[ix3+2]= (protInitialCA[i].z - protFinalCA[i].z);
	  
	  //Produce deltaX_S vector. 
	  deltaX_S[ix3]      = (protInitialCA[i].x - protInterCA[i].x);
	  deltaX_S[ix3+1]    = (protInitialCA[i].y - protInterCA[i].y);
	  deltaX_S[ix3+2]    = (protInitialCA[i].z - protInterCA[i].z);
	}
      else
	{
	  deltaX_S_obs[ix3]  =0.0;
	  deltaX_S_obs[ix3+1]=0.0; 
	  deltaX_S_obs[ix3+2]=0.0; 
	  
	  //Produce deltaX_S vector. 
	  deltaX_S[ix3]      =0.0; 
	  deltaX_S[ix3+1]    =0.0; 
	  deltaX_S[ix3+2]    =0.0; 
	}
      //Produce squared norm of deltaX_S_obs vector. 
      normSqrd+=( (deltaX_S_obs[ix3]*deltaX_S_obs[ix3]) + \
		  (deltaX_S_obs[ix3+1]*deltaX_S_obs[ix3+1]) + \
		  (deltaX_S_obs[ix3+2]*deltaX_S_obs[ix3+2]) );

    }

  //Apply equation 13 in Tekpinar, Zheng, DOI: 10.1002/prot.22755. 
  for(i=0; i<(3*N); i++)
    {
      RC+=(  (deltaX_S[i]*deltaX_S_obs[i]) / (normSqrd)   );
    }

  free(deltaX_S);
  free(deltaX_S_obs);
  return RC;
}

double fractional_RC(int N/*System size*/, CAcoord *protInitialCA, CAcoord *protInterCA, CAcoord *protFinalCA, int mode)
{
  //Purpose: Apply Equation 6 in Coarse-grained modeling of conformational
  //transitions underlying the processive stepping of myosin V dimer along
  //filamentous actin, Wenjun Zheng, Proteins, 2011. 
  // int mode: If mode==1, it produces reaction coordinate from 0 to 1.
  //           If mode==2, it produces reaction coordinate from -1 to 1.
  int i=0;
  double RC=0.0;
  double w[N];

  for(i=0; i<N; i++)
    {
      if(protInitialCA[i].selTag==1)
	w[i]=1.0;
      else
	w[i]=0.0;
    }
  double rmsd_s1=superimpose(N, protInitialCA, protInterCA, w, 0); 
  printf("RMSD(S,1): %.3lf\n", rmsd_s1);

  double rmsd_s2=superimpose(N, protInterCA, protFinalCA, w, 0); 
  printf("RMSD(S,2): %.3lf\n", rmsd_s2);

  double rmsd_s_obs=superimpose(N, protInitialCA, protFinalCA, w, 0); 
  printf("RMSD(S,Obs): %.3lf\n", rmsd_s_obs);

  if(mode==1)
    RC=0.5*(1.0 + ( ((rmsd_s1*rmsd_s1) - (rmsd_s2*rmsd_s2)) / (rmsd_s_obs*rmsd_s_obs) ) );
  else if(mode==2)
    RC=( ((rmsd_s1*rmsd_s1) - (rmsd_s2*rmsd_s2)) / (rmsd_s_obs*rmsd_s_obs) ) ;
  else
    {
      fprintf(stderr, "ERROR: Unknown mode has been selected. \nMode has to be either 1 or 2!\n");
      fprintf(stderr, "If mode==1, it produces reaction coordinate from 0 to 1.");
      fprintf(stderr, "If mode==2, it produces reaction coordinate from -1 to 1.");
      exit(EXIT_FAILURE);
    }
  return RC;
}
double modifiedFractional_RC(int N/*System size*/, CAcoord *protInitialCA, CAcoord *protInterCA, CAcoord *protFinalCA, int mode)
{
  //Purpose: Apply Equation 6 in Coarse-grained modeling of conformational
  //transitions underlying the processive stepping of myosin V dimer along
  //filamentous actin, Wenjun Zheng, Proteins, 2011. 
  //I modified this reaction coordinate so that initial state will be (-1.0) 
  //instead off 0.0 in original case. 
  // int mode: If mode==1, it produces reaction coordinate from 0 to 1.
  //           If mode==2, it produces reaction coordinate from -1 to 1.

  int i=0;
  double RC=0.0;
  double w[N];

  for(i=0; i<N; i++)
    {
      if(protInitialCA[i].selTag==1)
	w[i]=1.0;
      else
	w[i]=0.0;
    }
  double rmsd_s1=superimpose(N, protInitialCA, protInterCA, w, 0); 
  printf("RMSD(S,1): %.3lf\n", rmsd_s1);

  double rmsd_s2=superimpose(N, protInterCA, protFinalCA, w, 0); 
  printf("RMSD(S,2): %.3lf\n", rmsd_s2);

  double rmsd_s_obs=superimpose(N, protInitialCA, protFinalCA, w, 0); 
  printf("RMSD(S,Obs): %.3lf\n", rmsd_s_obs);

  if(mode==1)
    RC=0.5*(1.0 + ( (rmsd_s1 - rmsd_s2) / (rmsd_s_obs) ) );
  else if(mode==2)
    RC= ( (rmsd_s1 - rmsd_s2) / (rmsd_s_obs) ) ;
  else
    {
      fprintf(stderr, "ERROR: Unknown mode has been selected. \nMode has to be either 1 or 2!\n");
      fprintf(stderr, "If mode==1, it produces reaction coordinate from 0 to 1.");
      fprintf(stderr, "If mode==2, it produces reaction coordinate from -1 to 1.");
      exit(EXIT_FAILURE);
    }

  return RC;
}


