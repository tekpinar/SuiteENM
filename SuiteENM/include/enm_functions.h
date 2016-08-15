#ifndef  ENM_FUNCTIONS_H
#define  ENM_FUNCTIONS_H
int  isCoordValid(CAcoord *atom, int i);
int  isContact   (CAcoord *atom, int i, int j, double *dstnc_ij);
int  isContact_v2(CAcoord *atom, int i, int j, double *dstnc_ij, double R_cutoff);
void getPairInfo_new(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, int *map1);  
void getPairInfo_v2(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, int *map1);   
void getPairInfo_back(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, int *map1); 
void getAllPairInfo(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo);
void getAllPairInfo_v2(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo);
void getAllPairInfo_v3(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, double R_cutoff);
void getAllNonBondedPairsInfo(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo, double R_cutoff);
void getSSandBondedPairsInfo(char *ssArray, pairInfoNew *pairSS, int N, CAcoord *atomFixed, int *pairSSNo, double *weight, int is_ss_on);
void getAllPairInfo_v4(pairInfoNew *pair, int N/*Number of CA*/, CAcoord *protCA, int *pairNo, double exponent, double coeff, int printDetails/*1:yes, 0:No*/);
double **constructHessian_new(int conformationNo, int N/*Number of residues*/, double **hessian, CAcoord *atomMoving, int *pairNo, pairInfoNew *pair);
double gradientANDhessian_ENM(int conformationNo, int N/*Number of residues*/, 
			      double *grad_ENM, double **H_ENM, 
			      CAcoord *atomMoving, 
			      int *pairNo, pairInfoNew *pair,
			      int *pairInitialNo, pairInfoNew *pairInitial);
double gradientANDhessian_mENM(int conformationNo, int N/*Number of residues*/, 
			       double *grad_ENM, double **H_ENM, 
			       CAcoord *atomMoving, 
			       int *pairNo, pairInfoNew *pair,
			       int *pairInitialNo, pairInfoNew *pairInitial, double exponent);
double gradientANDhessian_Collision(int conformationNo, int N/*Number of residues*/, 
				    double *grad_col, double **H_col, 
				    CAcoord *atomMoving, 
				    int *pairNo, pairInfoNew *pair,
				    int *pairInitialNo, pairInfoNew *pairInitial, 
				    double R_collision);
  
#endif
