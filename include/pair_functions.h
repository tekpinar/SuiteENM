#ifndef PAIR_FUNCTIONS_H
#define PAIR_FUNCTIONS_H

//int isCoordValid(CAcoord *atom, int i);
//int isContact(CAcoord *atom, int i, int j, double *dstnc_ij);
void pairDistanceDistro(int N, int *pairNo);
void pairCorrPlot1(int N, int *pairNo, pairInfoNew* pair);
void pairCorrPlot2(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair);
void pairCorrPlot3(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair);
void pairCorrPlot4(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair);
double pairDistroGaussian1(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair, double binIndex/*Initially use 1.0*/);
double pairDistroGaussian2(int N, CAcoord *atomMoving, int *pairNo,  pairInfoNew *pair, double binIndex/*Initially use 1.0*/, double increment, int atomNo, int component);
double pairDistroGaussian3(int N, CAcoord *atom, int *pairNo,  pairInfoNew *pair, 
			   int *pairInitialNo,  pairInfoNew *pairInitial,
			   double binIndex/*Initially use 1.0*/);
double pairDistroGaussian3_I(int N, CAcoord *atomMoving, 
			     int *pairNo,  pairInfoNew *pair, 
			     int *pairInitialNo,  pairInfoNew *pairInitial, 
			     double binIndex, double binSize);

double pairDistroGaussian3_I_water(int N, CAcoord *atomMoving, 
			     int *pairNo,  pairInfoNew *pair, 
			     int *pairInitialNo,  pairInfoNew *pairInitial, 
				   double binIndex, double binSize, double N_pairSum);

double pairDistroGaussian4_I(int N, CAcoord *atomMoving, 
			     int *pairNo,  pairInfoNew *pair, 
			     int *pairInitialNo,  pairInfoNew *pairInitial, 
			     double binIndex, double binSize);

double pairDistroGaussian5_I(int N, CAcoord *atomMoving, 
			     int *pairNo,  pairInfoNew *pair, 
			     int *pairInitialNo,  pairInfoNew *pairInitial, 
			     double binIndex, double binSize);

void plotModelandTarget(int N, CAcoord *atomsFixed1, pairInfoNew *pair1, int *pair1No, 
			       CAcoord *atomsFixed2, pairInfoNew *pair2, int *pair2No, 
			       double increment, int atomNo, int component);

void gradient_P_model    (int N, double * gradP, CAcoord *atomMoving, pairInfoNew *pair, int *pairNo, double binIndex);
void gradient_P_numerical(int N, double * gradP, CAcoord *atomMoving,  pairInfoNew *pair1, int *pair1No, double binIndex);

double E_p(CAcoord *atomMoving, int atomCount1, int *pair1No,  pairInfoNew *pair1,
	   CAcoord *atom2Fixed, int atomCount2, int *pair2No,  pairInfoNew *pair2, double maxBigger);

void gradient_E_p    (int N, double *grad_E_p, CAcoord *atomMoving,  int *pair1No,  pairInfoNew *pair1,
		      CAcoord *atomsFixed2, int *pair2No,  pairInfoNew *pair2, double maxBigger);
void gradient_E_p_log(int N, double *grad_E_p, CAcoord *atomMoving,  int *pair1No,  pairInfoNew *pair1,
		      CAcoord *atomsFixed2, int *pair2No,  pairInfoNew *pair2, double maxBigger);

void gradient_E_p_numerical2(int N, double *grad_E_p,  CAcoord *atomMoving, int *pair1No,  pairInfoNew *pair1,
			     CAcoord *atom2Fixed, int *pair2No,  pairInfoNew *pair2, double increment, double maxBigger);

void hessian_P_model    (int N, double **H_p, CAcoord *atomMoving, int *pairNo, pairInfoNew *pair, double binSize);
void hessian_P_numerical(int N, double **H_p, CAcoord *atomMoving, int *pairNo, pairInfoNew *pair, double binSize, double increment);

void hessian_E_p_numerical (int N, double **H_E_p,
			    CAcoord *atomMoving,  CAcoord *atom2Fixed,
			    int *pair1No, pairInfoNew *pair1, 
			    int *pair2No, pairInfoNew *pair2, 
			    double increment, double maxBigger);
void hessian_E_p           (int N, double **H_E_p,
			    CAcoord *atomMoving, CAcoord *atomsFixed2,
			    int *pair1No, pairInfoNew *pair1, 
			    int *pair2No, pairInfoNew *pair2, 
			    double maxBigger);

void gradientANDhessian_P  (int N, double * grad, double **H_p, CAcoord *atomMoving, 
			    pairInfoNew *pair, int *pairNo, 
			    pairInfoNew *pairInitial, int *pairInitialNo, 
			    double binIndex);

void gradientANDhessian_E_p(int N, double **H_E_p, double *grad_E_p,
			    CAcoord *atomMoving,  CAcoord *atomsFixed2,
			    int *pair1No,         pairInfoNew *pair1, 
			    int *pair2No,         pairInfoNew *pair2, 
			    int *pairInitialNo,   pairInfoNew *pairInitial, 
			    pdf *target,
			    double maxBigger,     double scaler);

/*From now on, I will use integral form of pair distribution function described in Gorba, Tama.*/
/*All functions in integral form will have a '_I' suffix.                                      */
/* double pairDistroGaussian3_I(int N, CAcoord *atomMoving,  */
/* 			     int *pairNo,  pairInfoNew *pair,  */
/* 			     int *pairInitialNo,  pairInfoNew *pairInitial,  */
/* 			     double binIndex, double binSize); */
/* void gradientANDhessian_P_I(int N, double * grad, double **H_p, CAcoord *atomMoving, */
/* 			    pairInfoNew *pair, int *pairNo,  */
/* 			    pairInfoNew *pairInitial, int *pairInitialNo,  */
/* 			    double binIndex, double binSize); */
/* void gradientANDhessian_E_p_I(int N/\*System size*\/, double **H_E_p, double *grad_E_p, */
/* 			      CAcoord *atomMoving,  CAcoord *atomsFixed2, */
/* 			      int *pair1No,         pairInfoNew *pair1,  */
/* 			      int *pair2No,         pairInfoNew *pair2,  */
/* 			      int *pairInitialNo,   pairInfoNew *pairInitial,  */
/* 			      pdf *target, */
/* 			      double maxBigger,     double scaler); */
//void getAllPairInfo(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo);
void getAllPairInfo_v2(pairInfoNew *pair, int N, CAcoord *atomFixed, int *pairNo);
#endif
