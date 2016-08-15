#ifndef  RC_FUNCTIONS_H
#define  RC_FUNCTIONS_H
double reactionCoordinate(int N/*System size*/, CAcoord *protInitialCA, CAcoord *protInterCA, CAcoord *protFinalCA);
double fractional_RC(int N/*System size*/, CAcoord *protInitialCA, CAcoord *protInterCA, CAcoord *protFinalCA, int mode);
double modifiedFractional_RC(int N/*System size*/, CAcoord *protInitialCA, CAcoord *protInterCA, CAcoord *protFinalCA, int mode);
//void GROEL_RC(int frameNo, int *prot1Counts, CAcoord *prot1CA, CAcoord *protInterCA, int *prot2Counts, CAcoord *prot2CA, char *rc_output_file);
//void hemoglobin_RC(int frameNo, int *prot1Counts, CAcoord *prot1CA, CAcoord *protInterCA, int *prot2Counts, CAcoord *prot2CA, char *rc_output_file, char chainID);
//void hbClamshell_RC(int frameNo, int *prot1Counts, CAcoord *prot1CA, CAcoord *protInterCA, int *prot2Counts, CAcoord *prot2CA, char *rc_output_file);
//void hbClamshell_frRC(int frameNo, int *prot1Counts, CAcoord *prot1CA, CAcoord *protInterCA, int *prot2Counts, CAcoord *prot2CA, char *rc_output_file);
#endif
