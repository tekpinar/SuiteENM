#ifndef  WAXS_FUNCTIONS_H
#define  WAXS_FUNCTIONS_H
void q_x(int J_max, double q, double q_x_array[J_max]);
void q_y(int J_max, double q, double q_y_array[J_max]);
void q_z(int J_max, double q, double q_z_array[J_max]);

gsl_complex X1_q_tilde_v1(int N/*Number of atoms*/, pdb_v23 *info, int k/*A replacement for j.*/, int J_max, double q, double q_x_array[J_max], double q_y_array[J_max], double q_z_array[J_max], double *f_l);
void X1_q_tilde_allS_j_v1(int totalFrameNum, char framesList[][255], int J_max, double q, double waterWeight, double rho_s, int isSolventOn, double q_x_array[J_max], double q_y_array[J_max], double q_z_array[J_max], gsl_complex **X1_S);
gsl_complex X1_q_tilde_v2(int frameNo, int frm_beg, int frm_end, float *crdX, float *crdY, float *crdZ,	int k/*A replacement for j.*/, int J_max, double q, double q_x_array[J_max], double q_y_array[J_max], double q_z_array[J_max], double *f_l);
//void parseSWAXSexperimental(int n_exp, char *experimentalFile, double *exp_q, double *exp_I_q, double *sigma_q);
#endif
