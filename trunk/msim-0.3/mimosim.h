
/* $Id: mimosim.h,v 1.4 2004/07/31 00:05:40 miguel Exp $ */

/* main program functions */
int parser(void);
int source(void);
int coder(void);
int modulator(void);
int channel(void);
int chestimator(void);
int vblastreceiver(void);
int mlreceiver(void);
int allocate_mem( int );
void winddown(void);
/* found in debug.c */
void debug1(int iterations, double anorm);
void debug2(int a1);
void debug4(void);
void debug8(int a1, double d1, double d2);
void debug16(void);
void debug32(void);
void debug64(int a1, int a2);
void debug128(int a1, double *pW);
void debug256(int a1);
void debug1024(int a1, int column);
void debug2048(int a1);
void debug4096(int i, int Mr);
/* found in nr_functions.c */
double nr_ran2(long *idum);
double nr_ran2seed(long *idum);
double nr_ran2noise(long *idum);
double nr_ran2Hch(long *idum);
double nr_ran2Hrx(long *idum);
double nr_gasdevnoise(long *idum);
double nr_gasdevHch(long *idum);
double nr_gasdevHrx(long *idum);
double *nr_vector(int nl, int nh);
unsigned int *nr_ivector(int nl, int nh);
double **nr_matrix(int nrl, int nrh, int ncl, int nch);
void nr_free_vector(double *v, int nl);
void nr_free_ivector(unsigned int *v, int nl);
void nr_free_matrix(double **m, int nrl, int ncl);
void nr_free_imatrix(unsigned int **m, int nrl, int ncl);
void nr_svdcmp(double **a, int m, int n, double w[], double **v);
/* found in support.c */
void pinv_vblast_svd(int i);
void pinv_vblast_formula(int i);
void pinv_vblast_qr(int i);
void pinv_ml(double **H);
int argmin(double **G, int a, int b);
void demodulate(double x, double y, int infopointer);
double slice(double a);
/* found in support_ml.c */
void red(int k, int l);
void swap(int k, int kmax);
void lll(void);
void thinqr(void);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void inverse(double **A, double **Ai, int col);
void inverse2(void);
void decode(void);
/* defines */
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#ifdef NR_LICENSED
static double maxarg1,maxarg2;
#define DMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif
