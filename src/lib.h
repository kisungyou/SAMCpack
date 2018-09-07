#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.14159265
//#define RAND_MAX 0xffffffff
float *vector(int nl, int nh);
float **matrix(int nrl,int nrh,int ncl,int nch);
float **convert_matrix(float *a,int nrl,int nrh,int ncl,int nch);
double *dvector(int nl,int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
int *ivector(int nl,int nh);
int **imatrix(int nrl,int nrh,int ncl,int nch);
float **submatrix(float** a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);
void free_vector(float *v,int nl,int nh);
void free_dvector(double *v,int nl,int nh);
void free_ivector(int *v,int nl,int nh);
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch);
void free_convert_matrix(float **b,int nrl,int nrh,int ncl,int nch);
void nrerror(const char* error_text);
#define TINY 1.0e-20;
void ludcmp(double **a,int n, int *indx,  double *d);
#undef TINY
void lubksb(double **a, int n, int *indx, double *b);
double matrix_logdet(double **X, int n);
/* Y=inv(X), return d=log(det(X)) */ 
double matrix_inverse(double **X, double **Y, int n);
/* Y=inv(X), return d=log(det(X)) */
int matrix_inverse_diag(double **X, double **Y, double *diag, int n);
double matrix_trace(double **A,int p);
int matrix_sum(double **A,double **B,double **C,int n,int p);
/* Matrix: A: n by p; B: p by m;  C: n by m */
int matrix_multiply(double **A,double **B,double **C,int n,int p,int m);
int matrix_vector_prod(double **A,double *b,double *d,int n,int p);
double vector_matrix_vector(double *a,double **X, double *b,int m,int n);
void copy_vector(double *a,double *b,int p);                                                                                                                     
void copy_matrix(double **a,double **b,int n,int p);
int choldc(double **a,int n,double** D);
/* calculate log(Gamma(s))  */
double loggamma(double xx);
/* calculate log(k!) */
double logpum(int k);

/* generate the random variable form Gamma(a,b) */
double Rgamma(double a,double b);

/* Generate a random variable from Beta(1,k), where
  the first parameter is 1, the second parameter is b */
double Rbeta(double b);

/* Generate deviates from Dirichlet(a1,a2,\ldots,a_k) */
int RDirichlet(double *w,double *a,int k);

double gasdev();

double Rgasdev(double mean, double variance);

int RNORM(double *x,double *mu,double **Sigma,int p);

int Rwishart(double **B,double df,double **Sigma,int p);

/* calculated the log-density of  z~gamma(a,b) */
double dloggamma(double x,double a,double b);

double dloggauss(double z,double mean,double variance);

double DLOGGAUSS(double *z,double *mean,double **variance, int p);

double Dlogwishart(double **D,double df,double **Sigma,int p);

int uniform_direction(double *d, int n);

int dmaxclass(double *z,int n);
                                                                                                                        
int imaxclass(int *z,int n);

int binary_trans(int k,int l,int *d);
                                                                                                                                         
double logsum(double a,double b);

double maxvector(double *x,int n);

double minvector(double *x,int n);
                                                                                                                                        
double sample_variance(double *x,int n);

/* Return the value ln[Gamma(xx)] for xx>0 */
double gammln(double xx);

#define ITMAX 100
#define EPS 3.0e-7
void gser(double *gamser,double a,double x,double *gln);

#undef ITMAX
#undef EPS

#define ITMAX 100
#define EPS 3.0e-7
void gcf(double *gammcf,double a,double x,double* gln);
#undef ITMAX
#undef EPS


double gammp(double a,double x);

/* Return the CDF of the standard normal distribution */
double Gaussp(double x);

/* return Gamma'(z)/Gamma(z)   */
/* Refer to "mathematics handbook pp.287" */
double diGamma(double z);

double correlation(double *z1,double* z2,int p);

int permut_sample(int *sam, int n);
int random_order(int *x, int n);

/* Generate a subset sample of size M from the set 1:N */
int subset_sample(int *x,int M,int N);

void indexx(int n,double *arrin,int *indx);

