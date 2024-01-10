#ifndef __SOURCE_H__
#define __SOURCE_H__


/* History: Aug 16 2021 Initial coding
            Jan 19 2022 Add vector to return LRT in re-rand test
           
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define EQUAL_EPS 1e-6
#define ALMOST_ZERO 1e-16
#define MINLOGARG 1.0e-300
#define LOGLIKE_FUZZ 1.0e-6

#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}

#define IARG_DEBUG 0
#define IARG_N 1
#define IARG_J 2
#define IARG_num_rep_EM 3
#define IARG_num_rand 4
#define IARG_num_rep 5
#define IARG_n_failure 6
#define IARG_FUNCTION 7
#define IARG_n_trt 8
#define IARG_print 9
#define IARG_betaMat_nr 10

#define DARG_tstar 0
#define DARG_eps_EM 1

struct maxem {

  double *betaMat; /* betaMat_nr x J matrix of initial beta values stored as a vector. Points to beta passed in from R */
  int conv;
  int conv_niter;
  double loglik;
  double loglik_marg;
  double loglik_marg0;
  double *beta;
  double **W;

};
typedef struct maxem MAXEM;


struct mystr {
  int DEBUG;
  int FUNCTION;    /* 0=EM, 1=... */
  int J;
  int N;           /* number of obs */
  int n_failure;
  int num_rep_EM;
  int num_rand;
  int num_rep;
  int n_trt;
  int conv;
  int conv_niter;
  int print;
  int init_W;  /* flag for initializing W matrix in call to EM, default=1 */

  int *trt; /* treatment vector 0-1 */
  int *event_status;
  int *timeGEtstar_trt;
  int *oneMinusTimeGEtstar_trt;
  int sumTimeGEtstar_trt;
  int *timeGEtstar_trt_event;
  int *tallLEtstar; 
  int *delta_l;  /* passed in, length n_failure */
  int betaMat_nr; /* Number of rows (attempts) of initial beta values */

  double tstar;
  double eps_EM;

  double *time;        /* passed in */
  double *t_all;       /* passed in */
  double *lambda0;      /* passed in */
  double *lambda;      /* Length n_failure */
  double **W;          /* N x J returned from EM */
  double *beta0;        /* passed in */
  double *beta;  
  double *exp_beta;
  double *W_beta;      /* W %*% beta */
  double *W_exp_beta;  /* W %*% exp(beta) */
  double *P;
  double *logP;        /* log(P) */
  double **Mfn;        /* n_failure by N */
  double *Mstarfn;     /* n_failure */
  double *M1fn;        /* N */
  double Mstar1fn;
  double *M2fn;        /* N */
  double *Afn;        /* n_failure */
  double *Bfn;        /* n_failure */
  double *Cfn;        /* J-1 */
  double *Dfn;        /* J-1 */
  double *tmpvecN;     /* N */
  double loglik;
  double loglik_marg;
  double loglik_marg0;
  
  MAXEM *maxem;   /* Structure used for multiple sets of initial betas */
};
typedef struct mystr MYSTR;

struct myrand {

  /* Observed data */
  int *trt; 
  int *event_status;
  double *time;
  double **W; 

  /* Permuted data */
  int *trtRand;

  int sumNrand;
  double p_loglikm;  /* p-value from EM marginal likelihoods */
  double p_LRT;

  double *LRT_vec;   /* Vector to store LRT tests, passed in from R */

};
typedef struct myrand MYRAND;


/* printing */
void compare_dVec(double *, double *, int);
void compare_iVec(int *, int *, int);
void compare_iMat(int **, int **, int, int);
void print_dVec(double *, int, char *);
void print_iVec(int *, int, char *);
void print_dMat(double **, int, int, char *);
void print_iMat(int **, int, int, char *);

/* memory */
double * dVec_alloc(int, int, double);
int * iVec_alloc(int, int, int);
double ** dMat_alloc(int, int, int, double);
int ** iMat_alloc(int, int, int, int);
void matrix_free(void **, int);

/* matrix alg */
void copy_dVec(double *, double *, int);
void copy_iVec(int *, int *, int);
void copy_dMat(double **, double **, int, int);
double sum_dVec(double *, int);
int sum_iVec(int *, int);
void permute_iVec(int *, int, int *);
double dotProd(double *, double *, int);
double dotProd_di(double *, int *, int);
double dotProd_logdi(double *, int *, int);
void matTimesVec(double **, double *, int, int, double *);
void matTimesiVec(double **, int *, int, int, double *);
double sumLogdVec(double *, int);
void vecProd_di(double *, int *, int, double *);
void set_dVec(double *, int, double);
void matIntoVecByRow(double **, int, int, double *);

/*
void C_EM(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *,
          int *, double *, double *, double *, double *);
void C_LRT(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *,
           int *, double *, double *, double *, double *, double *);
void C_ReRandLRT(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *,
                int *, double *, double *, double *);
*/





#endif







