#include "./source.h"


/* Printing */
void print_dVec(double *vec, int n, char *name)
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
  }
  Rprintf("\n");
}
void print_iVec(int *vec, int n, char *name)
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %d ", vec[i]);
  }
  Rprintf("\n");
}

void print_dMat(double **mat, int nr, int nc, char *name)
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %g ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_iMat(int **mat, int nr, int nc, char *name)
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %d ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

/* Compare */
void compare_dVec(double *v1, double *v2, int n)
{
  int i;
  double maxdiff, diff;

  maxdiff = -1.0;
  for (i=0; i<n; i++) {
    diff = fabs(v1[i] - v2[i]);
    if (diff > maxdiff) maxdiff = diff;
  }
  Rprintf("ddiff=%g\n", maxdiff);

} 

void compare_iVec(int *v1, int *v2, int n)
{
  int i;
  int maxdiff, diff;

  maxdiff = -1;
  for (i=0; i<n; i++) {
    diff = abs(v1[i] - v2[i]);
    if (diff > maxdiff) maxdiff = diff;
  }
  Rprintf("idiff=%d\n", maxdiff);

} 

void compare_iMat(int **v1, int **v2, int nr, int nc)
{
  int i, j;
  int maxdiff, diff;

  maxdiff = -1;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      diff = abs(v1[i][j] - v2[i][j]);
      if (diff > maxdiff) maxdiff = diff;
    }
  }
  Rprintf("idiff=%d\n", maxdiff);

}

/* Memory */
double * dVec_alloc(int n, int initFlag, double initVal)
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

int * iVec_alloc(int n, int initFlag, int initVal)
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: iVec_alloc */

/* Function to allocate a double matrix */
double ** dMat_alloc(int nrow, int ncol, int initFlag, double initVal)
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  if (ncol > 0) {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);
  } else {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = NULL;
  }

  return(mat);

} /* END: dMat_alloc */

/* Function to allocate a integer matrix */
int ** iMat_alloc(int nrow, int ncol, int initFlag, int initVal)
{
  int i, **mat, **ptr;

  mat = (int **) malloc(nrow*sizeof(int *));
  CHECK_MEM(mat);
  if (ncol > 0) {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = iVec_alloc(ncol, initFlag, initVal);
  } else {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = NULL;
  }

  return(mat);

} /* END: iMat_alloc */

/* Function to free a matrix */
void matrix_free(void **x, int n)
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* matrix alg */
void copy_iVec(int *v1, int *v2, int n)
{
  int i, *p1, *p2;
  
  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) *p1 = *p2;

} /* END: copy_iVec */

void copy_dVec(double *v1, double *v2, int n)
{
  int i;
  double *p1, *p2;
  
  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) *p1 = *p2;

} /* END: copy_dVec */

void copy_dMat(double **m1, double **m2, int nr, int nc)
{
  int i, j;
  
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) m1[i][j] = m2[i][j];
  }

} /* END: copy_dMat */

void permute_iVec(int *vec, int n, int *ret)
{
  /*
  To shuffle an array a of n elements (indices 0..n-1):
  for i from n-1 downto 1 do
     j = random integer such that 0 <= j <= i
     exchange a[j] and a[i]
  */

  int i, j, tmp;

  copy_iVec(ret, vec, n);
  
  for (i=n-1; i>0; i--) {
    j = floor(runif(0.0, i+1.0));
    if (j > i) j = i;

    tmp    = ret[i];
    ret[i] = ret[j];
    ret[j] = tmp;
  }

} /* END: permute_iVec */

double dotProd(double *v1, double *v2, int n)
{
  int i;
  double *pd1, *pd2, sum=0.0;

  for (i=0, pd1=v1, pd2=v2; i<n; i++, pd1++, pd2++) sum += *pd1 * *pd2;

  return(sum);

} /* END: dotProd */

double dotProd_di(double *v1, int *v2, int n)
{
  int i, *pd2;
  double *pd1, sum=0.0;

  for (i=0, pd1=v1, pd2=v2; i<n; i++, pd1++, pd2++) sum += *pd1 * *pd2;

  return(sum);

} /* END: dotProd_di */

double dotProd_logdi(double *v1, int *v2, int n)
{
  int i, *pd2;
  double *pd1, sum=0.0;

  for (i=0, pd1=v1, pd2=v2; i<n; i++, pd1++, pd2++) sum += log(*pd1) * *pd2;

  return(sum);

} /* END: dotProd_logdi */


void matTimesVec(double **mat, double *vec, int nr, int nc, double *ret)
{
  int i;
  double **prow, *pret;

  for (i=0, pret=ret, prow=mat; i<nr; i++, pret++, prow++) *pret = dotProd(*prow, vec, nc);

}

void matTimesiVec(double **mat, int *vec, int nr, int nc, double *ret)
{
  int i;
  double **prow, *pret;

  for (i=0, pret=ret, prow=mat; i<nr; i++, pret++, prow++) *pret = dotProd_di(*prow, vec, nc);

}

double sum_dVec(double *v, int n)
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += v[i];
  return(sum);

} 

void set_dVec(double *v, int n, double val)
{
  int i;
  double *p;

  for (i=0, p=v; i<n; i++, p++) *p = val;

} 

int sum_iVec(int *v, int n)
{
  int i;
  int sum=0;

  for (i=0; i<n; i++) sum += v[i];
  return(sum);

} 

double sumLogdVec(double *vec, int n)
{
  int i;
  double sum=0.0, *p;

  for (i=0, p=vec; i<n; i++, p++) sum += log(*p);

  return(sum);

} /* END: sumLogdVec */

void vecProd_di(double *v1, int *v2, int n, double *ret)
{
  int i, *pd2;
  double *pd1, *pret;

  for (i=0, pd1=v1, pd2=v2, pret=ret; i<n; i++, pd1++, pd2++, pret++) *pret = *pd1 * *pd2;

} /* END: vecProd_di */

void matIntoVecByRow(double **mat, int nr, int nc, double *ret)
{
  int i, j;
  double *pd, *prow; 

  pd = ret;
  for (i=0; i<nr; i++) {
    for (j=0, prow=mat[i]; j<nc; j++, prow++) {
      *pd = *prow; 
      pd++;
    } 
  }

}




