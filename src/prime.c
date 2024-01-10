#include "./source.h"


static void M_fn(double *t, double *t_all, int N, int Nfail, double **ret)
{
  /*
  M=matrix(0, nrow=n.failure, ncol=N)
  for (l in 2:(n.failure+1))
    for (i in 1:N)
    {
      M[(l-1),i]=(min(t.all[l], t[i])-t.all[l-1])*ifelse(t[i]>=t.all[l-1], 1 ,0)
    } 
  */


  int i, k, km1; 
  double tmin, tk, tkm1, ti, test, **prow, *pd1, *pt;

  /* Return by row */
  for (k=1, prow=ret; k<Nfail+1; k++, prow++) {
    km1  = k - 1;
    tk   = t_all[k];
    tkm1 = t_all[km1];
    test = tkm1 - EQUAL_EPS;

    for (i=0, pt=t, pd1=*prow; i<N; i++, pt++, pd1++) {
      ti = *pt;
      if (ti >= test) {
        tmin = MIN(tk, ti);
        *pd1 = tmin - tkm1;
      } else {
        *pd1 = 0.0;
      }
    }
  }

} /* END: M_fn */


static void M1_fn(double **Mfn, double *lambda, int n, int nfail, double *ret)
{
  /* as.vector(t(lambda)%*%Mfn)
     Mfn is nfail by n, lambda is length nfail,
     ret length n */
  int i, j;
  double *pd1, *pret, sum;

  for (i=0, pret=ret; i<n; i++, pret++) {
    sum = 0.0;
    for (j=0, pd1=lambda; j<nfail; j++, pd1++) sum += *pd1 * Mfn[j][i]; 
    *pret = sum;
  }

} /* END: M1_fn */

static void Mstar_fn(double *tall, double tstar, int *tallLEtstar, int nfail, double *ret)
{
  /*
  tvec   <- t.all[1:n.failure]
  tmp1   <- tvec <= tstar
  tmp2   <- pmin(rep(tstar, n.failure), t.all[-1]) - tvec
  M      <- tmp1*tmp2
  dim(M) <- c(n.failure, 1)
  */
  int i, *pi1;
  double *pret;

  pi1  = tallLEtstar;
  pret = ret;
  for (i=0; i<nfail; i++, pi1++, pret++) {
    if (*pi1) {
      *pret = MIN(tstar, tall[i+1]) - tall[i];
    } else {
      *pret = 0.0;
    }
  }

} /* END: Mstar_fn */
 
static double Mstar1_fn(double *Mstarfn, double *lambda, int n)
{
  /* t(lambda)%*%Mstarfn */

  int i;
  double sum=0.0, *pd1, *pd2;

  for (i=0, pd1=Mstarfn, pd2=lambda; i<n; i++, pd1++, pd2++) sum += *pd1 * *pd2;

  return(sum);

} /* END: Mstar1_fn */ 

static void M2_fn(double *M1fn, int n, double Mstar1fn, double *ret)
{
  /* M1fn - rep(Mstar1fn, N) */
  int i;
  double *pd1, *pd2;

  pd1 = M1fn;
  pd2 = ret;
  for (i=0; i<n; i++, pd1++, pd2++) *pd2 = *pd1 - Mstar1fn;

} /* END: M2_fn */

static void A_fn(MYSTR *mystr, double *ret)
{
  /*
  M=rep(0, n.failure)
  for (l in 1:n.failure)
  {
    M[l]=sum(ifelse(t>=tstar, 1, 0)*trt*Mstarfn[l]) + t(Mfn[l,])%*%(1-trt*ifelse(t>=tstar, 1, 0))  
  } 
  M
  */

  int *timeGEtstar_trt=mystr->timeGEtstar_trt, sumTimeGEtstar_trt=mystr->sumTimeGEtstar_trt,
      n_failure=mystr->n_failure, N=mystr->N, *pi1, i;
  double *Mstarfn=mystr->Mstarfn, **Mfn=mystr->Mfn, *pret, *pd1;

  /* compute Mfn%*%(1-tmp) above and store in ret */
  matTimesiVec(Mfn, mystr->oneMinusTimeGEtstar_trt, n_failure, N, ret);

  for (i=0, pret=ret, pd1=Mstarfn, pi1=timeGEtstar_trt; i<n_failure; i++, pret++, pd1++, pi1++) {
    *pret = sumTimeGEtstar_trt * *pd1 + *pret; 
  }

} /* END: A_fn */

static void B_fn(MYSTR *mystr, double *ret)
{
  /*
  M=rep(0, n.failure)
  for (l in 1:n.failure)
  {
    M[l]=(t(Mfn[l,]-rep(Mstarfn[l], N))%*%(trt*ifelse(t>=tstar, 1, 0)*(W%*%exp(beta))))
  } 
  M
  */

  int i, j, N=mystr->N;
  double sum, **prow, *pret, *pd1, *pd2, *pd3, d2, *tmpvecN=mystr->tmpvecN;

  /* compute trt*ifelse(t>=tstar, 1, 0)*(W%*%exp(beta)) */
  vecProd_di(mystr->W_exp_beta, mystr->timeGEtstar_trt, N, tmpvecN);

  for (i=0, pret=ret, prow=mystr->Mfn, pd2=mystr->Mstarfn; i<mystr->n_failure; i++, pret++, prow++, pd2++) {
    sum = 0.0;
    d2  = *pd2;
    for (j=0, pd1=*prow, pd3=tmpvecN; j<N; j++, pd1++, pd3++) sum += (*pd1 - d2) * *pd3;
    *pret = sum;
  } 

} /* END: B_fn */

static void C_fn(MYSTR *mystr, double *ret)
{
  /*
    X    = rep(0, (J-1))
    for (j in 1:(J-1))
      {X[j] = sum(W[,j]*trt*ifelse(t>=tstar, 1, 0)*event_status)}
    X
  */
  int i, j, N=mystr->N, *timeGEtstar_trt_event=mystr->timeGEtstar_trt_event, *pi1;
  double **W=mystr->W, sum, *pret;

  for (j=0, pret=ret; j<mystr->J-1; j++, pret++) {
    sum = 0.0;
    for (i=0, pi1=timeGEtstar_trt_event; i<N; i++, pi1++) {
      sum += W[i][j]* *pi1;
    }
    *pret = sum;
  }

} /* END: C_fn */

static void D_fn(MYSTR *mystr, double *ret)
{
  /* 
  X    = rep(0, (J-1))
  for (j in 1:(J-1))
    {X[j] = sum(W[,j]*trt*ifelse(t>=tstar, 1, 0)*M2fn)}
  X
  */

  int i, j, N=mystr->N;
  double sum, *tmpvecN=mystr->tmpvecN, *pd, **W=mystr->W;

  /* compute trt*ifelse(t>=tstar, 1, 0)*M2fn) */
  vecProd_di(mystr->M2fn, mystr->timeGEtstar_trt, N, tmpvecN);

  for (j=0; j<mystr->J-1; j++) {
    sum = 0.0;
    for (i=0, pd=tmpvecN; i<N; i++, pd++) sum += W[i][j]* *pd; 
    ret[j] = sum;
  }

} /* END: D_fn */

static double loglik_fn(MYSTR *mystr, double *ll_marg, double *ll_marg0)
{
  /* 
  r1       = sum(delta.l%*%log(lambda))
  r2       = -sum(ifelse(dat$X >= tstar, 1, 0)*dat$trt)*Mstar1fn
  r3       = -t(M1fn)%*%(1-dat$trt*ifelse(dat$X>=tstar, 1, 0))
  r4       = sum((W%*%beta)*dat$trt*ifelse(dat$X>=tstar, 1, 0)*dat$event_status)
  r5       = -sum((W%*%(exp(beta)))*dat$trt*ifelse(dat$X>=tstar, 1, 0)*M2fn)
  r6       = sum(dat$trt*(W%*%log(P)))
  r7       = sum(dat$trt*replace(rowSums(W*log(W)), is.na(rowSums(W*log(W))), 0))
  
  loglik       = r1 + r2 + r3 + r4 + r5 + r6 
  loglik.marg  = r1 + r2 + r3 + r4 + r5 + r6 - r7 
  loglik.marg0 = r1 + r2 + r3 + r4 + r5      
  */

  int N=mystr->N, *timeGEtstar_trt_event=mystr->timeGEtstar_trt_event, *timeGEtstar_trt=mystr->timeGEtstar_trt,
      *trt=mystr->trt, J=mystr->J, trti, i, j;
  double r1, r2, r3, r4, r5, r6, r7, **W=mystr->W, *W_beta=mystr->W_beta, *W_exp_beta=mystr->W_exp_beta,
         *M2fn=mystr->M2fn, loglik, *logP=mystr->logP, tmp=0.0, *Wrow, Wij;

  /*r1  = sumLogdVec(mystr->lambda, mystr->n_failure);*/
  r1  = dotProd_logdi(mystr->lambda, mystr->delta_l, mystr->n_failure);
  r2  = -(mystr->Mstar1fn)*mystr->sumTimeGEtstar_trt;
  r3  = -dotProd_di(mystr->M1fn, mystr->oneMinusTimeGEtstar_trt, N);
  r4  = 0.0;
  r5  = 0.0;
  r6  = 0.0;
  r7  = 0.0;

  for (i=0; i<N; i++) {
    trti = trt[i];
    r4  += W_beta[i]*timeGEtstar_trt_event[i];
    r5  -= W_exp_beta[i]*timeGEtstar_trt[i]*M2fn[i];
    Wrow = W[i];
    r6  += dotProd(Wrow, logP, J)*trti;

    /* Compute rowSums(W*log(W)) */
    tmp = 0.0;
    for (j=0; j<J; j++) {
      Wij = Wrow[j];
      if (Wij > MINLOGARG) tmp += Wij*log(Wij);
    }
    r7  += tmp*trti;
  }

  loglik    = r1 + r2 + r3 + r4 + r5 + r6;
  *ll_marg  = r1 + r2 + r3 + r4 + r5 + r6 - r7;
  *ll_marg0 = r1 + r2 + r3 + r4 + r5;   
 
/*
Rprintf("%g %g %g %g %g %g %g\n", r1, r2, r3, r4, r5, r6, r7);
*/

  return(loglik);

} /* END: loglik_fn */

static void pdf_r_fn(MYSTR *mystr, double **ret)
{
  /*
  pdf.r    = matrix(0, nrow=N, ncol=(J))
  for (i in 1:N)
    for (j in 1:(J))
    {
      r1        = exp(beta[j]*dat$event_status[i]*ifelse(dat$X[i]>=tstar, 1, 0)*dat$trt[i])
      r2        = exp(-exp(beta[j])*ifelse(dat$X[i]>=tstar, 1, 0)*dat$trt[i]*M2fn[i])
      pdf.r[i,j]  = r1*r2
    }
  pdf.r
  */

  int i, j, J=mystr->J, *timeGEtstar_trt_event=mystr->timeGEtstar_trt_event, N=mystr->n_trt,
      *timeGEtstar_trt=mystr->timeGEtstar_trt, tmpi, *pi1, *pi2;
  double *beta=mystr->beta, *exp_beta=mystr->exp_beta, *M2fn=mystr->M2fn, r1, r2, 
          tmpd, *pd, *pbeta, *pebeta, *vec, *pvec;

  for (i=0, pi1=timeGEtstar_trt_event, pi2=timeGEtstar_trt, pd=M2fn; i<N; i++, pi1++, pi2++, pd++) {
    tmpi = *pi1;
    tmpd = *pd * *pi2;
    vec  = ret[i]; 
    for (j=0, pbeta=beta, pebeta=exp_beta, pvec=vec; j<J; j++, pbeta++, pebeta++, pvec++) {
      r1    = exp(*pbeta * tmpi);
      r2    = exp(-*pebeta * tmpd);
      *pvec = r1*r2;
    }
  }

} /* END: pdf_r_fn */

static void EM_init(MYSTR *mystr) 
{
  int i, n_trt=mystr->n_trt, J=mystr->J, N=mystr->N, Jm1;
  double **W=mystr->W, *P=mystr->P, *pd;

  mystr->loglik      = 0.0;
  mystr->loglik_marg = 0.0;
  mystr->conv        = 0;
  mystr->conv_niter  = -1;
  Jm1                = J - 1;

  /* treated subs are first */
  if (mystr->init_W) {
    for (i=0; i<N; i++) {
      pd = W[i]; 
      if (i < n_trt) {
        copy_dVec(pd, P, J);
      } else {
        set_dVec(pd, J, 0.0);
        pd[Jm1] = 1.0;   
      }
    }
  }
  copy_dVec(mystr->beta,   mystr->beta0,   J);
  copy_dVec(mystr->lambda, mystr->lambda0, mystr->n_failure);

}

static void update_beta_obj(MYSTR *mystr)
{

  double *beta=mystr->beta, **W=mystr->W, *exp_beta=mystr->exp_beta, *W_beta=mystr->W_beta,
         *W_exp_beta=mystr->W_exp_beta;
  int i, N=mystr->N, J=mystr->J;

  for (i=0; i<J; i++) exp_beta[i] = exp(beta[i]);
  matTimesVec(W, beta, N, J, W_beta);
  matTimesVec(W, exp_beta, N, J, W_exp_beta);

}

static int EM_main1(MYSTR *mystr)
{
  double loglik, loglik0, *time, *t_all, **Mfn, tstar, *Mstarfn, *lambda, *M1fn, *M2fn,
         **W, tmpd, *dvec, *Cfn, *Dfn, *beta, *Afn, *Bfn, ll_marg, ll_marg0, eps, *P;
  int iter, N, J, nfail, *tallLEtstar, i, j, n_trt, *delta_l, conv, print, badBetaFlag, DEBUG;

  DEBUG = mystr->DEBUG;
  EM_init(mystr); 

  N           = mystr->N;
  J           = mystr->J;
  nfail       = mystr->n_failure;
  time        = mystr->time;
  t_all       = mystr->t_all;
  Mfn         = mystr->Mfn;
  tstar       = mystr->tstar;
  tallLEtstar = mystr->tallLEtstar;
  Mstarfn     = mystr->Mstarfn;
  lambda      = mystr->lambda;
  M1fn        = mystr->M1fn;
  M2fn        = mystr->M2fn;
  iter        = 1;
  W           = mystr->W;
  n_trt       = mystr->n_trt;
  Cfn         = mystr->Cfn;
  Dfn         = mystr->Dfn;
  delta_l     = mystr->delta_l;
  beta        = mystr->beta;
  Afn         = mystr->Afn;
  Bfn         = mystr->Bfn;
  print       = mystr->print;
  eps         = mystr->eps_EM;
  P           = mystr->P;
  badBetaFlag = 0;
  loglik      = -9999.0;

  M_fn(time, t_all, N, nfail, Mfn);
  Mstar_fn(t_all, tstar, tallLEtstar, nfail, Mstarfn);
  M1_fn(Mfn, lambda, N, nfail, M1fn);
  mystr->Mstar1fn = Mstar1_fn(Mstarfn, lambda, nfail);
  M2_fn(M1fn, N, mystr->Mstar1fn, M2fn);

  if (DEBUG) {
    print_dVec(Mstarfn, nfail, "Mstarfn");
    print_dVec(M1fn, N, "M1fn");
    Rprintf("Mstar1fn = %g\n", mystr->Mstar1fn);
    print_dVec(M2fn, N, "M2fn");
  }

  beta[J-1]  = 0.0;
  update_beta_obj(mystr);

  loglik0  = loglik_fn(mystr, &ll_marg, &ll_marg0);
  conv     = 0;

  if (print) Rprintf("Initial loglike = %g\n", loglik0);
  if (DEBUG) Rprintf("ll.marg = %g, ll.marg.0 = %g\n", ll_marg, ll_marg0);

  while(iter <= mystr->num_rep_EM) {
    iter++;

    pdf_r_fn(mystr, W);

    for (i=0; i<n_trt; i++) {
      dvec = W[i];
      tmpd = dotProd(P, dvec, J);
      for (j=0; j<J; j++) {
        dvec[j] = P[j]*dvec[j]/tmpd;
      }
    }

    C_fn(mystr, Cfn);
    D_fn(mystr, Dfn);
    for (j=0; j<J-1; j++) {
      tmpd = log(Cfn[j]/Dfn[j]);
      if (!R_FINITE(tmpd)) {
        badBetaFlag = 1;
        break;
      }
      beta[j] = tmpd;
    }

    if (badBetaFlag) {
      if (print) Rprintf("Non-finite value for beta in iteration %d\n", iter);
      if (DEBUG) {
          print_dVec(Cfn, J-1, "Cfn");
          print_dVec(Dfn, J-1, "Dfn");
      }
      break;
    }

    update_beta_obj(mystr);
    A_fn(mystr, Afn);
    B_fn(mystr, Bfn);
    for (j=0; j<nfail; j++) lambda[j] = delta_l[j]/(Afn[j] + Bfn[j]);
    M1_fn(Mfn, lambda, N, nfail, M1fn);
    mystr->Mstar1fn = Mstar1_fn(Mstarfn, lambda, nfail);
    M2_fn(M1fn, N, mystr->Mstar1fn, M2fn);

    if (DEBUG) {
      print_dVec(Afn, nfail, "Afn");
      print_dVec(Bfn, nfail, "Bfn");
      print_dVec(lambda, nfail, "lambda");
      print_dVec(M1fn, N, "M1fn");
      Rprintf("Mstar1fn = %g\n", mystr->Mstar1fn);
      print_dVec(M2fn, N, "M2fn");
    }

    loglik = loglik_fn(mystr, &ll_marg, &ll_marg0);
    tmpd   = fabs(loglik - loglik0);
    if (tmpd <= eps) {  
      conv = 1;
      break;
    }

    if (print > 1) {
      if (iter > 1) {
        Rprintf("Iter=%d, loglike=%g, diff=%g\n", iter, loglik, tmpd);
      } else {
        Rprintf("Iter=%d, loglike=%g\n", iter, loglik);
      }     
    }  
    loglik0 = loglik;
  } 


  if (print) {
    if (conv) {
      Rprintf("EM algorithm converged in %d iteration(s), final loglike = %g\n", iter, loglik);
    } else {
      Rprintf("EM algorithm did not converge\n");
    }
  }   

  mystr->conv         = conv;
  mystr->conv_niter   = iter;
  mystr->loglik       = loglik;
  mystr->loglik_marg  = ll_marg;
  mystr->loglik_marg0 = ll_marg0;

  return(conv);

} /* END: EM_main1 */

static void EM_max_save(MYSTR *mystr)
{
  MAXEM* maxem=mystr->maxem;
  
  maxem->conv         = mystr->conv;
  maxem->conv_niter   = mystr->conv_niter;
  maxem->loglik       = mystr->loglik;
  maxem->loglik_marg  = mystr->loglik_marg;
  maxem->loglik_marg0 = mystr->loglik_marg0;
  copy_dVec(maxem->beta, mystr->beta, mystr->J);
  copy_dMat(maxem->W, mystr->W, mystr->N, mystr->J);

}

static void EM_max_copy(MYSTR *mystr)
{
  MAXEM* maxem=mystr->maxem;
  
  mystr->conv         = maxem->conv;
  mystr->conv_niter   = maxem->conv_niter;
  mystr->loglik       = maxem->loglik;
  mystr->loglik_marg  = maxem->loglik_marg;
  mystr->loglik_marg0 = maxem->loglik_marg0;
  copy_dVec(mystr->beta, maxem->beta, mystr->J);
  copy_dMat(mystr->W, maxem->W, mystr->N, mystr->J);

}

static void EM_max_init(MYSTR *mystr)
{
  MAXEM* maxem;

  /* mystr must be set first */
  if (!mystr->beta0) error("beta0 is not set in mystr");
  if (mystr->betaMat_nr < 2) error("mystr->betaMat_nr < 2");

  mystr->maxem = malloc(sizeof(MAXEM));
  CHECK_MEM(mystr->maxem);
  maxem = mystr->maxem;
  
  maxem->betaMat      = mystr->beta0; /* beta vector passed in from R */
  maxem->conv         = 0;
  maxem->conv_niter   = 0;
  maxem->loglik       = -1.0;
  maxem->loglik_marg  = -1.0;
  maxem->loglik_marg0 = -1.0;
  maxem->beta         = dVec_alloc(mystr->J, 0, 0.0);
  maxem->W            = dMat_alloc(mystr->N, mystr->J, 0, 0.0);

}

static void EM_max_free(MYSTR *mystr)
{
  MAXEM* maxem=mystr->maxem;
  
  if (maxem->beta) free(maxem->beta);
  if (maxem->W) matrix_free((void **) maxem->W, mystr->N);
  free(maxem);
  mystr->maxem = NULL;

}

static int EM_main(MYSTR *mystr)
{
  MAXEM *maxem;
  int i, ninit=mystr->betaMat_nr, conv=-1, J, firstFlag;
  double *betaMat;

  if (ninit < 2) {
    conv = EM_main1(mystr);
  } else {
    maxem     = mystr->maxem;
    betaMat   = maxem->betaMat;
    J         = mystr->J;
    firstFlag = 0;
    for (i=0; i<ninit; i++) {
      /* Copy initial betas into beta0 */
      copy_dVec(mystr->beta0, &betaMat[i*J], J);
      conv = EM_main1(mystr);
      if (conv) {
        if (!firstFlag) {
          /* save these results */
          firstFlag = 1;
          EM_max_save(mystr);
        } else {
          /* save if we have a larger likelihood */
          if (mystr->loglik > maxem->loglik + LOGLIKE_FUZZ) EM_max_save(mystr);
        }
      }
    }
    if (firstFlag) {
      conv = 1;
      /* Copy results at the max likelihood */
      EM_max_copy(mystr);
    } else {
      conv = 0;
    }
  }

  return(conv);

} /* END: EM_main */

static void mystr_null(MYSTR *mystr)
{
  mystr->DEBUG                   = -1;
  mystr->FUNCTION                = -1;    
  mystr->J                       = -1;                  
  mystr->N                       = -1;          
  mystr->n_failure               = -1;
  mystr->num_rep_EM              = -1;
  mystr->num_rand                = -1;
  mystr->num_rep                 = -1;
  mystr->n_trt                   = -1;
  mystr->conv                    = -1;
  mystr->conv_niter              = -1;
  mystr->trt                     = NULL; 
  mystr->event_status            = NULL; 
  mystr->timeGEtstar_trt         = NULL;
  mystr->oneMinusTimeGEtstar_trt = NULL;
  mystr->sumTimeGEtstar_trt      = -1;
  mystr->timeGEtstar_trt_event   = NULL;
  mystr->tallLEtstar             = NULL; 
  mystr->tstar                   = -1.0;
  mystr->eps_EM                  = -1.0;
  mystr->time                    = NULL;        
  mystr->t_all                   = NULL;       
  mystr->lambda0                 = NULL;
  mystr->lambda                  = NULL;
  mystr->W                       = NULL;          
  mystr->beta0                   = NULL;     
  mystr->beta                    = NULL; 
  mystr->exp_beta                = NULL;
  mystr->W_beta                  = NULL;      
  mystr->W_exp_beta              = NULL;  
  mystr->P                       = NULL;
  mystr->logP                    = NULL;      
  mystr->Mfn                     = NULL;        
  mystr->Mstarfn                 = NULL;  
  mystr->M1fn                    = NULL;  
  mystr->Mstar1fn                = -1.0;
  mystr->M2fn                    = NULL;     
  mystr->Afn                     = NULL;     
  mystr->Bfn                     = NULL;       
  mystr->Cfn                     = NULL;     
  mystr->Dfn                     = NULL;       
  mystr->tmpvecN                 = NULL;
  mystr->loglik                  = -1.0;
  mystr->loglik_marg             = -1.0;
  mystr->loglik_marg0            = -1.0;
  mystr->init_W                  = 1;     /* Initialize W matrix in call to EM by default */

  /* For multiple set of initial betas */
  mystr->betaMat_nr              = -1;
  mystr->maxem                   = NULL;
  
}

static void myrand_null(MYRAND *myrand)
{

  myrand->trt          = NULL; 
  myrand->event_status = NULL;
  myrand->time         = NULL;
  myrand->W            = NULL; 
  myrand->trtRand      = NULL;

  myrand->sumNrand     = 0;
  myrand->p_loglikm    = -9999.0;  
  myrand->p_LRT        = -9999.0;

}

static void set_trt_obj(MYSTR *mystr)
{
  int i, *ivec1, *ivec2, N=mystr->N, *event=mystr->event_status;
  double *dvec1, tstar;

  ivec1 = mystr->timeGEtstar_trt;
  ivec2 = mystr->trt;
  dvec1 = mystr->time;
  tstar = mystr->tstar;
  for (i=0; i<N; i++) {
    if (dvec1[i] >= tstar) {
      ivec1[i] = ivec2[i];
    } else {
      ivec1[i] = 0;
    }
  }
  ivec1 = mystr->timeGEtstar_trt;
  ivec2 = mystr->oneMinusTimeGEtstar_trt;
  for (i=0; i<N; i++) ivec2[i] = 1 - ivec1[i];
  ivec1 = mystr->timeGEtstar_trt;
  ivec2 = mystr->timeGEtstar_trt_event;
  for (i=0; i<N; i++) ivec2[i] = ivec1[i]*event[i];
  mystr->sumTimeGEtstar_trt = sum_iVec(mystr->timeGEtstar_trt, N);

}

static void mystr_init(MYSTR *mystr, int *iargs, double *dargs, double *time, int *event, 
      int *trt, double *lambda, double *t_all, double *beta, double *P, int *delta_l,
      double *ret_lambda, double *ret_beta)
{
  int debug, N, J, nfail, i, *ivec1;
  double tstar, *dvec1;

  mystr_null(mystr);
 
  debug               = iargs[IARG_DEBUG];
  if (debug) Rprintf("Begin: mystr_init\n");

  mystr->DEBUG        = debug;
  mystr->N            = iargs[IARG_N];
  mystr->J            = iargs[IARG_J];
  mystr->num_rep_EM   = iargs[IARG_num_rep_EM];
  mystr->num_rand     = iargs[IARG_num_rand];
  mystr->num_rep      = iargs[IARG_num_rep];
  mystr->n_failure    = iargs[IARG_n_failure];
  mystr->FUNCTION     = iargs[IARG_FUNCTION];
  mystr->n_trt        = iargs[IARG_n_trt];
  mystr->print        = iargs[IARG_print];
  mystr->betaMat_nr   = iargs[IARG_betaMat_nr];

  mystr->conv         = 0;
  mystr->conv_niter   = -1;
  mystr->trt          = trt;
  mystr->event_status = event;
  mystr->delta_l      = delta_l;

  mystr->tstar      = dargs[DARG_tstar];
  mystr->eps_EM     = dargs[DARG_eps_EM];
  mystr->time       = time;
  mystr->t_all      = t_all;
  mystr->lambda0    = lambda;
  mystr->beta0      = beta;
  mystr->lambda     = ret_lambda;
  mystr->beta       = ret_beta;
  mystr->P          = P;

  N                 = mystr->N;
  J                 = mystr->J;
  nfail             = mystr->n_failure;
  
  mystr->timeGEtstar_trt         = iVec_alloc(N, 0, 0);
  mystr->oneMinusTimeGEtstar_trt = iVec_alloc(N, 0, 0);
  mystr->timeGEtstar_trt_event   = iVec_alloc(N, 0, 0);
  mystr->tallLEtstar             = iVec_alloc(nfail+1, 0, 0);
  mystr->W                       = dMat_alloc(N, J, 1, 0.0);
  mystr->Mfn                     = dMat_alloc(nfail+1, N, 0, 0.0);
  mystr->exp_beta                = dVec_alloc(J, 0, 0);
  mystr->W_beta                  = dVec_alloc(N, 0, 0);
  mystr->W_exp_beta              = dVec_alloc(N, 0, 0);
  mystr->logP                    = dVec_alloc(J, 0, 0);
  mystr->Mstarfn                 = dVec_alloc(nfail, 0, 0);
  mystr->M1fn                    = dVec_alloc(N, 0, 0);
  mystr->M2fn                    = dVec_alloc(N, 0, 0);
  mystr->Afn                     = dVec_alloc(nfail, 0, 0);
  mystr->Bfn                     = dVec_alloc(nfail, 0, 0);
  mystr->Cfn                     = dVec_alloc(J-1, 0, 0);
  mystr->Dfn                     = dVec_alloc(J-1, 0, 0);
  mystr->tmpvecN                 = dVec_alloc(N, 0, 0);

  set_trt_obj(mystr);

  ivec1 = mystr->tallLEtstar;
  dvec1 = mystr->t_all;
  tstar = mystr->tstar;
  for (i=0; i<nfail+1; i++) {
    if (dvec1[i] <= tstar) {
      ivec1[i] = 1;
    } else {
      ivec1[i] = 0; 
    }
  }

  dvec1 = mystr->exp_beta;
  for (i=0; i<J; i++) dvec1[i] = exp(beta[i]);
  dvec1 = mystr->logP;
  for (i=0; i<J; i++) dvec1[i] = log(P[i]);

  /* For multiple sets of initial betas */
  if (mystr->betaMat_nr > 1) EM_max_init(mystr);

  if (debug) Rprintf("End: mystr_init\n");

}

static void myrand_init(MYSTR *mystr, MYRAND *myrand)
{
  int N=mystr->N, J=mystr->J;

  myrand_null(myrand);

  /* Observed data */
  myrand->trt          = iVec_alloc(N, 0, 0);
  myrand->event_status = iVec_alloc(N, 0, 0);
  myrand->time         = dVec_alloc(N, 0, 0);
  myrand->W            = dMat_alloc(N, J, 0, 0);

  /* Permuted treatment */
  myrand->trtRand      = iVec_alloc(N, 0, 0);

  copy_iVec(myrand->trt, mystr->trt, N);
  copy_iVec(myrand->event_status, mystr->event_status, N);
  copy_dVec(myrand->time, mystr->time, N);
  copy_dMat(myrand->W, mystr->W, N, J);

}

static void mystr_free(MYSTR *mystr)
{

  if (mystr->DEBUG) Rprintf("Begin: mystr_free\n");

  if (mystr->timeGEtstar_trt)         free(mystr->timeGEtstar_trt);    
  if (mystr->oneMinusTimeGEtstar_trt) free(mystr->oneMinusTimeGEtstar_trt);
  if (mystr->timeGEtstar_trt_event)   free(mystr->timeGEtstar_trt_event); 
  if (mystr->tallLEtstar)             free(mystr->tallLEtstar); 
  if (mystr->W)                       matrix_free((void **) mystr->W, mystr->N);
  if (mystr->Mfn)                     matrix_free((void **) mystr->Mfn, mystr->n_failure+1);
  if (mystr->exp_beta)                free(mystr->exp_beta);    
  if (mystr->W_beta)                  free(mystr->W_beta);
  if (mystr->W_exp_beta)              free(mystr->W_exp_beta);
  if (mystr->logP)                    free(mystr->logP);
  if (mystr->Mstarfn)                 free(mystr->Mstarfn);
  if (mystr->M1fn)                    free(mystr->M1fn);
  if (mystr->M2fn)                    free(mystr->M2fn);
  if (mystr->Afn)                     free(mystr->Afn);
  if (mystr->Bfn)                     free(mystr->Bfn);
  if (mystr->Cfn)                     free(mystr->Cfn);
  if (mystr->Dfn)                     free(mystr->Dfn);
  if (mystr->tmpvecN)                 free(mystr->tmpvecN);

  if (mystr->maxem)                   EM_max_free(mystr);

  if (mystr->DEBUG) Rprintf("End: mystr_free\n");

} 
 
static void myrand_free(MYSTR *mystr, MYRAND *myrand)
{
  if (mystr->DEBUG) Rprintf("Begin: myrand_free\n");

  if (myrand->trt) free(myrand->trt);
  if (myrand->event_status) free(myrand->event_status);
  if (myrand->time) free(myrand->time);
  if (myrand->W) matrix_free((void **) myrand->W, mystr->N);
  if (myrand->trtRand) free(myrand->trtRand);

  if (mystr->DEBUG) Rprintf("End: myrand_free\n");
}


void C_EM(int *iargs, double *dargs, double *time, int *event, int *trt, double *t_all, 
   double *lambda0, int *delta_l, double *beta0, double *P, 
   int *ret_conv, double *ret_lambda, double *ret_beta, double *ret_W, double *ret_loglik)
{
  MYSTR mystr;

  mystr_init(&mystr, iargs, dargs, time, event, trt, lambda0, t_all, beta0, P, delta_l,
             ret_lambda, ret_beta);

  *ret_conv     = EM_main(&mystr);
  ret_loglik[0] = mystr.loglik;
  ret_loglik[1] = mystr.loglik_marg;
  ret_loglik[2] = mystr.loglik_marg0;
  matIntoVecByRow(mystr.W, mystr.N, mystr.J, ret_W);  

  mystr_free(&mystr);

  return;

} /* END:  C_EM */

static void B_fn_marg(MYSTR *mystr, double *ret)
{
  /*
  W = beta = 0

  M=rep(0, n.failure)
  for (l in 1:n.failure)
  {
    M[l]=(t(Mfn[l,]-rep(Mstarfn[l], N))%*%(trt*ifelse(t>=tstar, 1, 0)*(W%*%exp(beta))))
  } 
  M
  */

  int i, j, N=mystr->N, *pi1;
  double sum, **prow, *pret, *pd1, *pd2, d2;

  for (i=0, pret=ret, prow=mystr->Mfn, pd2=mystr->Mstarfn; i<mystr->n_failure; i++, pret++, prow++, pd2++) {
    sum = 0.0;
    d2  = *pd2;
    for (j=0, pd1=*prow, pi1=mystr->timeGEtstar_trt; j<N; j++, pd1++, pi1++) sum += (*pd1 - d2) * *pi1;
    *pret = sum;
  } 

} /* END: B_fn_marg */

static double loglik_fn_marg(MYSTR *mystr, double *lambda)
{
  /* 
  beta = 0, W = 0, P = 0, only need loglik.marg0

  r1       = sum(delta.l%*%log(lambda))
  r2       = -sum(ifelse(dat$X >= tstar, 1, 0)*dat$trt)*Mstar1fn
  r3       = -t(M1fn)%*%(1-dat$trt*ifelse(dat$X>=tstar, 1, 0))
  r4       = sum((W%*%beta)*dat$trt*ifelse(dat$X>=tstar, 1, 0)*dat$event_status) = 0
  r5       = -sum((W%*%(exp(beta)))*dat$trt*ifelse(dat$X>=tstar, 1, 0)*M2fn)
  r6       = sum(dat$trt*(W%*%log(P)))

  loglik.marg0 = r1 + r2 + r3 + r4 + r5      
  */

  int N=mystr->N, *trt=mystr->trt, i;
  double r1, r2, r3, r5, loglik;

  /*r1  = sumLogdVec(lambda, mystr->n_failure);*/
  r1  = dotProd_logdi(lambda, mystr->delta_l, mystr->n_failure);
  r2  = -(mystr->Mstar1fn)*mystr->sumTimeGEtstar_trt;
  r3  = -dotProd_di(mystr->M1fn, mystr->oneMinusTimeGEtstar_trt, N);
  r5  = -dotProd_di(mystr->M2fn, mystr->timeGEtstar_trt, N);
  /*r6  = mystr->logP[mystr->J - 1]*mystr->n_trt;*/

  loglik = r1 + r2 + r3 + r5;   
 
  return(loglik);

} /* END: loglik_fn_marg */


static int LRT_main(MYSTR *mystr, double *LRT, double *ll_marg0_null)
{
  int conv, N, nfail, i, *delta_l, DEBUG=mystr->DEBUG;
  double *t_all, *Afn, *Bfn, *lam0, *M1fn, *Mstarfn, ll_marg, *lambda, ll_max, llmarg_max, llmarg0_max;

  if (DEBUG) {
    Rprintf("Begin LRT_main\n");
    Rprintf("Calling EM_main\n");
  }
  conv = EM_main(mystr);
  if (DEBUG) {
    Rprintf("Finished EM_main, converged = %d\n", conv);
  }
  if (!conv) return(conv);

  N       = mystr->N;
  nfail   = mystr->n_failure;
  t_all   = mystr->t_all;
  Afn     = mystr->Afn;
  Bfn     = mystr->Bfn;
  lam0    = mystr->tmpvecN;  /* Note: B_fn_marg does not use tmpvecN */
  M1fn    = mystr->M1fn;
  Mstarfn = mystr->Mstarfn;
  delta_l = mystr->delta_l;
  lambda  = mystr->lambda;

  if (DEBUG) Rprintf("Begin M_fn\n");
  M_fn(mystr->time, t_all, mystr->N, nfail, mystr->Mfn);
  if (DEBUG) Rprintf("Begin Mstar_fn\n");
  Mstar_fn(t_all, mystr->tstar, mystr->tallLEtstar, nfail, Mstarfn);

  if (DEBUG) Rprintf("Begin M1_fn\n");
  M1_fn(mystr->Mfn, lambda, N, nfail, M1fn);
  if (DEBUG) Rprintf("Begin Mstar1_fn\n");
  mystr->Mstar1fn = Mstar1_fn(Mstarfn, lambda, nfail);
  if (DEBUG) Rprintf("Begin M2_fn\n");
  M2_fn(M1fn, N, mystr->Mstar1fn, mystr->M2fn);

  if (DEBUG) Rprintf("Begin loglik_fn\n");
  ll_max = loglik_fn(mystr, &llmarg_max, &llmarg0_max);

  if (DEBUG) Rprintf("Begin A_fn\n");
  A_fn(mystr, Afn);
  if (DEBUG) Rprintf("Begin B_fn\n");
  B_fn_marg(mystr, Bfn);

  if (DEBUG) Rprintf("Compute lambda.0\n");
  for (i=0; i<nfail; i++) lam0[i] = ((double) delta_l[i])/(Afn[i] + Bfn[i]);

  if (DEBUG) Rprintf("Begin M1_fn\n");
  M1_fn(mystr->Mfn, lam0, N, nfail, M1fn);
  if (DEBUG) Rprintf("Begin Mstar1_fn\n");
  mystr->Mstar1fn = Mstar1_fn(Mstarfn, lam0, nfail);
  if (DEBUG) Rprintf("Begin M2_fn\n");
  M2_fn(M1fn, N, mystr->Mstar1fn, mystr->M2fn);

  if (DEBUG) Rprintf("Begin loglik_fn_marg\n");
  ll_marg        = loglik_fn_marg(mystr, lam0);
  *LRT           = 2.0*(llmarg_max - ll_marg);
  *ll_marg0_null = ll_marg;
 
  return(conv);

}

void C_LRT(int *iargs, double *dargs, double *time, int *event, int *trt, double *t_all,
  double *lambda0, int *delta_l, double *beta0, double *P, 
  int *ret_conv, double *ret_lambda, double *ret_beta, double *ret_W, 
  double *ret_loglik, double *ret_LRT)
{
  MYSTR mystr;
  double ll_marg0_null=-9999.0;

  mystr_init(&mystr, iargs, dargs, time, event, trt, lambda0, t_all, beta0, P, delta_l,
             ret_lambda, ret_beta);

  *ret_conv     = LRT_main(&mystr, ret_LRT, &ll_marg0_null);
  ret_loglik[0] = mystr.loglik;
  ret_loglik[1] = mystr.loglik_marg;
  ret_loglik[2] = ll_marg0_null;
  matIntoVecByRow(mystr.W, mystr.N, mystr.J, ret_W);  

  mystr_free(&mystr);

  return;

} 

static void orderData(MYSTR *mystr, MYRAND *myrand)
{
  int i, N, J, trtIndex, nonTrtIndex, *trtRand, *trt, *event, *event0;
  double *time, *time0, **W, **W0;

  trtIndex    = 0;
  nonTrtIndex = mystr->n_trt; 
  N           = mystr->N;
  J           = mystr->J;
  trtRand     = myrand->trtRand;
  trt         = mystr->trt;
  event       = mystr->event_status;
  time        = mystr->time;
  W           = mystr->W;
  event0      = myrand->event_status;
  time0       = myrand->time;
  W0          = myrand->W;

  for (i=0; i<N; i++) {
    if (trtRand[i]) {
      trt[trtIndex]      = 1;
      time[trtIndex]     = time0[i];
      event[trtIndex]    = event0[i];
      copy_dVec(W[trtIndex], W0[i], J);
      trtIndex++;
    } else {
      trt[nonTrtIndex]      = 0;
      time[nonTrtIndex]     = time0[i];
      event[nonTrtIndex]    = event0[i];
      copy_dVec(W[nonTrtIndex], W0[i], J);
      nonTrtIndex++;
    } 
  }

} 

static int ReRandLRT(MYSTR *mystr, MYRAND *myrand)
{
  int iter, N, sumNrand, conv, LRT_cnt, loglikm_cnt, DEBUG=mystr->DEBUG;
  double LRT, LRT_obs, loglikm, loglikm_obs, tmpd, *LRT_vec;

  if (DEBUG) Rprintf("Begin: ReRandLRT\n");
  N              = mystr->N;
  sumNrand       = 0;
  mystr->init_W  = 1;
  LRT_cnt        = 0;
  loglikm_cnt    = 0;
  LRT_vec        = myrand->LRT_vec;
 
  conv = LRT_main(mystr, &LRT_obs, &tmpd);
  if (!conv) return(sumNrand);
  loglikm_obs    = mystr->loglik_marg;
  
  /* Copy final estimates from observed data to use as initial estimates.
     After myrand_init, myrand contains observed trt, event_status, time, W. 
  */
  copy_dVec(mystr->beta0, mystr->beta, mystr->J);
  copy_dVec(mystr->lambda0, mystr->lambda, mystr->n_failure);
  copy_dMat(myrand->W, mystr->W, N, mystr->J);
  mystr->init_W = 0;      /* Do not initialize W matrix in EM alg, use final obs W */
  mystr->betaMat_nr = 0;  /* Only use one set of beta values (the MLEs from above) */

  for (iter=0; iter<mystr->num_rand; iter++) {
    permute_iVec(myrand->trt, N, myrand->trtRand);
    orderData(mystr, myrand);
    set_trt_obj(mystr);

    conv = LRT_main(mystr, &LRT, &tmpd);
    if (conv) {
      sumNrand++;
      loglikm = mystr->loglik_marg;
      if (LRT >= LRT_obs) LRT_cnt++;
      if (loglikm >= loglikm_obs) loglikm_cnt++;
      LRT_vec[iter] = LRT;
    }
  }

  myrand->sumNrand = sumNrand;
  if (sumNrand) {
    myrand->p_LRT     = ((double) LRT_cnt)/sumNrand;
    myrand->p_loglikm = ((double) loglikm_cnt)/sumNrand;
  }

  if (DEBUG) Rprintf("End: ReRandLRT, sumNrand = %d\n", sumNrand);

  return(sumNrand);
 
} 

void C_ReRandLRT(int *iargs, double *dargs, double *time, int *event, int *trt,
  double *t_all, double *lambda0, int *delta_l, double *beta0, double *P,
  int *ret_sumNrand, double *ret_p_LRT, double *ret_p_loglikm, double *ret_LRT_vec)
{
  MYSTR mystr;
  MYRAND myrand;
  double *ret_lambda, *ret_beta;
  int DEBUG=iargs[IARG_DEBUG];

  if (DEBUG) Rprintf("Begin: C_ReRandLRT\n");
  ret_lambda = dVec_alloc(iargs[IARG_n_failure], 0, 0.0);
  ret_beta   = dVec_alloc(iargs[IARG_J], 0, 0.0);

  mystr_init(&mystr, iargs, dargs, time, event, trt, lambda0, t_all, beta0, P, delta_l,
             ret_lambda, ret_beta);
  myrand_init(&mystr, &myrand);
  myrand.LRT_vec = ret_LRT_vec;

  GetRNGstate();

  *ret_sumNrand  = ReRandLRT(&mystr, &myrand);
  *ret_p_LRT     = myrand.p_LRT;
  *ret_p_loglikm = myrand.p_loglikm;

  PutRNGstate();  

  myrand_free(&mystr, &myrand);
  mystr_free(&mystr);
  free(ret_lambda);
  free(ret_beta);

  return;

} 











