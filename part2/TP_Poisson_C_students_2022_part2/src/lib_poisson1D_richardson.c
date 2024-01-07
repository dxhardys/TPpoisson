/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include <math.h>
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
  double max = sin(*la * M_PI_2 * (1.0 / (*la + 1)));
  return 4 * max * max;
}

double eigmin_poisson1D(int *la){
  double min = sin(M_PI_2 * (1.0 / (*la + 1)));
  return 4 * min * min;
}

double richardson_alpha_opt(int *la){
  double opt = eigmax_poisson1D(la) + eigmin_poisson1D(la);
  return (2/opt);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double *rhs_temp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la,RHS,1,rhs_temp,1);

  //  Initialiation de  r^0
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,rhs_temp,1);

  //  norme de rhs_temp
  double check = cblas_dnrm2(*la,rhs_temp,1);
  *nbite=0;
  resvec[*nbite] = check;

  while(check > (*tol) && *nbite < *maxit)
  {
    cblas_dcopy(*la,RHS,1,rhs_temp,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,rhs_temp,1);
    cblas_daxpy(*la,(*alpha_rich),rhs_temp,1,X,1);
    check = cblas_dnrm2(*la,rhs_temp,1);
    (*nbite)++;
    resvec[*nbite] = check;
  }
  free(rhs_temp);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
 int n = *la ;
 int m = *lab ;
 int k = *kv ;
 int ld;
  for(int i = 0; i<n; i++)
  {
    ld = i*m;
    MB[ld+k] = AB[ld+k];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int n = *la ;
  int m = *lab ;
  int k = *kv ;
  int ld;
  for(int i = 0; i<n; i++)
  {
    ld = i*m;
    MB[ld+k] = AB[ld+k];
    MB[ld+k+1] = AB[ld+k+1];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

  double *rhs_tmp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la,RHS,1,rhs_tmp,1);
  //  Initialisation de r^0
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,rhs_tmp,1);
  int info; int NHRS = 1;
  int *ipiv = (int *)calloc(*la, sizeof(int));
  int m = (*ku)-1;
  dgbtrf_(la,la,kl,&m,MB,lab,ipiv,&info);
  //  Norme de rhs_tmp
  double check = cblas_dnrm2(*la,rhs_tmp,1);
  *nbite=0;
  resvec[*nbite] = check;
  
  while((check > (*tol)) && ((*nbite) < (*maxit)))
  {
    cblas_dcopy(*la,RHS,1,rhs_tmp,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,rhs_tmp,1);
    // Résolution du système
    dgbtrs_("N",la,kl,&m,&NHRS,MB,lab,ipiv,rhs_tmp,la,&info,*la);
    // Mise à jour de la solution X
    cblas_daxpy(*la,1.0,rhs_tmp,1,X,1);
    //calcul de la nouvel norme 
    check = cblas_dnrm2(*la,rhs_tmp,1);
    (*nbite)++;
    resvec[*nbite] = check;
  }
}

