/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include <time.h>
#include "lib_poisson1D.h"

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X,*RHS_2, *RHS_3, *EX_RHS, *RHS_4;
  double **AAB;
  double *AB;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));


  // Using dgbmv
  RHS_2 =(double *) malloc(sizeof(double)*la);
  printf("DGBMV\n");
  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB+1,lab,EX_SOL,1,0.0,RHS_2,1);
  write_vec(RHS_2, &la, "RHS_2.dat");

  // Validation for dgbmv
  relres = relative_forward_error(RHS,RHS_2, &la);
  printf("The relative forward error for dgbmv is relres = %e\n",relres);




  /* LU Factorization */
  
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS_2,&la,&T0,&T1);
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "data/LU_solution_exact.dat");
  /* Validation of LU for tridiagonal matrix (dgbtrf/dgbtrs) */
  relres = relative_forward_error(EX_SOL,RHS_2, &la);
  printf("\nThe relative forward error for dgbtrf/dgbtrs is relres = %e\n",relres);



  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  if (IMPLEM == TRI || IMPLEM == TRF){
    /* Solution (Triangular) */
    if (info==0){
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS_2, &la, &info);
      write_vec(RHS_2, &la, "data/LU_solution.dat");
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}

    }else{
      printf("\n INFO = %d\n",info);
    }
  }
 
  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) {
    RHS_3=(double *) malloc(sizeof(double)*la);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(RHS_3,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv,RHS_3, &la, &info);
    write_xy(RHS_3, X, &la, "data/SOL_direct.dat");

    // Relative forward error for dgbsv
    relres = relative_forward_error(EX_SOL,RHS_3,&la);
    printf("\nThe relative forward error for dgbsv is relres = %e\n",relres);

  }
 

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS_2,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  /* Validation of our LU method for tridiagonal matrix */
  relres = relative_forward_error(EX_SOL,RHS_2, &la);
  printf("\nThe relative forward error for dgbtrftridiag is relres = %e\n",relres);







  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);


  clock_t top;
  
  FILE *direct = fopen("complx/direct.dat", "w");
  FILE *dgbtrf = fopen("complx/dgbtrf.dat", "w");

 

  for(int i = 100; i < 10000; i += 100)
  {
    // dgbtrftridiag only
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
    top = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    fprintf(dgbtrf, "%f ",((double)(clock() - top)/CLOCKS_PER_SEC));

    // dgbtrf only
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
    top = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    fprintf(dgbtrf, "%f\n",((double)(clock() - top)/CLOCKS_PER_SEC));

    // dgbtrftridiag with dgbtrs
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
    top = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info, la);    
    fprintf(direct, "%f ",((double)(clock() - top)/CLOCKS_PER_SEC));
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
   
    // dgbtrf with dgbtrs
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
    top = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info, la);    
    fprintf(direct, "%f ",((double)(clock() - top)/CLOCKS_PER_SEC));
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    // dgbsv only
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
    top = clock();
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info);
    fprintf(direct, "%f\n",((double)(clock() - top)/CLOCKS_PER_SEC));
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  }


  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
