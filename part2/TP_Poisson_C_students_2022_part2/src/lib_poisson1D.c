/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

double relative_forward_error(double *x, double *y, int* la){
  // norme de x
  double norm = cblas_dnrm2(*la,x,1);
  cblas_daxpy(*la,-1.0,y,1,x,1);

  // calcul de l'erreur avant 
  return (cblas_dnrm2(*la,x,1)/norm);
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  double pivot = 0;

  for (int i = 1; i < (*la); i++)
  {
    pivot = AB[i*(*lab)+1]*AB[(i-1)*(*lab)+3]; 
    pivot /= AB[(i-1)*(*lab)+2]; 
    AB[i*(*lab)+2] -= pivot; 
    AB[(i-1)*(*lab)+3] /= AB[(i-1)*(*lab)+2];
  }
  return *info;
}

void poisson_1d_csr(int n, int** row_ptr, int** col_idx, double** values) {
    int nnz = 3 * n - 2; 
    *row_ptr = (int*)malloc((n + 1) * sizeof(int));
    *col_idx = (int*)malloc(nnz * sizeof(int));
    *values = (double*)malloc(nnz * sizeof(double));

    // Initialisation des pointeurs
    int* row_ptr_array = *row_ptr;
    int* col_idx_array = *col_idx;
    double* values_array = *values;

    int index = 0;
    for (int i = 0; i < n; ++i) {
        row_ptr_array[i] = index;

        col_idx_array[index] = i;
        values_array[index] = 2.0;
        ++index;

        if (i > 0) {
            col_idx_array[index] = i - 1;
            values_array[index] = -1.0;
            ++index;
        }

        if (i < n - 1) {
            col_idx_array[index] = i + 1;
            values_array[index] = -1.0;
            ++index;
        }
    }

    row_ptr_array[n] = index;  
}

void poisson_1d_csc(int n, int** col_ptr, int** row_idx, double** values) {
    int nnz = 3 * n - 2; 
    *col_ptr = (int*)malloc((n + 1) * sizeof(int));
    *row_idx = (int*)malloc(nnz * sizeof(int));
    *values = (double*)malloc(nnz * sizeof(double));

    int* col_ptr_array = *col_ptr;
    int* row_idx_array = *row_idx;
    double* values_array = *values;

    // Remplissage des valeurs de la matrice de Poisson 1D
    int index = 0;
    for (int i = 0; i < n; ++i) {
        col_ptr_array[i] = index;

        row_idx_array[index] = i;
        values_array[index] = 2.0;
        ++index;

        if (i > 0) {
            row_idx_array[index] = i - 1;
            values_array[index] = -1.0;
            ++index;
        }

        if (i < n - 1) {
            row_idx_array[index] = i + 1;
            values_array[index] = -1.0;
            ++index;
        }
    }

    col_ptr_array[n] = index; 
}
