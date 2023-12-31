/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){


int k = 0;
int j = 0;
int m = *la ; // colonnes
int n = *lab ; // lignes
int v = *kv ; // index diagonale 

// On parcours les colonne de la matrice 
for(int i = 0; i<m;++i){
  k=i*n;
  if(v >= 0 ){
    for(int j = 0 ; j< v ; ++j){
      AB[k+j] =0.0;
    }
  }

  AB[k+v]=-1.0; // On met les élements de la surdiagonale à -1.0
  AB[k+v+1]= 2.0; // On met les éléments de la diagonale principale a 2.0
  AB[k+v+2]=-1.0; // On met les éléments de la sous diagonale a -1.0

}

// On remplit le reste de la matrice avec des 0
AB[0] = 0.0;
  if(v == 1){
    AB[1]=0 ;
  }
AB[n*m-1]=0.0;

    
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){

int k = 0;
int m = *la ;
int n = *lab ;
int v = *kv ;
  // On parcours les colonnes de la matrice 
  for(int i = 0 ; i<m ;++i){
    k=i*n ;

    if(v >=0){
      for(int j = 0; j < v ; j++){
        AB[k+j] = 0.0 ;
      }
    }
    AB[k+v]=0.0; // On met les élements de la surdiagonale a 0.0
    AB[k+v+1]=1.0; // On met les élement de la diagonale a 1.0
    AB[k+v+2]=0.0; //On met les élement de la sousdiagonale a 0.0
  }
  AB[1] = 0.0; // deuxieme élément de la matrice a 0
  AB[n*m-1] = 0.0 ; //dernier élement a 0
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  
  int m = *la ;
  // On met les élement (entre premier et dernier) du vecteur RHS a 0
  for (int i = 0 ; i < m ;++i){
    RHS[i] = 0.0 ;
  }
  RHS[0] = *BC0 ; //on met la premiere valeur du vecteur RHS a la valeur BC0
  RHS[m-1] = *BC1 ; // on met la derniere valeur du vecteur RHS a la valeur BC1

}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int m = *la ;
  double d = (*BC1)-(*BC0); // calcul du delta 
  //calcul de la solution analytique
  for(int i = 0 ; i < m ;++i){
    EX_SOL[i] = *BC0 + X[i] * d ;
  }
}  

void set_grid_points_1D(double* x, int* la){
  // cas d'erreur
  if (x == NULL || la == NULL || *la <= 0) {
    return;
  }
  int m = *la ;
  double a = 1.0 /(1.0 * m+1);
  for(int i = 0 ; i < m ; ++i){
    x[i] = (i+1)* a ;
  }

}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab){
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
