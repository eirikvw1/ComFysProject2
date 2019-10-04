#include <armadillo>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <tuple>
#include "time.h"
#define _USE_MATH_DEFINES
using namespace arma;
using namespace std;
#ifndef func
#define func

tuple<mat,mat> Jacobi_rotate ( mat,mat, int, int, int );
double offdiag(mat, int *, int *,int);
void Write_To_File(char *,vec, vec,vec, int,double,int);
tuple<mat,mat> Fyll_Matrise_mat( mat,mat,double,int);
double* Lager_vek(int );
vec eigen_rot_vek(mat,int);
vec GetEigenvalue(int, double);
void WriteVec(double *, int);
void Write_Vec_2(vec, vec,int);
bool Test_Orth(mat, int);
bool Test_Maks();
double Aoffdiag(mat, int);
void Write_Time_File(char *fil_navn, int n, double time, double time_arm,int inter, int iter2);


#endif



/*
Offdiag returns the largest non diagonal elements of the matrix 
*/

double offdiag(mat A, int * p, int * q, int n){
  double max=0;double aij=0;
  for (int i = 0; i < n; ++i){
     for ( int j = 0; j < n; ++j){
            aij = fabs(A(i,j));
            if ( aij > max && i!=j){
              max = aij;*p=i; *q=j;
            }
        }
      }
   return max;
}

/*
AOffdiag returns the largest non diagonal elements of the matrix 
*/
double Aoffdiag(mat A, int n){
  double max=0;double aij=0;
  for (int i = 0; i < n; ++i){
     for ( int j = 0; j < n; ++j){
            aij = fabs(A(i,j));
            if ( aij > max && i!=j){
              max = aij;
            }
        }
      }
   return max;
}

// --------------------------------------------------------------------------------------------------

tuple<mat, mat> Jacobi_rotate ( mat B,mat S, int k, int l, int n ){
  double s, c;
  if ( B(k,l) != 0.0 ) {
    double t, tau;
    tau = (B(l,l) - B(k,k))/(2*B(k,l));
    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }
    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il;double r_ik, r_il;
  a_kk = B(k,k);
  a_ll = B(l,l);
  B(k,k) = c*c*a_kk - 2.0*c*s*B(k,l) + s*s*a_ll;
  B(l,l) = s*s*a_kk + 2.0*c*s*B(k,l) + c*c*a_ll;
  B(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  B(l,k) = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = B(i,k);
      a_il = B(i,l);
      B(i,k) = c*a_ik - s*a_il;
      B(k,i) = B(i,k);
      B(i,l) = c*a_il + s*a_ik;
      B(l,i) = B(i,l);
    }
//And finally the new eigenvectors
  r_ik = S(i,k);
  r_il = S(i,l);
  S(i,k) = c*r_ik - s*r_il;
  S(i,l) = c*r_il + s*r_ik;
  }
  return make_tuple(B, S);
}

// --------------------------------------------------------------------------------------------------


tuple<mat,mat> Fyll_Matrise_mat( mat A ,mat R,double h, int n){
    double e=-1/(h*h),d=2/(h*h), v;

  for (int i = 0; i < n; ++i){
    v+=h;
    for(int j = 0; j<n; j++){
        if(i == j){
          A(i,j)=d+ (v*v);
          R(i,j)=1;
        }
        else if((i == j-1)||(i==j+1)){
          A(i,j)=e;
        }
        else{
          A(i,j)=0;
        }

    }
  }
  return make_tuple(A, R);
}

// --------------------------------------------------------------------------------------------------

vec GetEigenvalue(int n, double h){
  double g=h;
  double d= 2/(h*h);
  double a=-1/(h*h);
  vec vek(n);
  for(int i=0 ;i<n;i++){
     g+=h;
    vek(i)=d+g+ 2*a*cos((((double)i+1)*M_PI)/((double)n));
  }
  return vek;
}

// --------------------------------------------------------------------------------------------------

double* Lager_vek(int n){
  double* Vektor = new double [n];
  for(int i=0;i<n;i++){
    Vektor[i]= 0;
  }
  return Vektor;
}

// --------------------------------------------------------------------------------------------------

void WriteVec(double* Matrix, int n){
  for(int i=0;i < n;i++){
    cout<<Matrix[i]<<"  ";
  }
     cout << endl;
}


// --------------------------------------------------------------------------------------------------

void Write_Vec_2(vec Arma, vec rotate,int n){
  cout<<"Aarma"<<" || "<<" Rotate" <<endl;
  for(int i=0;i < n;i++){
    cout<<Arma(i)<<" || "<<rotate(i)<<endl;
  }
  return;
}
// --------------------------------------------------------------------------------------------------
vec eigen_rot_vek(mat matrise, int n){
  vec vektor(n);
  for(int i = 0;i<n;i++){
    for(int j = 0; j<n;j++){
      if(i==j){
        vektor(i)= matrise(i,j);
      }
    }
  }
  return vektor;
}

void Write_To_File(char *fil_navn, vec Arma, vec Roter,vec Ana,int n, double clk,int iterations){
  ofstream myfile;
  myfile.open(fil_navn);
  myfile.precision(14);

  myfile<<"rotation clock: " << clk << " us" <<endl;
  myfile<<n<<" "<<iterations<<endl;
  for(int j = 0; j<n; j++){
    myfile<<Arma(j)<<" "<<Roter(j)<<"  "<< Ana(j) <<n<<endl;
  }
  myfile.close();
  return;
}
void Write_Time_File(char *fil_navn, int n, double time, double time_arm,int inter, int inter2){
  ofstream myfile(fil_navn, ios_base::app);
  myfile.precision(8);
  myfile<<n<<" "<<time << " "<< time_arm <<" "<< inter << " "<<inter2<<endl;
  myfile.close();
  return;
}



//----- Unit tests --------------------------------------------------------------------

/*
Test_Orth:
Hvis ortonormal så er egenvektor ganget seg selv lik 1 og ganget med alle andre like 0.
*/

bool Test_Orth(mat eigenvec, int n){
  mat EigenT = eigenvec.t();
  double dot1 = 0;double dot2 = 0;
  for(int i = 0;i<n;i++){
    dot1+=EigenT(0,i)*EigenT(1,i);
    dot2+=EigenT(0,i)*EigenT(0,i);
  }  

  return (dot1<1e-9 && (dot2<1.0+1e-9)&&(dot2>1.0-1e-9)); // because fuck double comparison...
  }

/*
Test_Maks:
lager matrise med kjente verdier og sjekker om vi finner den største verdien som ikke er på diagonalen.
*/

bool Test_Maks(){
  int lengde= 4;
  mat Test_mat(lengde,lengde);
  //fylller matrise med stigene verdier +1
  double sum=0;
  for (int i = 0; i < lengde ; ++i){
      for ( int j = 0 ;j < lengde ; ++j){
        sum+= 1.0;
        Test_mat(i,j)=sum;
      }
  }
  double max; int pp; int qq;
  max = offdiag(Test_mat, &pp, &qq,lengde);
  return ((max>15.0-1E-9)&&(max<15.0+1E-9));
}