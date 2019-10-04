#include "func.h"
using namespace arma;
using namespace std;



int main(int argc, char** argv){
  if((argc <= 3)||(atoi(argv[1])<=0)){
    cout << "Need 3 arguments and the first argument must be a number larger than 0"<< endl;
    cout << "First argument is the size of dimensions of n x n matrix"<< endl;
    cout << "Third argument is the range for points"<< endl;
    return 1;
  }
  int n=atoi(argv[1]);
  int maxiter=atoi(argv[2]);
  double gernse=atoi(argv[3]);
  double h=(gernse)/((double)n);
  int iterations = 0;
  int iterations2 = 0;

  mat A(n, n); mat R(n,n);
  A.zeros();R.zeros();
  tie(A,R)=Fyll_Matrise_mat(A,R,h,n);
  mat D=A;

  int p; int q;
  clock_t start, finnish;
  double time = 0; double max_verdi = 1.0;
  double tol =  1.0e-9;

  cout << "-----------------------------------------------------------------------------" << endl;
  if (!Test_Maks()){
    cout << "Couldn't find the correct max nondiagonal values of matrix, program will be terminated. "<< endl;
    return 1;
  }else{
    cout << "Unit test 1: Find correct non diagonal max values of matrix passed." << endl;
  }
  cout << "-----------------------------------------------------------------------------" << endl;

  cout << "Running Jacobi rotational algorithm on "<<n<<" x "<<n<<" matrix" << endl;
  while ((iterations < maxiter) && (max_verdi >tol ) ){
    max_verdi = offdiag(A, &p, &q, n);
    start = clock();
    tie(A, R)=Jacobi_rotate(A,R,p, q, n);
    finnish = clock();
    time += (finnish-start)/(CLOCKS_PER_SEC/1000000);

    iterations++;
  }
  cout << "Algorithm completed, time used: "<< time<< " us, iterations used: " << iterations<<"/"<<maxiter << endl;
double time2 = 0;
clock_t start2, finnish2;
double Amax_verdi=1;
 cout << "-----------------------------------------------------------------------------" << endl;
 cout << "Running Armadillo method on "<<n<<" x "<<n<<" matrix" << endl;
  while ((iterations2 < maxiter) && (Amax_verdi >tol ) ){
    Amax_verdi =Aoffdiag( A, n);
    start2 = clock();
    mat D = D.t()*D;
    finnish2= clock();
    time2 += (finnish2-start2)/(CLOCKS_PER_SEC/1000000);

    iterations2++;
  }
  cout << "Algorithm completed, time used: "<< time2<< " us, iterations used: " << iterations2<<"/"<<maxiter << endl;
  cout << "-----------------------------------------------------------------------------" << endl;
  if(!Test_Orth(R,n)){
    cout << "The eigenvectors is not othrogonal "<< endl;
    return 1;
  }else{
    cout << "Unit test 2: Test for orthogonal eigenvectors passed." << endl;
  }
  cout << "-----------------------------------------------------------------------------" << endl;

  // Generate a symmetric matrix
  vec eigval(n);
  mat eigvec(n,n);
  vec egval_rot(n);
  vec analytisk=GetEigenvalue(n,h);
  egval_rot= eigen_rot_vek(A,n);
  egval_rot = sort(egval_rot);
  eig_sym(eigval, eigvec, D);

  char *filnavn ="Egenverdier.txt";
  Write_To_File(filnavn,eigval, egval_rot,analytisk,n, time,iterations);
  char *filnavn2 ="ArmaStats.txt";
  Write_Time_File(filnavn2, n, time, time2, iterations, iterations2);
  cout << "Running plotter.py" << endl;
  system("/usr/bin/python plotter.py");
  cout << "Run complete" << endl;
  cout << "-----------------------------------------------------------------------------" << endl;
  return 0;
}//main slutt
