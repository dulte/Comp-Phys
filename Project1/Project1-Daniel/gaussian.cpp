#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <armadillo>
#include "time.h"
#include <string>

/*
Denne filen brukes ved å gi 3 argumenter. Foerste argument bestemmer antall iterasjoner. Velger man 1 her, vil programmet
lage en fil med loesningen for N'en gitt i argument 2. Velger man iterasjon >1 vil programmet lage sine egene N'er definert under
og lage 2 filer med feilene for disse N'ene. Argument 2 er N. Denne betyr bare noe om iterasjon > 1. Argument 3 bestemmer
om programmet skal bruke den generelle eller spesielle (general/specific) algoritmen.

Datafilene lagres i en mappe ved navn "data".
*/

using namespace std;


double known_function(double x);
void gauss_elemination_general(int N, double *a, double* b, double* c, double* u, double* f);
void gauss_elemination_spesific(int N, double *a_inv, double* u, double* f);
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
double calculate_error(double,double);
void compute_solution(int N, double * u, double *exact, bool use_general);


//Main
int main(int argc, char *argv[])
{

  //Div variables
  int N;
  int iterations = 1; // 25; //If this variable is set to 1, you will get the numerical approx. with the N above, else you wil get the log error with different N's
  bool use_general_algorithm; //If we are going to use the general or specific algorithm
  bool give_error;
  ofstream outFile_u;


  //Checks if the arguments have the right
  if (argc < 4){
    cout << "You did not give enough arguments. It needs iterations, type of algorithm and error or plot" << endl;
    exit(1);
  }

  if (atof(argv[1]) < 1){
    cout << "To small a number for iterations, it has to be 1 or larger" << endl;
    exit(1);
  }
  else{
    iterations = atof(argv[1]);
    cout << "Running " << iterations << " iteration(s)" << endl;
  }

  if (strcmp(argv[2], "general") == 0){
    use_general_algorithm = true;
  }
  else if(strcmp(argv[2], "specific") == 0){
    use_general_algorithm = false;
  }
  else{
    cout << "Name of algorithm is not right" << endl;
    exit(0);
  }

  if (strcmp(argv[3], "error") == 0){
    give_error = true;
  }
  else if(strcmp(argv[3], "plot") == 0){
    give_error = false;
  }
  else{
    cout << "Say if you want error or plot" << endl;
    exit(0);
  }


  ofstream outFile_error("logerror.bin");
  ofstream outFile_h("logh.bin");

  //Opens files to write out error
  double* max_errors = new double[iterations];
  double* log_h = new double[iterations];

  for(int j = 0; j < iterations; j++){

      //If we are looking at error, this increase the N every iteration
      if (give_error){
          N = pow(10,1+j/4.0); //Increasing N this way to get better resolution
          log_h[j] = -(1+j/4.0);
          cout << "N: " << N << endl;
      }
      else{
        N = pow(10,1+j); //Increasing N this way to get better resolution
        log_h[j] = -(1+j);
        cout << "N: " << N << endl;
      }






      if (!give_error){
        string filename_numerical_solution = "numericalSolution_" + to_string(N) + ".bin";
        //Open files to write results
        outFile_u.open(filename_numerical_solution);
      }





      //Array that need updating depending on N
      double* u_numericalSolution = new double[N+2];
      double* exactSol = new double[N+2];
      double* error = new double[N+2];

      //Calculates the exact and numerical solution
      compute_solution(N,u_numericalSolution, exactSol, use_general_algorithm);


      for(int k = 0; k < N+1; k++){
          error[k] = calculate_error(u_numericalSolution[k],exactSol[k]);
      }

      max_errors[j] = *max_element(error + 2, error + N-1); //Cuts of the end points, because the the exact solution is 0 and the error is indefined


      if(!give_error){
          writeArrayToFile(outFile_u, u_numericalSolution, N+2); //saves the graph instead of the error
          outFile_u.close();
      }
      //Freeing memory for the next iteration
      delete [] u_numericalSolution;
      delete [] exactSol;
      delete [] error;
  }



  if (give_error){
      writeArrayToFile(outFile_error, max_errors, iterations);
      writeArrayToFile(outFile_h,log_h,iterations);
  }


  outFile_h.close();
  outFile_error.close();


  return 0;

}


//Function we want to integrate
double known_function(double x){
    return 100*exp(-10*x);
}


//Function for doing the calculation of the solution
void compute_solution(int N, double * u, double * exact, bool use_general){

    double dt = 1.0/(N+1);

    double* a = new double[N+2];
    double* b = new double[N+2];
    double* c = new double[N+2];

    double* f = new double[N+2];

    double* a_special_inv = new double[N+2]; //A special version calculated for the algoritem for this special kind of matrix

    //Loop to fill the array needed to solve the problem
    for (int i = 0; i  < N+1; i++){
        a[i] = 2;
        b[i] = -1;
        c[i] = -1;

        f[i] = dt*dt*known_function(dt*i);
        exact[i] = 1 - (1-exp(-10))*(dt*i) - exp(-10*dt*i);

        a_special_inv[i] = (double)i/(i+1);

    }


    //Start the clock when start the gaussian elemination
    clock_t start, finish;
    start = clock();

    if (use_general){
        gauss_elemination_general(N,a,b,c,u,f);
    }
    else{
        gauss_elemination_spesific(N,a_special_inv,u,f);
    }



    finish = clock();
    if (use_general){
    cout << "Time used for N =" << N << " with the general algorithm is: " << float((finish-start))/CLOCKS_PER_SEC << " sec." << endl;
  }
    else{
  cout << "Time used for N =" << N << " with the specific algorithm is: " << float((finish-start))/CLOCKS_PER_SEC << " sec." << endl;
}
    //Delets the arrays to free memory
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] f;
    delete [] a_special_inv;
}



//The general algorithm
void gauss_elemination_general(int N, double *a, double* b, double* c, double* u, double* f) {

    //Forward loop
    for(int i = 2; i < N+1; i++){

        double d = c[i-1]/a[i-1];
        a[i] -= b[i-1]*d;
        f[i] -= f[i-1]*d;
    }

    u[N] = f[N]/a[N];

    //Backward loop
    for(int i = N-1; i > 0; i--){
        u[i] = (f[i] - b[i]*u[i+1])/a[i];
    }

    u[N+1] = 0;
}


//The specific algorithm
void gauss_elemination_spesific(int N, double *a_inv, double* u, double* f){
    for(int i = 2; i < N+1; i++){
        f[i] += f[i-1]*a_inv[i-1];
    }
    u[N] = f[N]*a_inv[N];

    for(int i = N-1; i > 0; i--){
        u[i] = (f[i] + u[i+1])*a_inv[i];
    }
    u[N+1] = 0;
}

//Function for calculating the log of the error
double calculate_error(double computed,double exact){
    if (computed == 0  || exact == 0){
        return 0.0;
    }
    else{
        return log10(abs((computed-exact)/exact));
    }
}


//Writes arrays to Files
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
    outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
