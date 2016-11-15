#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>
#include <mpi.h>

using namespace std;
using namespace arma;

void output(ofstream & file, vec & expecVal, int Nspin, double temp, double MCsteps);
void metropolis(double Nspins,double MCSteps, double temp, mat & spinMatrix, vec & expecVal, ofstream & outFile,int waitSteps);

inline int periodic(int index, int limit){
    return (index + limit)%(limit);
}

int main(int argc, char *argv[])
{
    int numProcs,rank;
    ofstream outFile;
    ofstream outFileParameters;
    int Nspins,monteCarloSteps;
    double temp;
    mat spinMatrix;
    int waitSteps;



    //  MPI initializations
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //cout << rank << endl;
    //cout << "hei " << stod(argv[1]) << " " << "rank: " << rank << endl;


    if (rank == 0 && argc < 4){
        cout << "You did not give enough arguments. Number of spins, temperature and number of monte carlo iterations." << endl;
        //cout << "To few arguments given!" << endl;

        exit(1);
    }

    if (rank == 0){
        outFileParameters.open("parametersOneTemp.txt");
        outFile.open("dataOneTemp.txt");
        Nspins = stoi(argv[1]);
        monteCarloSteps = stoi(argv[3]);
        temp = stod(argv[2]);

    }

    if (rank == 0){

        outFileParameters << "temp " << temp << "\n";
        outFileParameters << "L " << Nspins << "\n";
        outFileParameters << "monteCarloSteps " << monteCarloSteps << "\n";
        outFileParameters << "coresUsed " << numProcs << "\n";
    }



    MPI_Bcast(&Nspins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&monteCarloSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    spinMatrix = ones(Nspins,Nspins);
    waitSteps = 0;



    double timeSinceLast = MPI_Wtime();

    //return 0;
    double timeStart = MPI_Wtime();
//    if (rank == 0){
//        clock_t begin = clock();
//    }



    // cout << "I am " << rank << " and will do T=" << temp << endl;
    vec localExpectationValues = zeros(5);

    // cout << "Doing some work: " << Nspins << " and " << monteCarloSteps << endl;
    metropolis(Nspins,monteCarloSteps, temp,spinMatrix,localExpectationValues,outFile, waitSteps);

    vec totalExpectationValues = zeros(5);


    for (int i = 0; i < totalExpectationValues.size(); i++){
        MPI_Reduce(&localExpectationValues[i], &totalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
//    if (rank == 0){
//        output(outFile, totalExpectationValues, Nspins, temp, monteCarloSteps*numProcs);
//    }



    if (rank == 0){
        outFile.close();

        double done = double(MPI_Wtime() - timeStart);
        cout << "Program done after: "<< done << " sec" << endl;
        outFileParameters << "time " << done << "\n";

        outFileParameters.close();

    }

    MPI_Finalize();

    return 0;
}


void metropolis(double Nspins,double MCSteps, double temp, mat & spinMatrix, vec & expecVal, ofstream & outFile, int waitSteps){

    random_device rd;
    int seed = rd();
    mt19937_64 gen(seed);
    //cout << rd() << endl;

    uniform_real_distribution<double> distribution(0.0,1.0);
    // cout << "I have an awesome random number: " << distribution(gen) << endl;


    //double dE[5] = {4,2,0,-2,-4};
    double thatFactor[17];
    double k_b = 1;//1.38064852e-23;
    double beta = -1.0/(k_b*temp);

    int writingFreq = 1;

    double energy = 0;
    double magneticMoment = 0;

    double acceptedFlips = 0;

    double spinNorm = 1.0/(Nspins*Nspins);


    //Init
    for (int i = 0; i < Nspins; i++){
        for (int j = 0; j<Nspins; j++){
            spinMatrix(i,j) = 1.0;
            //spinMatrix(i,j) = 2*round(distribution(gen)) - 1.0;


        }
    }

    for (int i = 0; i < Nspins; i++){
        for (int j = 0; j<Nspins; j++){
            energy -= spinMatrix(i,j)*(spinMatrix(periodic(i-1,Nspins),j)+ spinMatrix(i,periodic(j-1,Nspins)));

            magneticMoment += double(spinMatrix(i,j));
        }
    }

    for (double j = -8; j <= 8; j+= 4){
        thatFactor[int(j)+8] = exp(beta*j);
    }




    //MC loop
    for (int i = 0; i < MCSteps; i++){




        for (int j = 0; j<Nspins*Nspins; j++){

            int x_index = int(distribution(gen)*double(Nspins));
            int y_index = int(distribution(gen)*double(Nspins));


            int dEnergy = 2*spinMatrix(x_index,y_index)*( spinMatrix(periodic(x_index +1,Nspins),y_index)+ spinMatrix(periodic(x_index -1,Nspins),y_index)+
                                                                                                          spinMatrix(x_index,periodic(y_index +1,Nspins))+spinMatrix(x_index,periodic(y_index-1,Nspins)));

            //cout << thatFactor[int(dEnergy)+8] << endl;
            if (distribution(gen) <= thatFactor[int(dEnergy)+8]){

                spinMatrix(x_index,y_index) *= -1;
                energy += double(dEnergy);
                magneticMoment += double(2*spinMatrix(x_index,y_index));

                if (i > waitSteps){
                    acceptedFlips++;
                }

            }


        }

        if (i == 500){
            for (int k=0; k<Nspins; k++) {
                for (int l=0; l<Nspins; l++) {
                    if (spinMatrix(k,l) == 1.0) {
                        cout << "o ";
                    } else {
                        cout << "  ";
                    }
                }
                cout << endl;
            }
            cout << endl;
            //spinMatrix.print();
        }


        if (i > waitSteps){
            expecVal[0] += energy;
            expecVal[1] += energy*energy;
            expecVal[2] += (magneticMoment);
            expecVal[3] += magneticMoment*magneticMoment;
            expecVal[4] += fabs(magneticMoment);
        }




        if (i%writingFreq == 0 && i > waitSteps){
            outFile << (i+1- waitSteps) << " " << energy*spinNorm << " " <<expecVal[0]*spinNorm/(i+1 - waitSteps) << " " << fabs(magneticMoment)*spinNorm << " " << expecVal[4]*spinNorm/(i+1- waitSteps) << " " << acceptedFlips*spinNorm << " " << (expecVal[1]/(i+1- waitSteps) -  (expecVal[0]/(i+1- waitSteps))*(expecVal[0]/(i+1- waitSteps))) << " " << seed << endl;
        }



    }




}

void output(ofstream & file, vec & expecVal, int Nspins, double temp, double MCsteps){
    double factor = 1.0/MCsteps;
    double spinNorm = 1.0/(Nspins*Nspins);

    double expectE = expecVal[0]*factor;
    double expectESquared = expecVal[1]*factor;
    double expectM = expecVal[2]*factor;
    double expectMSquared = expecVal[3]*factor;
    double expectMFabs = expecVal[4]*factor;

    double varE = (expectESquared - expectE*expectE)*spinNorm;
    double varM = (expectMSquared - expectM*expectM)*spinNorm;


    double Cv = (expectESquared - expectE*expectE)/temp/temp*spinNorm;
    double sus = (expectMSquared - expectM*expectM)/temp*spinNorm;
    //cout << temp << " " << expectE*spinNorm << " " <<  expectESquared*spinNorm << " " <<  Cv << " " <<  expectM*spinNorm <<  " " << expectMSquared*spinNorm << " " << expectMFabs << " " <<  sus << "\n";

    file << temp << " " << expectE*spinNorm << " " <<  expectESquared*spinNorm << " " << varE << " " <<  Cv << " " <<  expectM*spinNorm <<  " " << expectMSquared*spinNorm << " " << expectMFabs << " " << varM << " " << sus << "\n";
}
