#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>
#include <mpi.h>

using namespace std;
using namespace arma;

void output(ofstream & file, vec & expecVal, int Nspin, double temp, double MCsteps, int waitNSteps);
void metropolis(double Nspins, double MCSteps, double temp, mat & spinMatrix, vec & expecVal, int waitNSteps);

inline int periodic(int index, int limit){
    return (index + limit)%(limit);
}

int main(int argc, char *argv[])
{
    int numProcs,rank;
    ofstream outFile;
    ofstream outFileParameters;
    int Nspins,monteCarloSteps,nTemp;
    double minTemp,maxTemp;
    mat spinMatrix;
    int waitNSteps;

    //  MPI initializations
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //cout << rank << endl;
    //cout << "hei " << stod(argv[1]) << " " << "rank: " << rank << endl;
    cout << "I'm workng"<<endl;


    if (rank == 0 && argc < 6){
        cout << "You did not give enough arguments. Number of spins, minimal temperature and max temperature, number of point and number of monte carlo iterations." << endl;
        //cout << "To few arguments given!" << endl;

        exit(1);
    }

    if (rank == 0){
        outFileParameters.open("parameters.txt");
        outFile.open("data.txt");
        Nspins = stoi(argv[1]);
        monteCarloSteps = stoi(argv[5]);
        nTemp = stoi(argv[4]);
        minTemp = stod(argv[2]);
        maxTemp = stod(argv[3]);
    }

    if (rank == 0){

        outFileParameters << "minimumTemp " << minTemp << "\n";
        outFileParameters << "maximumTemp " << maxTemp << "\n";
        outFileParameters << "tempSteps " << nTemp << "\n";
        outFileParameters << "L " << Nspins << "\n";
        outFileParameters << "monteCarloSteps " << monteCarloSteps << "\n";
        outFileParameters << "coresUsed " << numProcs << "\n";
        outFileParameters.close();
    }



    MPI_Bcast(&Nspins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&monteCarloSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&minTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    spinMatrix = ones(Nspins,Nspins);
    waitNSteps = 5000;


    double dt = (maxTemp-minTemp)/double(nTemp);

    double timeSinceLast = MPI_Wtime();

    //return 0;
    double timeStart = MPI_Wtime();
//    if (rank == 0){
//        clock_t begin = clock();
//    }

    cout << "I am " << rank << " and will do some cool statistical mechanics RIGHT NOW!" << endl;
    for (double temp = minTemp; temp <= maxTemp; temp+=dt){

        // cout << "I am " << rank << " and will do T=" << temp << endl;
        vec localExpectationValues = zeros(5);

        // cout << "Doing some work: " << Nspins << " and " << monteCarloSteps << endl;
        metropolis(Nspins,monteCarloSteps, temp,spinMatrix,localExpectationValues, waitNSteps);

        vec totalExpectationValues = zeros(5);


        for (int i = 0; i < totalExpectationValues.size(); i++){
            MPI_Reduce(&localExpectationValues[i], &totalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        if (rank == 0){
            output(outFile, totalExpectationValues, Nspins, temp, monteCarloSteps*numProcs,waitNSteps);
        }


        if (rank == 0){
            cout << "Done for temp " << temp << " " << "It took " << double(MPI_Wtime() - timeSinceLast)  << " sec." << endl;;
            timeSinceLast = MPI_Wtime();
        }

    }

    if (rank == 0){
        outFile.close();

        cout << "Program done after: "<< double(MPI_Wtime() - timeStart) << " sec" << endl;

    }

    MPI_Finalize();

    return 0;
}


void metropolis(double Nspins,double MCSteps, double temp, mat & spinMatrix, vec & expecVal, int waitNSteps){

    random_device rd;
    mt19937_64 gen(rd());

    uniform_real_distribution<double> distribution(0.0,1.0);
    // cout << "I have an awesome random number: " << distribution(gen) << endl;


    //double dE[5] = {4,2,0,-2,-4};
    double thatFactor[17];
    double k_b = 1;//1.38064852e-23;
    double beta = -1.0/(k_b*temp);

    double energy = 0;
    double magneticMoment = 0;

    for (int i = 0; i < Nspins; i++){
        for (int j = 0; j<Nspins; j++){
            spinMatrix(i,j) = 1.0;

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

            }


        }

        if (i > waitNSteps){
            expecVal[0] += energy;
            expecVal[1] += fabs(energy*energy);
            expecVal[2] += (magneticMoment);
            expecVal[3] += magneticMoment*magneticMoment;
            expecVal[4] += fabs(magneticMoment);
        }

//        expecVal[0] += energy;
//        expecVal[1] += energy*energy;
//        expecVal[2] += (magneticMoment);
//        expecVal[3] += magneticMoment*magneticMoment;
//        expecVal[4] += fabs(magneticMoment);




    }




}

void output(ofstream & file, vec & expecVal, int Nspins, double temp, double MCsteps, int waitNSteps){
    double factor = 1.0/(MCsteps- waitNSteps);
    double spinNorm = 1.0/(Nspins*Nspins);
    double expectE = expecVal[0];
    double expectESquared = expecVal[1];
    double expectM = expecVal[2]*factor;
    double expectMSquared = expecVal[3]*factor;
    double expectMFabs = expecVal[4]*factor;

    double varE = (expectESquared - expectE*expectE)
    double varM = (expectMSquared - expectM*expectM)*spinNorm;


    double Cv = (expectESquared - expectE*expectE)/temp/temp*spinNorm;
    double sus = (expectMSquared - expectM*expectM)/temp*spinNorm;
    //cout << temp << " " << expectE*spinNorm << " " <<  expectESquared*spinNorm << " " <<  Cv << " " <<  expectM*spinNorm <<  " " << expectMSquared*spinNorm << " " << expectMFabs << " " <<  sus << "\n";

    file << temp << " " << expectE*spinNorm << " " <<  expectESquared*spinNorm << " " << varE << " " <<  Cv << " " <<  expectM*spinNorm <<  " " << expectMSquared*spinNorm << " " << expectMFabs*spinNorm << " " << varM << " " << sus << "\n";
}
