#include <iostream>
#include <armadillo>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>


using namespace std;
using namespace arma;

void output(ofstream & file, vec & expecVal, int Nspin, double temp, double MCsteps);
void metropolis(double Nspins,double MCSteps, double temp, mat & spinMatrix, vec & expecVal);

inline int periodic(int index, int limit){
    return (index + limit)%(limit);
}

int main(int argc, char *argv[])
{

    if (argc < 6){
        cout << "You did not give enough arguments. Number of spins, minimal temperature and max temperature, number of point and number of monte carlo iterations." << endl;
        //cout << "To few arguments given!" << endl;

        exit(1);
    }


    ofstream outFile("data.txt");
    int Nspins = stoi(argv[1]);
    int monteCarloSteps = stoi(argv[5]);
    mat spinMatrix = ones(Nspins,Nspins);
    int nTemp = stoi(argv[4]);
    double minTemp = stod(argv[2]);
    double maxTemp = stod(argv[3]);



    double dt = (maxTemp-minTemp)/double(nTemp);
    cout << minTemp << " " << maxTemp << " " << dt << endl;

    double timeSinceLast = 0;

    //return 0;
    clock_t begin = clock();

    for (double temp = minTemp; temp <= maxTemp; temp+=dt){


        vec expectationValues = zeros(4);


        metropolis(Nspins,monteCarloSteps, temp,spinMatrix,expectationValues);
        output(outFile, expectationValues, Nspins, temp, monteCarloSteps);



        cout << "Done for temp " << temp << endl;
        cout << "It took " << double(clock() - timeSinceLast) / double(CLOCKS_PER_SEC) << " sec." << endl;
        timeSinceLast = clock();
    }

    outFile.close();

    cout << "Program done after: "<< double(clock() - begin)/ double(CLOCKS_PER_SEC) << " sec" << endl;
    return 0;
}


void metropolis(double Nspins,double MCSteps, double temp, mat & spinMatrix, vec & expecVal){


    random_device rd;
    mt19937_64 gen(rd());

    uniform_real_distribution<double> distribution(0.0,1.0);


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


        expecVal[0] += energy;
        expecVal[1] += energy*energy;
        expecVal[2] += fabs(magneticMoment);
        expecVal[3] += magneticMoment*magneticMoment;

    }




}

void output(ofstream & file, vec & expecVal, int Nspins, double temp, double MCsteps){
    double factor = 1.0/MCsteps;
    double spinNorm = 1.0/(Nspins*Nspins);

    double expectE = expecVal[0]*factor;
    double expectESquared = expecVal[1]*factor;
    double expectM = expecVal[2]*factor;
    double expectMSquared = expecVal[3]*factor;


    double Cv = (expectESquared - expectE*expectE)/temp/temp*spinNorm;
    double sub = (expectMSquared - expectM*expectM)/temp*spinNorm;
    cout << temp << " " << expecVal[0]*spinNorm << " " <<  expecVal[1]*spinNorm << " " <<  Cv << " " <<  expecVal[2]*spinNorm <<  " " << expecVal[3]*spinNorm << " " <<  sub << "\n";

    file << temp << " " << expecVal[0]*spinNorm << " " <<  expecVal[1]*spinNorm << " " <<  Cv << " " <<  expecVal[2]*spinNorm <<  " " << expecVal[3]*spinNorm << " " <<  sub << "\n";
}
