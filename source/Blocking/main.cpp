#include "iostream"
#include "fstream"
#include "vector"
#include "cmath"
#include "iomanip"

using namespace std;
ofstream outfile;

int main() {

    ifstream readFileEnergy("out4-4.d", ios::in);
    vector <double> energy;
    vector <double> energySquared;
    double tmp1;
    double tmp2;
    double tmp3;
    double tmp4;
    double average = 0;
    double averageSquared = 0;
    double std = 0;
    double averageSTD = 0;
    double mean = 0;
    int count = 0;
    int numberOfBlocks = 0;
    outfile.open("outputSTD.d");
    readFileEnergy.open("../outfiles/HeliumSimpleAnalytical_samples");
    while (!readFileEnergy.eof()) {
        readFileEnergy >> tmp1 >> tmp2 >> tmp3 >> tmp4;
        energy.push_back(tmp1);
        energySquared.push_back(tmp2);
        count++;
    }
    count--;
    energy.pop_back();
    energySquared.pop_back();
    readFileEnergy.close();
    for(int i = 100; i < count/5+1; i += 100) {
        mean = 0;
        averageSTD = 0;
        for(int j = 0; j < count; j += i) {
            numberOfBlocks = count / i;
            average = 0;
            averageSquared = 0;
            std = 0;
            for(int k = j; k < j+i; k++){
                average += energy[k];
                averageSquared += energySquared[k];
            }
            average /= (double)i;
            averageSquared /= (double)i;
            std = sqrt((averageSquared-average*average)/(double)i);
            mean += average;
            averageSTD += std;
        }
        mean /= (double)numberOfBlocks;
        averageSTD /= (double)numberOfBlocks;
        outfile << setw(15) << setprecision(8) << i << "\t" << averageSTD  << "\t" << mean << endl;
    }
outfile.close();
    return 0;
}
