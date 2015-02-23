#include "iostream"
#include "fstream"
#include "vector"
#include "cmath"
#include "iomanip"

using namespace std;
ofstream outfile;

int main() {

    ifstream readFileEnergy ("out4-4.d", ios::in);
    vector <double> energy;
    vector <double> energySquared;
    double tmp1;
    double tmp2;
    double average;
    double averageSquared;
    double std;
    int count;
    readFileEnergy >> count;
    outfile.open("outputSTD.d");
    for (int i = 0; i < count; i++) {
        readFileEnergy >> tmp1 >> tmp2;
        energy.push_back(tmp1);
        energySquared.push_back(tmp2);
    }
    for(int i = 1000; i < count/5+1; i += 1000) {
        for(int j= 0; j < count +1; j++) {
            average += energy[j];
            averageSquared += energySquared[j];
        }
        average /= (double)i;
        averageSquared /= (double)i;
        std = sqrt(averageSquared-average*average);
        outfile << setw(15) << setprecision(8) << i << "\t" << std  << "\t" << average << endl;
    }
outfile.close();
    return 0;
}
