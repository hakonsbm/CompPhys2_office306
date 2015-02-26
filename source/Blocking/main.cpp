#include "iostream"
#include "fstream"
#include "vector"
#include "cmath"
#include "iomanip"

using namespace std;

int main() {

    ifstream inFile;
    ofstream outFile;
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
    int dividend = 1000;
    int numberOfBlocks = 0;
    int maxBlockSize = 0;
    int minBlockSize = 0;
    outFile.open("outputSTD.dat");
    inFile.open("input.dat");
    while (!inFile.eof()) {
        inFile >> tmp1 >> tmp2 >> tmp3 >> tmp4;
        energy.push_back(tmp1);
        energySquared.push_back(tmp2);
        count++;
    }
    count--;
    energy.pop_back();
    energySquared.pop_back();
    inFile.close();
    maxBlockSize = count / 4;
    minBlockSize = count / dividend;
    if(minBlockSize < 1000) {
        minBlockSize = 1000;
        dividend = count / minBlockSize;
    }
    dividend--;
    for(int i = minBlockSize; i < maxBlockSize + 1; i += count / dividend - count / (dividend + 1)) {
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
        outFile << setw(15) << setprecision(8) << i << "\t" << averageSTD  << "\t" << mean << endl;
        dividend--;
    }
    outFile.close();
    return 0;
}
