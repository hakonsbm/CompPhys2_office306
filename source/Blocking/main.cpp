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
    double tmp1;
    double tmp2;
    double tmp3;
    double tmp4;
    double average = 0;
    double std = 0;
    double mean = 0;
    double meanSquared = 0;
    int count = 0;
    int numberOfBlocks = 0;
    int maxBlockSize = 0;
    int minBlockSize = 100;
    outFile.open("outputSTD.dat");
    inFile.open("input.dat");
    while (!inFile.eof()) {
        inFile >> tmp1 >> tmp2 >> tmp3 >> tmp4;
        energy.push_back(tmp1);
        count++;
    }
    count--;
    energy.pop_back();
    energySquared.pop_back();
    inFile.close();
    maxBlockSize = count / 4;
    for(int i = minBlockSize; i < maxBlockSize + 1; i++) {
        if(count % i != 0)continue;
        mean = 0;
        meanSquared = 0;
        numberOfBlocks = count / i;
        for(int j = 0; j < count; j += i) {
            average = 0;
            for(int k = j; k < j+i; k++){
                average += energy[k];
            }
            average /= (double)i;
            mean += average;
            meanSquared += average*average;
        }
        mean /= (double)numberOfBlocks;
        meanSquared /= (double)numberOfBlocks;
        std = sqrt((meanSquared-mean*mean)/(double)numberOfBlocks);
        outFile << setw(15) << setprecision(8) << i << "\t" << std << "\t" << mean << endl;
    }
    outFile.close();
    return 0;
}
