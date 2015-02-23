#include "blockinganalyzer.h"
#include "vmcsolver.h"
#include "iostream"
#include "fstream"
#include "vector"

using namespace std;

BlockingAnalyzer::BlockingAnalyzer() {}



void runStatisticalAnalysis() {
    
    ifstream readFileEnergy ("out4-4.d", ios::in);
    ifstream readFileEnergySquared ("out4-4.d", ios::in);
    vector <double> energy;
    vector <double> energySquared;
    double tmp1;
    double tmp2;
    double average;
    double averageSquared;
    double std;
    int count;
    readFileEnergy >> count;
    for (int i = 0; i < count; i++) {
        readFileEnergy >> tmp1;
        readFileEnergySquared >> tmp2;
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
    }
}
