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
    double tmp5;
    double tmp6;
    double tmp7;
    double tmp8;
    double tmp9;
    double tmp10;
    double tmp11;
    double tmp12;
    double tmp13;
    double tmp14;
    double tmp15;
    double tmp16;
    double tmp17;
    double tmp18;
    double tmp19;
    double tmp20;
    double tmp21;
    double tmp22;
    double tmp23;
    double tmp24;
    double tmp25;
    double tmp26;
    double tmp27;
    double tmp28;
    double tmp29;
    double tmp30;
    double tmp31;
    double tmp32;
    double tmp33;
    double tmp34;
    double average = 0;
    double std = 0;
    double mean = 0;
    double meanSquared = 0;
    int count = 0;
    int numberOfBlocks = 0;
    int maxBlockSize = 0;
    int minBlockSize = 10;
    outFile.open("outputSTD.dat");
    inFile.open("../source/outfiles/Helium_blockingSamples");
    cout << "reading file" <<endl;
    while (inFile && count<=1000000) {
        inFile >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7 >> tmp8 >> tmp9 >> tmp10;// >> tmp11 >> tmp12 >> tmp13 >> tmp14 >> tmp15 >> tmp16 >> tmp17 >> tmp18 >> tmp19 >> tmp20 >> tmp21 >> tmp22 >> tmp23 >> tmp24 >> tmp25 >> tmp26 >> tmp27 >> tmp28 >> tmp29 >> tmp30 >> tmp31 >> tmp32 >> tmp33 >> tmp34;
        energy.push_back(tmp1);
        count++;
    }
    count--;
    energy.pop_back();
    inFile.close();
    maxBlockSize = count / 4;
    cout << "blocking calculation" << minBlockSize << " " <<  maxBlockSize <<endl;
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
