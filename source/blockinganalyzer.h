#ifndef BLOCKINGANALYZER_H
#define BLOCKINGANALYZER_H

//class VMCSolver {};

class BlockingAnalyzer
{
public:
    BlockingAnalyzer();
    void runStatisticalAnalysis();
    double getMean() {return m_mean;}
    double getStd() {return m_std;}

private:
    char const *m_fileName = "out4-4.d";
    double m_std;
    double m_mean;
    int m_totalEntries;
    int m_blockSize;
};

#endif // BLOCKINGANALYZER_H
