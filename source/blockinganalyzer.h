#ifndef BLOCKINGANALYZER_H
#define BLOCKINGANALYZER_H

class BlockingAnalyzer
{
public:
    BlockingAnalyzer();
    void runStatisticalAnalysis();
    double getMean() {return m_average;}
    double getStd() {return m_std;}

private:
    double m_std;
    double m_average;
    double m_averageSquared;
};

#endif // BLOCKINGANALYZER_H
