#ifndef NDELAYCONSTDIST_HH
#define NDELAYCONSTDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "include/ElementNames.hh"
#include <iomanip>

class NDelayConstDist
{
    public:
        NDelayConstDist(int tableEndPos);
        virtual ~NDelayConstDist();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

        int tabEndPos;
        vector <int> numRegs, numIncEner;
        vector <double> delConst;
        vector <int*> regEndPos, intScheme;
        vector <double*> incEner, prob;
    protected:
    private:
};

#endif // NDELAYCONSTDIST_HH
