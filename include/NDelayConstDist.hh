#ifndef NDELAYCONSTDIST_HH
#define NDELAYCONSTDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;
/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the delayed neutron average lifetime distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class NDelayConstDist
{
    public:
        NDelayConstDist(int numDNFam);
        virtual ~NDelayConstDist();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);

        int numDNPrecursorFam;
        vector <int> numRegs, numIncEner;
        vector <double> delConst;
        vector <int*> regEndPos, intScheme;
        vector <double*> incEner, prob;
    protected:
    private:
};

#endif // NDELAYCONSTDIST_HH
