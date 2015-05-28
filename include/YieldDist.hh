#ifndef YIELDDIST_HH
#define YIELDDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "include/ElementNames.hh"
#include <iomanip>

class YieldDist
{
    public:
        YieldDist();
        virtual ~YieldDist();
        virtual void ExtractMCNPData(stringstream stream, int &count)=0;
        virtual void WriteG4NDLData(stringstream data)=0;
    protected:
    private:
};

#endif // YIELDDIST_HH
