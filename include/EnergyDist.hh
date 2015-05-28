#ifndef ENERGYDIST_HH
#define ENERGYDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "zlib.h"
#include <dirent.h>
#include "include/ElementNames.hh"
#include <iomanip>


class EnergyDist
{
    public:
        EnergyDist();
        virtual ~EnergyDist();
        virtual void ExtractMCNPData(stringstream data, int count&)=0;
        virtual void WriteG4NDLData(stringstream data)=0;
        virtual string IdentifyYourSelf()=0;
        virtual double GetAverageOutEnergy()=0;
    protected:
    private:
};

#endif // ENERGYDIST_HH
