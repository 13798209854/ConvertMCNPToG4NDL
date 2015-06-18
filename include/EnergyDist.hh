#ifndef ENERGYDIST_HH
#define ENERGYDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include <iomanip>
using namespace std;

/*
Created By: Wesley Ford June 17, 2015

This class is the mother class of all the classes responsible for the extraction of the out-going energy distribution data from MCNP
*/

class EnergyDist
{
    public:
        EnergyDist();
        virtual ~EnergyDist();
        virtual void ExtractMCNPData(stringstream &stream, int &count)=0;
        virtual void WriteG4NDLData(stringstream &stream)=0;
        virtual string IdentifyYourSelf()=0;
        virtual double GetAverageOutEnergy()
        {
            return 0.;
        }
    protected:
    private:
};

#endif // ENERGYDIST_HH
