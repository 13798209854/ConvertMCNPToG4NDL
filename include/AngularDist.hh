#ifndef ANGULARDIST_HH
#define ANGULARDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "include/ElementNames.hh"
#include <iomanip>


class AngularDist
{
    public:
        AngularDist();
        virtual ~AngularDist();
        virtual void ExtractMCNPData(stringstream data)=0;
        virtual void WriteG4NDLData(stringstream data)=0;
        virtual void SetPoint(stringstream data, int &count, double incNEner)=0;
        SetTemperature(double temp) {temperature=temp;}
        vector<double> incNEnerVec;
        double temperature=0;
    protected:
    private:
};

#endif // ANGULARDIST_HH
