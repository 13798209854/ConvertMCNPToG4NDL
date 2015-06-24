#ifndef ANGULARDIST_HH
#define ANGULARDIST_HH

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

This class is the mother class for all the classes responsible for the extraction of the energy independant out-going angular distribution data from MCNP
*/

class AngularDist
{
    public:
        AngularDist();
        virtual ~AngularDist();
        virtual void ExtractMCNPData(stringstream &stream, int &count)=0;
        virtual void WriteG4NDLData(stringstream &stream)=0;
        virtual void SetPoint(stringstream &stream, int &count, double incNEner)=0;
        void SetTemperature(double temp) {temperature=temp;}
        vector<double> incNEnerVec;
        double temperature=0;
    protected:
    private:
};

#endif // ANGULARDIST_HH
