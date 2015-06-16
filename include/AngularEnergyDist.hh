#ifndef ANGULARENERGYDIST_HH
#define ANGULARENERGYDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "include/ElementNames.hh"
#include <iomanip>

class AngularEnergyDist
{
    public:
        AngularEnergyDist();
        virtual ~AngularEnergyDist();
        virtual void ExtractMCNPData(stringstream data)=0;
        virtual void WriteG4NDLData(stringstream data)=0;
        virtual void SetPoint(stringstream data, int &count, double incNEner)=0;
        SetTemperature(double temp) {temperature=temp;}

        double Interpolate(int aScheme, double x, double x1, double x2, double y1, double y2) const;
        double Histogram(double , double , double , double y1, double ) const;
        double LinearLinear(double x, double x1, double x2, double y1, double y2) const;
        double LinearLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLinear(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double Random(double , double , double , double y1, double y2) const;

        vector<double> incNEnerVec;
        double temperature=0;
    protected:
    private:
};

#endif // ANGULARENERGYDIST_HH
