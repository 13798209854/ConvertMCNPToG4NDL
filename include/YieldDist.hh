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
        virtual void SubtractPrompt(YieldDist* promptYieldDist)=0;
        virtual void SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield)=0;
        virtual void ConvertToLinDist(int *regEndPos, int &numIncEner, double *incEner, double *yield)=0;

        double Interpolate(int aScheme,
            double x, double x1, double x2, double y1, double y2) const;
        double Histogram(double , double , double , double y1, double ) const;
        double LinearLinear(double x, double x1, double x2, double y1, double y2) const;
        double LinearLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLinear(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double Random(double , double , double , double y1, double y2) const;
        virtual string IdentifyYourSelf()=0;

    protected:
    private:
};

#endif // YIELDDIST_HH
