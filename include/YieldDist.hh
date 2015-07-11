#ifndef YIELDDIST_HH
#define YIELDDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "../include/CSDist.hh"

using namespace std;
/*
Created By: Wesley Ford June 17, 2015

This class is the mother class to all the classes responsible for the extraction of the out-going neutron yield distribution data from MCNP
*/

class YieldDist
{
    public:
        YieldDist();
        virtual ~YieldDist();
        virtual void ExtractMCNPData(stringstream &stream, int &count)=0;
        virtual void WriteG4NDLData(stringstream &stream)=0;
        virtual void SubtractPrompt(YieldDist* &promptYieldDist)=0;
        virtual void SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield)=0;
        virtual void SubtractPrompt(double *&totalCoeff, int &totalNumCoeff)=0;
        virtual void ConvertToLinDist(int *regEndPos, int &numIncEner, double *&incEner, double *&yield)=0;
        virtual void AddData(YieldDist* nYield, CSDist** nCSDist, int index, int numProc, int numProc2)=0;
        virtual void AddData(int &numRegsSum, int &numIncEnerSum, int* &regEndPosSum, int* &intSchemeSum, double* &incEnerSum, double* &yieldSum, CSDist** nCSDist, int index, int numProc, int numProc2)=0;
        virtual void SetYieldData(int &numRegsSum, int &numIncEnerSum, int* &regEndPosSum, int* &intSchemeSum, double* &incEnerSum, double* &yieldSum, CSDist** nCSDist, int index, int numProc, int numProc2)=0;
        virtual string IdentifyYourSelf()=0;

        double Interpolate(int aScheme, double x, double x1, double x2, double y1, double y2) const;
        double Histogram(double , double , double , double y1, double ) const;
        double LinearLinear(double x, double x1, double x2, double y1, double y2) const;
        double LinearLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLinear(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2) const;
    protected:
    private:
};

#endif // YIELDDIST_HH
