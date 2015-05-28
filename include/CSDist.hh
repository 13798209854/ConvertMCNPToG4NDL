#ifndef CSDIST_HH
#define CSDIST_HH

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

class CSDist
{
    public:
        CSDist();
        virtual ~CSDist();
        virtual void ExtractMCNPData(stringstream data)=0;
        virtual void WriteG4NDLYieldData(stringstream data)=0;
        virtual void WriteG4NDLCSData(stringstream data)=0;
        virtual void SetCSData(double *enerVec, double *csVec)=0;
        virtual double* GetEnergyVec()=0;
        virtual double* GetCSVec()=0;
        virtual int GetCSSize()=0;
        virtual string IdentifyYourSelf()=0;
    protected:
    private:
};

#endif // CSDIST_HH
