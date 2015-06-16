#ifndef CSDIST1DTAB_HH
#define CSDIST1DTAB_HH

#include "CSDist.hh"

class CSDist1DTab: public CSDist
{
    public:
        CSDist1DTab(double *energyVec);
        CSDist1DTab(double *energyVec, int numCSEner);
        virtual ~CSDist1DTab();
        void ExtractMCNPData(stringstream data, int count&);
        void WriteG4NDLCSData(stringstream data);
        double Interpolate(double x);
        void SetNCSData(double* enerCSVec, int &csEnerStart, double* csVec, int &csSize)
        {
            enerCSVec=enerVec;
            csEnerStart=startEner;
            csVec=CSVec;
            csSize=CSVecSize;
        }
        void IdentifyYourSelf();
        double* GetEnergyVec();
        {
            return (enerVec+startEner-1);
        }
        double* GetCSVec();
        {
            return CSVec;
        }
        int GetCSSize();
        {
            return CSVecSize;
        }

        int startEner, CSVecSize;
        double *CSVec, *enerVec;
        bool elastic=false;
    protected:
    private:
};

#endif // CSDIST1DTAB_HH
