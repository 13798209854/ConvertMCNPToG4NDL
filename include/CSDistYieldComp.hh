#ifndef CSDISTYIELDCOMP_HH
#define CSDISTYIELDCOMP_HH

#include "CSDist.hh"


class CSDistYieldComp : public CSDist
{
    public:
        CSDistYieldComp();
        virtual ~CSDistYieldComp();
        void ExtractMCNPData(stringstream data, int count&);
        void WriteG4NDLYieldData(stringstream data);
        void WriteG4NDLCSData(stringstream data);
        string IdentifyYourSelf();
        double Interpolate(double x);
        void SetCSData(CSDist* nCSData)
        {
            nCSData->SetNCSData(enerCSVec, csEnerStart, csVec, csSize);
        }

        int mtMulti, numRegs, mtNum, numIncEner, csSize, csEnerStart;
        int *regEndPos, *intScheme
        double *yieldVec, *enerVec, *enerCSVec, *csVec;
    protected:
    private:
};

#endif // CSDISTYIELDCOMP_HH
