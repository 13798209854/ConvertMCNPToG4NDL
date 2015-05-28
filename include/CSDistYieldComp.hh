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
        void SetCSData(double *csEnerVec, double *cSVec, int cSSize)
        {
            enerCSVec = csEnerVec;
            csVec = cSVec;
            csSize= cSSize;
        }

        int mtMulti, numRegs, numIncEner, csSize;
        int *regEndPos, *intScheme
        double *yieldVec, *enerVec, *enerCSVec, *csVec;
    protected:
    private:
};

#endif // CSDISTYIELDCOMP_HH
