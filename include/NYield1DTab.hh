#ifndef NYield1DTab_HH
#define NYield1DTab_HH

#include "YieldDist.hh"

class NYield1DTab: public YieldDist
{
    public:
        NYield1DTab();
        virtual ~NYield1DTab();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

        int numRegs, numIncEner;
        int *regEndPos, *intScheme;
        double *incEner; // contains incoming neutron energy
        double *yield;
    protected:
    private:
};

#endif // NYield1DTab_HH
