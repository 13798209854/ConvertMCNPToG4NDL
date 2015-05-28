#ifndef ENERDISTEVAPSPEC_HH
#define ENERDISTEVAPSPEC_HH

#include "EnergyDist.hh"


class EnerDistEvapSpec : public EnergyDist
{
    public:
        EnerDistEvapSpec();
        virtual ~EnerDistEvapSpec();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        string IdentifyYourSelf()
        {
            return "EnerDistEvapSpec";
        }

        int numRegs, numIncEner;
        int *regEndPos, *intScheme;
        double *incEner, *tValues, restrictEner;
    protected:
    private:
};

#endif // ENERDISTEVAPSPEC_HH
