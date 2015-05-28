#ifndef ENERDISTEQPENERBINS_HH
#define ENERDISTEQPENERBINS_HH

#include "EnergyDist.hh"

class EnerDistEqPEnerBins: public EnergyDist
{
    public:
    EnerDistEqPEnerBins();
    virtual ~EnerDistEqPEnerBins();
    void ExtractMCNPData(stringstream stream, int &count);
    void WriteG4NDLData(stringstream data);
    string IdentifyYourSelf()
    {
        return "EnerDistEqPEnerBins";
    }

    int numReg, numIncEner, numOutEnerPerInc;
    int *regEndPos, *intScheme;
    double *incEner;
    double **outEner;

    protected:
    private:
};

#endif // ENERDISTEQPENERBINS_HH
