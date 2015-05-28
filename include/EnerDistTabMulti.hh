#ifndef ENERDISTTABMULTI_HH
#define ENERDISTTABMULTI_HH

#include "EnergyDist.hh"


class EnerDistTabMulti : public EnergyDist
{
    public:
        EnerDistTabMulti();
        virtual ~EnerDistTabMulti();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        string IdentifyYourSelf()
        {
            return "EnerDistTabMulti";
        }

        int numRegs, numIncEner, numOutEnerPerIn;
        int *regEndPos, *intScheme;
        double *incEner, **tValues
    protected:
    private:
};

#endif // ENERDISTTABMULTI_HH
