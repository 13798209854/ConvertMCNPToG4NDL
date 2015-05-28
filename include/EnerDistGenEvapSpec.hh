#ifndef ENERDISTGENEVAPSPEC_HH
#define ENERDISTGENEVAPSPEC_HH

#include "EnergyDist.hh"

class EnerDistGenEvapSpec: public EnergyDist
{
    public:
        EnerDistGenEvapSpec();
        virtual ~EnerDistGenEvapSpec();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        string IdentifyYourSelf()
        {
            return "EnerDistGenEvapSpec";
        }

        int numRegs, numIncEner, numNormOutEnerDistPoints;
        int *regEndPos, *intScheme;
        double *incEner, *outEnerMulti, *normOutEnerDist;

    protected:
    private:
};

#endif // ENERDISTGENEVAPSPEC_HH
