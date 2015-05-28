#ifndef ENERDISTMAXWELLFISSPEC_HH
#define ENERDISTMAXWELLFISSPEC_HH

#include "EnergyDist.hh"

class EnerDistMaxwellFisSpec: public EnergyDist
{
    public:
        EnerDistMaxwellFisSpec();
        virtual ~EnerDistMaxwellFisSpec();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        string IdentifyYourSelf()
        {
            return "EnerDistMaxwellFisSpec";
        }

        int numRegs, numIncEner;
        int *regEndPos, *intScheme;
        double *incEner, *tValues, restrictEner;
    protected:
    private:
};

#endif // ENERDISTMAXWELLFISSPEC_HH
