#ifndef ENERDISTWATTSPEC_HH
#define ENERDISTWATTSPEC_HH

#include "EnergyDist.hh"


class EnerDistWattSpec : public EnergyDist
{
    public:
        EnerDistWattSpec();
        virtual ~EnerDistWattSpec();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        string IdentifyYourSelf()
        {
            return "EnerDistWattSpec";
        }

        int numRegsA, numRegsb, numIncEnerA, numIncEnerB;
        //numRegs number of interpolation regions
        //numIncEner number of incoming neutron energy
        int *regEndPosA, *regEndPosB, *intSchemeA, *intSchemeB;
        //regEndPos containes the end positon for each region
        //intScheme is the scheme for interpolating between incEner points
        double *incEnerA, *incEnerB, *aValues, *bValues;
        double rejectEner;
    protected:
    private:
};

#endif // ENERDISTWATTSPEC_HH
