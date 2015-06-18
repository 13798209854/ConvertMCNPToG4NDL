#ifndef ENERDISTWATTSPEC_HH
#define ENERDISTWATTSPEC_HH

#include "EnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the Watt spectrum out-going neutron energy distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class EnerDistWattSpec : public EnergyDist
{
    public:
        EnerDistWattSpec();
        virtual ~EnerDistWattSpec();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        string IdentifyYourSelf()
        {
            return "EnerDistWattSpec";
        }

        int numRegsA, numRegsB, numIncEnerA, numIncEnerB;
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
