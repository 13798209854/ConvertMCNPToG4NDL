#ifndef ENERDISTLEVSCAT_HH
#define ENERDISTLEVSCAT_HH

#include "EnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the level scattering out-going neutron energy distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class EnerDistLevScat: public EnergyDist
{
    public:
        EnerDistLevScat(double enerRangeStart, double enerRangeEnd);
        virtual ~EnerDistLevScat();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        string IdentifyYourSelf()
        {
            return "EnerDistLevScat";
        }

        double firstHalfEq, secondHalfEq, enerStart, enerEnd;
    protected:
    private:
};

#endif // ENERDISTLEVSCAT_HH
