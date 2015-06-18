#ifndef ENERDIST1PHENER_HH
#define ENERDIST1PHENER_HH

#include "EnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the discrete out-going photon energy distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class EnerDist1PhEner: public EnergyDist
{
    public:
        EnerDist1PhEner(double AWR, double incEnerLow, double incEnerHigh);
        virtual ~EnerDist1PhEner();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        double GetAverageOutEnergy()
        {
            if(photonType!=2)
            {
                return photonEn;
            }
            else
            {
                // assume average incoming neutron energy is 1eV
                cout << "### Use of large energy distribution approximation" << endl;
                return (photonEn+1*awr/(awr+1));
            }
        }
        string IdentifyYourSelf()
        {
            return "EnerDist1PhEner";
        }

        double photonType; // states whether the photn is a primary or non primary photon,
        //which also tells us whether the energy is the photon energy (0,1) or the binding energy (2)
        double photonEn, awr, incEnLow, incEnHigh;
    protected:
    private:
};

#endif // ENERDIST1PHENER_HH
