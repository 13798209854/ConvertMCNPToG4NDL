#ifndef ANGENDISTNBODY_HH
#define ANGENDISTNBODY_HH

#include "AngularEnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy dependant N-body phase space out-going angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/


class AngEnDistNBody : public AngularEnergyDist
{
    public:
        AngEnDistNBody();
        virtual ~AngEnDistNBody();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        void ConvertToEnerAndAngDist(EnergyDist **enDist, AngularDist **angDist, int &numAngEner)
        {
            cout << "Error tried to convert AngEnDistNBody to an energy distributin and an Angular distribution, this is not possible" << endl;
        }

        int numBodies;
        double particleMassRatio;
    protected:
    private:
};

#endif // ANGENDISTNBODY_HH
