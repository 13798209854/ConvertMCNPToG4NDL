#ifndef ANGDISTISOP_HH
#define ANGDISTISOP_HH

#include "AngularDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy independant isotropic out-going photon angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class AngDistIsoP: public AngularDist
{
    public:
        AngDistIsoP();
        virtual ~AngDistIsoP();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        void SetPoint(stringstream &stream, int &count, double incNEner)
        {
            incNEnerVec.push_back(incNEner);
        }
    protected:
    private:
};

#endif // ANGDISTISOP_HH
