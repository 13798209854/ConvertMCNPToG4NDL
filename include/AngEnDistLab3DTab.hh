#ifndef AngEnDistLab3DTab_HH
#define AngEnDistLab3DTab_HH

#include "AngularEnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy independant tabular labratory out-going angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/


class AngEnDistLab3DTab : public AngularEnergyDist
{
    public:
        AngEnDistLab3DTab();
        virtual ~AngEnDistLab3DTab();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);

        int numRegs, numIncEner /*startEnerDist*/;
        int *regEndPos, *intScheme1, *intScheme2, *outAngDistPos, *numPAngPoints;
        double *incEner; // contains incoming neutron energy
        double **outAng;
        int **outEnerDistPos, **intScheme3, **numPEnerPoints;
        double ***outEner, ***outEnProb, ***outEnSumProb;
    protected:
    private:
};

#endif // AngEnDistLab3DTab_HH
