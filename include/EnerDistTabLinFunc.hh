#ifndef ENERDISTTABLINFUNC_HH
#define ENERDISTTABLINFUNC_HH

#include "EnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the tabular linear function out-going neutron energy distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class EnerDistTabLinFunc : public EnergyDist
{
    public:
        EnerDistTabLinFunc(/*int EnerDistStart*/);
        virtual ~EnerDistTabLinFunc();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        string IdentifyYourSelf()
        {
            return "EnerDistTabLinFunc";
        }

        int numRegs, numIncEner, startEnerDist;
        //numRegs number of interpolation regions
        //numIncEner number of incoming neutron energy
        int *regEndPos, *intScheme, *distPos, *numPEnerPoints;
        //regEndPos containes the end positon for each region
        //intScheme is the scheme for interpolating between incEner points
        //distPos contains the positions of the outgoing energy-probability data for each incoming energy
        //numPEnerPoints is the total number of probability energy points
        double *incEner; // contains incoming neutron energy
        double **outProb, **outEnerOffset, **outEnerMulti;
        //outSumProb cumulative sum of outProb up to this incoming neutron energy
    protected:
    private:
};

#endif // ENERDISTTABLINFUNC_HH
