#ifndef ENERDISTCONTAB_Hh
#define ENERDISTCONTAB_Hh

#include "EnergyDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the continous tabular out-going energy distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class EnerDistConTab: public EnergyDist
{
    public:
        EnerDistConTab(/*int EnerDistStart*/);
        virtual ~EnerDistConTab();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        double GetAverageOutEnergy();
        string IdentifyYourSelf()
        {
            return "EnerDistConTab";
        }

        int numRegs, numIncEner, startEnerDist;
        //numRegs number of interpolation regions
        //numIncEner number of incoming neutron energy
        int *regEndPos, *intScheme1, *intScheme2, *distPos, *numPEnerPoints, *numDiscreteEnerPoints;
        //regEndPos containes the end positon for each region
        //intScheme1 is the scheme for interpolating between incEner points
        //intScheme2 is the scheme for interpolating between outEner points
        //distPos contains the positions of the outgoing energy-probability data for each incoming energy
        //numPEnerPoints is the total number of probability energy points
        //numDiscreteEnerPoints is the number of probability energy points that are discrete, numDiscreteEnerPoints<=numPEnerPoints
        double *incEner; // contains incoming neutron energy
        double **outEner, **outProb, **outSumProb;
        //outEner outgoing energy
        //outProb probability of this happening
        //outSumProb cumulative sum of outProb up to this incoming neutron energy

    protected:
    private:
};

#endif // ENERDISTCONTAB_Hh
