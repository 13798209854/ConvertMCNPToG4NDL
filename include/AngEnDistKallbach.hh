#ifndef AngEnDistKallbach_HH
#define AngEnDistKallbach_HH

#include "AngularEnergyDist.hh"

#define numDistSample 10

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy dependant Kallbach out-going angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class AngEnDistKallbach : public AngularEnergyDist
{
    public:
        AngEnDistKallbach(/*int EnerDistStart*/);
        virtual ~AngEnDistKallbach();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        void ConvertToEnerAndAngDist(EnergyDist **enDist, AngularDist **angDist, int &numAngEner);

        int numRegs, numIncEner /*startEnerDist*/;
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
        double outAng[numDistSample];
        double **outEner, **outEnerProb, **outEnerSumProb, **rFraction, **angDistSlope;
        double ***outAngProb;
        //outEner outgoing energy
        //outEnerProb probability of this happening
        //outEnerSumProb cumulative sum of outEnerProb up to this incoming neutron energy
    protected:
    private:
};

#endif // AngEnDistKallbach_HH
