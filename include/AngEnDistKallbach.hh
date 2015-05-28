#ifndef AngEnDistKallbach_HH
#define AngEnDistKallbach_HH

#include "AngularDist.hh"


class AngEnDistKallbach : public AngularDist
{
    public:
        AngEnDistKallbach(int EnerDistStart);
        virtual ~AngEnDistKallbach();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

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
        double **outEner, **outProb, **outSumProb, **rFraction, **angDistSlope;
        //outEner outgoing energy
        //outProb probability of this happening
        //outSumProb cumulative sum of outProb up to this incoming neutron energy
    protected:
    private:
};

#endif // AngEnDistKallbach_HH
