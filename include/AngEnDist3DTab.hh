#ifndef ANGENDIST3DTAB_HH
#define ANGENDIST3DTAB_HH

#include "AngularEnergyDist.hh"


class AngEnDist3DTab : public AngularEnergyDist
{
    public:
        AngEnDist3DTab(int EnerDistStart);
        virtual ~AngEnDist3DTab();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

        int numRegs, numIncEner, startEnerDist;
        //numRegs number of interpolation regions
        //numIncEner number of incoming neutron energy
        int *regEndPos, *intScheme1, *intScheme2, *outEnDistPos, *numPEnerPoints, *numDiscreteEnerPoints;
        //regEndPos containes the end positon for each region
        //intScheme1 is the scheme for interpolating between incEner points
        //intScheme2 is the scheme for interpolating between outEner points
        //distPos contains the positions of the outgoing energy-probability data for each incoming energy
        //numPEnerPoints is the total number of probability energy points
        //numDiscreteEnerPoints is the number of probability energy points that are discrete, numDiscreteEnerPoints<=numPEnerPoints
        double *incEner; // contains incoming neutron energy
        double **outEner, **outEnerProb, **outEnerSumProb;
        //outEner outgoing energy
        //outProb probability of this happening
        //outSumProb cumulative sum of outProb up to this incoming neutron energy
        int **outAngDistPos, **intScheme3, **numPAngPoints;
        bool **isoDist;
        double ***outAng, ***outAngProb, ***outAngSumProb;

    protected:
    private:
};

#endif // ANGENDIST3DTAB_HH
