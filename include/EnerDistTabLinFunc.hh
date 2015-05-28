#ifndef ENERDISTTABLINFUNC_HH
#define ENERDISTTABLINFUNC_HH

#include "EnergyDist.hh"


class EnerDistTabLinFunc : public EnergyDist
{
    public:
        EnerDistTabLinFunc(int EnerDistStart);
        virtual ~EnerDistTabLinFunc();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
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
