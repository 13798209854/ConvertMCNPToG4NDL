#ifndef ANGENDISTNBODY_HH
#define ANGENDISTNBODY_HH

#include "AngularEnergyDist.hh"


class AngEnDistNBody : public AngularEnergyDist
{
    public:
        AngEnDistNBody();
        virtual ~AngEnDistNBody();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

        int numBodies;
        double particleMassRatio;
    protected:
    private:
};

#endif // ANGENDISTNBODY_HH
