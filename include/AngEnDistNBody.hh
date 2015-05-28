#ifndef ANGENDISTNBODY_HH
#define ANGENDISTNBODY_HH

#include "AngularDist.hh"


class AngEnDistNBody : public AngularDist
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
