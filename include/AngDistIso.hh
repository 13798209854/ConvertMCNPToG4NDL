#ifndef ANGDISTISO_HH
#define ANGDISTISO_HH

#include "AngularDist.hh"

class AngDistIso: public AngularDist
{
    public:
        AngDistIso();
        virtual ~AngDistIso();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        void SetPoint(stringstream stream, int &count, double incNEner)
        {
            incNEnerVec.push_back(incNEner);
        }
    protected:
    private:
};

#endif // ANGDISTISO_HH
