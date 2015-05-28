#ifndef ANGDISTISOP_HH
#define ANGDISTISOP_HH

#include "AngularDist.hh"

class AngDistIsoP: public AngularDist
{
    public:
        AngDistIsoP();
        virtual ~AngDistIsoP();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        void SetPoint(stringstream stream, int &count, double incNEner)
        {
            incNEnerVec.push_back(incNEner);
        }
    protected:
    private:
};

#endif // ANGDISTISOP_HH
