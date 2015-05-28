#ifndef ENERDISTLEVSCAT_HH
#define ENERDISTLEVSCAT_HH

#include "EnergyDist.hh"

class EnerDistLevScat: public EnergyDist
{
    public:
        EnerDistLevScat();
        virtual ~EnerDistLevScat();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        string IdentifyYourSelf()
        {
            return "EnerDistLevScat";
        }

        double firstHalfEq, secondHalfEq, enerStart, enerEnd;
    protected:
    private:
};

#endif // ENERDISTLEVSCAT_HH
