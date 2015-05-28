#ifndef ENERDIST1PHENER_HH
#define ENERDIST1PHENER_HH

#include "EnergyDist.hh"

class EnerDist1PhEner: public EnergyDist
{
    public:
        EnerDist1PhEner(double AWR);
        virtual ~EnerDist1PhEner();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        double GetAverageOutEnergy()
        {
            if(photonType!=2)
            {
                return photonEN;
            }
            else
            {
                // assume average incoming neutron energy is 1eV
                return (photonEN+1*awr/(awr+1));
            }
        }
        string IdentifyYourSelf()
        {
            return "EnerDist1PhEner";
        }

        double photonType // states whether the photn is a primary or non primary photon,
        //which also tells us whether the energy is the photon energy (0,1) or the binding energy (2)
        double photonEN, awr;
    protected:
    private:
};

#endif // ENERDIST1PHENER_HH
