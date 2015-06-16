#ifndef AngEnDistLab3DTab_HH
#define AngEnDistLab3DTab_HH

#include "AngularEnergyDist.hh"


class AngEnDistLab3DTab : public AngularEnergyDist
{
    public:
        AngEnDistLab3DTab();
        virtual ~AngEnDistLab3DTab();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

        int numRegs, numIncEner, startEnerDist;
        int *regEndPos, *intScheme1, *intScheme2, *outAngDistPos, *numPAngPoints;
        double *incEner; // contains incoming neutron energy
        double **outAng;
        int **outEnerDistPos, **intScheme3, **numPEnerPoints;
        double ***outEner, ***outAngProb, ***outAngSumProb;
    protected:
    private:
};

#endif // AngEnDistLab3DTab_HH
