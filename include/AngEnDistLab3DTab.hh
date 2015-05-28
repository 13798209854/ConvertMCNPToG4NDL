#ifndef AngEnDistLab3DTab_HH
#define AngEnDistLab3DTab_HH

#include "AngularDist.hh"


class AngEnDistLab3DTab : public AngularDist
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
        double ***outEner, ***outProb, ***outSumProb;
    protected:
    private:
};

#endif // AngEnDistLab3DTab_HH
