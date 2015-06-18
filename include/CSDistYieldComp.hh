#ifndef CSDISTYIELDCOMP_HH
#define CSDISTYIELDCOMP_HH

#include "CSDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy dependant tabular out-going photon multiplicity/cross-section distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class CSDistYieldComp : public CSDist
{
    public:
        CSDistYieldComp();
        virtual ~CSDistYieldComp();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLYieldData(stringstream &stream);
        void WriteG4NDLCSData(stringstream &stream);
        string IdentifyYourSelf();
        double Interpolate(double x);
        void SetCSData(CSDist* nCSData)
        {
            nCSData->SetCSData(enerCSVec, csEnerStart, csVec, csSize);
        }
        void SetCSData(double* enerCSVec, int &csEnerStart, double* csVec, int &csSize)
        {

        }
        int mtMulti, numRegs, mtNum, numIncEner, csSize, csEnerStart;
        int *regEndPos, *intScheme;
        double *yieldVec, *enerVec, *enerCSVec, *csVec;
    protected:
    private:
};

#endif // CSDISTYIELDCOMP_HH
