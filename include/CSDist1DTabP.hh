#ifndef CSDIST1DTABP_HH
#define CSDIST1DTABP_HH

#include "CSDist.hh"
/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy dependant tabular out-going photon cross-section distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class CSDist1DTabP: public CSDist
{
    public:
        CSDist1DTabP();
        virtual ~CSDist1DTabP();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLCSData(stringstream &stream);
        void WriteG4NDLYieldData(stringstream &stream);
        string IdentifyYourSelf();
        double Interpolate(double x);
        double GetAvgCS();
        double GetAvgCS(double ener);
        void SetCSData(CSDist* nCSData)
        {
            nCSData->SetCSData(enerCSVec, csEnerStart, csVec, csSize);
        }
        void SetCSData(double* &enerCSVecSet, int &csEnerStartSet, double* &csVecSet, int &csSizeSet)
        {
            cout << "this function has not been implemented" << endl;
        }
        void AddData(CSDist *secDist)
        {
            cout << "this function has not been implemented" << endl;
        }
        void AddData(int &csEnerStartSet, double* &csVecSet, int &csSizeSet)
        {
            cout << "this function has not been implemented" << endl;
        }
        int startEner, CSVecSize, csSize, mtNum, csEnerStart;
        double *CSVec, *csVec, *enerCSVec;
        bool elastic=false;
    protected:
    private:
};

#endif // CSDIST1DTABP_HH
