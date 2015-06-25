#ifndef CSDIST1DTAB_HH
#define CSDIST1DTAB_HH

#include "CSDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy dependant tabular out-going neutron cross-section distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class CSDist1DTab: public CSDist
{
    public:
        CSDist1DTab(double *energyVec);
        CSDist1DTab(double *energyVec, int numCSEner);
        virtual ~CSDist1DTab();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLCSData(stringstream &stream);
        void WriteG4NDLYieldData(stringstream &stream)
        {

        }
        void SetCSData(CSDist* nCSData)
        {

        }
        double Interpolate(double x);
        void SetCSData(double* &enerCSVecSet, int &csEnerStartSet, double* &csVecSet, int &csSizeSet)
        {
            enerCSVecSet=enerVec;
            csEnerStartSet=startEner;
            csVecSet=CSVec;
            csSizeSet=CSVecSize;
        }
        string IdentifyYourSelf();
        /*double* GetEnergyVec()
        {
            return (enerVec+startEner-1);
        }
        double* GetCSVec()
        {
            return CSVec;
        }
        int GetCSSize()
        {
            return CSVecSize;
        }*/

        int startEner, CSVecSize;
        double *CSVec, *enerVec;
        bool elastic=false;
    protected:
    private:
};

#endif // CSDIST1DTAB_HH
