#ifndef NYield1DTab_HH
#define NYield1DTab_HH

#include "YieldDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the tabular out-going neutron yield distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class NYield1DTab: public YieldDist
{
    public:
        NYield1DTab();
        NYield1DTab(YieldDist *totalYield);
        virtual ~NYield1DTab();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        void SubtractPrompt(YieldDist* promptYieldDist);
        void SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield);
        void SubtractPrompt(double *totalCoeff, int &totalNumCoeff)
        {

        }
        void ConvertToLinDist(int *regEndPos, int &numIncEner, double *incEner, double *yield)
        {

        }
        string IdentifyYourSelf()
        {
            return "NYield1DTab";
        }

        int numRegs, numIncEner;
        int *regEndPos, *intScheme;
        double *incEner; // contains incoming neutron energy
        double *yield;
    protected:
    private:
};

#endif // NYield1DTab_HH
