#ifndef NYIELDPOLYFUNC_HH
#define NYIELDPOLYFUNC_HH

#include "YieldDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the polynomial expansion of the out-going neutron yield distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class NYieldPolyFunc: public YieldDist
{
    public:
        NYieldPolyFunc();
        virtual ~NYieldPolyFunc();
        void WriteG4NDLData(stringstream &stream);
        void ExtractMCNPData(stringstream &stream, int &count);
        void SubtractPrompt(YieldDist* &promptYieldDist);
        void SubtractPrompt(double *&totalCoeff, int &totalNumCoeff);
        void SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield);
        void ConvertToLinDist(int *regEndPos, int &numIncEner, double *&incEner, double *&yield);

        string IdentifyYourSelf()
        {
            return "NYieldPolyFunc";
        }
        int numCoeff;
        double *coeff;
    protected:
    private:
};

#endif // NYIELDPOLYFUNC_HH
