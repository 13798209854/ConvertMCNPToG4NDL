#ifndef NYIELDPOLYFUNC_HH
#define NYIELDPOLYFUNC_HH

#include "YieldDist.hh"

class NYieldPolyFunc: public YieldDist
{
    public:
        NYieldPolyFunc();
        virtual ~NYieldPolyFunc();
        void WriteG4NDLData(stringstream data);
        void ExtractMCNPData(stringstream stream, int &count);
        void SubtractPrompt(YieldDist* promptYieldDist);
        void SubtractPrompt(double *totalCoeff, int &totalNumCoeff);
        void SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield);
        void ConvertToLinDist(int *regEndPos, int &numIncEner, double *incEner, double *yield);

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
