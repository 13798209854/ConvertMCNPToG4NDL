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
        int numCoeff;
        double *coeff;
    protected:
    private:
};

#endif // NYIELDPOLYFUNC_HH
