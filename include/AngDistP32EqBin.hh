#ifndef ANGDISTP32EQBIN_Hh
#define ANGDISTP32EQBIN_Hh

#include "AngularDist.hh"

class AngDistP32EqBin: public AngularDist
{
    public:
        AngDistP32EqBin(int startAngularPDist);
        virtual ~AngDistP32EqBin();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        void SetPoint(stringstream stream, int &count, double incNEner)
        {
            incNEnerVec.push_back(incNEner);
            double temp;

                //the angular distribution is represented by a probability density function where each of the 33 points, 32 bins, are spaced along the entire cosine of the scattering angle axis
            //such that the integral of the probability density with the cosine of the scattering angle inbetween them is kept constant, meaning that the bins are in a equilprobable arrangement
            // the first bin probably corresponds to a cosine value of -1 and and the last bin being some where around 1 depending on the spacing

            angVec.push_back(new double [33]);
            for(int k=0; k<33; k++, count++)
            {
                stream >> temp;
                angVec.back()[k]=temp;
            }
        }

        /*int startAngPDist;
        int numIncEner;
        double *incEner;
        int *probDistPos;*/
        vector<double*> angVec;
    protected:
    private:
};

#endif // ANGDISTP32EQBIN_Hh
