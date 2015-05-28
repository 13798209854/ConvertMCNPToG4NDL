#ifndef AngDist2DTabular.cc_HH
#define AngDist2DTabular.cc_HH

#include "AngularDist.hh"

class AngDist2DTabular.cc: public AngularDist
{
    public:
        AngDist2DTabular.cc();
        virtual ~AngDist2DTabular.cc();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);
        void SetPoint(stringstream stream, int &count, double incNEner)
        {
            //the angular distribution is represented as a table of cosine and prob
            int intTemp;
            double temp;

            incNEnerVec.push_back(incNEner);
            stream >> intTemp;
            count++;
            intSchemeAng.push_back(intTemp);

            stream >> intTemp;
            count++;
            numAngProb.push_back(intTemp);

            angVec.push_back(new double[intTemp]);
            angProbVec.push_back(new double[intTemp]);

            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> temp;
                angVec.back()[k]=temp;
            }
            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> temp;
                angProbVec.back()[k]=temp;
            }
            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> dummy;
            }
        }
        vector <double*> angVec, angProbVec;
        vector<int> intSchemeAng, numAngProb;
    protected:
    private:
};

#endif // AngDist2DTabular.cc_HH
