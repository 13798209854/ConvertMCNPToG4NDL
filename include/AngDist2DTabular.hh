#ifndef AngDist2DTabular_HH
#define AngDist2DTabular_HH

#include "AngularDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy independant tabular neutron out-going angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class AngDist2DTabular: public AngularDist
{
    public:
        AngDist2DTabular();
        virtual ~AngDist2DTabular();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        void SetPoint(stringstream &stream, int &count, double incNEner)
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
        double dummy;
    protected:
    private:
};

#endif // AngDist2DTabular_HH
