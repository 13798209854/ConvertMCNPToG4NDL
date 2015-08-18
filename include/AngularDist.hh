#ifndef ANGULARDIST_HH
#define ANGULARDIST_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "../include/CSDist.hh"

using namespace std;
/*
Created By: Wesley Ford June 17, 2015

This class is the mother class for all the classes responsible for the extraction of the energy independant out-going angular distribution data from MCNP
*/

class AngularDist
{
    public:
        AngularDist();
        virtual ~AngularDist();
        virtual void ExtractMCNPData(stringstream &stream, int &count)=0;
        virtual void WriteG4NDLData(stringstream &stream)=0;
        virtual void SetPoint(stringstream &stream, int &count, double incNEner)=0;
        virtual void AddEnergyVec(vector<double> &incNEnerVecSum)
        {
            int i=0, j=0;
            while((i<int(incNEnerVec.size()))&&(j<int(incNEnerVecSum.size())))
            {
                if(incNEnerVec[i]<incNEnerVecSum[j])
                {
                    incNEnerVecSum.insert(incNEnerVecSum.begin()+j, incNEnerVec[i]);
                    i++; j++;
                }
                else if(incNEnerVec[i]>incNEnerVecSum[j])
                {
                    j++;
                }
                else
                {
                    i++; j++;
                }
            }

            for(;i<int(incNEnerVec.size());i++)
            {
                incNEnerVecSum.push_back(incNEnerVec[i]);
            }
        }
        virtual void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)=0;
        virtual void SumAngularData(vector<AngularDist*> &angDistList, vector<CSDist*> &pCSDistList, int &numAngEner)=0;
        virtual void AddAngleVec(vector<double> &temp, double incNEner)=0;
        virtual double GetAngleProb(double incNEner, double angle)=0;
        virtual bool CheckData()=0;
        /*
        virtual void AddData(AngularDist *secDist)=0;
        virtual void AddData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2)=0;
        */
        virtual void SetData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2, double &temp)=0;
        virtual string IdentifyYourSelf()=0;
        void SetTemperature(double temp) {temperature=temp;}

        double Interpolate(int aScheme, double x, double x1, double x2, double y1, double y2) const;
        double Histogram(double , double , double , double y1, double ) const;
        double LinearLinear(double x, double x1, double x2, double y1, double y2) const;
        double LinearLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLinear(double x, double x1, double x2, double y1, double y2) const;
        double LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2) const;
        double Random(double , double , double , double y1, double y2) const;

        vector<double> incNEnerVec;
        double temperature=0;
    protected:
    private:
};

#endif // ANGULARDIST_HH
