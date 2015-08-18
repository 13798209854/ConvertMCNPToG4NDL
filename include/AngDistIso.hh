#ifndef ANGDISTISO_HH
#define ANGDISTISO_HH

#include "AngularDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the out-going energy independant isotropic neutron angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class AngDistIso: public AngularDist
{
    public:
        AngDistIso();
        virtual ~AngDistIso();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        bool CheckData()
        {
            return true;
        }
        string IdentifyYourSelf()
        {
            return "AngDistIso";
        }
        void SetPoint(stringstream &stream, int &count, double incNEner)
        {
            incNEnerVec.push_back(incNEner);
        }
        void AddEnergyVec(vector<double> &incNEnerVecSum)
        {
            if(incNEnerVec.size()==0)
            {
                incNEnerVec.push_back(0.);
                incNEnerVec.push_back(20);
            }
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
        void AddAngleVec(vector<double> &temp, double incNEner)
        {
            if(incNEnerVec.size()<1)
                return;

            int low=0;
            for(; low<int(incNEnerVec.size()-1); low++)
            {
                if(incNEnerVec[low]>incNEner)
                {
                    break;
                }
            }
            if(low>0)
                low--;

            int i, j;
            double angVec[2]={-1,1};
            i=0; j=0;
            while((i<2)&&(j<int(temp.size())))
            {
                if(angVec[i]<temp[j])
                {
                    temp.insert(temp.begin()+j, angVec[i]);
                    i++; j++;
                }
                else if(angVec[i]>temp[j])
                {
                    j++;
                }
                else
                {
                    i++; j++;
                }
            }

            for(;i<2;i++)
            {
                temp.push_back(angVec[i]);
            }
        }

        double GetAngleProb(double incNEner, double angle)
        {
            int low=0;
            if(incNEnerVec.size()<1)
                return 0.;
            else if(incNEnerVec.size()==1)
            {
                if(incNEnerVec[0]!=incNEner)
                {
                    return 0.;
                }
            }
            else if(incNEnerVec.size()==2)
            {
                if((incNEnerVec[0]>incNEner)||(incNEnerVec[1]<incNEner))
                {
                    return 0.;
                }
            }
            else
            {
                if((incNEnerVec[0]>incNEner)||(incNEnerVec[incNEnerVec.size()-1]<incNEner))
                {
                    return 0.;
                }
                for(; low<int(incNEnerVec.size()-1); low++)
                {
                    if(incNEnerVec[low]>incNEner)
                    {
                        break;
                    }
                }
                if(low>0)
                    low--;
            }

            return 1.0;
        }

        void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)
        {
            cout << "this function has not been implemented" << endl;
        }

        void SumAngularData(vector<AngularDist*> &angDistList, vector<CSDist*> &pCSDistList, int &numAngEner)
        {
            cout << "this function has not been implemented" << endl;
        }

        void SetData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2, double &temp)
        {
            temp = temperature;
            enerVec=incNEnerVec;
            numAngProb2.assign(incNEnerVec.size(),2);
            intSchemeAng2.assign(incNEnerVec.size(),1);
            angVec2.assign(incNEnerVec.size(), new double [2]);
            angProbVec2.assign(incNEnerVec.size(), new double [2]);

            for(int i=0; i<int(incNEnerVec.size()); i++)
            {
                angVec2[i][0]=-1.0;
                angProbVec2[i][0]=1.0;

                angVec2[i][1]=1.0;
                angProbVec2[i][1]=1.0;
            }
        }

    protected:
    private:
};

#endif // ANGDISTISO_HH
