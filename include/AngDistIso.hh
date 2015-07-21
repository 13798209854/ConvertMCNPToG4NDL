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
        string IdentifyYourSelf()
        {
            return "AngDistIso";
        }
        void SetPoint(stringstream &stream, int &count, double incNEner)
        {
            incNEnerVec.push_back(incNEner);
        }
        /*
        void AddData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2)
        {
            for(int i=0; i<int(incNEnerVec.size()); i++)
            {
                for(int j=0; j<int(enerVec.size()); j++)
                {
                    if(incNEnerVec[i]<enerVec[j])
                    {
                        enerVec.insert(enerVec.begin()+j, incNEnerVec[i]);
                        numAngProb2.insert(numAngProb2.begin()+j, 2);
                        intSchemeAng2.insert(intSchemeAng2.begin()+j, 1);
                        angVec2.insert(angVec2.begin()+j, new double [2]);
                        angProbVec2.insert(angProbVec2.begin()+j, new double [2]);
                        angVec2[j][0]=-1;
                        angProbVec2[j][0]=1;

                        angVec2[j][1]=1;
                        angProbVec2[j][1]=1;
                        break;
                    }
                    else if(incNEnerVec[i]==enerVec[j])
                    {
                        double sum2=0.;
                        for(int m=0; m<numAngProb2[j]; m++)
                        {
                            sum2 += angProbVec2[j][m];
                        }

                        for(int k=0; k<numAngProb2[j]; k++)
                        {
                            angProbVec2[j][k]/=sum2;
                            angProbVec2[j][k]+=1;
                        }
                        break;
                    }
                    else if(j==int(enerVec.size()))
                    {
                        j++;
                        enerVec.insert(enerVec.begin()+j, incNEnerVec[i]);
                        numAngProb2.insert(numAngProb2.begin()+j, 2);
                        intSchemeAng2.insert(intSchemeAng2.begin()+j, 1);
                        angVec2.insert(angVec2.begin()+j, new double [2]);
                        angProbVec2.insert(angProbVec2.begin()+j, new double [2]);
                        angVec2[j][0]=-1;
                        angProbVec2[j][0]=1;

                        angVec2[j][1]=1;
                        angProbVec2[j][1]=1;
                        break;
                    }
                }
            }
        }
        void AddData(AngularDist *secDist)
        {

        }
        */
        void AddAngleVec(vector<double> &temp, double incNEner)
        {
            int low=0;
            for(; low<int(incNEnerVec.size()-1); low++)
            {
                if(incNEnerVec[low]>incNEner)
                {
                    break;
                }
            }
            if((low==0)||(low==int(incNEnerVec.size()-1)))
                return;
            else
                low--;

            int i, j, cond=low+2;
            double angVec[2]={-1,1};
            for(; low<cond; low++)
            {
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
        }

        double GetAngleProb(double incNEner, double angle)
        {
            int low=0;
            for(; low<int(incNEnerVec.size()); low++)
            {
                if(incNEnerVec[low]>incNEner)
                {
                    break;
                }
            }
            if((low==0)||(low==int(incNEnerVec.size())))
                return 0.;

            return 1.0;
        }

        void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)
        {

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
                angVec2[i][0]=-1;
                angProbVec2[i][0]=1;

                angVec2[i][1]=1;
                angProbVec2[i][1]=1;
            }
        }

    protected:
    private:
};

#endif // ANGDISTISO_HH
