#ifndef AngDist2DTabularP_HH
#define AngDist2DTabularP_HH

#include "AngularDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy independant tabular neutron out-going angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class AngDist2DTabularP: public AngularDist
{
    public:
        AngDist2DTabularP();
        AngDist2DTabularP(AngularDist* angDist);
        virtual ~AngDist2DTabularP();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        bool CheckData()
        {
            double sum;
            if(incNEnerVec.size()==0)
            {
                return false;
            }
            else
            {
                for(int i=0; i<int(incNEnerVec.size()); i++)
                {
                    sum=0.;
                    for(int j=0; j<numAngProb[i]; j++)
                    {
                        sum+=angProbVec[i][j];
                    }
                    if(sum==0.)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        string IdentifyYourSelf()
        {
            return "AngDist2DTabularP";
        }
        void SetPoint(stringstream &stream, int &count, double incNEner)
        {
            //the angular distribution is represented as a table of cosine and prob
            int intTemp;
            double temp, sum;

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
            sum=0.;
            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> temp;
                angProbVec.back()[k]=temp;
                sum+=temp;
            }
            if(sum==0.)
            {
                cout << "Error in the collection of the angular data AngDist2DTabularP.hh:60" << endl;
            }
            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> dummy;
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

            int i, j, cond;
            if(incNEnerVec.size()==1)
                cond=low+1;
            else
                cond=low+2;

            for(; low<cond; low++)
            {
                i=0; j=0;
                while((i<numAngProb[low])&&(j<int(temp.size())))
                {
                    if(angVec[low][i]<temp[j])
                    {
                        temp.insert(temp.begin()+j, angVec[low][i]);
                        i++; j++;
                    }
                    else if(angVec[low][i]>temp[j])
                    {
                        j++;
                    }
                    else
                    {
                        i++; j++;
                    }
                }

                for(;i<numAngProb[low];i++)
                {
                    temp.push_back(angVec[low][i]);
                }
            }
        }
        double GetAngleProb(double incNEner, double angle)
        {
            int low=0, lowAng;
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

            int count=0, cond;
            if(incNEnerVec.size()==1)
                cond=low+1;
            else
                cond=low+2;
            double sumProb;
            double prob[2];
            for(; low<cond; low++)
            {
                sumProb=0.;
                lowAng=0;
                for(; lowAng<numAngProb[low]-1; lowAng++)
                {
                    if(angVec[low][lowAng]>angle)
                    {
                        break;
                    }
                }
                if(lowAng!=0)
                    lowAng--;

                for(int i=0; i<numAngProb[low]; i++)
                {
                    sumProb+=angProbVec[low][i];
                }

                if(sumProb!=0.)
                    prob[count]=AngularDist::Interpolate(intSchemeAng[low], angle, angVec[low][lowAng], angVec[low][lowAng+1], angProbVec[low][lowAng], angProbVec[low][lowAng+1])/sumProb;
                else
                    prob[count]=0.;
                count++;
            }
            return AngularDist::Interpolate(intSchemeAng[low-2], incNEner, incNEnerVec[low-2], incNEnerVec[low-1], prob[0], prob[1]);
        }
        void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)
        {
            cout << "this function has not been implemented" << endl;
        }
        void SumAngularData(vector<AngularDist*> &angDistList, vector<CSDist*> &pCSDistList, int &numAngEner)
        {
            double sumCS, curCS, sumCheck;

            intSchemeAng.clear(); numAngProb.clear(); incNEnerVec.clear();
            for(int i=0; i<int(angVec.size()); i++)
            {
                if(angVec[i])
                    delete [] angVec[i];
                angVec[i]=NULL;
            }
            for(int i=0; i<int(angProbVec.size()); i++)
            {
                if(angProbVec[i])
                    delete [] angProbVec[i];
                angProbVec[i]=NULL;
            }

            vector<double> temp;
            for(int i=0; i<int(angDistList.size()); i++)
            {
                if(!(angDistList[i]->CheckData()))
                {
                    cout << "Error in angular data AngDist2DTabularP.hh:235" << endl;
                }
                angDistList[i]->AddEnergyVec(incNEnerVec);
            }
            numAngEner=incNEnerVec.size();

            //we set the interpolation to always be linear
            intSchemeAng.assign(numAngEner,2);

            for(int m=0; m<int(incNEnerVec.size()); m++)
            {
                sumCS=0.;
                for(int i=0; i<int(angDistList.size()); i++)
                {
                    if(pCSDistList[i])
                    {
                        sumCS+=max(0.,pCSDistList[i]->GetAvgCS());
                        angDistList[i]->AddAngleVec(temp, incNEnerVec[m]);
                    }
                }
                numAngProb.push_back(temp.size());
                angVec.push_back(new double [temp.size()]);
                angProbVec.push_back(new double [temp.size()]);
                for(int i=0; i<int(temp.size()); i++)
                {
                    angVec.back()[i]=temp[i];
                    angProbVec.back()[i]=0.;
                }
                temp.clear();
                if(sumCS!=0.)
                {
                    sumCheck=0.;
                    for(int k=0; k<numAngProb.back(); k++)
                    {
                        for(int i=0; i<int(angDistList.size()); i++)
                        {
                            if(pCSDistList[i])
                            {
                                curCS = max(0.,pCSDistList[i]->GetAvgCS());
                                angProbVec.back()[k]+=curCS*angDistList[i]->GetAngleProb(incNEnerVec[m], angVec.back()[k])/sumCS;
                            }
                        }
                    }
                    for(int k=0; k<numAngProb.back(); k++)
                    {
                        sumCheck += angProbVec.back()[k];
                    }
                    if(sumCheck==0.)
                    {
                        cout << "Error in the summation of the angular probability data AngDist2DTabularP.hh:253" << endl;
                    }
                }
                else
                {
                    for(int i=0; i<numAngProb.back(); i++)
                    {
                        angProbVec.back()[i]=1.0/numAngProb.back();
                    }
                }
            }
        }
        void SetData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2, double &temp)
        {
            cout << "Error This funtion has not been implemented yet" << endl;
        }
        vector <double*> angVec, angProbVec;
        vector<int> intSchemeAng, numAngProb;
        double dummy;
    protected:
    private:
};

#endif // AngDist2DTabularP_HH
