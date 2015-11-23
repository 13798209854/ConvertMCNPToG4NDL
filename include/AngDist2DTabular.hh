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
        AngDist2DTabular(int numIncEnerTemp, double *incEnerTemp, int *intSchemeTemp, vector<double> *outAngTemp, vector<double> *outAngProbTemp);
        AngDist2DTabular(AngularDist* angDist);
        virtual ~AngDist2DTabular();
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
                    if(sum!=sum)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        string IdentifyYourSelf()
        {
            return "AngDist2DTabular";
        }
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
            angProbSumVec.push_back(new double[intTemp]);

            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> temp;
                angVec.back()[k]=temp;
            }
            double sum=0.;
            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> temp;
                angProbVec.back()[k]=temp;
                sum+=angProbVec.back()[k];
            }
            if(sum==0.)
            {
                cout << "Error in the collection of the angular data : AngDist2DTabular.hh:61" << endl;
            }
            for(int k=0; k<numAngProb.back(); k++, count++)
            {
                stream >> temp;
                angProbSumVec.back()[k]=temp;
            }
            for(int k=0; k<numAngProb.back(); k++)
            {
                // here we set correct the angular prob so that it is integrated over its angular regime
                if(angVec.back()[k]==0.0)
                    angVec.back()[k]=1.0e-12;

                if(k==0)
                {
                    if(numAngProb.back()==2)
                        angProbVec.back()[k] = angProbSumVec.back()[k]-angProbSumVec.back()[k-1];
                    else
                        angProbVec.back()[k] = angProbSumVec.back()[k];
                }
                else
                    angProbVec.back()[k]=angProbSumVec.back()[k]-angProbSumVec.back()[k-1];
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
            if(incNEnerVec.size()==1)
                return max(0.,prob[0]);
            else
                return max(0.,AngularDist::Interpolate(2, incNEner, incNEnerVec[cond-2], incNEnerVec[cond-1], prob[0], prob[1]));
        }
        void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)
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
            for(int i=startList; i<endList; i++)
            {
                for(int j=0; j<int(angDistList[i].size()); j++)
                {
                    if(!(angDistList[i][j]->CheckData()))
                    {
                        cout << "Error in angular data AngDist2DTabular.hh:234" << endl;
                    }
                    angDistList[i][j]->AddEnergyVec(incNEnerVec);
                }
            }
            numAngEner+=incNEnerVec.size();

            //we set the interpolation to always be linear
            intSchemeAng.assign(incNEnerVec.size(),2);

            for(int m=0; m<int(incNEnerVec.size()); m++)
            {
                sumCS=0.;
                for(int i=startList; i<endList; i++)
                {
                    if(nCSDistList[i])
                    {
                        sumCS+=max(0.,nCSDistList[i]->GetAvgCS());
                        for(int j=0; j<int(angDistList[i].size()); j++)
                        {
                            angDistList[i][j]->AddAngleVec(temp, incNEnerVec[m]);
                        }
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
                        for(int i=startList; i<endList; i++)
                        {
                            if(nCSDistList[i])
                            {
                                curCS = max(0.,nCSDistList[i]->GetAvgCS());
                                for(int j=0; j<int(angDistList[i].size()); j++)
                                {
                                    angProbVec.back()[k]+=curCS*angDistList[i][j]->GetAngleProb(incNEnerVec[m], angVec.back()[k])/sumCS;
                                }
                            }
                        }
                    }
                    for(int k=0; k<numAngProb.back(); k++)
                    {
                        sumCheck += angProbVec.back()[k];
                    }
                    if(sumCheck<=0.)
                    {
                        cout << "Error in the summation of the angular probability data : AngDist2DTabular.hh:259" << endl;
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
        void SumAngularData(vector<AngularDist*> &angDistList, vector<CSDist*> &pCSDistList, int &numAngEner)
        {
            cout << "this function has not been implemented" << endl;
        }
        void SetData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2, double &temp)
        {
            cout << "this function has not been implemented" << endl;
        }
        vector <double*> angVec, angProbVec, angProbSumVec;
        vector<int> intSchemeAng, numAngProb;
        double dummy;
    protected:
    private:
};

#endif // AngDist2DTabular_HH
