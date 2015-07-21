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
        AngDist2DTabular(AngularDist* angDist);
        virtual ~AngDist2DTabular();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
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

        /*
        void AddData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2)
        {
            // change algorythm to include CS weigthing and add the distributions together rather than just adding on the points
            for(int i=0; i<int(incNEnerVec.size()); i++)
            {
                for(int j=0; j<int(enerVec.size()); j++)
                {
                    if(incNEnerVec[i]<enerVec[j])
                    {
                        enerVec.insert(enerVec.begin()+j, incNEnerVec[i]);
                        numAngProb2.insert(numAngProb2.begin()+j, numAngProb[i]);
                        intSchemeAng2.insert(intSchemeAng2.begin()+j, intSchemeAng[i]);
                        angVec2.insert(angVec2.begin()+j, angVec[i]);
                        angProbVec2.insert(angProbVec2.begin()+j, angProbVec[i]);
                        break;
                    }
                    else if(incNEnerVec[i]==enerVec[j])
                    {
                        vector<double> tempAng, tempAngProb;
                        double *newAng, *newProb, sum=0., sum2=0.;
                        int k=0,l=0;

                        for(int m=0; m<numAngProb[i]; m++)
                        {
                            sum += angProbVec[i][m];
                        }

                        for(int m=0; m<numAngProb2[j]; m++)
                        {
                            sum2 += angProbVec2[j][m];
                        }

                        while((k<numAngProb2[j])||(l<numAngProb[i]))
                        {
                            if(k==numAngProb2[j])
                            {
                                tempAng.push_back(angVec[i][l]);
                                tempAngProb.push_back(angProbVec[i][l]/sum);
                                tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng2[j], tempAng.back(), angVec2[j][k-2], angVec2[j][k-1], angProbVec2[j][k-2], angProbVec2[j][k-1])/sum2;
                                l++;
                            }
                            else if(l==numAngProb[i])
                            {
                                tempAng.push_back(angVec2[j][k]);
                                tempAngProb.push_back(angProbVec2[j][k]/sum2);
                                tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng[i], tempAng.back(), angVec[i][l-2], angVec[i][l-1], angProbVec[i][l-2], angProbVec[i][l-1])/sum;
                                k++;
                            }
                            else if(angVec2[j][k]>angVec[i][l])
                            {
                                tempAng.push_back(angVec[i][l]);
                                tempAngProb.push_back(angProbVec[i][l]/sum);
                                if(k!=0)
                                    tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng2[j], tempAng.back(), angVec2[j][k-1], angVec2[j][k], angProbVec2[j][k-1], angProbVec2[j][k])/sum2;
                                else
                                    tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng2[j], tempAng.back(), angVec2[j][k], angVec2[j][k+1], angProbVec2[j][k], angProbVec2[j][k+1])/sum2;
                                l++;
                            }
                            else if(angVec2[j][k]<angVec[i][l])
                            {
                                tempAng.push_back(angVec2[j][k]);
                                tempAngProb.push_back(angProbVec2[j][k]/sum2);
                                if(l!=0)
                                    tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng[i], tempAng.back(), angVec[i][l-1], angVec[i][l], angProbVec[i][l-1], angProbVec[i][l])/sum;
                                else
                                    tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng[i], tempAng.back(), angVec[i][l], angVec[i][l+1], angProbVec[i][l], angProbVec[i][l+1])/sum;
                                k++;
                            }
                            else
                            {
                                tempAng.push_back(angVec2[j][k]/sum2);
                                tempAngProb.push_back(angProbVec2[j][k]/sum2+angProbVec[i][l]/sum);
                                k++; l++;
                            }
                        }
                        newAng = new double [tempAng.size()];
                        newProb = new double [tempAng.size()];
                        numAngProb2[i]=tempAng.size();
                        for(int h=0; h<int(tempAng.size()); h++)
                        {
                            newAng[h]=tempAng[h];
                            newProb[h]=tempAngProb[h];
                        }
                        delete [] angVec2[j];
                        delete [] angProbVec2[j];
                        angVec2[j]=newAng;
                        angProbVec2[j]=newProb;
                        break;
                    }
                    else if(j==int(enerVec.size()))
                    {
                        j++;
                        enerVec.insert(enerVec.begin()+j, incNEnerVec[i]);
                        numAngProb2.insert(numAngProb2.begin()+j, numAngProb[i]);
                        intSchemeAng2.insert(intSchemeAng2.begin()+j, intSchemeAng[i]);
                        angVec2.insert(angVec2.begin()+j, angVec[i]);
                        angProbVec2.insert(angProbVec2.begin()+j, angProbVec[i]);
                        break;
                    }
                }
            }
        }
        void AddData(AngularDist *secDist)
        {
            secDist->AddData(incNEnerVec, angVec, angProbVec, intSchemeAng, numAngProb);
        }
        */
        void AddAngleVec(vector<double> &temp, double incNEner)
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
                return;
            else
                low--;

            int i, j, cond=low+2;
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
            for(; low<int(incNEnerVec.size()); low++)
            {
                if(incNEnerVec[low]>incNEner)
                {
                    break;
                }
            }
            if((low==0)||(low==int(incNEnerVec.size())))
                return 0.;
            else
                low--;

            int count=0, cond=low+2;
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

                prob[count]=AngularDist::Interpolate(intSchemeAng[low], angle, angVec[low][lowAng], angVec[low][lowAng+1], angProbVec[low][lowAng], angProbVec[low][lowAng+1])/sumProb;
                count++;
            }
            return AngularDist::Interpolate(intSchemeAng[low-2], incNEner, incNEnerVec[low-2], incNEnerVec[low-1], prob[0], prob[1]);
        }

        void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)
        {
            double sumCS, curCS;

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
                    angDistList[i][j]->AddEnergyVec(incNEnerVec);
                }
            }
            numAngEner=incNEnerVec.size();

            //we set the interpolation to always be linear
            intSchemeAng.assign(numAngEner,2);

            for(int m=0; m<int(incNEnerVec.size()); m++)
            {
                sumCS=0.;
                for(int i=startList; i<endList; i++)
                {
                    if(nCSDistList[i])
                    {
                        sumCS+=max(0.,nCSDistList[i]->Interpolate(incNEnerVec[m]));
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
                for(int k=0; k<numAngProb.back(); k++)
                {
                    for(int i=startList; i<endList; i++)
                    {
                        if(nCSDistList[i])
                        {
                            curCS = max(0.,nCSDistList[i]->Interpolate(incNEnerVec[m]));
                            for(int j=0; j<int(angDistList[i].size()); j++)
                            {
                                angProbVec.back()[k]+=curCS*angDistList[i][j]->GetAngleProb(incNEnerVec[m], angVec.back()[k])/sumCS;
                            }
                        }
                    }
                }
            }
        }
        void SetData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2, double &temp)
        {

        }
        vector <double*> angVec, angProbVec;
        vector<int> intSchemeAng, numAngProb;
        double dummy;
    protected:
    private:
};

#endif // AngDist2DTabular_HH
