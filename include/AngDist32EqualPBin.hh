#ifndef AngDist32EqualPBin_HH
#define AngDist32EqualPBin_HH

#include "AngularDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy independant 32 equal probability bin out-going angular distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class AngDist32EqualPBin: public AngularDist
{
    public:
        AngDist32EqualPBin();
        virtual ~AngDist32EqualPBin();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        string IdentifyYourSelf()
        {
            return "AngDist32EqualPBin";
        }
        void SetPoint(stringstream &stream, int &count, double incNEner)
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
                        numAngProb2.insert(numAngProb2.begin()+j, 32);
                        intSchemeAng2.insert(intSchemeAng2.begin()+j, 1);
                        angVec2.insert(angVec2.begin()+j, angVec[i]);
                        angProbVec2.insert(angProbVec2.begin()+j, new double [32]);
                        for(int k=0; k<32; k++)
                        {
                            angProbVec2[j][k]=1/(32*(angVec[i][k+1]-angVec[i][k]));
                        }
                        break;
                    }
                    else if(incNEnerVec[i]==enerVec[j])
                    {
                        vector<double> tempAng, tempAngProb;
                        double *newAng, *newProb, sum2=0.;
                        int k=0,l=0;

                        for(int m=0; m<numAngProb2[j]; m++)
                        {
                            sum2 += angProbVec2[j][m];
                        }

                        while((k<numAngProb2[j])||(l<32))
                        {
                            if(k==numAngProb2[j])
                            {
                                tempAng.push_back(angVec[i][l]);
                                tempAngProb.push_back(1/(32*(angVec[i][l+1]-angVec[i][l])));
                                tempAngProb.back()+=AngularDist::Interpolate(intSchemeAng2[j], tempAng.back(), angVec2[j][k-2], angVec2[j][k-1], angProbVec2[j][k-2], angProbVec2[j][k-1])/sum2;
                                l++;
                            }
                            else if(l==32)
                            {
                                tempAng.push_back(angVec2[j][k]);
                                tempAngProb.push_back(angProbVec2[j][k]/sum2);
                                tempAngProb.back()+=AngularDist::Interpolate(2, tempAng.back(), angVec[i][l-2], angVec[i][l-1], 1/(32*(angVec[i][l-1]-angVec[i][l-2])), 1/(32*(angVec[i][l]-angVec[i][l-1])));
                                k++;
                            }
                            else if(angVec2[j][k]>angVec[i][l])
                            {
                                tempAng.push_back(angVec[i][l]);
                                tempAngProb.push_back(1/(32*(angVec[i][l+1]-angVec[i][l])));
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
                                    tempAngProb.back()+=AngularDist::Interpolate(2, tempAng.back(), angVec[i][l-1], angVec[i][l], 1/(32*(angVec[i][l]-angVec[i][l-1])), 1/(32*(angVec[i][l+1]-angVec[i][l])));
                                else
                                    tempAngProb.back()+=AngularDist::Interpolate(2, tempAng.back(), angVec[i][l], angVec[i][l+1], 1/(32*(angVec[i][l+1]-angVec[i][l])), 1/(32*(angVec[i][l+2]-angVec[i][l+1])));
                                k++;
                            }
                            else
                            {
                                tempAng.push_back(angVec2[j][k]);
                                tempAngProb.push_back(angProbVec2[j][k]/sum2+1/(32*(angVec[i][l+1]-angVec[i][l])));
                                k++; l++;
                            }
                        }
                        newAng = new double [tempAng.size()];
                        newProb = new double [tempAng.size()];
                        numAngProb2[j]=tempAng.size();
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
                        numAngProb2.insert(numAngProb2.begin()+j, 32);
                        intSchemeAng2.insert(intSchemeAng2.begin()+j, 1);
                        angVec2.insert(angVec2.begin()+j, angVec[i]);
                        angProbVec2.insert(angProbVec2.begin()+j, new double [32]);
                        for(int k=0; k<32; k++)
                        {
                            angProbVec2[j][k]=1/(32*(angVec[i][k+1]-angVec[i][k]));
                        }
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
                while((i<32)&&(j<int(temp.size())))
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

                for(;i<32;i++)
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
                for(; lowAng<32; lowAng++)
                {
                    if(angVec[low][lowAng]>angle)
                    {
                        break;
                    }
                }
                if(lowAng!=0)
                    lowAng--;

                for(int i=0; i<32; i++)
                {
                    sumProb+=1/(32*(angVec[low][i+1]-angVec[low][i]));
                }

                prob[count]=1/(32*(angVec[low][lowAng+1]-angVec[low][lowAng]))/sumProb;
                count++;
            }
            return AngularDist::Interpolate(2, incNEner, incNEnerVec[low], incNEnerVec[low+1], prob[0], prob[1]);
        }

        void SumAngularData(vector<AngularDist*> *angDistList, CSDist **nCSDistList, int startList, int endList, int &numAngEner)
        {

        }

        void SetData(vector<double> &enerVec, vector <double*> &angVec2, vector <double*> &angProbVec2, vector<int> &intSchemeAng2, vector<int> &numAngProb2, double &temp)
        {
            temp=temperature;
            enerVec=incNEnerVec;
            numAngProb2.assign(incNEnerVec.size(),32);
            intSchemeAng2.assign(incNEnerVec.size(),1);
            angVec2.assign(incNEnerVec.size(), new double [32]);
            angProbVec2.assign(incNEnerVec.size(), new double [32]);

            for(int i=0; i<int(incNEnerVec.size()); i++)
            {
                for(int j=0; j<32; j++)
                {
                    angVec2[i][j]=angVec[i][j];
                    angProbVec2[i][j]=1/(32*(angVec[i][j+1]-angVec[i][j]));
                }
            }
        }

    vector<double*> angVec;
    protected:
    private:
};

#endif // AngDist32EqualPBin_HH
