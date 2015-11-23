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
        bool CheckData()
        {
            if(incNEnerVec.size()>0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
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
                    sumProb+=1.0/(32*(angVec[low][i+1]-angVec[low][i]));
                }

                prob[count]=1.0/(32*(angVec[low][lowAng+1]-angVec[low][lowAng]))/sumProb;
                count++;
            }
            if(incNEnerVec.size()==1)
                return max(0.,prob[0]);
            else
                return max(0.,AngularDist::Interpolate(2, incNEner, incNEnerVec[cond-2], incNEnerVec[cond-1], prob[0], prob[1]));
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
                    angProbVec2[i][j]=1.0/(32*(angVec[i][j+1]-angVec[i][j]));
                }
            }
        }

    vector<double*> angVec;
    protected:
    private:
};

#endif // AngDist32EqualPBin_HH
