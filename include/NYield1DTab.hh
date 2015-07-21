#ifndef NYield1DTab_HH
#define NYield1DTab_HH

#include "YieldDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the tabular out-going neutron yield distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class NYield1DTab: public YieldDist
{
    public:
        NYield1DTab();
        NYield1DTab(YieldDist *totalYield);
        NYield1DTab(int TYR, double ener1, double ener2, CSDist** nCSDist, int index, int numProc, int numProc2)
        {
            numRegs=1;
            numIncEner=2;
            regEndPos = new int [numRegs];
            intScheme = new int [numRegs];
            regEndPos[0]=2;
            intScheme[0]=2;
            incEner = new double [numIncEner];
            yield = new double [numIncEner];
            incEner[0]=ener1;
            incEner[1]=ener2;

            double sumCS=0.;
            for(int j=numProc2; j<numProc; j++)
            {
                sumCS+=max(0.,nCSDist[j]->Interpolate(ener1));
            }
            yield[0]=TYR*nCSDist[index]->Interpolate(ener1)/sumCS;
            sumCS=0.;
            for(int j=numProc2; j<numProc; j++)
            {
                sumCS+=max(0.,nCSDist[j]->Interpolate(ener2));
            }
            yield[1]=TYR*nCSDist[index]->Interpolate(ener2)/sumCS;
        }
        NYield1DTab(YieldDist* nYield, CSDist** nCSDist, int index, int numProc, int numProc2)
        {
            SetYieldData(numRegs, numIncEner, regEndPos, intScheme, incEner, yield, nCSDist, index, numProc, numProc2);
        }

        virtual ~NYield1DTab();
        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLData(stringstream &stream);
        void SubtractPrompt(YieldDist* &promptYieldDist);
        void SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield);
        void SubtractPrompt(double *&totalCoeff, int &totalNumCoeff)
        {

        }
        void ConvertToLinDist(int *regEndPos, int &numIncEner, double *&incEner, double *&yield)
        {

        }
        void AddData(YieldDist* nYield, CSDist** nCSDist, int index, int numProc, int numProc2)
        {
            if(nYield->IdentifyYourSelf()!="NYield1DTab")
            {
                YieldDist* temp = new NYield1DTab(nYield);
                delete nYield;
                nYield=temp;
            }

            nYield->AddData(numRegs, numIncEner, regEndPos, intScheme, incEner, yield, nCSDist, index, numProc, numProc2);

        }
        void AddData(int &numRegsSum, int &numIncEnerSum, int* &regEndPosSum, int* &intSchemeSum, double* &incEnerSum, double* &yieldSum, CSDist** nCSDist, int index, int numProc, int numProc2)
        {
            int reg, reg1, low;
            vector<int> regVec;
            vector<double> aboveEn, aboveY, belowEn, belowY, cenEn, cenY;
            double sumCS;
            reg1=0;
            for(int i=0; i<numIncEner; i++)
            {
                while(regEndPos[reg1]<=i)
                    reg1++;
                if(incEner[i]>incEnerSum[numIncEnerSum-1])
                {
                    sumCS=0.;
                    for(int j=numProc2; j<numProc; j++)
                    {
                        sumCS+=max(0.,nCSDist[j]->Interpolate(incEner[i]));
                    }
                    aboveEn.push_back(incEner[i]);
                    aboveY.push_back(nCSDist[index]->Interpolate(incEner[i])*yield[i]/sumCS+max(0.,Interpolate( intScheme[numRegsSum-1], incEner[i], incEnerSum[numIncEnerSum-2], incEnerSum[numIncEnerSum-1], yieldSum[numIncEnerSum-2], yieldSum[numIncEnerSum-1])));
                }
                else if(incEner[i]<incEnerSum[0])
                {
                    sumCS=0.;
                    for(int j=numProc2; j<numProc; j++)
                    {
                        sumCS+=max(0.,nCSDist[j]->Interpolate(incEner[i]));
                    }
                    belowEn.push_back(incEner[i]);
                    belowY.push_back(nCSDist[index]->Interpolate(incEner[i])*yield[i]/sumCS+max(0.,Interpolate( intScheme[0], incEner[i], incEnerSum[0], incEnerSum[1], yieldSum[0], yieldSum[1])));
                }
                else
                {
                    reg=0;
                    for(int j=0; j<numIncEnerSum; j++)
                    {
                        while((regEndPosSum[reg]<=j)&&(numRegs-1>reg))
                            reg++;
                        if(incEnerSum[j]==incEner[i])
                        {
                            break;
                        }
                        else if(incEnerSum[j]>incEner[i])
                        {
                            sumCS=0.;
                            for(int j=numProc2; j<numProc; j++)
                            {
                                sumCS+=max(0.,nCSDist[j]->Interpolate(incEner[i]));
                            }
                            cenEn.push_back(incEner[i]);
                            cenY.push_back(nCSDist[index]->Interpolate(incEner[i])*yield[i]/sumCS+max(0.,Interpolate( intScheme[reg], incEner[i], incEnerSum[j-1], incEnerSum[j], yieldSum[j-1], yieldSum[j])));
                            regVec.push_back(reg1);
                        }
                    }
                }
            }


            for(int i=0; i<numIncEnerSum; i++)
            {
                reg=0;
                for(low=0; low<numIncEner-1; low++)
                {
                    if(incEner[low]>incEnerSum[i])
                    {
                        break;
                    }
                    while((regEndPos[reg]<=low)&&(numRegs-1>reg))
                        reg++;
                }
                if(low!=0)
                    low--;

                sumCS=0.;
                for(int j=numProc2; j<numProc; j++)
                {
                    sumCS+=max(0.,nCSDist[j]->Interpolate(incEnerSum[i]));
                }

                yieldSum[i]+=max(0.,nCSDist[index]->Interpolate(incEnerSum[i])*Interpolate( intScheme[reg], incEnerSum[i], incEner[low], incEner[low+1], yield[low], yield[low+1])/sumCS);
            }

            numIncEnerSum+=belowEn.size()+aboveEn.size()+cenEn.size();
            regEndPosSum[0]+=belowEn.size();
            regEndPosSum[numRegsSum-1]+=+belowEn.size()+aboveEn.size()+cenEn.size();
            double *tempEner = new double [numIncEnerSum];
            double *tempYield = new double [numIncEnerSum];
            int count=0;

            for(int i=0; i<numIncEnerSum; i++)
            {
                if(i<int(belowEn.size()))
                {
                    tempEner[i]=belowEn[i];
                    tempYield[i]=belowY[i];
                }
                else if(i>=(numIncEnerSum-int(aboveEn.size())))
                {
                    tempEner[i]=aboveEn[i-numIncEnerSum+aboveEn.size()];
                    tempYield[i]=aboveY[i-numIncEnerSum+aboveEn.size()];
                }
                else
                {
                    while((count!=int(cenEn.size()))&&(tempEner[i]>cenEn[count]))
                    {
                        for(int j=regVec[count]; j<numRegsSum; j++)
                        {
                            regEndPosSum[j]++;
                        }
                        tempEner[i]=cenEn[count];
                        tempYield[i]=cenY[count];
                        count++;
                    }
                    tempEner[i]=incEnerSum[i-belowEn.size()];
                    tempYield[i]=yieldSum[i-belowEn.size()];
                }

            }
            delete [] incEnerSum;
            delete [] yieldSum;

            incEnerSum=tempEner;
            yieldSum=tempYield;

        }

        void SetYieldData(int &numRegsSum, int &numIncEnerSum, int* &regEndPosSum, int* &intSchemeSum, double* &incEnerSum, double* &yieldSum, CSDist** nCSDist, int index, int numProc, int numProc2)
        {
            numRegsSum=numRegs;
            numIncEnerSum=numIncEner;
            regEndPosSum=new int [numRegsSum];
            intSchemeSum=new int [numRegsSum];
            for(int i=0; i<numRegs; i++)
            {
                regEndPosSum[i]=regEndPos[i];
                intSchemeSum[i]=intScheme[i];
            }

            double sumCS;
            incEnerSum = new double[numIncEner];
            yieldSum = new double[numIncEner];
            for(int i=0; i<numIncEner; i++)
            {
                sumCS=0.;
                for(int j=numProc2; j<numProc; j++)
                {
                    sumCS+=max(0.,nCSDist[j]->Interpolate(incEner[i]));
                }
                incEnerSum[i]=incEner[i];
                yieldSum[i]=max(0.,nCSDist[index]->Interpolate(incEner[i])*yield[i]/sumCS);
            }
        }

        string IdentifyYourSelf()
        {
            return "NYield1DTab";
        }

        int numRegs, numIncEner;
        int *regEndPos, *intScheme;
        double *incEner; // contains incoming neutron energy
        double *yield;
    protected:
    private:
};

#endif // NYield1DTab_HH
