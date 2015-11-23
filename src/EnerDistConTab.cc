#include "../include/EnerDistConTab.hh"

EnerDistConTab::EnerDistConTab(/*int EnerDistStart*/)
{
    /*startEnerDist =  EnerDistStart;*/
}

EnerDistConTab::EnerDistConTab(int numRegsTemp, int *regEndPosTemp, int *intScheme1Temp, int numIncEnerTemp, double *incEnerTemp, int *intScheme2Temp, vector<double> *outSumEnTemp, vector<double> *outSumEnProbTemp)
{
    distPos = NULL;

    numRegs = numRegsTemp;
    regEndPos = new int[numRegsTemp];
    intScheme1 = new int[numRegsTemp];
    for(int i=0; i<numRegsTemp; i++)
    {
        regEndPos[i] = regEndPosTemp[i];
        intScheme1[i] = intScheme1Temp[i];
    }

    numIncEner = numIncEnerTemp;
    incEner = new double[numIncEner];
    intScheme2  = new int[numIncEner];
    numPEnerPoints  = new int[numIncEner];
    numDiscreteEnerPoints  = new int[numIncEner];
    outEner = new double* [numIncEner];
    outProb = new double* [numIncEner];
    outSumProb = new double* [numIncEner];
    for(int i=0; i<numIncEnerTemp; i++)
    {
        incEner[i] = incEnerTemp[i];
        intScheme2[i] = intScheme2Temp[i];
        numPEnerPoints[i] = outSumEnTemp[i].size();
        numDiscreteEnerPoints[i] = 0;
        outEner[i] = new double [numPEnerPoints[i]];
        outProb[i] = new double [numPEnerPoints[i]];
        outSumProb[i] = new double [numPEnerPoints[i]];
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            outEner[i][j]=outSumEnTemp[i][j];
            outProb[i][j]=outSumEnProbTemp[i][j];
            if(j>0)
                outSumProb[i][j] = outSumProb[i][j-1]+outProb[i][j];
            else
                outSumProb[i][j] = outProb[i][j];
        }
    }
}

EnerDistConTab::~EnerDistConTab()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme1)
        delete [] intScheme1;
    if(intScheme2)
        delete [] intScheme2;
    if(distPos)
        delete [] distPos;
    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(incEner)
        delete [] incEner;
    if(numDiscreteEnerPoints)
        delete [] numDiscreteEnerPoints;

    for(int i=0; i<numIncEner; i++)
    {
        if(outEner[i])
            delete [] outEner[i];
        if(outProb[i])
            delete [] outProb[i];
        if(outSumProb[i])
            delete [] outSumProb[i];
    }

    if(outEner)
        delete [] outEner;
    if(outProb)
        delete [] outProb;
    if(outSumProb)
        delete [] outSumProb;
}

void EnerDistConTab::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;

    stream >> numRegs; count++;
    if(numRegs==0)
    {
        numRegs=1;
        regEndPos = new int[numRegs];
        intScheme1 = new int[numRegs];

        stream >> numIncEner; count++;
        regEndPos[0]=numIncEner;
        intScheme1[0]=2;
    }
    else
    {
        regEndPos = new int[numRegs];
        intScheme1 = new int[numRegs];

        for(int i=0; i<numRegs; i++, count++)
        {
            stream >> intTemp;
            regEndPos[i]=intTemp;
        }

        for(int i=0; i<numRegs; i++, count++)
        {
            stream >> intTemp;
            intScheme1[i]=intTemp;
        }
        stream >> numIncEner; count++;
    }

    incEner = new double[numIncEner];
    distPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPEnerPoints = new int[numIncEner];
    numDiscreteEnerPoints = new int[numIncEner];
    outEner = new double *[numIncEner];
    outProb = new double *[numIncEner];
    outSumProb = new double *[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        distPos[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++)
    {
        //This part is not needed and potentially erroneous
        /*
        for(;count<(startEnerDist+distPos[i]-1); count++)
        {
            stream >> dummy;
        }
        */

        stream >> intTemp; count++;
        if(intTemp>10)
        {
            numDiscreteEnerPoints[i] = int(intTemp/10);
            intScheme2[i]=intTemp-numDiscreteEnerPoints[i]*10;
        }
        else
        {
            numDiscreteEnerPoints[i] = 0;
            intScheme2[i]=intTemp;
        }
        if(intScheme2[i]>5)
            intScheme2[i]=2;

        stream >> intTemp; count++;
        numPEnerPoints[i]=intTemp;
        outEner[i] = new double [numPEnerPoints[i]];
        outProb[i] = new double [numPEnerPoints[i]];
        outSumProb[i] = new double [numPEnerPoints[i]];

        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEner[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outProb[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outSumProb[i][j] = temp;
        }
    }
}

//for Capture data, Fission
void EnerDistConTab::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 4
    //convert this to G4NDL theRepresentationType=1

    // here we set correct the energy prob so that it is integrated over its energy regime
    for(int i=0; i<numIncEner; i++)
    {
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            if(outEner[i][j]==0.)
                    outEner[i][j]=1.0e-12;
            if(j==numPEnerPoints[i]-1)
            {
                outProb[i][j] = outSumProb[i][j];
            }
            else
                outProb[i][j] = outSumProb[i][j+1]-outSumProb[i][j];
        }
    }

    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    double sum, enerRange;
    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;

        if(numPEnerPoints[i]!=1)
        {
            stream << std::setw(14) << std::right << numPEnerPoints[i];
            // assume linear interpolation
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

            sum=0.;
            enerRange=0.;
            for(int j=0; j<numPEnerPoints[i]; j++)
            {
                stream << std::setw(14) << std::right << outEner[i][j]*1000000;
                if(j>0)
                {
                    enerRange += outEner[i][j]-outEner[i][j-1];
                }
                stream << std::setw(14) << std::right << outProb[i][j];
                sum+=outProb[i][j];
                if(((j+1)%3==0)||(j==numPEnerPoints[i]-1))
                    stream << '\n';
            }
            if(sum == 0.)
            {
                cout << "Error in EnerDistConTab::WriteG4NDLData" << '\n';
            }
            if(enerRange==0.)
            {
                cout << "Error in EnerDistConTab::WriteG4NDLData" << '\n';
            }
        }
        else
        {
            stream << std::setw(14) << std::right << 2;
            // assume linear interpolation
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << intScheme2[i] << '\n';

            for(int j=0; j<2; j++)
            {
                if(outEner[i][0]==0.)
                    outEner[i][0]=1.0e-12;
                stream << std::setw(14) << std::right << (0.999+double(j)*0.002)*outEner[i][0]*1000000;
                stream << std::setw(14) << std::right << 0.5;
                if(((j+1)%3==0)||(j==1))
                    stream << '\n';
            }
        }
    }

}

double EnerDistConTab::GetAverageOutEnergy()
{
    int i;
    double probSum=0, avgEner1=0, avgEner2=0;
    for(i=0; i<numIncEner-1; i++)
    {
        // assume average incoming neutron energy is 1eV
        if(incEner[i]>1.0e-06)
            break;
    }
    if(i>0)
        i--;

    for(int j=0; j<numPEnerPoints[i]-1; j++)
    {
        probSum += (outProb[i][j+1]-outProb[i][j])/2;
        avgEner1 += (outEner[i][j+1]-outEner[i][j])*(outProb[i][j+1]-outProb[i][j])/2;
    }
    if(probSum==0)
        return incEner[i]*1000000;
    avgEner1/=probSum;

    if(numIncEner>1)
    {
        probSum=0;

        for(int j=0; j<numPEnerPoints[i+1]-1; j++)
        {
            probSum += (outProb[i+1][j+1]-outProb[i+1][j])/2;
            avgEner2 += (outEner[i+1][j+1]-outEner[i+1][j])*(outProb[i+1][j+1]-outProb[i+1][j])/2;
        }
        if(probSum==0)
            return incEner[i]*1000000;
        avgEner2/=probSum;

        if((incEner[i+1]-incEner[i])==0)
            return incEner[i]*1000000;
        return ((1.0-incEner[i])*(avgEner2-avgEner1)/(incEner[i+1]-incEner[i])+avgEner1)*1000000;
    }
    else
    {
        return avgEner1*1000000;
    }
}
