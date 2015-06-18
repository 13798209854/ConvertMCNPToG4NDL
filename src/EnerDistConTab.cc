#include "../include/EnerDistConTab.hh"

EnerDistConTab::EnerDistConTab(/*int EnerDistStart*/)
{
    /*startEnerDist =  EnerDistStart;*/
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
    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << numPEnerPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outEner[i][j]*1000000;
            stream << std::setw(14) << std::right << outProb[i][j];
            if(j%3==0)
                stream << '\n';
        }
    }

}

double EnerDistConTab::GetAverageOutEnergy()
{
    int i;
    double probSum=0, avgEner1=0, avgEner2=0;
    for(i=0; i<numIncEner; i++)
    {
        // assume average incoming neutron energy is 1eV
        if(incEner[i]>1.0)
            break;
    }
    if(i==0)
        i++;

    for(int j=0; j<numPEnerPoints[i-1]-1; j++)
    {
        probSum += (outProb[i-1][j+1]-outProb[i-1][j])/2;
        avgEner1 += (outEner[i-1][j+1]-outEner[i-1][j])*(outProb[i-1][j+1]-outProb[i-1][j])/2;
    }
    avgEner1/=probSum;
    probSum=0;

    for(int j=0; j<numPEnerPoints[i]-1; j++)
    {
        probSum += (outProb[i][j+1]-outProb[i][j])/2;
        avgEner2 += (outEner[i][j+1]-outEner[i][j])*(outProb[i][j+1]-outProb[i][j])/2;
    }
    avgEner2/=probSum;

    return ((1.0-incEner[i-1])*(avgEner2-avgEner1)/(incEner[i]-incEner[i-1])+avgEner1)*1000000;

}
