#include "AngEnDistKallbach.hh"

AngEnDistKallbach::AngEnDistKallbach(int EnerDistStart)
{
    startEnerDist =  EnerDistStart;
}

AngEnDistKallbach::~AngEnDistKallbach()
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
        if(rFraction[i])
            delete [] rFraction[i];
        if(angDistSlope[i])
            delete [] angDistSlope[i];
    }

    if(outEner)
        delete [] outEner;
    if(outProb)
        delete [] outProb;
    if(outSumProb)
        delete [] outSumProb;
    if(rFraction)
        delete [] rFraction;
    if(angDistSlope)
        delete [] angDistSlope;
}

void AngEnDistKallbach::ExtractMCNPData(stringstream stream, int &count)
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
    rFraction = new double *[numIncEner];
    angDistSlope = new double *[numIncEner];

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
        for(;count<(startEnerDist+distPos[i]-1); count++)
        {
            stream >> dummy;
        }

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
        rFraction[i] = new double [numPEnerPoints[i]];
        angDistSlope[i] = new double [numPEnerPoints[i]];

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
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            rFraction[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            angDistSlope[i][j] = temp;
        }
    }
}

void AngEnDistKallbach::WriteG4NDLData(stringstream data)
{


}
