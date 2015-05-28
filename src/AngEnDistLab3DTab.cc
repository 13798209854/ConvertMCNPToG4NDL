#include "AngEnDistLab3DTab.hh"

AngEnDistLab3DTab::AngEnDistLab3DTab(int EnerDistStart)
{
    startEnerDist =  EnerDistStart;
}

AngEnDistLab3DTab::~AngEnDistLab3DTab()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme1)
        delete [] intScheme1;
    if(intScheme2)
        delete [] intScheme2;
    if(outAngDistPos)
        delete [] outAngDistPos;
    if(incEner)
        delete [] incEner;

    for(int i=0; i<numIncEner; i++)
    {
        if(outAng[i])
            delete [] outAng[i];
        if(outEnerDistPos[i])
            delete [] outEnerDistPos[i];
        if(intScheme3[i])
            delete [] intScheme3[i];
        if(numPEnerPoints[i])
            delete [] numPEnerPoints[i];
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            if(outEner[i][j])
                delete [] outEner[i][j];
            if(outProb[i][j])
                delete [] outProb[i][j];
            if(outSumProb[i][j])
                delete [] outSumProb[i][j];
        }
        if(outEner[i])
            delete [] outEner[i];
        if(outProb[i])
            delete [] outProb[i];
        if(outAngSumProb[i])
            delete [] outSumProb[i];
    }
    if(numPAngPoints)
        delete [] numPAngPoints;
    if(outAng)
        delete [] outAng;
    if(outEnerDistPos)
        delete [] outEnerDistPos;
    if(intScheme3)
        delete [] intScheme3;
    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(outEner)
        delete [] outEner;
    if(outProb)
        delete [] outProb;
    if(outSumProb)
        delete [] outSumProb;
}

void AngEnDistLab3DTab::ExtractMCNPData(stringstream stream, int &count)
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
    outAngDistPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPAngPoints = new int[numIncEner];
    outAng = new double *[numIncEner];
    outEnerDistPos = new int *[numIncEner];
    intScheme3 = new int *[numIncEner];
    numPEnerPoints = new int *[numIncEner];
    outEner = new double **[numIncEner];
    outProb = new double **[numIncEner];
    outSumProb = new double **[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        outAngDistPos[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++)
    {
        //Check to make sure the -1 is needed
        for(;count<(startEnerDist+outAngDistPos[i]-1); count++)
        {
            stream >> dummy;
        }

        stream >> intTemp; count++;
        intScheme2[i]=intTemp;

        stream >> intTemp; count++;
        numPAngPoints[i]=intTemp;
        outAng[i] = new double [numPAngPoints[i]];
        outEnerDistPos[i] = new int [numPAngPoints[i]];
        intScheme3[i] = new int [numPAngPoints[i]];
        numPEnerPoints[i] = new int [numPAngPoints[i]];
        outEner[i] = new double *[numPAngPoints[i]];
        outProb[i] = new double *[numPAngPoints[i]];
        outSumProb[i] = new double *[numPAngPoints[i]];

        for(int j=0; j<numPAngPoints[i]; j++, count++)
        {
            stream >> temp;
            outAng[i][j] = temp;
        }
        for(int j=0; j<numPAngPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerDistPos[i][j] = temp;
        }
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            //Check to make sure the -1 is needed
            for(;count<(startEnerDist+outEnerDistPos[i]-1); count++)
            {
                stream >> dummy;
            }
            stream >> intTemp; count++;
            intScheme3[i][j] = intTemp;

            stream >> intTemp; count++;
            numPEnerPoints[i][j] = intTemp;

            outEner = new double [numPEnerPoints[i][j]];
            outAngProb = new double [numPEnerPoints[i][j]];
            outAngSumProb = new double [numPEnerPoints[i][j]];

            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEner[i][j][k] = temp;
            }
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outProb[i][j][k] = temp;
            }
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outSumProb[i][j][k] = temp;
            }
        }
    }
}

void AngEnDistLab3DTab::WriteG4NDLData(stringstream data)
{


}
