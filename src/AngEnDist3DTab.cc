#include "AngEnDist3DTab.hh"

AngEnDist3DTab::AngEnDist3DTab(int EnerDistStart)
{
    startEnerDist =  EnerDistStart;
}

AngEnDist3DTab::~AngEnDist3DTab()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme1)
        delete [] intScheme1;
    if(intScheme2)
        delete [] intScheme2;
    if(outEnDistPos)
        delete [] outEnDistPos;
    if(incEner)
        delete [] incEner;
    if(numDiscreteEnerPoints)
        delete [] numDiscreteEnerPoints;

    for(int i=0; i<numIncEner; i++)
    {
        if(outEner[i])
            delete [] outEner[i];
        if(outEnerProb[i])
            delete [] outEnerProb[i];
        if(outEnerSumProb[i])
            delete [] outEnerSumProb[i];
        if(isoDist[i])
            delete [] isoDist[i];
        if(outAngDistPos[i])
            delete [] outAngDistPos[i];
        if(intScheme3[i])
            delete [] intScheme3[i];
        if(numPAngPoints[i])
            delete [] numPAngPoints[i];
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            if(outAng[i][j])
                delete [] outAng[i][j];
            if(outAngProb[i][j])
                delete [] outAngProb[i][j];
            if(outAngSumProb[i][j])
                delete [] outAngSumProb[i][j];
        }
        if(outAng[i])
            delete [] outAng[i];
        if(outAngProb[i])
            delete [] outAngProb[i];
        if(outAngSumProb[i])
            delete [] outAngSumProb[i];
    }
    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(outEner)
        delete [] outEner;
    if(outEnerProb)
        delete [] outEnerProb;
    if(outEnerSumProb)
        delete [] outEnerSumProb;
    if(isoDist)
        delete [] isoDist;
    if(outAngDistPos)
        delete [] outAngDistPos;
    if(intScheme3)
        delete [] intScheme3;
    if(numPAngPoints)
        delete [] numPAngPoints;
    if(outAng)
        delete [] outAng;
    if(outAngProb)
        delete [] outAngProb;
    if(outAngSumProb)
        delete [] outAngSumProb;
}

void AngEnDist3DTab::ExtractMCNPData(stringstream stream, int &count)
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
    outEnDistPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPEnerPoints = new int[numIncEner];
    numDiscreteEnerPoints = new int[numIncEner];
    outEner = new double *[numIncEner];
    outEnerProb = new double *[numIncEner];
    outEnerSumProb = new double *[numIncEner];
    isoDist = new double *[numIncEner];
    outAngDistPos = new int *[numIncEner];
    intScheme3 = new int *[numIncEner];
    numPAngPoints = new int *[numIncEner];
    outAng = new double **[numIncEner];
    outAngProb = new double **[numIncEner];
    outAngSumProb = new double **[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        outEnDistPos[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++)
    {
        //Check to make sure the -1 is needed
        for(;count<(startEnerDist+outEnDistPos[i]-1); count++)
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
        outEnerProb[i] = new double [numPEnerPoints[i]];
        outEnerSumProb[i] = new double [numPEnerPoints[i]];
        isoDist[i] = new bool [numPEnerPoints[i]];
        outAngDistPos = new int [numPEnerPoints[i]];
        intScheme3 = new int [numPEnerPoints[i]];
        numPAngPoints = new int [numPEnerPoints[i]];
        outAng = new double *[numPEnerPoints[i]];
        outAngProb = new double *[numPEnerPoints[i]];
        outAngSumProb = new double *[numPEnerPoints[i]];

        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEner[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerProb[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerSumProb[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outAngDistPos[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            if(outAngDistPos[i][j]>0)
            {
                for(;count<(startEnerDist+outEnDistPos[i]); count++)
                {
                    stream >> dummy;
                }
            }
            else
            {
                isoDist[i][j]=true;
                continue;
            }
            stream >> intTemp; count++;
            intScheme3[i][j] = intTemp;

            stream >> intTemp; count++;
            numPAngPoints[i][j] = intTemp;

            outAng = new double [numPAngPoints[i][j]];
            outAngProb = new double [numPAngPoints[i][j]];
            outAngSumProb = new double [numPAngPoints[i][j]];

            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outAng[i][j][k] = temp;
            }
            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outAngProb[i][j][k] = temp;
            }
            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outAngSumProb[i][j][k] = temp;
            }
        }
    }
}

void AngEnDist3DTab::WriteG4NDLData(stringstream data)
{


}
