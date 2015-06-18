#include "../include/AngEnDistLab3DTab.hh"

AngEnDistLab3DTab::AngEnDistLab3DTab(/*int EnerDistStart*/)
{
    //startEnerDist =  EnerDistStart;
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
            if(outEnProb[i][j])
                delete [] outEnProb[i][j];
            if(outEnSumProb[i][j])
                delete [] outEnSumProb[i][j];
        }
        if(outEner[i])
            delete [] outEner[i];
        if(outEnProb[i])
            delete [] outEnProb[i];
        if(outEnSumProb[i])
            delete [] outEnSumProb[i];
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
    if(outEnProb)
        delete [] outEnProb;
    if(outEnSumProb)
        delete [] outEnSumProb;
}

void AngEnDistLab3DTab::ExtractMCNPData(stringstream &stream, int &count)
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
    outEnProb = new double **[numIncEner];
    outEnSumProb = new double **[numIncEner];

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
        //not needed, and potentially eroneous
        /*
        for(;count<(startEnerDist+outAngDistPos[i]-1); count++)
        {
            stream >> dummy;
        }
        */

        stream >> intTemp; count++;
        intScheme2[i]=intTemp;

        stream >> intTemp; count++;
        numPAngPoints[i]=intTemp;
        outAng[i] = new double [numPAngPoints[i]];
        outEnerDistPos[i] = new int [numPAngPoints[i]];
        intScheme3[i] = new int [numPAngPoints[i]];
        numPEnerPoints[i] = new int [numPAngPoints[i]];
        outEner[i] = new double *[numPAngPoints[i]];
        outEnProb[i] = new double *[numPAngPoints[i]];
        outEnSumProb[i] = new double *[numPAngPoints[i]];

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
            //not needed, and potentially erroneous
            /*
            for(;count<(startEnerDist+outEnerDistPos[i]-1); count++)
            {
                stream >> dummy;
            }
            */
            stream >> intTemp; count++;
            intScheme3[i][j] = intTemp;

            stream >> intTemp; count++;
            numPEnerPoints[i][j] = intTemp;

            outEner[i][j] = new double [numPEnerPoints[i][j]];
            outEnProb[i][j] = new double [numPEnerPoints[i][j]];
            outEnSumProb[i][j] = new double [numPEnerPoints[i][j]];

            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEner[i][j][k] = temp;
            }
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEnProb[i][j][k] = temp;
            }
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEnSumProb[i][j][k] = temp;
            }
        }
    }
}

void AngEnDistLab3DTab::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 67
    //convert this to G4NDL DistLaw=7

    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << numPAngPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numPAngPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

        for(int j=0; j<numPAngPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outAng[i][j];
            stream << std::setw(14) << std::right << numPEnerPoints[i][j];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i][j] << std::setw(14) << std::right << intScheme3[i][j] << '\n';

            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream << std::setw(14) << std::right << outEner[i][j][k]*1000000;
                stream << std::setw(14) << std::right << outEnProb[i][j][k];
                if(k%3==0)
                    stream << '\n';
            }
        }
    }
}
