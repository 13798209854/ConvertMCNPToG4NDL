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
        if(outEnerProb[i])
            delete [] outEnerProb[i];
        if(outEnerSumProb[i])
            delete [] outEnerSumProb[i];
        if(rFraction[i])
            delete [] rFraction[i];
        if(angDistSlope[i])
            delete [] angDistSlope[i];
    }

    if(outEner)
        delete [] outEner;
    if(outEnerProb)
        delete [] outEnerProb;
    if(outEnerSumProb)
        delete [] outEnerSumProb;
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

    for(int i=0; i<numDistSample; i++)
    {
        outAng[i] = -1+2*i/numDistSample;
    }

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
    outEnerProb = new double *[numIncEner];
    outEnerSumProb = new double *[numIncEner];
    rFraction = new double *[numIncEner];
    angDistSlope = new double *[numIncEner];
    outAngProb = new double **[numIncEner];

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
        //this is not needed, and potentially erroneous
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
        outEnerProb[i] = new double [numPEnerPoints[i]];
        outEnerSumProb[i] = new double [numPEnerPoints[i]];
        rFraction[i] = new double [numPEnerPoints[i]];
        angDistSlope[i] = new double [numPEnerPoints[i]];
        outAngProb[i] = new double *[numIncEner];

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
            rFraction[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            angDistSlope[i][j] = temp;
        }
        // create the angle probability distribution from the given function for this out-going energy
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            outAngProb[i][j] = new double [numDistSample];
            for(int k=0; k<numDistSample; k++)
            {
                outAngProb[i][j][k] = 0.5*angDistSlope[i][j]*(std::cosh(angDistSlope[i][j]*outAng[k])+rFraction[i][j]*std::sinh(angDistSlope[i][j]*outAng[k]))/(std::sinh(angDistSlope[i][j]));
            }
        }
    }

}

void AngEnDistKallbach::WriteG4NDLData(stringstream stream)
{
    //this is MCNP Law 44
    //convert this to G4NDL DistLaw=7
    // may be able to exactly convert it to DistLaw=1

     stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n'

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << numDistSample;
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numDistSample << std::setw(14) << std::right << 2 << '\n';

        for(int j=0; j<numDistSample; j++)
        {
            stream << std::setw(14) << std::right << outAng[j];
            stream << std::setw(14) << std::right << numPEnerPoints[i];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';
            for(int k=0; k<numPEnerPoints[i]; k++)
            {
                stream << std::setw(14) << std::right << outEner[i][k]*1000000;
                stream << std::setw(14) << std::right << outEnerProb[i][k]*outAngProb[i][k][j] << '\n';
            }
        }
    }
}
