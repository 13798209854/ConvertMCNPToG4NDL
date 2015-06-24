#include "../include/AngEnDist3DTab.hh"

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

void AngEnDist3DTab::ExtractMCNPData(stringstream &stream, int &count)
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
    outEnDistPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPEnerPoints = new int[numIncEner];
    numDiscreteEnerPoints = new int[numIncEner];
    outEner = new double *[numIncEner];
    outEnerProb = new double *[numIncEner];
    outEnerSumProb = new double *[numIncEner];
    isoDist = new bool *[numIncEner];
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
        //Not needed and dangerous since we don't know exactly where outEnDistPos[i] points
        /*
        for(;count<(startEnerDist+outEnDistPos[i]-1); count++)
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
        isoDist[i] = new bool [numPEnerPoints[i]];
        outAngDistPos[i] = new int [numPEnerPoints[i]];
        intScheme3[i] = new int [numPEnerPoints[i]];
        numPAngPoints[i] = new int [numPEnerPoints[i]];
        outAng[i] = new double *[numPEnerPoints[i]];
        outAngProb[i] = new double *[numPEnerPoints[i]];
        outAngSumProb[i] = new double *[numPEnerPoints[i]];

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
                for(;count<(startEnerDist+outAngDistPos[i][j]-1); count++)
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
            if(intTemp==0)
                intScheme3[i][j] = 1;
            else
                intScheme3[i][j] = 2;

            stream >> intTemp; count++;
            numPAngPoints[i][j] = intTemp;

            outAng[i][j] = new double [numPAngPoints[i][j]];
            outAngProb[i][j] = new double [numPAngPoints[i][j]];
            outAngSumProb[i][j] = new double [numPAngPoints[i][j]];

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

void AngEnDist3DTab::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 61
    //convert this to G4NDL DistLaw=7
    //some approximations are used in the transformation

    // this creates a new out-going angle array, for each incoming energy. each of the angles in the new angle arrays
    // has a new probability dist which is dependant on the out-going energy. We did this to bring the 3D table format
    // (order of properties) inline with that used by G4NDL DistLaw=7
    int *sumAngPoints = new int [numIncEner], *newNumAngPoints = new int [numIncEner];
    double *angMax= new double [numIncEner], *angMin = new double [numIncEner];
    double **outAngNew = new double* [numIncEner];
    double ***outEnProbNew = new double** [numIncEner];

    for(int i=0; i<numIncEner; i++)
    {
        sumAngPoints[i]=0;
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            sumAngPoints[i]+=numPAngPoints[i][j];
            angMax[i]=0;
            angMin[i]=0;
            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                if(outAng[i][j][k]>angMax[i])
                    angMax[i]=outAng[i][j][k];
                if((outAng[i][j][k]<angMin[i])||(angMin[i]==0))
                    angMin[i]=outAng[i][j][k];
            }
        }
    }

    for(int i=0; i<numIncEner; i++)
    {
        newNumAngPoints[i] = floor(5*sumAngPoints[i]/numPEnerPoints[i]);
        outAngNew[i] = new double [newNumAngPoints[i]];
        outEnProbNew[i] = new double* [newNumAngPoints[i]];

        for(int j=0; j<newNumAngPoints[i]; j++)
        {
            outAngNew[i][j] = (angMax[i]-angMin[i])*j/newNumAngPoints[i]+angMin[i]+(angMax[i]-angMin[i])*0.5/newNumAngPoints[i];
            outEnProbNew[i][j] = new double [numPEnerPoints[i]];

            for(int k=0; k<numPEnerPoints[i]; k++)
            {
                int l;
                for(l=0; l<numPAngPoints[i][k]; l++)
                {
                    if(outAng[i][k][l]>outAngNew[i][j])
                    {
                        l--;
                        break;
                    }
                }
                if(l<0)
                    l=0;
                outEnProbNew[i][j][k]=outEnerProb[i][k]*Interpolate(intScheme3[i][k], outAngNew[i][j], outAng[i][k][l], outAng[i][k][l+1], outAngProb[i][k][l], outAngProb[i][k][l+1]);
            }
        }
    }

    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << newNumAngPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << newNumAngPoints[i] << std::setw(14) << std::right << 2 << '\n';

        for(int j=0; j<newNumAngPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outAngNew[i][j];
            stream << std::setw(14) << std::right << numPEnerPoints[i];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

            for(int k=0; k<numPEnerPoints[i]; k++)
            {
                stream << std::setw(14) << std::right << outEner[i][k]*1000000;
                stream << std::setw(14) << std::right << outEnProbNew[i][j][k];
                if(k%3==0)
                    stream << '\n';
            }
        }
    }

    for(int i=0; i<numIncEner; i++)
    {
        delete [] outAngNew[i];
        for(int j=0; j<newNumAngPoints[i]; j++)
        {
            delete [] outEnProbNew[i][j];
        }
        delete [] outEnProbNew[i];
    }

    delete [] outEnProbNew;
    delete [] outAngNew;
    delete [] angMax;
    delete [] angMin;
    delete [] sumAngPoints;
    delete [] newNumAngPoints;
}
