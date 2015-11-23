#include "../include/EnerDistTabLinFunc.hh"

EnerDistTabLinFunc::EnerDistTabLinFunc(/*int EnerDistStart*/)
{
    /*startEnerDist =  EnerDistStart;*/
}

EnerDistTabLinFunc::~EnerDistTabLinFunc()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme)
        delete [] intScheme;
    if(distPos)
        delete [] distPos;
    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(incEner)
        delete [] incEner;

    for(int i=0; i<numIncEner; i++)
    {
        if(outEnerOffset[i])
            delete [] outEnerOffset[i];
        if(outEnerMulti[i])
            delete [] outEnerMulti[i];
        if(outProb[i])
            delete [] outProb[i];
    }

    if(outEnerOffset)
        delete [] outEnerOffset;
    if(outEnerMulti)
        delete [] outEnerMulti;
    if(outProb)
        delete [] outProb;
}

void EnerDistTabLinFunc::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;

    stream >> numRegs; count++;
    if(numRegs==0)
    {
        numRegs=1;
        regEndPos = new int[numRegs];
        intScheme = new int[numRegs];

        stream >> numIncEner; count++;
        regEndPos[0]=numIncEner;
        intScheme[0]=2;
    }
    else
    {
        regEndPos = new int[numRegs];
        intScheme = new int[numRegs];

        for(int i=0; i<numRegs; i++, count++)
        {
            stream >> intTemp;
            regEndPos[i]=intTemp;
        }

        for(int i=0; i<numRegs; i++, count++)
        {
            stream >> intTemp;
            intScheme[i]=intTemp;
        }
        stream >> numIncEner; count++;
    }

    incEner = new double[numIncEner];
    distPos = new int[numIncEner];
    numPEnerPoints = new int[numIncEner];
    outProb = new double *[numIncEner];
    outEnerOffset = new double *[numIncEner];
    outEnerMulti = new double *[numIncEner];

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
        //Not needed, and potentially erroneous
        /*
        for(;count<(startEnerDist+distPos[i]-1); count++)
        {
            stream >> dummy;
        }
        */

        stream >> intTemp; count++;
        numPEnerPoints[i]=intTemp;
        outProb[i] = new double [numPEnerPoints[i]];
        outEnerOffset[i] = new double [numPEnerPoints[i]];
        outEnerMulti[i] = new double [numPEnerPoints[i]];

        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            if(j!=0)
                outProb[i][j] = temp-outProb[i][j-1];
            else
                outProb[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerOffset[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerMulti[i][j] = temp;
        }
    }
}

//For Fission
void EnerDistTabLinFunc::WriteG4NDLData(stringstream &stream)
{
//this is MCNP law 22
//there is no direct translation for this law in G4NDL but it can be made to fit theRepresentationType=1 with some assumptions
//use the average outgoing energy from each linear function and assign it the probability for the function and interpolate between them
//or to be more precise, create a fine grid of out-going energies, and a corresponding probability vector, then input the incoming energies boundaries
//for which the linear functions are valid, into each linear function, and add the probability of the lines to each of the out going probability elements whose
//coresponding out-going energy falls with the range of the linear function.
//then normalize the created out-going energy probability distribution

    double *outEner, *outProbDist, *outEnerLow, *outEnerHigh, low, high, emax=-1, emin=-1, tempLow, tempHigh, sum=0., enerRange;

    //we ignore the given interpolation scheme since MCNP ignores it and uses a histogram scheme instead
    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << numIncEner;
    stream << std::setw(14) << std::right << 1 << '\n';

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << 2*numPEnerPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2*numPEnerPoints[i] << std::setw(14) << std::right << 2 << '\n';

        outEner=new double [2*numPEnerPoints[i]];
        outProbDist=new double [2*numPEnerPoints[i]];
        outEnerLow=new double [2*numPEnerPoints[i]];
        outEnerHigh=new double [2*numPEnerPoints[i]];

        if(i==numIncEner-1)
        {
            low=incEner[i]-0.5*(incEner[i]-incEner[i-1]);
            high=incEner[i]+0.5*(incEner[i]-incEner[i-1]);
        }
        else
        {
            low=incEner[i]-0.5*(incEner[i+1]-incEner[i]);
            high=incEner[i]+0.5*(incEner[i+1]-incEner[i]);
        }

        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            tempLow=outEnerMulti[i][j]*(low-outEnerOffset[i][j]);
            tempHigh=outEnerMulti[i][j]*(high-outEnerOffset[i][j]);
            if((emin>tempLow)||(emin==-1))
            {
                emin=tempLow;
            }
            if((emin>tempHigh))
            {
                emin=tempHigh;
            }
            if(emax<tempHigh)
            {
                emax=tempHigh;
            }
            if(emax<tempLow)
            {
                emax=tempLow;
            }
        }

        for(int j=0; j<2*numPEnerPoints[i]; j++)
        {
            outEner[j]=(emax-emin)*j/(2*numPEnerPoints[i])+(emax-emin)/(4*numPEnerPoints[i])+emin;
            outEnerLow[j]=(emax-emin)*j/(2*numPEnerPoints[i])+emin;
            outEnerHigh[j]=(emax-emin)*(j+1)/(2*numPEnerPoints[i])+emin;
            outProbDist[j]=0;
        }

        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            tempLow=outEnerMulti[i][j]*(low-outEnerOffset[i][j]);
            tempHigh=outEnerMulti[i][j]*(high-outEnerOffset[i][j]);
            for(int k=0; k<2*numPEnerPoints[i]; k++)
            {
                if((outEnerLow[k]<tempLow)&&(outEnerHigh[k]>tempLow))
                {
                    outProbDist[k]+=outProb[i][j];
                    sum+=outProb[i][j];
                }
                if((outEnerLow[k]<tempHigh)&&(outEnerHigh[k]>tempHigh))
                {
                    outProbDist[k]+=outProb[i][j];
                    sum+=outProb[i][j];
                }
                if((tempLow<outEnerLow[k])&&(tempHigh>outEnerLow[k]))
                {
                    outProbDist[k]+=outProb[i][j];
                    sum+=outProb[i][j];
                }
                if((tempLow<outEnerHigh[k])&&(tempHigh>outEnerHigh[k]))
                {
                    outProbDist[k]+=outProb[i][j];
                    sum+=outProb[i][j];
                }
            }
        }

        for(int k=0; k<2*numPEnerPoints[i]; k++)
        {
            outProbDist[k]/=sum;
        }
        sum=0.;
        enerRange=0.;

        for(int j=0; j<2*numPEnerPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outEner[j]*1000000;
            if(j>0)
            {
                enerRange += outEner[j]-outEner[j-1];
            }
            stream << std::setw(14) << std::right << outProbDist[j];
            sum += outProbDist[j];
            if(((j+1)%3==0)||(j==2*numPEnerPoints[i]-1))
                stream << '\n';
        }
        if(sum==0.)
        {
            cout << "Error with the energy probability data" << endl;
        }
        if(enerRange==0.)
        {
            cout << "Error in EnerDistTabLinFunc::WriteG4NDLData" << '\n';
        }
        sum=0.;

        delete[] outEner;
        delete[] outProbDist;
        delete[] outEnerLow;
        delete[] outEnerHigh;
    }

}

