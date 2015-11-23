#include "../include/CSDistYieldComp.hh"

CSDistYieldComp::CSDistYieldComp()
{
    regEndPos=NULL; intScheme=NULL;
    yieldVec=NULL; enerVec=NULL; enerCSVec=NULL; csVec=NULL;
}

CSDistYieldComp::~CSDistYieldComp()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme)
        delete [] intScheme;
    if(enerVec)
        delete [] enerVec;
    if(yieldVec)
        delete [] yieldVec;
}

void CSDistYieldComp::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;

    stream >> mtNum; count++;
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

    enerVec = new double[numIncEner];
    yieldVec = new double[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        enerVec[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        yieldVec[i]=temp;
    }
}

//set up for Capture data
void CSDistYieldComp::WriteG4NDLCSData(stringstream &stream )
{
    double csNum;
    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << endl;

    int i,csIndex;
    for(i=0; i<numRegs; i++)
    {
        stream << regEndPos[i] << intScheme[i] << endl;
    }
    stream << '\n';
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int j=0; j<numIncEner; j++)
    {
        stream << std::setw(14) << std::right << enerVec[j]*1000000;
        for(i=csEnerStart-1; i<csSize+csEnerStart-2; i++)
        {
            if(enerCSVec[i]>enerVec[j])
            {
                break;
            }
        }
        i--;

        csIndex=i-csEnerStart+1;
        if(csIndex<0)
        {
            csIndex=0;
            i=csEnerStart-1;
        }

        csNum = (enerVec[i]-enerCSVec[i])*(csVec[csIndex+1]-csVec[csIndex])/(enerCSVec[i+1]-enerCSVec[i])+csVec[csIndex];
        stream << std::setw(14) << std::right << (yieldVec[j]*csNum);
        if((j+1)%3==0)
            stream << '\n';
    }
}

void CSDistYieldComp::WriteG4NDLYieldData(stringstream &stream )
{
    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << endl;

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i] << std::setw(14) << std::right << intScheme[i] << endl;
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << enerVec[i]*1000000;
        stream << std::setw(14) << std::right << yieldVec[i];
        if(((i+1)%3==0)||(i==numIncEner-1))
            stream << '\n';
    }
}

string CSDistYieldComp::IdentifyYourSelf()
{
    return "CSDistYieldComp";
}

double CSDistYieldComp::Interpolate(double x)
{
    int i, low, reg=0, csIndex;
    double csNum1, csNum2;

    for(i=0; i<numIncEner-1; i++)
    {
        while((i>=regEndPos[reg])&&(numRegs-1>reg))
            reg++;
        if(enerVec[i]>x)
        {
            break;
        }
    }
    i--;
    if(i<0)
        i=0;

    for(low=csEnerStart-1; low<csSize+csEnerStart-2; low++)
    {
        if(enerCSVec[low]>enerVec[i])
        {
            break;
        }
    }

    low--;
    csIndex=low-csEnerStart+1;
    if(csIndex<0)
    {
        csIndex=0;
        low=csEnerStart-1;
    }

    if(csSize>1)
        csNum1 = max(0.,(enerVec[i]-enerCSVec[low])*(csVec[csIndex+1]-csVec[csIndex])/(enerCSVec[low+1]-enerCSVec[low])+csVec[csIndex]);
    else
        csNum1 = csVec[0];

    if(numIncEner>1)
    {
        for(low=csEnerStart-1; low<csSize+csEnerStart-2; low++)
        {
            if(enerCSVec[low]>enerVec[i+1])
            {
                break;
            }
        }
        low--;
        csIndex=low-csEnerStart+1;
        if(csIndex<0)
        {
            csIndex=0;
            low=csEnerStart-1;
        }

        if(csSize>1)
            csNum2 = max(0.,(enerVec[i+1]-enerCSVec[low])*(csVec[csIndex+1]-csVec[csIndex])/(enerCSVec[low+1]-enerCSVec[low])+csVec[csIndex]);
        else
            csNum2 = csVec[0];
    }
    else
    {
        csNum2 = csNum1;
    }

    if(numIncEner>1)
        return max(0.,CSDist::Interpolate(intScheme[reg], x, enerVec[i], enerVec[i+1], yieldVec[i]*csNum1, yieldVec[i+1]*csNum2));
    else
        return yieldVec[0]*csNum1;
}

double CSDistYieldComp::GetAvgCS()
{
    double csSum=0., csNum1, csNum2;
    int i = csEnerStart-1, csIndex=0;

    if(!setAvgCS)
    {
        for(int j=0; j<numIncEner-1; j++)
        {
            for(; i<csSize+csEnerStart-2; i++)
            {
                if(enerCSVec[i]>enerVec[j])
                {
                    break;
                }
            }
            i--;

            csIndex=i-csEnerStart+1;
            if(csIndex<0)
            {
                csIndex=0;
                i=csEnerStart-1;
            }

            // we use the average CS of the closest two points instead of the interpolated cs to avoid the avgCS being set to zero in the case that both of the intepolated CS being zero (ie the boundaries)
            csNum1 = ((csVec[csIndex+1]+csVec[csIndex])/2)*yieldVec[j];

            for(; i<csSize+csEnerStart-2; i++)
            {
                if(enerCSVec[i]>enerVec[j+1])
                {
                    break;
                }
            }
            i--;

            csIndex=i-csEnerStart+1;
            if(csIndex<0)
            {
                csIndex=0;
                i=csEnerStart-1;
            }

            // we use the average CS of the closest two points instead of the interpolated cs to avoid the avgCS being set to zero in the case that both of the intepolated CS being zero (ie the boundaries)
            csNum2 = ((csVec[csIndex+1]+csVec[csIndex])/2)*yieldVec[j+1];
            csSum += abs((csNum2+csNum1)*(enerVec[j+1]-enerVec[j])/2);
        }
        avgCS = csSum/(enerVec[numIncEner-1]-enerVec[0]);
        setAvgCS = true;
    }

    return avgCS;
}

double CSDistYieldComp::GetAvgCS(double ener)
{
    int i, low, reg=0, csIndex;
    double csNum1, csNum2;

    for(i=0; i<numIncEner-1; i++)
    {
        while((i>=regEndPos[reg])&&(numRegs-1>reg))
            reg++;
        if(enerVec[i]>ener)
        {
            break;
        }
    }
    i--;
    if(i<0)
        i=0;

    for(low=csEnerStart-1; low<csSize+csEnerStart-2; low++)
    {
        if(enerCSVec[low]>enerVec[i])
        {
            break;
        }
    }

    low--;
    csIndex=low-csEnerStart+1;
    if(csIndex<0)
    {
        csIndex=0;
        low=csEnerStart-1;
    }

    if(csSize>1)
        csNum1 = max(0.,(enerVec[i]-enerCSVec[low])*(csVec[csIndex+1]-csVec[csIndex])/(enerCSVec[low+1]-enerCSVec[low])+csVec[csIndex]);
    else
        csNum1 = csVec[0];

    if(numIncEner>1)
    {
        for(low=csEnerStart-1; low<csSize+csEnerStart-2; low++)
        {
            if(enerCSVec[low]>enerVec[i+1])
            {
                break;
            }
        }
        low--;
        csIndex=low-csEnerStart+1;
        if(csIndex<0)
        {
            csIndex=0;
            low=csEnerStart-1;
        }

        if(csSize>1)
            csNum2 = max(0.,(enerVec[i+1]-enerCSVec[low])*(csVec[csIndex+1]-csVec[csIndex])/(enerCSVec[low+1]-enerCSVec[low])+csVec[csIndex]);
        else
            csNum2 = csVec[0];
    }
    else
    {
        csNum2 = csNum1;
    }

    if(numIncEner>1)
        return (yieldVec[i]*csNum1+yieldVec[i+1]*csNum2)/2;
    else
        return yieldVec[0]*csNum1;
}
