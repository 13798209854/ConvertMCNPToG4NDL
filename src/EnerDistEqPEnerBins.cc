#include "../include/EnerDistEqPEnerBins.hh"

EnerDistEqPEnerBins::EnerDistEqPEnerBins()
{
    //ctor
}

EnerDistEqPEnerBins::~EnerDistEqPEnerBins()
{
    if(regEndPos)
        delete [] regEndPos;

    if(intScheme)
        delete [] intScheme;

    if(incEner)
        delete [] incEner;

    for(int i=0; i<numIncEner; i++)
    {
        if(outEner[i])
            delete [] outEner[i];
    }
    if(outEner)
        delete [] outEner;
}

void EnerDistEqPEnerBins::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;

    stream >> numReg; count++;
    if(numReg==0)
    {
        numReg=1;
        regEndPos = new int[numReg];
        intScheme = new int[numReg];

        stream >> numIncEner; count++;
        regEndPos[0]=numIncEner;
        intScheme[0]=2;
    }
    else
    {
        regEndPos = new int[numReg];
        intScheme = new int[numReg];

        for(int i=0; i<numReg; i++, count++)
        {
            stream >> intTemp;
            regEndPos[i]=intTemp;
        }

        for(int i=0; i<numReg; i++, count++)
        {
            stream >> intTemp;
            intScheme[i]=intTemp;
        }
        stream >> numIncEner; count++;
    }

    incEner = new double[numIncEner];
    outEner = new double *[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }

    stream >> numOutEnerPerInc; count++;

    for(int i=0; i<numIncEner; i++)
    {
        outEner[i]=new double[numOutEnerPerInc];

        for(int j=0; j<numOutEnerPerInc; j++, count++)
        {
            stream >> temp;
            outEner[i][j] = temp;
        }
    }
}

//For Fission
void EnerDistEqPEnerBins::WriteG4NDLData(stringstream &stream)
{
    // this is MCNP law 1
    // convert this to G4NDL theRepresentationType=1
    double enerRange;
    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numReg << '\n';
    for(int i=0; i<numReg; i++)
    {
        stream << std::setw(14) << regEndPos[i] << std::setw(14) << std::right << intScheme[i];
    }
    stream << '\n';
    for(int i=0; i<int(numIncEner); i++)
    {
        // note the histogram scheme is right biased, we checked
        if(numOutEnerPerInc>2)
        {
            stream << std::setw(14) << std::right << incEner[i]*1000000 << std::setw(14) << std::right << numOutEnerPerInc << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numOutEnerPerInc << std::setw(14) << std::right << 1 << '\n';

            enerRange=0.;
            for(int j=0; j<int(numOutEnerPerInc-1); j++)
            {
                stream << std::setw(14) << std::right << outEner[i][j]*1000000 << std::setw(14) << std::right << 1.0/(numOutEnerPerInc*(outEner[i][j+1]-outEner[i][j]));
                if(j>0)
                {
                    enerRange += outEner[i][j]-outEner[i][j-1];
                }
                if((j%3==0)||(j==numOutEnerPerInc-2))
                    stream << '\n';
            }
            if(enerRange==0.)
            {
                cout << "break here" << '\n';
            }
        }
        else if(numOutEnerPerInc<2)
        {
            stream << std::setw(14) << std::right << incEner[i]*1000000 << std::setw(14) << std::right << 0 << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 1 << '\n';
        }
        else
        {
            stream << std::setw(14) << std::right << incEner[i]*1000000 << std::setw(14) << std::right << 2 << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';

            for(int j=0; j<int(2); j++)
            {
                stream << std::setw(14) << std::right << outEner[i][j]*1000000 << std::setw(14) << std::right << 1.0/2;
                if((j>0)&&(outEner[i][j]==outEner[i][j-1]))
                {
                    cout << "break here" << '\n';
                }
                if((j%3==0)||(j==1))
                    stream << '\n';
            }
        }
    }
}
