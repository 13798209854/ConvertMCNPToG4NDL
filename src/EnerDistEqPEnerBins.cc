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

    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numReg << '\n';
    for(int i=0; i<numReg; i++)
    {
        stream << std::setw(14) << regEndPos[i] << std::setw(14) << std::right << intScheme[i];
    }
    stream << '\n';
    for(int i=0; i<int(numIncEner); i++)
    {
        // note the histogram scheme is right biased, we checked
        stream << std::setw(14) << std::right << incEner[i]*1000000 << std::setw(14) << std::right << numOutEnerPerInc << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numOutEnerPerInc << std::setw(14) << std::right << 1 << '\n';

        for(int j=0; j<int(numOutEnerPerInc); j++)
        {
            stream << std::setw(14) << std::right << outEner[i][j]*1000000 << std::setw(14) << std::right << 1/(numOutEnerPerInc*(outEner[i][j+1]-outEner[i][j]));
            if(j%3==0)
                stream << '\n';
        }
    }
}
