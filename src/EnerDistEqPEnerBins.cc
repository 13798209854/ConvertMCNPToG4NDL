#include "EnerDistEqPEnerBins.hh"

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

void EnerDistEqPEnerBins::ExtractMCNPData(stringstream stream, int &count)
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
    outEner = new double *[numIncEner]

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
void EnerDistEqPEnerBins::WriteG4NDLData(stringstream data)
{
    // this is MCNP law 1
    // convert this to G4NDL theRepresentationType=1

    stream << std::setw(14) << std::right << incEner.size() << numReg << '\n';
    for(int i=0; i<numReg; i++)
    {
        stream << std::setw(14) << regEndPos[i] << intScheme[i];
    }
    stream << '\n';
    for(int i=0; i<int(incEner.size()); i++)
    {
        // note the histogram scheme is right biased, we checked
        stream << std::setw(14) << std::right << incEner[i]*1000000 << outEner[i].size() << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << outEner[i].size() << 1 << '\n';

        for(int j=0; j<int(outEner[i].size()); j++)
        {
            stream << std::setw(14) << std::right << outEner[i][j]*1000000 << 1/(outEner[i].size()*(outEner[i][j+1]*1000000-outEner[i][j]*1000000)) << '\n';
        }
    }
}
