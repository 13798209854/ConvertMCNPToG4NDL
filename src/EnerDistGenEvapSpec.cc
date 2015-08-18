#include "../include/EnerDistGenEvapSpec.hh"

EnerDistGenEvapSpec::EnerDistGenEvapSpec()
{
    //ctor
}

EnerDistGenEvapSpec::~EnerDistGenEvapSpec()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme)
        delete [] intScheme;
    if(incEner)
        delete [] incEner;
    if(outEnerMulti)
        delete [] outEnerMulti;
    if(normOutEnerDist)
        delete [] normOutEnerDist;
}

void EnerDistGenEvapSpec::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;

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
    normOutEnerDist = new double[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        normOutEnerDist[i]=temp;
    }

    stream >> numNormOutEnerDistPoints; count++;
    outEnerMulti = new double[numNormOutEnerDistPoints];

    for(int i=0; i<numNormOutEnerDistPoints; i++, count++)
    {
        stream >> temp;
        outEnerMulti[i]=temp;
    }

}

//For Fission
void EnerDistGenEvapSpec::WriteG4NDLData(stringstream &stream)
{
    //this MCNP energy dist law 5
    //Convert to G4NDL theRepresentationType 5
    // since probability values are not associated with the given multipliers in the
    // mcnp data, assume that the multipliers are uniformly probable
    //or figure out what kind of a probability distribution MCNP uses for the evaporation scheme

    stream << std::setw(14) << std::right << numIncEner << '\n';
    stream << std::setw(14) << std::right << numRegs << '\n';
    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i] << std::setw(14) << std::right << intScheme[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000 << std::setw(14) << std::right << normOutEnerDist[i]*1000000 << '\n';
    }

    stream << std::setw(14) << std::right << numNormOutEnerDistPoints << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << numNormOutEnerDistPoints << std::setw(14) << std::right << 2 << '\n';

    //we assume that the outEnerMulti are randomly sample equaly but a different distribution could be used, try and find out
    for(int i=0; i<numNormOutEnerDistPoints; i++)
    {
        stream << std::setw(14) << std::right << outEnerMulti[i] << std::setw(14) << std::right << 1.0/(numNormOutEnerDistPoints);
        if((i%3==0)||(i==numNormOutEnerDistPoints-1))
            stream << '\n';
    }
}
