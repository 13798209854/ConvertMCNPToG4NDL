#include "EnerDistGenEvapSpec.hh"

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

void EnerDistGenEvapSpec::ExtractMCNPData(stringstream stream, int &count)
{
    int intTemp;
    double temp;

    stream >> numRegs; count++;
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
void EnerDistGenEvapSpec::WriteG4NDLData(stringstream stream)
{
    //this MCNP energy dist law 5
    //Convert to G4NDL theRepresentationType 5
    // since probability values are not associated with the given multipliers in the
    // mcnp data, assume that the multipliers are uniformly probable
    //or figure out what kind of a probability distribution MCNP uses for the evaporation scheme

    stream << numIncEner << '\n';
    stream << numRegs << '\n';
    for(int i=0; i<numRegs; i++)
    {
        stream << regEndPos[i] << intScheme[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << incEner[i]*1000000 << normOutEnerDist[i] << '\n';
    }

    stream << numNormOutEnerDistPoints << '\n';
    stream << 1 << '\n';
    stream << numNormOutEnerDistPoints << 2 << '\n';

    for(int i=0; i<numNormOutEnerDistPoints; i++)
    {
        stream << outEnerMulti[i] << 1/(outEnerMulti.size()) << '\n';
    }
}
