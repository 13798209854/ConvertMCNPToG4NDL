#include "EnerDistEvapSpec.hh"

EnerDistEvapSpec::EnerDistEvapSpec()
{
    //ctor
}

EnerDistEvapSpec::~EnerDistGenEvapSpec()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme)
        delete [] intScheme;
    if(incEner)
        delete [] incEner;
    if(tValues)
        delete [] tValues;
}

void EnerDistEvapSpec::ExtractMCNPData(stringstream stream, int &count)
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
    tValues = new double[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        tValues[i]=temp;
    }

    stream >> restrictEner; count++;
}

//For Fission
void EnerDistEvapSpec::WriteG4NDLData(stringstream data)
{
//this is MCNP law 9
//convert this to G4NDL theRepresentationType=9
//check the physics to make sure this is equivalent, appears to be the same except for the added coefficient in the MCNP version
//we ignore the provided restriction energy since G4NDL does not ask for it

    stream << numIncEner << '\n';
    stream << numRegs << '\n';
    for(int i=0; i<numRegs; i++)
    {
        stream << regEndPos[i] << intScheme[i] << '\n';
    }

    for(int i=0; i<numIncEner; i++)
    {
        stream << incEner[i] << tValues[i] << '\n';
    }

}
