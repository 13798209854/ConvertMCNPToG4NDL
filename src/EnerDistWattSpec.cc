#include "EnerDistWattSpec.hh"

EnerDistWattSpec::EnerDistWattSpec()
{
    //ctor
}

EnerDistWattSpec::~EnerDistWattSpec()
{
    if(regEndPosA)
        delete [] regEndPosA;
    if(regEndPosB)
        delete [] regEndPosB;
    if(intSchemeA)
        delete [] intSchemeA;
    if(intSchemeB)
        delete [] intSchemeB;
    if(incEnerA)
        delete [] incEnerA;
    if(incEnerB)
        delete [] incEnerB;
    if(aValues)
        delete [] aValues;
    if(bValues)
        delete [] bValues;

}

void EnerDistWattSpec::ExtractMCNPData(stringstream stream, int &count)
{
    int intTemp;
    double temp;

    stream >> numRegsA; count++;
    regEndPosA = new int[numRegsA];
    intSchemeA = new int[numRegsA];

    for(int i=0; i<numRegsA; i++, count++)
    {
        stream >> intTemp;
        regEndPosA[i]=intTemp;
    }

    for(int i=0; i<numRegsA; i++, count++)
    {
        stream >> intTemp;
        intSchemeA[i]=intTemp;
    }

    stream >> numIncEnerA; count++;
    incEnerA = new double[numIncEnerA];
    aValues = new double [numIncEnerA];

    for(int i=0; i<numIncEnerA; i++, count++)
    {
        stream >> temp;
        incEnerA[i]=temp;
    }
    for(int i=0; i<numIncEnerA; i++, count++)
    {
        stream >> intTemp;
        aValues[i]=intTemp;
    }

    stream >> numRegsB; count++;
    regEndPosB = new int[numRegsB];
    intSchemeB = new int[numRegsB];

    for(int i=0; i<numRegsB; i++, count++)
    {
        stream >> intTemp;
        regEndPosB[i]=intTemp;
    }

    for(int i=0; i<numRegsB; i++, count++)
    {
        stream >> intTemp;
        intSchemeB[i]=intTemp;
    }

    stream >> numIncEnerB; count++;
    incEnerB = new double[numIncEnerB];
    bValues = new double [numIncEnerB];

    for(int i=0; i<numIncEnerB; i++, count++)
    {
        stream >> temp;
        incEnerB[i]=temp;
    }
    for(int i=0; i<numIncEnerB; i++, count++)
    {
        stream >> intTemp;
        bValues[i]=intTemp;
    }

    stream >> rejectEner; count++;
}

//For Fission
void EnerDistWattSpec::WriteG4NDLData(stringstream data)
{
//this is MCNP energy distribution law 11
//convert this law to G4NDL law 11
//we ignore the rejection energy since G4NDL does not use it
//check the physics to make sure this is equivalent, appears to be the same except for the added coefficient in the MCNP version

    stream << numIncEnerA << '\n';
    stream << numRegsA << '\n';
    for(int i=0; i<numRegsA; i++)
    {
        stream << regEndPosA[i] << intSchemeA[i] << '\n';
    }

    for(int i=0; i<numIncEnerA; i++)
    {
        stream << incEnerA[i]*1000000 << aValues[i] << '\n';
    }

    stream << numIncEnerB << '\n';
    stream << numRegsB << '\n';
    for(int i=0; i<numRegsB; i++)
    {
        stream << regEndPosB[i] << intSchemeB[i] << '\n';
    }

    for(int i=0; i<numIncEnerB; i++)
    {
        stream << incEnerB[i]*1000000 << bValues[i] << '\n';
    }
}
