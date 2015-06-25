#include "../include/EnerDistWattSpec.hh"

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

void EnerDistWattSpec::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;

    stream >> numRegsA; count++;
    if(numRegsA==0)
    {
        numRegsA=1;
        regEndPosA = new int[numRegsA];
        intSchemeA = new int[numRegsA];

        stream >> numIncEnerA; count++;
        regEndPosA[0]=numIncEnerA;
        intSchemeA[0]=2;
    }
    else
    {
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
    }

    incEnerA = new double[numIncEnerA];
    aValues = new double [numIncEnerA];

    for(int i=0; i<numIncEnerA; i++, count++)
    {
        stream >> temp;
        incEnerA[i]=temp;
    }
    for(int i=0; i<numIncEnerA; i++, count++)
    {
        stream >> temp;
        aValues[i]=temp;
    }

    stream >> numRegsB; count++;
    if(numRegsB==0)
    {
        numRegsB=1;
        regEndPosB = new int[numRegsB];
        intSchemeB = new int[numRegsB];

        stream >> numIncEnerB; count++;
        regEndPosB[0]=numIncEnerB;
        intSchemeB[0]=2;
    }
    else
    {
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
    }

    incEnerB = new double[numIncEnerB];
    bValues = new double [numIncEnerB];

    for(int i=0; i<numIncEnerB; i++, count++)
    {
        stream >> temp;
        incEnerB[i]=temp;
    }
    for(int i=0; i<numIncEnerB; i++, count++)
    {
        stream >> temp;
        bValues[i]=temp;
    }

    stream >> rejectEner; count++;
}

//For Fission
void EnerDistWattSpec::WriteG4NDLData(stringstream &stream)
{
//this is MCNP energy distribution law 11
//convert this law to G4NDL law 11
//we ignore the rejection energy since G4NDL does not use it
//check the physics to make sure this is equivalent, appears to be the same except for the added coefficient in the MCNP version

    stream << std::setw(14) << std::right << numIncEnerA << '\n';
    stream << std::setw(14) << std::right << numRegsA << '\n';
    for(int i=0; i<numRegsA; i++)
    {
        stream << std::setw(14) << std::right << regEndPosA[i] << std::setw(14) << std::right << intSchemeA[i] << '\n';
    }

    for(int i=0; i<numIncEnerA; i++)
    {
        stream << std::setw(14) << std::right << incEnerA[i]*1000000 << std::setw(14) << std::right << aValues[i] << '\n';
    }

    stream << std::setw(14) << std::right << numIncEnerB << '\n';
    stream << std::setw(14) << std::right << numRegsB << '\n';
    for(int i=0; i<numRegsB; i++)
    {
        stream << std::setw(14) << std::right << regEndPosB[i] << std::setw(14) << std::right << intSchemeB[i] << '\n';
    }

    for(int i=0; i<numIncEnerB; i++)
    {
        stream << std::setw(14) << std::right << incEnerB[i]*1000000 << std::setw(14) << std::right << bValues[i];
        if(i%3==0)
            stream << '\n';
    }
}
