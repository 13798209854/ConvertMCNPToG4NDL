#include "EnDepNYieldDist.hh"

EnDepNYieldDist::EnDepNYieldDist()
{
    //ctor
}

EnDepNYieldDist::~EnDepNYieldDist()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme)
        delete [] intScheme;
    if(incEner)
        delete [] incEner;
    if(yield)
        delete [] yield;
}

void EnDepNYieldDist::ExtractMCNPData(stringstream stream, int &count)
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
    yield = new double[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        incEner[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        yield[i]=intTemp;
    }
}

void EnDepNYieldDist::WriteG4NDLData(stringstream data)
{

}
