#include "NYield1DTab.hh"

NYield1DTab::NYield1DTab()
{
    //ctor
}

NYield1DTab::~NYield1DTab()
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

void NYield1DTab::ExtractMCNPData(stringstream stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;

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
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        yield[i]=temp;
    }
}

void NYield1DTab::WriteG4NDLData(stringstream data)
{


}
