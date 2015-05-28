#include "AngDist2DTabular.cc.hh"

AngDist2DTabular.cc::AngDist2DTabular.cc()
{
    //ctor
}

AngDist2DTabular.cc::~AngDist2DTabular.cc()
{
    for(int i=0; i<int(incNEnerVec.size()), i++)
    {
        if(angVec[i])
            delete [] angVec[i];
        if(angProbVec[i])
            delete [] angProbVec[i];
    }
    if(angVec)
        delete [] angVec;
    if(angProbVec)
        delete [] angProbVec;
}

void AngDist2DTabular.cc::ExtractMCNPData(stringstream stream, int &count)
{

}

//set up for elastic files
void AngDist2DTabular.cc::WriteG4NDLData(stringstream data)
{
    for(int i=0; i<int(incNEnerVec.size()); i++)
    {
        stream << std::setw(14) << std::right << temperature << energyAngVec[i] << 0 << numAngProb[i] << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numAngProb[i] << intSchemeAng[i];
        for(int j=0; j<numAngProb[i] j++)
        {
            stream << std::setw(14) << std::right << angVec[i][j] << angProbVec[i][j] << '\n';
        }
    }
}
