#include "../include/AngDist2DTabular.hh"

AngDist2DTabular::AngDist2DTabular()
{
    //ctor
}

AngDist2DTabular::~AngDist2DTabular()
{
    for(int i=0; i<int(angVec.size()); i++)
    {
        if(angVec[i])
            delete [] angVec[i];
        if(angProbVec[i])
            delete [] angProbVec[i];
    }
}

void AngDist2DTabular::ExtractMCNPData(stringstream &stream, int &count)
{

}

//set up for elastic files
void AngDist2DTabular::WriteG4NDLData(stringstream &stream)
{
    for(int i=0; i<int(incNEnerVec.size()); i++)
    {
        stream << std::setw(14) << std::right << temperature << std::setw(14) << std::right << (incNEnerVec[i]*1000000) << std::setw(14) << std::right
                << 0 << std::setw(14) << std::right << numAngProb[i] << '\n';

        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numAngProb[i] << std::setw(14) << std::right << intSchemeAng[i] << '\n';
        for(int j=0; j<numAngProb[i]; j++)
        {
            stream << std::setw(14) << std::right << angVec[i][j] << std::setw(14) << std::right << angProbVec[i][j];
            if(j%3==0)
                stream << '\n';
        }
    }
}
