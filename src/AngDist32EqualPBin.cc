#include "AngDist32EqualPBin.hh"

AngDist32EqualPBin::AngDist32EqualPBin()
{
    //ctor
}

AngDist32EqualPBin::~AngDist32EqualPBin()
{
    for(int i=0; i<int(incNEnerVec.size()), i++)
    {
        if(angVec[i])
            delete [] angVec[i];
    }
}

void AngDist32EqualPBin::ExtractMCNPData(stringstream stream, int &count)
{
    //so far not needed
}

//set up for Elastic files
void AngDist32EqualPBin::WriteG4NDLData(stringstream stream)
{
    for(int i=0; i<int(incNEnerVec.size()); i++)
    {
        stream << std::setw(14) << std::right << temperature << incNEnerVec[i]*1000000 << 0 << 32 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 32 << 1 << '\n';
        for(int j=0; j<32; j++)
        {
            // note the histogram scheme is right biased, we checked
            stream << std::setw(14) << std::right << angVec[i][j] << 1/(32*(angVec[i][j+1]-angVec[i][j])) << '\n';
        }
    }
}
