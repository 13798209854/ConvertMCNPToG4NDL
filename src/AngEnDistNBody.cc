#include "../include/AngEnDistNBody.hh"

AngEnDistNBody::AngEnDistNBody()
{
    //ctor
}

AngEnDistNBody::~AngEnDistNBody()
{
    //dtor
}

void AngEnDistNBody::ExtractMCNPData(stringstream &stream, int &count)
{
    stream >> numBodies >> particleMassRatio; count = count+2;
}

void AngEnDistNBody::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 66
    //convert this to G4NDL DistLaw=6

    stream << std::setw(14) << std::right << particleMassRatio << std::setw(14) << std::right << numBodies << '\n';
}

