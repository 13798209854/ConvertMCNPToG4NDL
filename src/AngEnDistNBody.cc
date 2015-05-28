#include "AngEnDistNBody.hh"

AngEnDistNBody::AngEnDistNBody()
{
    //ctor
}

AngEnDistNBody::~AngEnDistNBody()
{
    //dtor
}

void AngEnDistNBody::ExtractMCNPData(stringstream stream, int &count)
{
    stream >> numBodies >> particleMassRatio; count = count+2;
}

void AngEnDistNBody::WriteG4NDLData(stringstream data)
{

}

