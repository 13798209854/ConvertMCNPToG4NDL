#include "EnerDist1PhEner.hh"

EnerDist1PhEner::EnerDist1PhEner(double AWR)
{
    awr=AWR;
}

EnerDist1PhEner::~EnerDist1PhEner()
{
    //dtor
}

void EnerDist1PhEner::ExtractMCNPData(stringstream stream, int &count)
{
    stream >> photonType >> photonEn; count=count+2;
}
void EnerDist1PhEner::WriteG4NDLData(stringstream data)
{

}
