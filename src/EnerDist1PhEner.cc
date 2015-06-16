#include "EnerDist1PhEner.hh"

EnerDist1PhEner::EnerDist1PhEner(double AWR, double, incEnerLow, double incEnerHigh)
{
    awr=AWR;
    incEnLow=incEnerLow;
    incEnHigh=incEnerHigh;
}

EnerDist1PhEner::~EnerDist1PhEner()
{
    //dtor
}

void EnerDist1PhEner::ExtractMCNPData(stringstream stream, int &count)
{
    stream >> photonType >> photonEn; count=count+2;
}
void EnerDist1PhEner::WriteG4NDLData(stringstream stream)
{
    stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n'
    stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2 << '\n';

    stream << std::setw(14) << std::right << incEnLow*1000000;
    stream << std::setw(14) << std::right << 1;
    // assume linear interpolation
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 1 << '\n';

    stream << std::setw(14) << std::right << photonEn*1000000;
    stream << std::setw(14) << std::right << 0.5 << '\n';

    stream << std::setw(14) << std::right << incEnerHigh*1000000;
    stream << std::setw(14) << std::right << 1;
    // assume linear interpolation
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 1 << '\n';

    stream << std::setw(14) << std::right << photonEn*1000000;
    stream << std::setw(14) << std::right << 0.5 << '\n';

}
