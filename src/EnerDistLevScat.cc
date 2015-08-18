#include "../include/EnerDistLevScat.hh"

EnerDistLevScat::EnerDistLevScat(double enerRangeStart, double enerRangeEnd)
{
    enerStart=enerRangeStart;
    enerEnd=enerRangeEnd;
}

EnerDistLevScat::~EnerDistLevScat()
{
    //dtor
}

void EnerDistLevScat::ExtractMCNPData(stringstream &stream, int &count)
{
    stream >> firstHalfEq >> secondHalfEq; count=count+2;
}

//For Fission
void EnerDistLevScat::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 3
//there is no direct translation for this law in G4NDL but it can be made to fit theRepresentationType=1 with some assumptions
//using two different incoming energies use the collected linear function to find the corresponding out going energies then
// enter these 2 incoming energy and outgoing energy points with equal probability and set the interpolation scheme to linear between the incoming points

    //assume that G4NDL wants the out-going energy in CM frame
    double eOut1 = secondHalfEq*(enerStart-firstHalfEq), eOut2 = secondHalfEq*(enerEnd-firstHalfEq);

    stream << std::setw(14) << std::right << 2 << '\n' << 1 << '\n';
    stream << std::setw(14) << 2 << std::setw(14) << std::right << 2 << '\n';

    // note the histogram scheme is right biased, we checked
    stream << std::setw(14) << std::right << enerStart*1000000 << std::setw(14) << std::right << 2 << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << eOut1*1000000*0.999 << std::setw(14) << std::right << 1.0 << '\n';
    stream << std::setw(14) << std::right << eOut1*1000000*1.001+1e-6 << std::setw(14) << std::right << 1.0 << '\n';

    stream << std::setw(14) << std::right << enerEnd*1000000 << std::setw(14) << std::right << 2 << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << eOut2*1000000*0.999 << std::setw(14) << std::right << 1.0 << '\n';
    stream << std::setw(14) << std::right << eOut2*1000000*1.001+1e-6 << std::setw(14) << std::right << 1.0 << '\n';
    stream << '\n';

}
