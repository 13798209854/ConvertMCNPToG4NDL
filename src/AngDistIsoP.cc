#include "AngDistIsoP.hh"

AngDistIsoP::AngDistIsoP()
{
    //ctor
}

AngDistIsoP::~AngDistIsoP()
{
    //dtor
}

void AngDistIsoP::ExtractMCNPData(stringstream stream, int &count)
{

}

//set up for capture files
void AngDistIsoP::WriteG4NDLData(stringstream data)
{
    if(incNEnerVec.size()>0)
    {
        for(int i=0; i<int(incNEnerVec.size()); i++)
        {
            // creates a normalized isotropic distribution
            stream << std::setw(14) << std::right << incNEnerVec[i]*1000000 << 2 << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << 2 << 2;
            stream << std::setw(14) << std::right << -1 << 0.5 << '\n';
            stream << std::setw(14) << std::right << 1 << 0.5 << '\n';
        }
    }
    else
    {
        // creates a normalized isotropic distribution
        stream << std::setw(14) << std::right << 0.0 << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << 2;
        stream << std::setw(14) << std::right << -1 << 0.5 << '\n';
        stream << std::setw(14) << std::right << 1 << 0.5 << '\n';

        stream << std::setw(14) << std::right << 2.0e+07 << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << 2;
        stream << std::setw(14) << std::right << -1 << 0.5 << '\n';
        stream << std::setw(14) << std::right << 1 << 0.5 << '\n';
    }

}
