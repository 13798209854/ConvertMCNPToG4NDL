#include "AngDistIso.hh"

AngDistIso::AngDistIso()
{
    //ctor
}

AngDistIso::~AngDistIso()
{
    //dtor
}

void AngDistIso::ExtractMCNPData(stringstream stream, int &count)
{

}

//set up for elastic files
void AngDistIso::WriteG4NDLData(stringstream data)
{
    if(incNEnerVec.size()>0)
    {
        for(int i=0; i<int(incNEnerVec.size()); i++)
        {
            // creates a normalized isotropic distribution
            stream << std::setw(14) << std::right << temperature << incNEnerVec[i] << 0 << 2 << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << 2 << 2;
            stream << std::setw(14) << std::right << -1 << 0.5 << '\n';
            stream << std::setw(14) << std::right << 1 << 0.5 << '\n';
        }
    }
    else
    {
        // creates a normalized isotropic distribution
        stream << std::setw(14) << std::right << temperature << 0.0 << 0 << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << 2;
        stream << std::setw(14) << std::right << -1 << 0.5 << '\n';
        stream << std::setw(14) << std::right << 1 << 0.5 << '\n';

        stream << std::setw(14) << std::right << temperature << 2.0e+07 << 0 << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << 2;
        stream << std::setw(14) << std::right << -1 << 0.5 << '\n';
        stream << std::setw(14) << std::right << 1 << 0.5 << '\n';
    }

}
