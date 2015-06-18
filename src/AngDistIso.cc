#include "../include/AngDistIso.hh"

AngDistIso::AngDistIso()
{
    //ctor
}

AngDistIso::~AngDistIso()
{
    //dtor
}

void AngDistIso::ExtractMCNPData(stringstream &stream, int &count)
{

}

//set up for elastic files
void AngDistIso::WriteG4NDLData(stringstream &stream)
{
    if(incNEnerVec.size()>0)
    {
        for(int i=0; i<int(incNEnerVec.size()); i++)
        {
            // creates a normalized isotropic distribution
            stream << std::setw(14) << std::right << temperature << std::setw(14) << std::right << incNEnerVec[i]*1000000
                    << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2;
            stream << std::setw(14) << std::right << -1 << std::setw(14) << std::right << 0.5 << '\n';
            stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 0.5 << '\n';
        }
    }
    else
    {
        // creates a normalized isotropic distribution
        stream << std::setw(14) << std::right << temperature << std::setw(14) << std::right << 0.0 << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2;
        stream << std::setw(14) << std::right << -1 << std::setw(14) << std::right << 0.5 << '\n';
        stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 0.5 << '\n';

        stream << std::setw(14) << std::right << temperature << std::setw(14) << std::right << 2.0e+07 << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2;
        stream << std::setw(14) << std::right << -1 << std::setw(14) << std::right << 0.5 << '\n';
        stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 0.5 << '\n';
    }

}
