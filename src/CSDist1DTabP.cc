#include "CSDist1DTabP.hh"

CSDist1DTabP::~CSDist1DTabP()
{
    delete [] CSVec;
}

CSDist1DTabP::void ExtractMCNPData(stringstream data, int count&)
{
    int intTemp;
    double temp;

    if(!elastic)
    {
        stream >> intTemp; count++;
        startEner = intTemp;
        stream >> intTemp; count++;
        CSVecSize = intTemp;
    }

    CSVec = new double[CSVecSize];

    for(int j=0; j<CSVecSize; j++, count++)
    {
        stream >> temp;
        CSVec[j] = temp;
    }
}

CSDist1DTabP::void WriteG4NDLCSData(stringstream stream )
{
    // may have to supplement G4NDL data for this section, hopefully this format is not needed
    // ask Dr. Buijs about how to convert this section with out knowing the photon energy
    stream << std::setw(14) << std::right << CSVecSize << std::setw(14) << std::right << 1 << endl;
    stream << CSVecSize << 2 << endl;

    stream << '\n';
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int j=0; j<CSVecSize; j++)
    {
        stream << std::setw(14) << std::right << enerCSVec[j+startEner]*1000000;
        stream << std::setw(14) << std::right << CSVec[j];
    }
}

CSDist1DTabP::void WriteG4NDLYieldData(stringstream stream )
{
    double csNum;
    stream << std::setw(14) << std::right << CSVecSize << std::setw(14) << std::right << 1 << endl;
    stream << CSVecSize << 2 << endl;

    for(int i=0; i<CSVecSize; i++)
    {
        stream << std::setw(14) << std::right << enerCSVec[i+startEner]*1000000;
        stream << std::setw(14) << std::right << CSVec[i]/csVec[i+startEner-csEnerStart];
    }
}

CSDist1DTabP::void IdentifyYourSelf()
{
    return "CSDist1DTabP";
}

double CSDist1DTabP::Interpolate(double x)
{
    int i;
    for(i=0; i<CSVecSize; i++)
    {
        if(enerVec[i+startEner-1]>x)
        {
            i--;
            break;
        }
    }
    if(i<0)
        i=0;

    return Interpolate(2, x, enerVec[i+startEner-1], enerVec[i+startEner], CSVec[i+startEner-1], CSVec[i+startEner]);
}
