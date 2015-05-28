#include "CSDist1DTabP.hh"

CSDist1DTabP::CSDist1DTabP(double *energyVec)
{
    enerVec=energyVec;
}

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
}

CSDist1DTabP::void WriteG4NDLYieldData(stringstream stream )
{
    // may have to supplement G4NDL data for this section, hopefully this format is not needed
    // ask Dr. Buijs about how to convert this section with out knowing the photon energy
}

CSDist1DTabP::void IdentifyYourSelf()
{
    return "CSDist1DTabP";
}
