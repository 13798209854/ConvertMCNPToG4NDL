#include "CSDist1DTab.hh"

CSDist1DTab::CSDist1DTab(double *energyVec)
{
    enerVec=energyVec;
}

CSDist1DTab::CSDist1DTab(double *energyVec, int numCSEner, int startEnerTable)
{
    startEner = startEnerTable;
    CSVecSize = numCSEner;
    elastic = true;
}

CSDist1DTab::~CSDist1DTab()
{
    delete [] CSVec;
}

CSDist1DTab::void ExtractMCNPData(stringstream data, int count&)
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

//Set up for Cross-section data
CSDist1DTab::void WriteG4NDLCSData(stringstream stream )
{
    stream << std::setw(14) << std::right << '0' << '\n';
    stream << std::setw(14) << std::right << CSVecSize;
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int i=0; i<CSVecSize; i++)
    {
        if(floor(i/3)==ceil(double(i)/3))
        {
            stream << '\n';
        }
        stream << std::setw(14) << std::right << enerVec[i+startEner-1];
        stream << std::setw(14) << std::right << CSVec[i];
    }
    stream << '\n';
}

CSDist1DTab::void IdentifyYourSelf()
{
    return "CSDist1DTab";
}
