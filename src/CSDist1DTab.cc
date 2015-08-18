#include "../include/CSDist1DTab.hh"

CSDist1DTab::CSDist1DTab(double *energyVec)
{
    enerVec=energyVec;
}

CSDist1DTab::CSDist1DTab(double *energyVec, int numCSEner)
{
    startEner = 1;
    CSVecSize = numCSEner;
    enerVec= energyVec;
    elastic = true;
}

CSDist1DTab::~CSDist1DTab()
{
    if(CSVec)
        delete [] CSVec;
}

void CSDist1DTab::ExtractMCNPData(stringstream &stream, int &count)
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
void CSDist1DTab::WriteG4NDLCSData(stringstream &stream )
{
    stream << std::setw(14) << std::right << CSVecSize;
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int i=0; i<CSVecSize; i++)
    {
        if((i%3==0)||(i==CSVecSize-1))
            stream << '\n';
        stream << std::setw(14) << std::right << enerVec[i+startEner-1]*1000000;
        stream << std::setw(14) << std::right << CSVec[i];
    }
    stream << '\n';
}

string CSDist1DTab::IdentifyYourSelf()
{
    return "CSDist1DTab";
}

double CSDist1DTab::Interpolate(double x)
{
    int i;
    for(i=0; i<CSVecSize-1; i++)
    {
        if(enerVec[i+startEner-1]>x)
        {
            break;
        }
    }
    i--;
    if(i<0)
        i=0;

    return max(0.,CSDist::Interpolate(2, x, enerVec[i+startEner-1], enerVec[i+startEner], CSVec[i], CSVec[i+1]));
}

double CSDist1DTab::GetAvgCS()
{
    if(!setAvgCS)
    {
        double csSum=0.;
        for(int i=0; i<(CSVecSize-1); i++)
        {
            csSum+=abs((CSVec[i]+CSVec[i+1])*(enerVec[i+startEner]-enerVec[i+startEner-1])/2);
        }
        avgCS = csSum/(enerVec[startEner-1]-enerVec[startEner+CSVecSize-2]);
        setAvgCS = true;
    }

    return avgCS;
}

double CSDist1DTab::GetAvgCS(double ener)
{
    int i;
    for(i=0; i<CSVecSize-1; i++)
    {
        if(enerVec[i+startEner-1]>ener)
        {
            break;
        }
    }
    i--;
    if(i<0)
        i=0;

    return (CSVec[i]+CSVec[i+1])/2;
}
