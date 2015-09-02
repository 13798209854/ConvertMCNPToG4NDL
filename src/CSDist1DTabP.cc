#include "../include/CSDist1DTabP.hh"

CSDist1DTabP::CSDist1DTabP()
{

}

CSDist1DTabP::~CSDist1DTabP()
{
    delete [] CSVec;
}

void CSDist1DTabP::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;

    stream >> intTemp; count++;
    startEner = intTemp;
    stream >> intTemp; count++;
    CSVecSize = intTemp;

    CSVec = new double[CSVecSize];

    for(int j=0; j<CSVecSize; j++, count++)
    {
        stream >> temp;
        CSVec[j] = temp;
    }
}

void CSDist1DTabP::WriteG4NDLCSData(stringstream &stream )
{
    stream << std::setw(14) << std::right << CSVecSize << std::setw(14) << std::right << 1 << endl;
    stream << CSVecSize << 2 << endl;

    stream << '\n';
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int j=0; j<CSVecSize; j++)
    {
        stream << std::setw(14) << std::right << enerCSVec[j+startEner-1]*1000000;
        stream << std::setw(14) << std::right << CSVec[j];
        if(j%3==0)
            stream << '\n';
    }
}

void CSDist1DTabP::WriteG4NDLYieldData(stringstream &stream )
{
    stream << std::setw(14) << std::right << CSVecSize << std::setw(14) << std::right << 1 << endl;
    stream << std::setw(14) << std::right << CSVecSize << std::setw(14) << std::right << 2 << endl;

    for(int i=0; i<CSVecSize; i++)
    {
        stream << std::setw(14) << std::right << enerCSVec[i+startEner-1]*1000000;
        if((i+startEner-csEnerStart)>(csSize-1))
        {
            if(csVec[csSize-1]!=0.)
                stream << std::setw(14) << std::right << CSVec[i]/csVec[csSize-1];
            else
                stream << std::setw(14) << std::right << 0.;
        }
        else if(0>(i+startEner-csEnerStart))
        {
            if(csVec[0]!=0.)
                stream << std::setw(14) << std::right << CSVec[i]/csVec[0];
            else
                stream << std::setw(14) << std::right << 0.;
        }
        else
        {
            if(csVec[i+startEner-csEnerStart]!=0.)
                stream << std::setw(14) << std::right << CSVec[i]/csVec[i+startEner-csEnerStart];
            else
                stream << std::setw(14) << std::right << 0.;
        }
        if(((i+1)%3==0)||(i==CSVecSize-1))
            stream << '\n';
    }
}

string CSDist1DTabP::IdentifyYourSelf()
{
    return "CSDist1DTabP";
}

double CSDist1DTabP::Interpolate(double x)
{
    int i;
    for(i=0; i<CSVecSize-1; i++)
    {
        if(enerCSVec[i+startEner-1]>x)
        {
            break;
        }
    }
    i--;
    if(i<0)
        i=0;

    return max(0.,CSDist::Interpolate(2, x, enerCSVec[i+startEner-1], enerCSVec[i+startEner], CSVec[i], CSVec[i+1]));
}

double CSDist1DTabP::GetAvgCS()
{
    if(!setAvgCS)
    {
        double csSum=0.;
        for(int i=0; i<CSVecSize-1; i++)
        {
            csSum+=abs((CSVec[i]+CSVec[i+1])*(enerCSVec[i+startEner]-enerCSVec[i+startEner-1])/2);
        }
        avgCS = csSum/(enerCSVec[startEner-1]-enerCSVec[startEner+CSVecSize-2]);
        setAvgCS = true;
    }

    return avgCS;
}

double CSDist1DTabP::GetAvgCS(double ener)
{
    int i;
    for(i=0; i<CSVecSize-1; i++)
    {
        if(enerCSVec[i+startEner-1]>ener)
        {
            break;
        }
    }
    i--;
    if(i<0)
        i=0;

    return (CSVec[i]+CSVec[i+1])/2;
}
