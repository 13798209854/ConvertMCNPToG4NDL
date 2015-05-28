#include "CSDistYieldComp.hh"

CSDistYieldComp::CSDistYieldComp()
{
    //ctor
}

CSDistYieldComp::~CSDistYieldComp()
{
    //dtor
}

CSDistYieldComp::void ExtractMCNPData(stringstream data, int count&)
{
    int intTemp;
    double temp;

    stream >> numRegs; count++;
    regEndPos = new int[numRegs];
    intScheme = new int[numRegs];

    for(int i=0; i<numRegs; i++, count++)
    {
        stream >> intTemp;
        regEndPos[i]=intTemp;
    }

    for(int i=0; i<numRegs; i++, count++)
    {
        stream >> intTemp;
        intScheme[i]=intTemp;
    }

    stream >> numIncEner; count++;
    enerVec = new double[numIncEner];
    yieldVec = new double[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        enerVec[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        yieldVec[i]=intTemp;
    }
}

//set up for Capture data
CSDistYieldComp::void WriteG4NDLCSData(stringstream stream )
{
    double csNum;
    stream << std::setw(14) << std::right << enerVec.size() << std::setw(14) << std::right << numRegs << endl;

    for(int i=0; i<numRegs; i++)
    {
        stream << regEndPos[i] << intScheme[i] << endl;
    }
    stream << '\n';
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int j=0; j<enerVec.size(); j++)
    {
        stream << std::setw(14) << std::right << enerVec[j];
        for(i=0; i<cSSize; i++)
        {
            // assume average incoming neutron energy is 1eV
            if(enerCSVec[i]>enerVec[j])
                break;
        }
        if(i==0)
            i++;
        csNum = (enerVec[i]-enerCSVec[i-1])*(csVec[i]-csVec[i-1])/(enerCSVec[i]-enerCSVec[i-1])+csVec[i-1];
        stream << std::setw(14) << std::right << (yieldVec[j]*csNum);
    }
}

CSDistYieldComp::void WriteG4NDLYieldData(stringstream stream )
{
    stream << std::setw(14) << std::right << enerVec.size() << std::setw(14) << std::right << numRegs << endl;

    for(int i=0; i<numRegs; i++)
    {
        stream << regEndPos[i] << intScheme[i] << endl;
    }
    stream << '\n';
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    for(int i=0; i<enerVec.size(); i++)
    {
        stream << std::setw(14) << std::right << enerVec[i];
        stream << std::setw(14) << std::right << yieldVec[i];
    }
}

string CSDistYieldComp::IdentifyYourSelf()
{
    return "CSDistYieldComp";
}
