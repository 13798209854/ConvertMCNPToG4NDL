#include "AngDistP32EqBin.hh"

AngDistP32EqBin::AngDistP32EqBin(int startAngularPDist)
{
    startAngPDist = startAngularPDist;
}

AngDistP32EqBin::~AngDistP32EqBin()
{
    /*
    if(incNEner)
        delete [] incNEner;
    if(probDistPos)
        delete [] probDistPos;
    for(int i=0; i<numIncEner, i++)
    {
        if(angVec[i])
            delete [] angVec[i];
    }
    if(angVec)
        delete [] angVec;
    */

    for(int i=0; i<int(incNEnerVec.size()), i++)
    {
        if(angVec[i])
            delete [] angVec[i];
    }
}

void AngDistP32EqBin::ExtractMCNPData(stringstream stream, int &count)
{
    /*double temp;
    int intTemp;
    string dummy;

    stream >> numIncEner; count++;
    incNEner = new double[numIncEner];
    probDistPos = new int[numIncEner];
    angVec = new double *[numIncEner];
    for (int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incNEner[i] = temp;
    }
    for (int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        probDistPos[i] = intTemp;
    }

    for(int i=0; i<33; i++)
    {
        for(int j=0; j<(startAngPDist+probDistPos[i]-1); j++, count++)
        {
            stream >> dummy;
        }
        angVec[i] = new double[33];
        for(int k=0; k<33; k++, count++)
        {
            stream >> temp;
            angVec.back()[k]=temp;
        }
    }
    */
}

void AngDistP32EqBin::WriteG4NDLData(stringstream data)
{
    for(int i=0; i<int(incNEnerVec.size()); i++)
    {
        stream << std::setw(14) << std::right << incNEnerVec[i] << 32 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 32 << 1 << '\n';
        for(int j=0; j<32; j++)
        {
            // note the histogram scheme is right biased, we checked
            stream << std::setw(14) << std::right << angVec[i][j] << 1/(32*(angVec[i][j+1]-angVec[i][j])) << '\n';
        }
    }
}
