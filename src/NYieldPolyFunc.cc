#include "NYieldPolyFunc.hh"

NYieldPolyFunc::NYieldPolyFunc()
{
    //ctor
}

NYieldPolyFunc::~NYieldPolyFunc()
{
    if (coeff)
        delete [] coeff;
}

void NYieldPolyFunc::ExtractMCNPData(stringstream stream, int &count)
{
    double temp;
    stream >> numCoeff; count++;
    coeff = new double[numCoeff];
    for(int i=0; i< numCoeff; i++)
    {
        stream >> temp;
        coeff[i] = temp;
    }
}

void NYieldPolyFunc::WriteG4NDLData(stringstream data)
{

}
