#include "../include/AngDist2DTabularP.hh"

AngDist2DTabularP::AngDist2DTabularP()
{
    //ctor
}

AngDist2DTabularP::AngDist2DTabularP(AngularDist *angDist)
{
    angDist->SetData(incNEnerVec, angVec, angProbVec, intSchemeAng, numAngProb, temperature);
}

AngDist2DTabularP::~AngDist2DTabularP()
{
    for(int i=0; i<int(angVec.size()); i++)
    {
        if(angVec[i])
            delete [] angVec[i];
        if(angProbVec[i])
            delete [] angProbVec[i];
    }
}

void AngDist2DTabularP::ExtractMCNPData(stringstream &stream, int &count)
{

}

//set up for elastic files
void AngDist2DTabularP::WriteG4NDLData(stringstream &stream)
{
    double sum;
    for(int i=0; i<int(incNEnerVec.size()); i++)
    {
        stream  << std::setw(14) << std::right << (incNEnerVec[i]*1000000) << std::setw(14) << std::right << numAngProb[i] << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numAngProb[i] << std::setw(14) << std::right << intSchemeAng[i] << '\n';
        sum=0.;
        for(int j=0; j<numAngProb[i]; j++)
        {
            stream << std::setw(14) << std::right << angVec[i][j] << std::setw(14) << std::right << angProbVec[i][j];
            sum += angProbVec[i][j];
            if(((j+1)%3==0)||(j==numAngProb[i]-1))
                stream << '\n';
        }
        if(sum==0.)
        {
            cout << "Error with angular probability data AngDist2DTabularP.cc:48" << endl;
        }
    }
}
