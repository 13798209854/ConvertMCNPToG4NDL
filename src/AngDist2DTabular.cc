#include "../include/AngDist2DTabular.hh"

AngDist2DTabular::AngDist2DTabular()
{
    //ctor
}

AngDist2DTabular::AngDist2DTabular(int numIncEnerTemp, double *incEnerTemp, int *intSchemeTemp, vector<double> *outAngTemp, vector<double> *outAngProbTemp)
{
    double sum;
    for(int i=0; i<numIncEnerTemp; i++)
    {
        incNEnerVec.push_back(incEnerTemp[i]);
        intSchemeAng.push_back(intSchemeTemp[i]);
        numAngProb.push_back(outAngTemp[i].size());
        angVec.push_back(new double [numAngProb[i]]);
        angProbVec.push_back(new double [numAngProb[i]]);
        sum = 0.;
        for(int j=0; j<numAngProb[i]; j++)
        {
            angVec[i][j]=outAngTemp[i][j];
            angProbVec[i][j]=outAngProbTemp[i][j];
            sum += angProbVec[i][j];
        }
        if(sum==0.)
        {
            cout << "Error with angular probability data AngDist2DTabular.cc:27" << endl;
        }
    }
}

AngDist2DTabular::AngDist2DTabular(AngularDist *angDist)
{
    angDist->SetData(incNEnerVec, angVec, angProbVec, intSchemeAng, numAngProb, temperature);
}

AngDist2DTabular::~AngDist2DTabular()
{
    for(int i=0; i<int(angVec.size()); i++)
    {
        if(angVec[i])
            delete [] angVec[i];
        if(angProbVec[i])
            delete [] angProbVec[i];
    }
}

void AngDist2DTabular::ExtractMCNPData(stringstream &stream, int &count)
{

}

//set up for elastic files
void AngDist2DTabular::WriteG4NDLData(stringstream &stream)
{
    double sum;
    for(int i=0; i<int(incNEnerVec.size()); i++)
    {
        stream << std::setw(14) << std::right << temperature << std::setw(14) << std::right << (incNEnerVec[i]*1000000) << std::setw(14) << std::right
                << 0 << std::setw(14) << std::right << numAngProb[i] << '\n';

        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numAngProb[i] << std::setw(14) << std::right << intSchemeAng[i] << '\n';
        sum=0.;
        for(int j=0; j<numAngProb[i]; j++)
        {
            stream << std::setw(14) << std::right << angVec[i][j] << std::setw(14) << std::right << angProbVec[i][j];
            sum+=angProbVec[i][j];
            if(((j+1)%3==0)||(j==numAngProb[i]-1))
                stream << '\n';
        }
        if(sum==0.)
        {
            cout << "Error with angular probability data AngDist2DTabular.cc:74" << endl;
        }
    }
}
