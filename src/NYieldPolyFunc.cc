#include "../include/NYieldPolyFunc.hh"

NYieldPolyFunc::NYieldPolyFunc()
{
    //ctor
}

NYieldPolyFunc::~NYieldPolyFunc()
{
    if (coeff)
        delete [] coeff;
}

void NYieldPolyFunc::ExtractMCNPData(stringstream &stream, int &count)
{
    double temp;
    stream >> numCoeff; count++;
    coeff = new double[numCoeff];
    for(int i=0; i< numCoeff; i++, count++)
    {
        stream >> temp;
        coeff[i] = temp;
    }
}

void NYieldPolyFunc::WriteG4NDLData(stringstream &stream)
{
    stream << std::setw(14) << std::right << numCoeff << '\n';
    for(int i=0; i< numCoeff; i++)
    {
        stream << std::setw(14) << std::right << coeff[i];
        if(((i+1)%6==0)||(i==numCoeff-1))
            stream << '\n';
    }
}

void NYieldPolyFunc::SubtractPrompt(YieldDist* &promptYieldDist)
{
    if(promptYieldDist->IdentifyYourSelf()=="NYieldPolyFunc")
    {
        promptYieldDist->SubtractPrompt(coeff, numCoeff);
    }
    else
    {
        cout << "### Warning: approximation caused by the mixing of two schemes" << endl;
    }
}

void NYieldPolyFunc::SubtractPrompt(double *&totalCoeff, int &totalNumCoeff)
{
    if(totalNumCoeff<numCoeff)
    {
        double *temp = new double [numCoeff];
        for(int i=0; i<totalNumCoeff; i++)
        {
            temp[i]=totalCoeff[i];
        }
        for(int i=totalNumCoeff; i<numCoeff; i++)
        {
            temp[i]=0;
        }

        delete [] totalCoeff;
        totalCoeff = temp;
        totalNumCoeff=numCoeff;
    }
    for(int i=0; i<numCoeff; i++)
    {
        totalCoeff[i]-=coeff[i];
    }
}

void NYieldPolyFunc::SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield)
{
    //this is an approximation since there is no way to exactly merge the two schemes into one
    cout << "### Warning: approximation caused by the mixing of two schemes" << endl;
    int j;
    double sum, energy;
    for(int i=0; i<totalNumIncEner; i++)
    {
        sum=0.;
        energy=1.;
        for(j=0; j<numCoeff; j++)
        {
            sum+=coeff[j]*energy;
            energy*=totalIncEner[i];
        }

        totalYield[i]-=sum;
    }

}

void NYieldPolyFunc::ConvertToLinDist(int *regEndPos, int &numIncEner, double *&incEner, double *&yield)
{
    int j;
    double sum, energy;

    numIncEner=(numCoeff-2)*100+2;
    if(numIncEner<1)
        numIncEner=1;
    regEndPos[0]=numIncEner;

    if(incEner)
        delete incEner;
    if(yield)
        delete yield;

    incEner = new double [numIncEner];
    yield = new double [numIncEner];

    for(int i=0; i<numIncEner; i++)
    {
        incEner[i]=20*(i/numIncEner);
        sum=0.;
        energy=1.;
        for(j=0; j<numCoeff; j++)
        {
            sum+=coeff[j]*energy;
            energy*=incEner[i];
        }
        yield[i]=sum;
    }
}
