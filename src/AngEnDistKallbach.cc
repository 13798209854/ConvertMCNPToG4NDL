#include "../include/AngEnDistKallbach.hh"

#ifndef NeutHPMod_Use
#define NeutHPMod_Use 0
#endif

AngEnDistKallbach::AngEnDistKallbach(/*int EnerDistStart*/)
{
    /*startEnerDist =  EnerDistStart;*/
}

AngEnDistKallbach::~AngEnDistKallbach()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme1)
        delete [] intScheme1;
    if(intScheme2)
        delete [] intScheme2;
    if(distPos)
        delete [] distPos;
    if(incEner)
        delete [] incEner;
    if(numDiscreteEnerPoints)
        delete [] numDiscreteEnerPoints;

    for(int i=0; i<numIncEner; i++)
    {
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            if(outAngProb[i][j])
                delete [] outAngProb[i][j];
        }
        if(outAngProb[i])
            delete [] outAngProb[i];
        if(outEner[i])
            delete [] outEner[i];
        if(outEnerProb[i])
            delete [] outEnerProb[i];
        if(outEnerSumProb[i])
            delete [] outEnerSumProb[i];
        if(rFraction[i])
            delete [] rFraction[i];
        if(angDistSlope[i])
            delete [] angDistSlope[i];
    }

    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(outAngProb)
        delete [] outAngProb;
    if(outEner)
        delete [] outEner;
    if(outEnerProb)
        delete [] outEnerProb;
    if(outEnerSumProb)
        delete [] outEnerSumProb;
    if(rFraction)
        delete [] rFraction;
    if(angDistSlope)
        delete [] angDistSlope;
}

void AngEnDistKallbach::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;

    for(int i=0; i<numDistSample; i++)
    {
        outAng[i] = -1+double(2*i)/numDistSample;
    }

    stream >> numRegs; count++;
    if(numRegs==0)
    {
        numRegs=1;
        regEndPos = new int[numRegs];
        intScheme1 = new int[numRegs];

        stream >> numIncEner; count++;
        regEndPos[0]=numIncEner;
        intScheme1[0]=2;
    }
    else
    {
        regEndPos = new int[numRegs];
        intScheme1 = new int[numRegs];

        for(int i=0; i<numRegs; i++, count++)
        {
            stream >> intTemp;
            regEndPos[i]=intTemp;
        }

        for(int i=0; i<numRegs; i++, count++)
        {
            stream >> intTemp;
            intScheme1[i]=intTemp;
        }
        stream >> numIncEner; count++;
    }

    incEner = new double[numIncEner];
    distPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPEnerPoints = new int[numIncEner];
    numDiscreteEnerPoints = new int[numIncEner];
    outEner = new double *[numIncEner];
    outEnerProb = new double *[numIncEner];
    outEnerSumProb = new double *[numIncEner];
    rFraction = new double *[numIncEner];
    angDistSlope = new double *[numIncEner];
    outAngProb = new double **[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        distPos[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++)
    {
        //this is not needed, and potentially erroneous
        /*
        for(;count<(startEnerDist+distPos[i]-1); count++)
        {
            stream >> dummy;
        }
        */

        stream >> intTemp; count++;
        if(intTemp>10)
        {
            numDiscreteEnerPoints[i] = int(intTemp/10);
            intScheme2[i]=intTemp-numDiscreteEnerPoints[i]*10;
        }
        else
        {
            numDiscreteEnerPoints[i] = 0;
            intScheme2[i]=intTemp;
        }

        stream >> intTemp; count++;
        numPEnerPoints[i]=intTemp;
        outEner[i] = new double [numPEnerPoints[i]];
        outEnerProb[i] = new double [numPEnerPoints[i]];
        outEnerSumProb[i] = new double [numPEnerPoints[i]];
        rFraction[i] = new double [numPEnerPoints[i]];
        angDistSlope[i] = new double [numPEnerPoints[i]];
        outAngProb[i] = new double *[numPEnerPoints[i]];

        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEner[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerProb[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerSumProb[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            rFraction[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            angDistSlope[i][j] = temp;
        }
        // create the angle probability distribution from the given function for this out-going energy
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            outAngProb[i][j] = new double [numDistSample];
            for(int k=0; k<numDistSample; k++)
            {
                outAngProb[i][j][k] = 0.5*angDistSlope[i][j]*(std::cosh(angDistSlope[i][j]*outAng[k])+rFraction[i][j]*std::sinh(angDistSlope[i][j]*outAng[k]))/(std::sinh(angDistSlope[i][j]));
            }
        }
    }

}

#if NeutHPMod_Use
void AngEnDistKallbach::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 61
    //convert this to G4NDL mod DistLaw=61
    //some approximations are used in the transformation

    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }
    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << numPEnerPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outEner[i][j]*1000000;
            stream << std::setw(14) << std::right << outEnerSumProb[i][j];
            stream << std::setw(14) << std::right << numDistSample;
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numDistSample << std::setw(14) << std::right << 2 << '\n';

            for(int k=0; k<numDistSample; k++)
            {
                stream << std::setw(14) << std::right << outAng[k];
                stream << std::setw(14) << std::right << outAngProb[i][j][k];
                if(((k+1)%3==0)||(k==numPEnerPoints[i]-1))
                    stream << '\n';
            }
        }
    }
}
#else
void AngEnDistKallbach::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 44
    //convert this to G4NDL DistLaw=7
    // may be able to exactly convert it to DistLaw=1

    // here we set correct the energy prob so that it is integrated over its energy regime
//    for(int i=0; i<numIncEner; i++)
//    {
//        for(int j=0; j<numPEnerPoints[i]; j++)
//        {
//            if(j==0)
//            {
//                if(outEner[i][j]==0.)
//                    outEner[i][j]=1.0e-12;
//
//                if(numPEnerPoints[i]==2)
//                    outEnerProb[i][j] = outEnerSumProb[i][j+1];
//                else
//                    outEnerProb[i][j] = outEnerSumProb[i][j];
//            }
//            else
//                outEnerProb[i][j] = outEnerSumProb[i][j]-outEnerSumProb[i][j-1];
//        }
//    }

     stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    double sum1, sum2;
    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << numDistSample;
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numDistSample << std::setw(14) << std::right << 2 << '\n';

//        sum1=0.;
//        for(int j=0; j<numDistSample; j++)
//        {
//            sum2 = 0.;
//            for(int k=0; k<numPEnerPoints[i]; k++)
//            {
//                sum2 += outEnerProb[i][k]*outAngProb[i][k][j];
//            }
//            sum1+=sum2;
//            if(sum2<=0.)
//            {
//                cout << "Error with angular energy probability data" << endl;
//            }
//        }

        for(int j=0; j<numDistSample; j++)
        {
            stream << std::setw(14) << std::right << outAng[j];
            stream << std::setw(14) << std::right << numPEnerPoints[i];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

            for(int k=0; k<numPEnerPoints[i]; k++)
            {
                stream << std::setw(14) << std::right << outEner[i][k]*1000000;
                stream << std::setw(14) << std::right << outEnerProb[i][k]*outAngProb[i][k][j];// /sum1*numDistSample;
                if(((k+1)%3==0)||(k==numPEnerPoints[i]-1))
                    stream << '\n';
            }
        }
    }
}
#endif
void AngEnDistKallbach::ConvertToEnerAndAngDist(EnergyDist **enDist, AngularDist **angDist, int &numAngEner)
{
    if(enDist[0])
        delete enDist[0];
    if(angDist[0])
        delete angDist[0];

    vector<double> *outAngConv = new vector<double> [numIncEner];
    vector<double> *outAngProbConv = new vector<double> [numIncEner];
    double sum;
    for(int i=0; i<numIncEner; i++)
    {
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            int k=0, l=0;
            while((k<numDistSample)&&(l<int(outAngConv[i].size())))
            {
                if(outAng[k]<outAngConv[i][l])
                {
                    outAngConv[i].insert(outAngConv[i].begin()+l, outAng[k]);
                    k++; l++;
                }
                else if(outAng[k]>outAngConv[i][l])
                {
                    l++;
                }
                else
                {
                    k++; l++;
                }
            }

            for(;k<numDistSample;k++)
            {
                outAngConv[i].push_back(outAng[k]);
            }
        }

        sum=0.;
        outAngProbConv[i].assign(outAngConv[i].size(), 0.);
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            int low=0;
            for(int l=0; l<int(outAngConv[i].size()); l++)
            {
                for(; low<numDistSample-1; low++)
                {
                    if(outAng[low]>outAngConv[i][l])
                    {
                        break;
                    }
                }
                if(low>0)
                    low--;

                if(numDistSample>1)
                    outAngProbConv[i][l] += Interpolate(2, outAngConv[i][l], outAng[low], outAng[low+1], outAngProb[i][j][low], outAngProb[i][j][low+1])*outEnerProb[i][j];
                else
                    outAngProbConv[i][l] += outAngProb[i][j][low]*outEnerProb[i][j];
                sum += outAngProbConv[i][l];
            }
        }
        if(sum!=0.)
        {
            for(int l=0; l<int(outAngConv[i].size()); l++)
            {
                outAngProbConv[i][l] /= sum;
            }
        }
        else
        {
            cout << "Angular probability is zero for all angles at this incoming neutron energy" << endl;
            for(int l=0; l<int(outAngConv[i].size()); l++)
            {
                outAngProbConv[i][l] = 1.0/outAngConv[i].size();
            }
        }
    }

    numAngEner+=numIncEner;
    angDist[0] = new AngDist2DTabular(numIncEner, incEner, intScheme2, outAngConv, outAngProbConv);


    //convert outEner and outEnerProb into vector<double> *
    vector<double> *outEnerConv = new vector<double> [numIncEner];
    vector<double> *outEnerProbConv = new vector<double> [numIncEner];
    double sumEn;
    for(int i=0; i<numIncEner; i++)
    {
        sumEn=0.;
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            outEnerConv[i].push_back(outEner[i][j]);
            outEnerProbConv[i].push_back(outEnerProb[i][j]);
            sumEn+=outEnerProb[i][j];
        }
        if(sumEn<=0.)
        {
            cout << "Error in AngEnDistKallbach::ConvertToEnerAndAngDist" << endl;
        }
    }
    enDist[0] = new EnerDistConTab(numRegs, regEndPos, intScheme1, numIncEner, incEner, intScheme2, outEnerConv, outEnerProbConv);

    delete [] outAngConv;
    delete [] outAngProbConv;
    delete [] outEnerConv;
    delete [] outEnerProbConv;
}
