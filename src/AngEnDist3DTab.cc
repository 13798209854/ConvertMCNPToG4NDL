#include "../include/AngEnDist3DTab.hh"

#ifndef NeutHPMod_Use
#define NeutHPMod_Use 0
#endif

AngEnDist3DTab::AngEnDist3DTab(int EnerDistStart)
{
    startEnerDist =  EnerDistStart;
}

AngEnDist3DTab::~AngEnDist3DTab()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme1)
        delete [] intScheme1;
    if(intScheme2)
        delete [] intScheme2;
    if(outEnDistPos)
        delete [] outEnDistPos;
    if(incEner)
        delete [] incEner;
    if(numDiscreteEnerPoints)
        delete [] numDiscreteEnerPoints;

    for(int i=0; i<numIncEner; i++)
    {
        if(outEner[i])
            delete [] outEner[i];
        if(outEnerProb[i])
            delete [] outEnerProb[i];
        if(outEnerSumProb[i])
            delete [] outEnerSumProb[i];
        if(isoDist[i])
            delete [] isoDist[i];
        if(outAngDistPos[i])
            delete [] outAngDistPos[i];
        if(intScheme3[i])
            delete [] intScheme3[i];
        if(numPAngPoints[i])
            delete [] numPAngPoints[i];
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            if(outAng[i][j])
                delete [] outAng[i][j];
            if(outAngProb[i][j])
                delete [] outAngProb[i][j];
            if(outAngSumProb[i][j])
                delete [] outAngSumProb[i][j];
        }
        if(outAng[i])
            delete [] outAng[i];
        if(outAngProb[i])
            delete [] outAngProb[i];
        if(outAngSumProb[i])
            delete [] outAngSumProb[i];
    }
    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(outEner)
        delete [] outEner;
    if(outEnerProb)
        delete [] outEnerProb;
    if(outEnerSumProb)
        delete [] outEnerSumProb;
    if(isoDist)
        delete [] isoDist;
    if(outAngDistPos)
        delete [] outAngDistPos;
    if(intScheme3)
        delete [] intScheme3;
    if(numPAngPoints)
        delete [] numPAngPoints;
    if(outAng)
        delete [] outAng;
    if(outAngProb)
        delete [] outAngProb;
    if(outAngSumProb)
        delete [] outAngSumProb;
}

void AngEnDist3DTab::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;
    double outEnerPropSumTemp;

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
    outEnDistPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPEnerPoints = new int[numIncEner];
    numDiscreteEnerPoints = new int[numIncEner];
    outEner = new double *[numIncEner];
    outEnerProb = new double *[numIncEner];
    outEnerSumProb = new double *[numIncEner];
    isoDist = new bool *[numIncEner];
    outAngDistPos = new int *[numIncEner];
    intScheme3 = new int *[numIncEner];
    numPAngPoints = new int *[numIncEner];
    outAng = new double **[numIncEner];
    outAngProb = new double **[numIncEner];
    outAngSumProb = new double **[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        outEnDistPos[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++)
    {
        //Not needed and dangerous since we don't know exactly where outEnDistPos[i] points
        /*
        for(;count<(startEnerDist+outEnDistPos[i]-1); count++)
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
        isoDist[i] = new bool [numPEnerPoints[i]];
        outAngDistPos[i] = new int [numPEnerPoints[i]];
        intScheme3[i] = new int [numPEnerPoints[i]];
        numPAngPoints[i] = new int [numPEnerPoints[i]];
        outAng[i] = new double *[numPEnerPoints[i]];
        outAngProb[i] = new double *[numPEnerPoints[i]];
        outAngSumProb[i] = new double *[numPEnerPoints[i]];

        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEner[i][j] = temp;
        }
        outEnerPropSumTemp=0.;
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerProb[i][j] = temp;
            outEnerPropSumTemp+=temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerSumProb[i][j] = temp;
        }
//        if(outEnerPropSumTemp!=0.)
//        {
//            for(int j=0; j<numPEnerPoints[i]; j++)
//            {
//                outEnerProb[i][j] /= outEnerPropSumTemp;
//            }
//        }
        for(int j=0; j<numPEnerPoints[i]; j++, count++)
        {
            stream >> temp;
            outAngDistPos[i][j] = temp;
        }
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            if(outAngDistPos[i][j]>0)
            {
                for(;count<(startEnerDist+outAngDistPos[i][j]-1); count++)
                {
                    stream >> dummy;
                }
            }
            else
            {
                isoDist[i][j]=true;
                continue;
            }

            stream >> intTemp; count++;
            if(intTemp==0)
                intScheme3[i][j] = 1;
            else
                intScheme3[i][j] = 2;

            stream >> intTemp; count++;
            numPAngPoints[i][j] = intTemp;

            outAng[i][j] = new double [numPAngPoints[i][j]];
            outAngProb[i][j] = new double [numPAngPoints[i][j]];
            outAngSumProb[i][j] = new double [numPAngPoints[i][j]];

            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outAng[i][j][k] = temp;
            }
            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outAngProb[i][j][k] = temp;
            }
            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outAngSumProb[i][j][k] = temp;
            }
        }
    }
}

#if NeutHPMod_Use
void AngEnDist3DTab::WriteG4NDLData(stringstream &stream)
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
            stream << std::setw(14) << std::right << numPAngPoints[i][j];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPAngPoints[i][j] << std::setw(14) << std::right << intScheme3[i][j] << '\n';

            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                stream << std::setw(14) << std::right << outAng[i][j][k];
                stream << std::setw(14) << std::right << outAngProb[i][j][k];
                if(((k+1)%3==0)||(k==numPEnerPoints[i]-1))
                    stream << '\n';
            }
        }
    }
}
#else
void AngEnDist3DTab::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 61
    //convert this to G4NDL DistLaw=7
    //some approximations are used in the transformation

    // this creates a new out-going angle array, for each incoming energy. each of the angles in the new angle arrays
    // has a new probability dist which is dependant on the out-going energy. We did this to bring the 3D table format
    // (order of properties) inline with that used by G4NDL DistLaw=7
    int *sumAngPoints = new int [numIncEner], *newNumAngPoints = new int [numIncEner];
    vector<double> *outAngNew = new vector<double> [numIncEner];
    double ***outEnProbNew = new double** [numIncEner];

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
    // here we set correct the energy prob so that it is integrated over its energy regime
    for(int i=0; i<numIncEner; i++)
    {
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            for(int k=0; k<numPAngPoints[i][j]; k++)
            {
                if(k==0)
                    outAngProb[i][j][k] = outAngSumProb[i][j][k];
                else
                    outAngProb[i][j][k] = outAngSumProb[i][j][k]-outAngSumProb[i][j][k-1];
            }
        }
    }

    for(int i=0; i<numIncEner; i++)
    {
        sumAngPoints[i]=0;
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            sumAngPoints[i]+=numPAngPoints[i][j];
            //take this out later when we have fixed the G4NeutronHpVector so that it doesn't make multiple copies of 0 every time it merges two data sets

            int k=0, l=0;
            while((k<numPAngPoints[i][j])&&(l<int(outAngNew[i].size())))
            {
                if(outAng[i][j][k]<outAngNew[i][l])
                {
                    outAngNew[i].insert(outAngNew[i].begin()+l, outAng[i][j][k]);
                    k++; l++;
                }
                else if(outAng[i][j][k]>outAngNew[i][l])
                {
                    l++;
                }
                else
                {
                    k++; l++;
                }
            }

            for(;k<numPAngPoints[i][j];k++)
            {
                outAngNew[i].push_back(outAng[i][j][k]);
            }
        }
    }

    for(int i=0; i<numIncEner; i++)
    {
        newNumAngPoints[i] = outAngNew[i].size();
        outEnProbNew[i] = new double* [newNumAngPoints[i]];

        for(int j=0; j<newNumAngPoints[i]; j++)
        {
            outEnProbNew[i][j] = new double [numPEnerPoints[i]];

            for(int k=0; k<numPEnerPoints[i]; k++)
            {
                int l;
                for(l=0; l<numPAngPoints[i][k]-1; l++)
                {
                    if(outAng[i][k][l]>outAngNew[i][j])
                    {
                        break;
                    }
                }
                l--;
                if(l<0)
                    l=0;
                //added on angle width so that the data would be properly interpereted by G4NeutronHPLabAngularEnergy.cc
                outEnProbNew[i][j][k]=outEnerProb[i][k]*max(0.,Interpolate(intScheme3[i][k], outAngNew[i][j], outAng[i][k][l], outAng[i][k][l+1], outAngProb[i][k][l], outAngProb[i][k][l+1]));
            }
        }
    }

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
        stream << std::setw(14) << std::right << newNumAngPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << newNumAngPoints[i] << std::setw(14) << std::right << 2 << '\n';

//        sum1=0.;
//        for(int j=0; j<newNumAngPoints[i]; j++)
//        {
//            sum2 = 0.;
//            for(int k=0; k<numPEnerPoints[i]; k++)
//            {
//                sum2 += outEnProbNew[i][j][k];
//            }
//            sum1+=sum2;
//            if(sum2<=0.)
//            {
//                //cout << "Error with angular energy probability data" << endl;
//            }
//        }
//        if(sum1<=0)
//        {
//            cout << "Error with angular energy probability data" << endl;
//        }

        for(int j=0; j<newNumAngPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outAngNew[i][j];
            stream << std::setw(14) << std::right << numPEnerPoints[i];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

            for(int k=0; k<numPEnerPoints[i]; k++)
            {
                stream << std::setw(14) << std::right << outEner[i][k]*1000000;
                stream << std::setw(14) << std::right << outEnProbNew[i][j][k];// /sum1;
                if(((k+1)%3==0)||(k==numPEnerPoints[i]-1))
                    stream << '\n';
            }
        }
    }

    for(int i=0; i<numIncEner; i++)
    {
        for(int j=0; j<newNumAngPoints[i]; j++)
        {
            delete [] outEnProbNew[i][j];
        }
        delete [] outEnProbNew[i];
    }

    delete [] outEnProbNew;
    delete [] outAngNew;
    delete [] sumAngPoints;
    delete [] newNumAngPoints;
}
#endif

void AngEnDist3DTab::ConvertToEnerAndAngDist(EnergyDist **enDist, AngularDist **angDist, int &numAngEner)
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
            while((k<numPAngPoints[i][j])&&(l<int(outAngConv[i].size())))
            {
                if(outAng[i][j][k]<outAngConv[i][l])
                {
                    outAngConv[i].insert(outAngConv[i].begin()+l, outAng[i][j][k]);
                    k++; l++;
                }
                else if(outAng[i][j][k]>outAngConv[i][l])
                {
                    l++;
                }
                else
                {
                    k++; l++;
                }
            }

            for(;k<numPAngPoints[i][j];k++)
            {
                outAngConv[i].push_back(outAng[i][j][k]);
            }
        }

        sum=0.;
        outAngProbConv[i].assign(outAngConv[i].size(), 0.);
        for(int j=0; j<numPEnerPoints[i]; j++)
        {
            int low=0;
            for(int l=0; l<int(outAngConv[i].size()); l++)
            {
                for(; low<numPAngPoints[i][j]-1; low++)
                {
                    if(outAng[i][j][low]>outAngConv[i][l])
                    {
                        break;
                    }
                }
                if(low>0)
                    low--;

                if(numPAngPoints[i][j]>1)
                    outAngProbConv[i][l] += Interpolate(2, outAngConv[i][l], outAng[i][j][low], outAng[i][j][low+1], outAngProb[i][j][low], outAngProb[i][j][low+1])*outEnerProb[i][j];
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

    if(!(angDist[0]->CheckData()))
    {
        cout << "Error in angular-energy data AngEnDist3DTab.cc:454" << endl;
    }

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
        if(sumEn==0.)
        {
            cout << "Error in AngEnDist3DTab::ConvertToEnerAndAngDist" << endl;
        }
    }
    enDist[0] = new EnerDistConTab(numRegs, regEndPos, intScheme1, numIncEner, incEner, intScheme2, outEnerConv, outEnerProbConv);

    delete [] outAngConv;
    delete [] outAngProbConv;
    delete [] outEnerConv;
    delete [] outEnerProbConv;
}
