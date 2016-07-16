#include "../include/AngEnDistLab3DTab.hh"

AngEnDistLab3DTab::AngEnDistLab3DTab(/*int EnerDistStart*/)
{
    //startEnerDist =  EnerDistStart;
}

AngEnDistLab3DTab::~AngEnDistLab3DTab()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme1)
        delete [] intScheme1;
    if(intScheme2)
        delete [] intScheme2;
    if(outAngDistPos)
        delete [] outAngDistPos;
    if(incEner)
        delete [] incEner;

    for(int i=0; i<numIncEner; i++)
    {
        if(outAng[i])
            delete [] outAng[i];
        if(outEnerDistPos[i])
            delete [] outEnerDistPos[i];
        if(intScheme3[i])
            delete [] intScheme3[i];
        if(numPEnerPoints[i])
            delete [] numPEnerPoints[i];
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            if(outEner[i][j])
                delete [] outEner[i][j];
            if(outEnProb[i][j])
                delete [] outEnProb[i][j];
            if(outEnSumProb[i][j])
                delete [] outEnSumProb[i][j];
        }
        if(outEner[i])
            delete [] outEner[i];
        if(outEnProb[i])
            delete [] outEnProb[i];
        if(outEnSumProb[i])
            delete [] outEnSumProb[i];
    }
    if(numPAngPoints)
        delete [] numPAngPoints;
    if(outAng)
        delete [] outAng;
    if(outEnerDistPos)
        delete [] outEnerDistPos;
    if(intScheme3)
        delete [] intScheme3;
    if(numPEnerPoints)
        delete [] numPEnerPoints;
    if(outEner)
        delete [] outEner;
    if(outEnProb)
        delete [] outEnProb;
    if(outEnSumProb)
        delete [] outEnSumProb;
}

void AngEnDistLab3DTab::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;

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
    outAngDistPos = new int[numIncEner];
    intScheme2 = new int[numIncEner];
    numPAngPoints = new int[numIncEner];
    outAng = new double *[numIncEner];
    outEnerDistPos = new int *[numIncEner];
    intScheme3 = new int *[numIncEner];
    numPEnerPoints = new int *[numIncEner];
    outEner = new double **[numIncEner];
    outEnProb = new double **[numIncEner];
    outEnSumProb = new double **[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> intTemp;
        outAngDistPos[i]=intTemp;
    }
    for(int i=0; i<numIncEner; i++)
    {
        //not needed, and potentially eroneous
        /*
        for(;count<(startEnerDist+outAngDistPos[i]-1); count++)
        {
            stream >> dummy;
        }
        */

        stream >> intTemp; count++;
        intScheme2[i]=intTemp;

        stream >> intTemp; count++;
        numPAngPoints[i]=intTemp;
        outAng[i] = new double [numPAngPoints[i]];
        outEnerDistPos[i] = new int [numPAngPoints[i]];
        intScheme3[i] = new int [numPAngPoints[i]];
        numPEnerPoints[i] = new int [numPAngPoints[i]];
        outEner[i] = new double *[numPAngPoints[i]];
        outEnProb[i] = new double *[numPAngPoints[i]];
        outEnSumProb[i] = new double *[numPAngPoints[i]];

        for(int j=0; j<numPAngPoints[i]; j++, count++)
        {
            stream >> temp;
            outAng[i][j] = temp;
        }
        for(int j=0; j<numPAngPoints[i]; j++, count++)
        {
            stream >> temp;
            outEnerDistPos[i][j] = temp;
        }
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            //not needed, and potentially erroneous
            /*
            for(;count<(startEnerDist+outEnerDistPos[i]-1); count++)
            {
                stream >> dummy;
            }
            */
            stream >> intTemp; count++;
            intScheme3[i][j] = intTemp;

            stream >> intTemp; count++;
            numPEnerPoints[i][j] = intTemp;

            outEner[i][j] = new double [numPEnerPoints[i][j]];
            outEnProb[i][j] = new double [numPEnerPoints[i][j]];
            outEnSumProb[i][j] = new double [numPEnerPoints[i][j]];

            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEner[i][j][k] = temp;
            }
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEnProb[i][j][k] = temp;
            }
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream >> temp; count++;
                outEnSumProb[i][j][k] = temp;
            }
        }
    }
}

void AngEnDistLab3DTab::WriteG4NDLData(stringstream &stream)
{
    //this is MCNP Law 67
    //convert this to G4NDL DistLaw=7

    stream << std::setw(14) << std::right << numIncEner << std::setw(14) << std::right << numRegs << '\n';

    for(int i=0; i<numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme1[i] << '\n';
    }

    double sum;
    for(int i=0; i<numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << numPAngPoints[i];
        // assume linear interpolation
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << numPAngPoints[i] << std::setw(14) << std::right << intScheme2[i] << '\n';

        for(int j=0; j<numPAngPoints[i]; j++)
        {
            stream << std::setw(14) << std::right << outAng[i][j];
            stream << std::setw(14) << std::right << numPEnerPoints[i][j];
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numPEnerPoints[i][j] << std::setw(14) << std::right << intScheme3[i][j] << '\n';

            sum=0.;
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                stream << std::setw(14) << std::right << outEner[i][j][k]*1000000;
                stream << std::setw(14) << std::right << outEnProb[i][j][k];
                sum += outEnProb[i][j][k];
                if(((k+1)%3==0)||(k==numPEnerPoints[i][j]-1))
                    stream << '\n';
            }
            if(sum<=0.)
            {
                cout << "Error with angular energy probability data" << endl;
            }
        }
    }
}

void AngEnDistLab3DTab::ConvertToEnerAndAngDist(EnergyDist **enDist, AngularDist **angDist, int &numAngEner)
{
    if(enDist[0])
        delete enDist[0];
    if(angDist[0])
        delete angDist[0];

    vector<double> *outAngConv = new vector<double> [numIncEner];
    vector<double> *outAngProbConv = new vector<double> [numIncEner];
    double angSum;
    for(int i=0; i<numIncEner; i++)
    {
        angSum=0.;
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            outAngConv[i].push_back(outAng[i][j]);
            outAngProbConv[i].push_back(0.);
            for(int k=0; k<numPEnerPoints[i][j]; k++)
            {
                outAngProbConv[i][j]+=outEnProb[i][j][k];
                angSum+=outEnProb[i][j][k];
            }
        }
        if(angSum!=0.)
        {
            for(int j=0; j<int(outAngProbConv[i].size()); j++)
            {
                outAngProbConv[i][j]/=angSum;
            }
        }
        else
        {
            cout << "Angular probability is zero for all angles at this incoming neutron energy" << endl;
            for(int l=0; l<int(outAngProbConv[i].size()); l++)
            {
                outAngProbConv[i][l] = 1.0/outAngProbConv[i].size();
            }
        }
    }

    numAngEner+=numIncEner;
    angDist[0] = new AngDist2DTabular(numIncEner, incEner, intScheme2, outAngConv, outAngProbConv);

    vector<double>* outSumEn = new vector<double> [numIncEner];
    vector<double>* outSumEnProb = new vector<double> [numIncEner];
    double sum;
    for(int i=0; i<numIncEner; i++)
    {
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            int k=0, l=0;
            while((k<numPEnerPoints[i][j])&&(l<int(outSumEn[i].size())))
            {
                if(outEner[i][j][k]<outSumEn[i][l])
                {
                    outSumEn[i].insert(outSumEn[i].begin()+l, outEner[i][j][k]);
                    k++; l++;
                }
                else if(outEner[i][j][k]>outSumEn[i][l])
                {
                    l++;
                }
                else
                {
                    k++; l++;
                }
            }

            for(;k<numPEnerPoints[i][j];k++)
            {
                outSumEn[i].push_back(outEner[i][j][k]);
            }
        }

        sum=0.;
        outSumEnProb[i].assign(outSumEn[i].size(), 0.);
        for(int j=0; j<numPAngPoints[i]; j++)
        {
            int low=0;
            for(int l=0; l<int(outSumEn[i].size()); l++)
            {
                for(; low<numPEnerPoints[i][j]-1; low++)
                {
                    if(outEner[i][j][low]>outSumEn[i][l])
                    {
                        break;
                    }
                }
                if(low>0)
                    low--;

                if(numPEnerPoints[i][j]>1)
                    outSumEnProb[i][l] += Interpolate(2, outSumEn[i][l], outEner[i][j][low], outEner[i][j][low+1], outEnProb[i][j][low], outEnProb[i][j][low+1]);
                else
                    outSumEnProb[i][l] += outEnProb[i][j][low];
                sum += outSumEnProb[i][l];
            }
        }
        if(sum!=0.)
        {
            for(int l=0; l<int(outSumEn[i].size()); l++)
            {
                outSumEnProb[i][l] /= sum;
            }
        }
        else
        {
            cout << "out-going energy dist prob zero for all out going energies at this incoming neutron energy" << endl;
            for(int l=0; l<int(outSumEn[i].size()); l++)
            {
                outSumEnProb[i][l] = 1.0/outSumEn[i].size();
            }
        }
    }

    enDist[0] = new EnerDistConTab(numRegs, regEndPos, intScheme1, numIncEner, incEner, intScheme2, outSumEn, outSumEnProb);

    delete [] outAngConv;
    delete [] outAngProbConv;
    delete [] outSumEn;
    delete [] outSumEnProb;
}
