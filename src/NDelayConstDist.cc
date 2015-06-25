#include "../include/NDelayConstDist.hh"

NDelayConstDist::NDelayConstDist(int numDNFam)
{
    numDNPrecursorFam = numDNFam;
}

NDelayConstDist::~NDelayConstDist()
{
    for(int i=0; i<int(regEndPos.size()); i++)
    {
        if(regEndPos[i])
            delete [] regEndPos[i];
    }
    for(int i=0; i<int(intScheme.size()); i++)
    {
        if(intScheme[i])
            delete [] intScheme[i];
    }
    for(int i=0; i<int(incEner.size()); i++)
    {
        if(incEner[i])
            delete [] incEner[i];
    }
    for(int i=0; i<int(prob.size()); i++)
    {
        if(prob[i])
            delete [] prob[i];
    }
}

void NDelayConstDist::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;

    for(int i=0; i<numDNPrecursorFam; i++)
    {
        stream >> temp >> intTemp; count=count+2;

        delConst.push_back(temp);
        numRegs.push_back(intTemp);
        if(numRegs.back()==0)
        {
            numRegs.back()=1;
            regEndPos.push_back(new int[numRegs.back()]);
            intScheme.push_back(new int[numRegs.back()]);

            stream >> intTemp; count++;
            regEndPos.back()[0]=intTemp;
            intScheme.back()[0]=2;
        }
        else
        {
            regEndPos.push_back(new int[numRegs.back()]);
            intScheme.push_back(new int[numRegs.back()]);

            for(int i=0; i<numRegs.back(); i++, count++)
            {
                stream >> intTemp;
                regEndPos.back()[i]=intTemp;
            }

            for(int i=0; i<numRegs.back(); i++, count++)
            {
                stream >> intTemp;
                intScheme.back()[i]=intTemp;
            }
            stream >> intTemp; count++;
        }

        numIncEner.push_back(intTemp);
        incEner.push_back(new double[numIncEner.back()]);
        prob.push_back(new double[numIncEner.back()]);

        for(int i=0; i<numIncEner.back(); i++, count++)
        {
            stream >> temp;
            incEner.back()[i]=temp;
        }
        for(int i=0; i<numIncEner.back(); i++, count++)
        {
            stream >> temp;
            prob.back()[i]=temp;
        }
    }
}

void NDelayConstDist::WriteG4NDLData(stringstream &stream)
{
    stream << std::setw(14) << std::right << numDNPrecursorFam;
    for(int i=0; i<numDNPrecursorFam; i++)
    {
        if(i%6==0)
            stream << endl;
        stream << std::setw(14) << std::right << delConst[i];
    }

}
