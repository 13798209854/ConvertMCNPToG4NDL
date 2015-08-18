#include "../include/NYield1DTab.hh"

NYield1DTab::NYield1DTab()
{
    //ctor
}

NYield1DTab::NYield1DTab(YieldDist *totalYield)
{
    if(totalYield->IdentifyYourSelf()!="NYield1DTab")
    {
        cout << "### Warning: conversion from polynomial to segmented linear distribution" << endl;

        numRegs = 1;
        regEndPos = new int[numRegs];
        intScheme = new int[numRegs];
        intScheme[0]=2;
        totalYield->ConvertToLinDist(regEndPos, numIncEner, incEner, yield);
    }
    else
    {
        cout << "### Error: NYield1DTab(YieldDist *totalYield) is only meant for conversion from polynomial to segmented linear distribution, \n" <<
                "not for copying NYield1DTab objects" << endl;
    }
}

NYield1DTab::~NYield1DTab()
{
    if(regEndPos)
        delete [] regEndPos;
    if(intScheme)
        delete [] intScheme;
    if(incEner)
        delete [] incEner;
    if(yield)
        delete [] yield;
}

void NYield1DTab::ExtractMCNPData(stringstream &stream, int &count)
{
    int intTemp;
    double temp;
    string dummy;

    stream >> numRegs; count++;
    if(numRegs==0)
    {
        numRegs=1;
        regEndPos = new int[numRegs];
        intScheme = new int[numRegs];

        stream >> numIncEner; count++;
        regEndPos[0]=numIncEner;
        intScheme[0]=2;
    }
    else
    {
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
    }

    incEner = new double[numIncEner];
    yield = new double[numIncEner];

    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        incEner[i]=temp;
    }
    for(int i=0; i<numIncEner; i++, count++)
    {
        stream >> temp;
        yield[i]=temp;
    }
}

void NYield1DTab::WriteG4NDLData(stringstream &stream)
{
    stream << std::setw(14) << std::right << numIncEner << '\n';
    stream << std::setw(14) << std::right << numRegs << '\n';
    for(int i=0; i < numRegs; i++)
    {
        stream << std::setw(14) << std::right << regEndPos[i];
        stream << std::setw(14) << std::right << intScheme[i];
    }

    for(int i=0; i < numIncEner; i++)
    {
        stream << std::setw(14) << std::right << incEner[i]*1000000;
        stream << std::setw(14) << std::right << yield[i];
        if((i%3==0)||(i==numIncEner-1))
            stream << '\n';
    }
}

void NYield1DTab::SubtractPrompt(YieldDist* &promptYieldDist)
{
    promptYieldDist->SubtractPrompt(numIncEner, incEner, yield);
}

void NYield1DTab::SubtractPrompt(int totalNumIncEner, double *totalIncEner, double *totalYield)
{
    //this is an approximation, two different interpolation schemes cannot be exactly represented by one so a rough fit is used
    cout << "### Warning: approximation caused by the mixing of two schemes" << endl;

    int reg, low;
    for(int i=0; i<totalNumIncEner; i++)
    {
        reg=0;
        for(low=0; low<numIncEner-1; low++)
        {
            if(incEner[low]>totalIncEner[i])
            {
                break;
            }
            while((regEndPos[reg]<=low)&&(numRegs-1>reg))
                reg++;
        }
        if(low!=0)
            low--;

        totalYield[i]-=max(0.,Interpolate( intScheme[reg], totalIncEner[i], incEner[low], incEner[low+1], yield[low], yield[low+1]));
    }

}
