#ifndef CSDIST1DTAB_HH
#define CSDIST1DTAB_HH

#include "CSDist.hh"

/*
Created By: Wesley Ford June 17, 2015

This class is responsible for the extraction of the energy dependant tabular out-going neutron cross-section distribution data from MCNP
and the writing of this data into the applicable G4NDL files.

To better understand the MCNP format that this class is built to extract from please refer to MCNP5 Manual Vol III
To better understand the G4NDL format that this class is built to write to, please refer to G4NDL Final State Decryption
*/

class CSDist1DTab: public CSDist
{
    public:
        CSDist1DTab(double *energyVec);
        CSDist1DTab(double *energyVec, int numCSEner);
        CSDist1DTab(CSDist1DTab *csDist)
        {
            startEner=csDist->startEner;
            CSVecSize=csDist->CSVecSize;
            enerVec=csDist->enerVec;
            CSVec= new double[CSVecSize];
            for(int i=0; i<CSVecSize; i++)
            {
                CSVec[i] = csDist->CSVec[i];
            }
        }

        virtual ~CSDist1DTab();

        void ExtractMCNPData(stringstream &stream, int &count);
        void WriteG4NDLCSData(stringstream &stream);
        void WriteG4NDLYieldData(stringstream &stream)
        {

        }
        void SetCSData(CSDist* nCSData)
        {

        }
        double Interpolate(double x);
        void SetCSData(double* &enerCSVecSet, int &csEnerStartSet, double* &csVecSet, int &csSizeSet)
        {
            enerCSVecSet=enerVec;
            csEnerStartSet=startEner;
            csVecSet=CSVec;
            csSizeSet=CSVecSize;
        }
        void AddData(CSDist *secDist)
        {
            secDist->AddData(startEner, CSVec, CSVecSize);
        }
        void AddData(int &csEnerStartSet, double* &csVecSet, int &csSizeSet)
        {
            int tempStartEner = min(startEner, csEnerStartSet);
            int tempCSVecSize = max(CSVecSize+startEner, csSizeSet+csEnerStartSet)-tempStartEner;
            double *tempCSVec = new double [tempCSVecSize];

            for(int i=0; i<tempCSVecSize; i++)
            {
                tempCSVec[i]=0;
                if((i+tempStartEner>=startEner)&&(i+tempStartEner<CSVecSize+startEner))
                {
                    tempCSVec[i]+=CSVec[i-startEner+tempStartEner];
                }
                else if(i+tempStartEner>=CSVecSize+startEner)
                {
                    tempCSVec[i]+=CSVec[CSVecSize+startEner-1];
                }
                else
                {
                    tempCSVec[i]+=max(0.,CSDist::Interpolate(1, enerVec[i+tempStartEner], enerVec[startEner], enerVec[startEner+1], CSVec[0], CSVec[1]));
                }
                if((i+tempStartEner>=csEnerStartSet)&&(i+tempStartEner<csSizeSet+csEnerStartSet))
                {
                    tempCSVec[i]+=csVecSet[i-csEnerStartSet+tempStartEner];
                }
                else if(i+tempStartEner>=csSizeSet+csEnerStartSet)
                {
                    tempCSVec[i]+=csVecSet[csSizeSet+csEnerStartSet-1];
                }
                else
                {
                    tempCSVec[i]+=max(0.,CSDist::Interpolate(1, enerVec[i+tempStartEner], enerVec[csEnerStartSet], enerVec[csEnerStartSet+1], csVecSet[0], csVecSet[1]));
                }
            }

            csEnerStartSet=tempStartEner;
            csSizeSet=tempCSVecSize;
            delete [] csVecSet;
            csVecSet=tempCSVec;

        }
        string IdentifyYourSelf();
        /*double* GetEnergyVec()
        {
            return (enerVec+startEner-1);
        }
        double* GetCSVec()
        {
            return CSVec;
        }
        int GetCSSize()
        {
            return CSVecSize;
        }*/

        int startEner, CSVecSize;
        double *CSVec, *enerVec;
        bool elastic=false;
    protected:
    private:
};

#endif // CSDIST1DTAB_HH
