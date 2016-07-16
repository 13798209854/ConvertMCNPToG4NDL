#ifndef ISOTOPEMASS_HH
#define ISOTOPEMASS_HH

#include <string>
#include <iostream>

using namespace std;

class IsotopeMass
{
    public:
        IsotopeMass();
        virtual ~IsotopeMass();
        static void ClearStore();
        static void SetIsotopeMass();
        static double GetIsotopeMass(int Z,int A)
        {
            int iso=0;
            if(elemBaseA[Z]>0)
            {
                if((A==0)&&(Z>0))
                {
                    double sum=0.;
                    for(int i=0; i<elemNumIso[Z]; i++)
                    {
                        sum+=isotopeMass[Z][i]*isoNatAbun[Z][i];
                    }
                    if(sum==0.)
                    {
                        cout << "\nisotope Z=" << Z << " A=" << A << " is not in the mass list" << endl;
                    }
                    return sum;
                }
                else
                    iso = A-elemBaseA[Z];
            }
            if((iso<elemNumIso[Z])&&(iso>=0))
            {
                if(isotopeMass[Z][iso]==0.)
                {
                    cout << "\nisotope Z=" << Z << " A=" << A << " is not in mass list" << endl;
                }
                return isotopeMass[Z][iso];
            }
            else
            {
                cout << "\nError: isotope " << Z << " " << A << " is beyond the scope of this container" << endl;
                return 0.;
            }

        }
        static double EstimateIsotopeMass(int Z,int A)
        {
            int iso=0;
            if(elemBaseA[Z]>0)
            {
                if((A==0)&&(Z>0))
                {
                    double sum=0.;
                    for(int i=0; i<elemNumIso[Z]; i++)
                    {
                        sum+=isotopeMass[Z][i]*isoNatAbun[Z][i];
                    }
                    if(sum==0.)
                    {
                        if(elemNumIso[Z]>0)
                        {
                            return isotopeMass[Z][int(elemNumIso[Z]/2)];
                        }
                        else
                        {
                            cout << "\nError: isotope " << Z << " " << A << " is beyond the scope of this container" << endl;
                            return 0.;
                        }
                    }
                    return sum;
                }
                else
                    iso = A-elemBaseA[Z];
            }
            if((iso<elemNumIso[Z])&&(iso>=0))
            {
                if(isotopeMass[Z][iso]==0.)
                {
                    cout << "\nisotope Z=" << Z << " A=" << A << " is not in mass list" << endl;
                }
                return isotopeMass[Z][iso];
            }
            else
            {
                if(elemNumIso[Z]>1)
                {
                    if(A>=(elemBaseA[Z]+elemNumIso[Z]-1))
                    {
                        return (isotopeMass[Z][elemNumIso[Z]-1]-isotopeMass[Z][elemNumIso[Z]-2])*(A-elemBaseA[Z]-elemNumIso[Z]-2)+isotopeMass[Z][elemNumIso[Z]-2];
                    }
                    else if(A<=elemBaseA[Z])
                    {
                        return (isotopeMass[Z][1]-isotopeMass[Z][0])*(A-elemBaseA[Z])+isotopeMass[Z][elemNumIso[Z]-2];
                    }
                }
                else if(elemNumIso[Z]==1)
                {
                    return isotopeMass[Z][0];
                }
                else
                {
                    cout << "\nError: isotope " << Z << " " << A << " is beyond the scope of this container" << endl;
                    return 0.;
                }
                cout << "\nError: isotope " << Z << " " << A << " is beyond the scope of this container" << endl;
                return 0.;
            }

        }
        static double GetIsotopeMassEnergy(int Z,int A)
        {
            return GetIsotopeMass(Z,A)*931.494061;
        }
        static double EstimateIsotopeMassEnergy(int Z,int A)
        {
            return EstimateIsotopeMass(Z,A)*931.494061;
        }
        static void GetNaturalAbundanceVec(int Z, double* &natAbun, int &vecSize, int &baseA)
        {
            natAbun=isoNatAbun[Z];
            vecSize=elemNumIso[Z];
            baseA=elemBaseA[Z];
        }
        static void GetMostNatAbunIso(int Z, int &A)
        {
            double maxAbun=0.;
            int index=0;
            for(int i=0; i<elemNumIso[Z]; i++)
            {
                if(maxAbun<isoNatAbun[Z][i])
                {
                    maxAbun=isoNatAbun[Z][i];
                    index=i;
                }
            }
            A=elemBaseA[Z]+index;
        }
        static double **isotopeMass;
        static double **isoNatAbun;
        static int *elemNumIso;
        static int *elemBaseA;
    protected:
    private:

};
#endif // ISOTOPEMASS_HH
