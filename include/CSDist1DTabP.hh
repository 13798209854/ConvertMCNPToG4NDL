#ifndef CSDIST1DTABP_HH
#define CSDIST1DTABP_HH


class CSDist1DTabP
{
    public:
        CSDist1DTabP();
        virtual ~CSDist1DTabP();
        void ExtractMCNPData(stringstream data, int count&);
        void WriteG4NDLCSData(stringstream data);
        void WriteG4NDLYieldData(stringstream data);
        void IdentifyYourSelf();
        double Interpolate(double x);
        void SetCSData(CSDist* nCSData)
        {
            nCSData->SetNCSData(enerCSVec, csEnerStart, csVec, csSize);
        }
        int startEner, CSVecSize, csSize, mtNum, csEnerStart;
        double *CSVec, *enerVec, *csVec, *enerCSVec;
        bool elastic=false;
    protected:
    private:
};

#endif // CSDIST1DTABP_HH
