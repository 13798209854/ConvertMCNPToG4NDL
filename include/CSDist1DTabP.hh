#ifndef CSDIST1DTABP_HH
#define CSDIST1DTABP_HH


class CSDist1DTabP
{
    public:
        CSDist1DTabP(double *energyVec);
        virtual ~CSDist1DTabP();
        void ExtractMCNPData(stringstream data, int count&);
        void WriteG4NDLCSData(stringstream data);
        void WriteG4NDLYieldData(stringstream data);
        void IdentifyYourSelf();

        int startEner, CSVecSize;
        double *CSVec, *enerVec;
        bool elastic=false;
    protected:
    private:
};

#endif // CSDIST1DTABP_HH
