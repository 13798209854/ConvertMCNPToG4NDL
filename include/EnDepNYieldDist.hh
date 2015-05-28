#ifndef ENDEPNYIELDDIST_HH
#define ENDEPNYIELDDIST_HH


class EnDepNYieldDist
{
    public:
        EnDepNYieldDist();
        virtual ~EnDepNYieldDist();
        void ExtractMCNPData(stringstream stream, int &count);
        void WriteG4NDLData(stringstream data);

        int numRegs, numIncEner;
        int *regEndPos, *intScheme;
        double *incEner, *yield; // contains incoming neutron energy
    protected:
    private:
};

#endif // ENDEPNYIELDDIST_HH
