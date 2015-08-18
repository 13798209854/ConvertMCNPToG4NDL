using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "zlib.h"
#include <dirent.h>
#include "include/ElementNames.hh"
#include <iomanip>

//classes responsible for the extraction of the energy independant out-going angular distribution data
#include "include/AngularDist.hh"
#include "include/AngDist2DTabular.hh"
#include "include/AngDist32EqualPBin.hh"
#include "include/AngDistIso.hh"
#include "include/AngDist2DTabularP.hh"
#include "include/AngDistIsoP.hh"
#include "include/AngDist32EqPBinP.hh"

//classes responsible for the extraction of the energy dependant out-going angular distribution data
#include "include/AngularEnergyDist.hh"
#include "include/AngEnDist3DTab.hh"
#include "include/AngEnDistKallbach.hh"
#include "include/AngEnDistLab3DTab.hh"
#include "include/AngEnDistNBody.hh"

//classes responsible for the extraction of the energy dependant cross-section data
#include "include/CSDist.hh"
#include "include/CSDist1DTab.hh"
#include "include/CSDist1DTabP.hh"
#include "include/CSDistYieldComp.hh"

//classes responsible for the extraction of the out-going energy distribution data
#include "include/EnergyDist.hh"
#include "include/EnerDist1PhEner.hh"
#include "include/EnerDistConTab.hh"
#include "include/EnerDistEqPEnerBins.hh"
#include "include/EnerDistEvapSpec.hh"
#include "include/EnerDistGenEvapSpec.hh"
#include "include/EnerDistLevScat.hh"
#include "include/EnerDistMaxwellFisSpec.hh"
#include "include/EnerDistTabLinFunc.hh"
#include "include/EnerDistTabMulti.hh"
#include "include/EnerDistWattSpec.hh"

//class responsible for the extraction of the delayed neutron average lifetime distribution data
#include "include/NDelayConstDist.hh"

//classes responsible for the extraction of the out-going neutron yield distribution data
#include "include/YieldDist.hh"
#include "include/NYield1DTab.hh"
#include "include/NYieldPolyFunc.hh"

//total number of processes
#define numProcess 78
//total number of processes not including those that make up the MT4 reaction
#define numProcess2 38

int CreateIsoCSData(stringstream &stream, string outDirName, bool ascii, double setTemperature, bool limitTemp, bool onlyCS);

bool DirectoryExists( const char* pzPath );

void MakeReactionListFile(string outDirName, string isoName, int MTRList[], int MTRListPos[]);

void MakeCSDataFile(string outDirName, string isoName, int MTRNum[], CSDist **csDist, bool ascii);

void MakeElasFSFile(string fileName, string isoName, double isoMass, int numAngEner[], vector<AngularDist*> *angDist, bool ascii);

void MakeCaptureFSFile(string outDirName, string isoName, double isoMass, double temperature, int TYRList[], CSDist **pCSVec, vector<AngularDist*> *angDistP,
                        vector<int> *enDisLawP, vector<EnergyDist*> *enerDistP, vector<int> MTRPList, int *extractOrderP, int *numAngEnerP,
                        vector<int> *enDisNumLawApplNRegP, vector<int> *enDisNumLawApplNEnP, vector<int*> *enDisSchemeVecP, vector<int*> *enDisRangeVecP,
                        vector<double*> *enDisEnApplVecP, vector<double*> *enDisProbApplVecP, double **energyAngVecP, CSDist *nCSVec, bool ascii);

void CreateMT2(int* MTRListPos, string outDirName, string isoName, int isoNum, double isoMass, double temperature, int *MTRList, CSDist* nCSVec[],
                        vector<int>* enDisLaw, vector<int>* enDisNumLawApplNReg, vector<int>* enDisNumLawApplNEn, vector<int*>* enDisSchemeVec,
                        vector<int*>* enDisRangeVec, vector<double*>* enDisEnApplVec, vector<double*>* enDisProbApplVec, vector<EnergyDist*>* enerDist,
                        YieldDist* nYieldReac[], double* reacQValue, int* numAngEner, vector<AngularDist*>* angDist, bool* angDistInEnDistFlag,
                        vector<AngularEnergyDist*>* angEnDist, CSDist*** pCSVec, int *TYRList, int** numAngEnerP, vector<AngularDist*>** angDistP, vector<int>** enDisLawP,
                        vector<int>** enDisNumLawApplNRegP, vector<int>** enDisNumLawApplNEnP, vector<int*>** enDisSchemeVecP, vector<int*>** enDisRangeVecP,
                        vector<double*>** enDisEnApplVecP, vector<double*>** enDisProbApplVecP,vector<EnergyDist*>** enerDistP, vector<int> &MTRPList,
                        double*** energyAngVecP, bool ascii);

void MakeFissionFSFile(int *MTRList, int *MTRListPos, string outDirName, string isoName, double isoMass, double temperature, CSDist **nCSVec, vector<int> *enDisLaw, vector<int> *enDisNumLawApplNReg, vector<int> *enDisNumLawApplNEn,
                        vector<int*> *enDisSchemeVec, vector<int*> *enDisRangeVec, vector<double*> *enDisEnApplVec, vector<double*> *enDisProbApplVec, vector<EnergyDist*> *enerDist,  vector<AngularEnergyDist*> *angEnDist,
                        vector<int> enDisLawND, vector<int> enDisNumLawApplNRegND, vector<int> enDisNumLawApplNEnND, vector<int*> enDisSchemeVecND, vector<int*> enDisRangeVecND,
                        vector<double*> enDisEnApplVecND, vector<double*> enDisProbApplVecND, vector<EnergyDist*> enerDistND, double *reacQValue, int *numAngEner, vector<AngularDist*> *angDist,
                        bool *angDistInEnDistFlag, CSDist **pCSVec, int *TYRList, bool promptYieldFlag, bool totalYieldFlag, YieldDist *dNPromptYieldDist, YieldDist *dNTotalYieldDist,
                        YieldDist *dNYield, NDelayConstDist *nDelConst, int *numAngEnerP, vector<AngularDist*> *angDistP, vector<int> *enDisLawP, vector<int> *enDisNumLawApplNRegP, vector<int> *enDisNumLawApplNEnP,
                        vector<int*> *enDisSchemeVecP, vector<int*> *enDisRangeVecP, vector<double*> *enDisEnApplVecP, vector<double*> *enDisProbApplVecP, vector<EnergyDist*> *enerDistP,
                        vector<int> MTRPList, double **energyAngVecP, bool ascii);

void CreateMT4(int* MTRListPos, string outDirName, string isoName, int isoNum, double isoMass, double temperature, int *MTRList, CSDist* nCSVec[],
                        vector<int>* enDisLaw, vector<int>* enDisNumLawApplNReg, vector<int>* enDisNumLawApplNEn, vector<int*>* enDisSchemeVec,
                        vector<int*>* enDisRangeVec, vector<double*>* enDisEnApplVec, vector<double*>* enDisProbApplVec, vector<EnergyDist*>* enerDist,
                        YieldDist* nYieldReac[], double* reacQValue, int* numAngEner, vector<AngularDist*>* angDist, bool* angDistInEnDistFlag,
                        vector<AngularEnergyDist*>* angEnDist, CSDist*** pCSVec, int *TYRList, int** numAngEnerP, vector<AngularDist*>** angDistP, vector<int>** enDisLawP,
                        vector<int>** enDisNumLawApplNRegP, vector<int>** enDisNumLawApplNEnP, vector<int*>** enDisSchemeVecP, vector<int*>** enDisRangeVecP,
                        vector<double*>** enDisEnApplVecP, vector<double*>** enDisProbApplVecP,vector<EnergyDist*>** enerDistP, vector<int> &MTRPList,
                        double*** energyAngVecP, bool ascii);

void MakeInElasticFSFile(int *MTRListPos, string outDirName, string isoName, int isoNum, double isoMass, double temperature, int *MTRList, CSDist *nCSVec[],
                        vector<int> *enDisLaw, vector<int> *enDisNumLawApplNReg, vector<int> *enDisNumLawApplNEn, vector<int*> *enDisSchemeVec,
                        vector<int*> *enDisRangeVec, vector<double*> *enDisEnApplVec, vector<double*> *enDisProbApplVec, vector<EnergyDist*> *enerDist,
                        YieldDist *nYieldReac[], double *reacQValue, int *numAngEner, vector<AngularDist*> *angDist, bool *angDistInEnDistFlag,
                        vector<AngularEnergyDist*> *angEnDist, CSDist **pCSVec, int *TYRList, int *numAngEnerP, vector<AngularDist*> *angDistP, vector<int> *enDisLawP,
                        vector<int> *enDisNumLawApplNRegP, vector<int> *enDisNumLawApplNEnP, vector<int*> *enDisSchemeVecP, vector<int*> *enDisRangeVecP,
                        vector<double*> *enDisEnApplVecP, vector<double*> *enDisProbApplVecP,vector<EnergyDist*> *enerDistP, vector<int> MTRPList,
                        double **energyAngVecP, bool ascii);

double Interpolate(int aScheme, double x, double x1, double x2, double y1, double y2);
double Histogram(double , double , double , double y1, double );
double LinearLinear(double x, double x1, double x2, double y1, double y2);
double LinearLogarithmic(double x, double x1, double x2, double y1, double y2);
double LogarithmicLinear(double x, double x1, double x2, double y1, double y2);
double LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2);
double Random(double , double , double , double y1, double y2);
void GetDataStream( string, std::stringstream&);
void SetDataStream( string, std::stringstream&, bool, bool overWrite=true);


// X Note: range values (NBT) seem to be absolute positions in the array, test to make sure this is true
// X create code interpreter using the functions from ExtractMaterialComposition to make sure that every dynamic variable is deleted
// X the energy in MCNP seems to be MeV by default, check that this is true and convert from MeV to eV
// X all distinct photons share the same energy distribution for the same incoming neutron reaction, whereas in MCNP each photon reaction has its own energy distribution
// X we can get around this by combining all of the energy distributions for every photon reaction into one big list with the the probability of each energy dist weighed
//  by its photon reactions cross-section or we can replace this section of MCNP data with the original G4NDL data
// X we don't use the out-going photons in our simulations thus far any ways so it doesn't matter too much if this is accurate
// X the chronological ordering of the process extraction is probably unnessary since they seem to already be in order
// X Check photon angular distribution
// X give the option to only extract CS data
// X we will have to check to make sure that the distributions are being used in the right energy regime, they may have to be listed in the G4NDL file in order of increasing energy
// X Maybe an error in the extraction of inelastic, only a few of the isotopes seem to have the identifier it is looking for
// X in AngEnDit3D classes check to make sure that the locators are being used properly, we assumed a -1 add on based off other examples
// X go through the Extract functions and check that they match the format of MCNP
// X go through the write functions and check that they match the format of G4NDL
// X make sure that the count variable is being properly kept up to date
// X Check for dynamic memory issues
// X make sure that region positions are set properly (they should be the vector index + 1)
// X add comments to major sections of the code
// X Format the output so that it is properly spaced
// X add warnings to the code
// X Check NU block extraction, unclear what XSS(JXS(2)) means versus JXS(2)
// X Since there is no inelastic scattering for 1001 MCNP doesn't bother giving a cross-section file which causes G4STORK to not use the MCNP data for 1001,
//  we fixed this issue by forcing MCNP to output a cross-section file for the reaction which sets the probability of the reaction to zero for all incoming neutron energies

// Add /Inelastic/Gammas/
// Compare liams converted data to ours
// We use the average CS when we weight reactions for mixing energy/ angular distributions. a more exact way would be completely mix the cross-section data with the dist probability data, preserving
//  the accuracy of the data by using an energy vector that is a combination of both schemes, this would heavily slow down the simulation for probably little accuracy gain
// We currently interpolate the CS data when we weigh the yied data, this causes errors because the boundary cs data is always zero so if there is only two yield points they will always be set to zero
//  to fix this error we need to add combione the energy vector of the CS data and the yied data
// Check photon enegy dist extraction in G4NDL and make sure that our approximation of mixing them is right
// Check the delayed neutron energy distribution weight vector and see if we need to multiply it by the delayed group weight
// Inelastic/F36/ is not created because the reaction does not appear in the MCNP data, check to make sure this is correct
// parrallelize code

//Takes in a directory of MCNP cross-section libraries, converts the data into the G4NDL format and then outputs the information in a given directory

int main(int argc, char **argv)
{
    ElementNames elementNames;
    elementNames.SetElementNames();

    string outputType="ascii";
    bool ascii=true, onlyCS=false, limitTemp=false;
    string word, libName="endf";
    char lib, version='7', check1, check2, check3;
    string inFileName, outDirName, fileName;
    int result=0;
    double temperature;

    stringstream stream;

    //Extracts user Inputs
    if(argc==7)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5] << ' ' << argv[6];
        stream >> inFileName >> outDirName >> outputType >> version >> onlyCS >> temperature;
        limitTemp=true;
    }
    else if(argc==6)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5];
        stream >> inFileName >> outDirName >> outputType >> version >> onlyCS;
    }
    else if(argc==5)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
        stream >> inFileName >> outDirName >> outputType >> version;
    }
    else if(argc==4)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3];
        stream >> inFileName >> outDirName >> outputType;
    }
    else if(argc==3)
    {
        stream << argv[1] << ' ' << argv[2];
        stream >> inFileName >> outDirName;
    }
    else
    {
        cout << "Incorrect number of inputs; give the MCNP file to be converted and the output directory for the created G4NDL files" << endl;
        elementNames.ClearStore();
        return 1;
    }

    if(outputType == "compressed"||outputType == "compress"||outputType == "Compressed"||outputType == "Compress"||outputType == "Zipped"||outputType == "Zip"||outputType == "zipped"||outputType == "zip" )
        ascii=false;

    stream.clear();
    stream.str("");

    if(inFileName.back()=='/')
    {
        libName.push_back(version);
        DIR *dir;
        struct dirent *ent;

        //goes through the given directory and converts the ENDF libraries that match the given vversion
        if ((dir = opendir (inFileName.c_str())) != NULL)
        {
            while ((ent = readdir (dir)) != NULL)
            {
                if(string(ent->d_name).substr(0, 5)==libName)
                {
                    fileName=inFileName+ent->d_name;
                    cout << "Opening file: " << fileName << endl;
                    // Gets data from the file and stores it into a data stream
                    GetDataStream(fileName, stream);

                    stream >> word;
                    while(stream)
                    {
                        check1=word[int(word.find_last_of('.')+1)];
                        check2=word[int(word.find_last_of('.')+2)];
                        check3=word[int(word.find_last_of('.')+3)];

                        //checks whether the word matches the beggining of an isotope data set identifier
                        if((check1==version)&&((check2>='0')&&(check2<='9'))/*&&((check3=='c')||(check3=='d'))*/)
                        {
                            if((check3=='c')||(check3=='d'))
                            {
                                // gets the elastic, inelastic, fission and capture CS data for the isotope
                                result += CreateIsoCSData(stream, outDirName, ascii, temperature, limitTemp, onlyCS);
                            }
                        }
                        stream >> word;
                    }
                    stream.str("");
                    stream.clear();
                }
            }
            closedir(dir);
        }
    }
    else
    {
        // Gets data from the file and stores it into a data stream
        cout << "Opening file: " << inFileName << endl;
        GetDataStream(inFileName, stream);
        lib = inFileName[int(inFileName.length()-3)];

        stream >> word;
        while(stream)
        {
            check1=word[int(word.find_last_of('.')+1)];
            check2=word[int(word.find_last_of('.')+2)];
            check3=word[int(word.find_last_of('.')+3)];

            //checks whether the word matches the beggining of an isotope data set identifier
            if((check1==lib)&&((check2>='0')&&(check2<='9'))/*&&((check3=='c')||(check3=='d'))*/)
            {
                if((check3=='c')||(check3=='d'))
                {
                    // gets the elastic, inelastic, fission and capture CS data for the isotope
                    result += CreateIsoCSData(stream, outDirName, ascii, temperature, limitTemp, onlyCS);
                }
            }
            stream >> word;
        }
    }

    system( ("chmod -R  777 "+outDirName).c_str());
    elementNames.ClearStore();

    if(result>1)
        result=1;
    return result;
}

//Extraxts the CS data from the MCNP file and stroes it in a G4NDL formatted file
int CreateIsoCSData(stringstream &stream, string outDirName, bool ascii, double setTemperature, bool limitTemp, bool onlyCS)
{
    //file name data
    ElementNames *elementNames;
    string temperature, isoName;
    string dummy, Z, A;
    stringstream numConv;
    int Znum, isoNum;
    double isoMass, dataTemperature;
    int yieldDistType;
    bool promptYieldFlag=false, totalYieldFlag=false, createMT4Flag=false, createMT2Flag=false;

    //fixed location data
    int numCSEner=0, numReactions=0, numReacSecN=0, numPReactions=0, numDNPrecursorFam=0;
    int startEnerTable, startElasticBlock, startMTRBlock;
    int startNUBlock, startQValBlock, startTYRBlock;
    int startLSIGBlock, startCSBlock;
    int startLANDBlock, startANDBlock;
    int startLDLWBlock, startDLWBlock;
    int startMTRPBlock, startLSIGPBlock, startCSPBlock;
    int startLANDPBlock, startANDPBlock;
    int startLDLWPBlock, startDLWPBlock;
    int startYPBlock;
    int startDelayedNUBlock, startDNEDLBlock, startDNEDBlock;
    int startFisBlock;
    int lastMTRPos=-1;

    //temperary data
    double temp;
    int intTemp;
    int index, count=1;
    char line[256];

    //fixed data arrays
    double *energyCSVec=NULL; // incoming neutron kinetic energy for CS data

    //for the following arrays, the indicies 0,1, and 2 correspond to elastic, capture, and total fission respectively
    //location arrays are set to -1 to show that the data has not been set
    // this array sets the order of the processes, must be manually enetered
    int MTRList[numProcess]={2,102,18,19,20,21,38,4,5,11,16,17,22,23,24,25,28,29,32,33,34,37,41,42,44,45,103,104,105,106,107,108,111,112,113,115,116,117,
                                51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,78,79,80,81,82,83,84,85,86,87,88,89,90,91};
    int MTRListPos[numProcess]; //position of the reaction MT in MT list
    for(int i=0; i<numProcess; i++)
    {
        MTRListPos[i]=-1;
    }
    double reacQValue[numProcess]; // Q-Value of the reaction
    for(int i=0; i<numProcess; i++)
    {
        reacQValue[i]=0.;
    }
    int TYRList[numProcess]; // the average number of neutrons released from the reaction,
    // the sign determines the the reference frame that the angular and energy dist where measured from, - means CM and + means Lab
    for(int i=0; i<numProcess; i++)
    {
        TYRList[i]=0;
    }
    int LSIGList[numProcess]; //starting point of the reaction cross-section data in the data file
    for(int i=0; i<numProcess; i++)
    {
        LSIGList[i]=-1;
    }
    CSDist* nCSVec[numProcess]; //cross section array for each reaction
    for(int i=0; i<numProcess; i++)
    {
        nCSVec[i]=NULL;
    }
    int LANDList[numProcess]; //starting point of the reaction's out going neutron angular distribution data
    for(int i=0; i<numProcess; i++)
    {
        LANDList[i]=-1;
    }
    int numAngEner[numProcess]; // the number of incoming neutron energy points for the angular distribution
    for(int i=0; i<numProcess; i++)
    {
        numAngEner[i]=-1;
    }
    int LDLWList[numProcess]; //starting point of the reaction's out going neutron energy distribution data
    for(int i=0; i<numProcess; i++)
    {
        LDLWList[i]=-1;
    }
    int extractOrder[numProcess]; // the order in which the data for each reaction will be extracted, pos [0] will be the first to be extracted,
    //where the number it contains corresponds to the index in the above arrays that the reaction belongs to
    for(int i=0; i<numProcess; i++)
    {
        extractOrder[i]=i;
    }

    int MTROrder[numProcess]; // specifies the order in which the MTR reactions occur in the dataset, each element of MTROrder is an index of the MTRList array,
    // pointing to a specific MTR reaction, MTROrder[0] points to the first MTR reaction to be extracted
    for(int i=0; i<numProcess; i++)
    {
        MTROrder[i]=-1;
    }

    // first 12 lines contain miscellaneous data and the NXS and JXS arrays page 118 of the MCNP manual
    //### in this section we extract the isotope naming, temperature and mass data

    stream >> isoMass >> dataTemperature;
    dataTemperature=dataTemperature*1000000/(8.6173324*(pow(10, -5)));
    numConv << dataTemperature;
    numConv >> temperature;

    if(limitTemp&&(dataTemperature!=setTemperature))
    {
        return 0;
    }

    numConv.clear();
    numConv.str("");

    if(outDirName.back()!='/')
        outDirName+='/';

    outDirName = outDirName+temperature+"K/";

    for(int i=0; i<6; i++)
    {
        stream.getline(line,256);
    }
    // Extracts the isotope name
    //start of NXS block
    stream >> dummy >> isoNum;
    Znum = int(isoNum/1000);
    numConv.clear();
    numConv.str("");
    numConv << Znum;
    numConv >> Z;
    numConv.clear();
    numConv.str("");
    numConv << (isoNum-Znum*1000);
    numConv >> A;
    numConv.clear();
    numConv.str("");

    if(A=="0")
        A="nat";

    isoName = Z+'_'+A+'_'+elementNames->GetName(Znum);

    //### in this section we extract the data block locations for the following sections

    stream >> numCSEner >> numReactions >> numReacSecN >> numPReactions >> dummy >> numDNPrecursorFam;

    energyCSVec = new double[numCSEner];

    for(int i=0; i<8; i++)
    {
        stream >> dummy;
    }

    // this is the extraction of the JSX array
    stream >> startEnerTable >> startNUBlock >> startMTRBlock >> startQValBlock >> startTYRBlock >> startLSIGBlock >> startCSBlock >> startLANDBlock >>
                startANDBlock >> startLDLWBlock >> startDLWBlock >> dummy >> startMTRPBlock >> startLSIGPBlock >> startCSPBlock >> startLANDPBlock >>
                startANDPBlock >> startLDLWPBlock >> startDLWPBlock >> startYPBlock >> startFisBlock >> dummy >> dummy >> startDelayedNUBlock >> dummy
                >> startDNEDLBlock >> startDNEDBlock;

    startElasticBlock = startEnerTable+3*numCSEner;

    for(int i=0; i<5; i++)
    {
        stream >> dummy;
    }

    //### in this section we extract incoming neutron energy, cross-section, neutron yield, and reaction Q-value data

    //Extract the ESZ Block to get the main incoming neutron energy grid and elastic cross-section
    for(;count<startEnerTable; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numCSEner; i++, count++)
    {
        stream >> temp;
        energyCSVec[i]=temp;
    }

    for(;count<startElasticBlock; count++)
    {
        stream >> dummy;
    }

    nCSVec[0] = new CSDist1DTab(energyCSVec, numCSEner);
    nCSVec[0]->ExtractMCNPData(stream, count);

    //Extract NU block here for fission yield
    //have not iterated through this section yet
    YieldDist *dNPromptYieldDist=NULL;
    YieldDist *dNTotalYieldDist=NULL;
    int startTotalYield=0;
    if(startNUBlock!=0)
    {
        if((count+numCSEner)!=startNUBlock)
        {
            cout << "Error: counter is off at the start of the NU Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be equal to " << (startNUBlock-numCSEner) << endl;
        }
        for(;count<startNUBlock; count++)
        {
            stream >> dummy;
        }
        stream >> startTotalYield; count++;
        if(startTotalYield>0)
        {
            yieldDistType=startTotalYield;
            //###Check this: if the delayed neutron yield exists, assume the NU block contains the prompt yield
            (startDelayedNUBlock!=0) ? promptYieldFlag=true : totalYieldFlag=true;

            if(promptYieldFlag)
            {
                if(yieldDistType==1)
                {
                    dNPromptYieldDist = new NYieldPolyFunc;
                }
                else
                {
                    dNPromptYieldDist = new NYield1DTab;
                }

                dNPromptYieldDist->ExtractMCNPData(stream, count);
            }
            else
            {
                if(yieldDistType==1)
                {
                    dNTotalYieldDist = new NYieldPolyFunc;
                }
                else
                {
                    dNTotalYieldDist = new NYield1DTab;
                }

                dNTotalYieldDist->ExtractMCNPData(stream, count);
            }

        }
        else if(startTotalYield<0)
        {
            promptYieldFlag=true;

            stream >> yieldDistType; count++;
            if(yieldDistType==1)
            {
                dNPromptYieldDist = new NYieldPolyFunc;
            }
            else
            {
                dNPromptYieldDist = new NYield1DTab;
            }
            dNPromptYieldDist->ExtractMCNPData(stream, count);

            if(count!=(startNUBlock-startTotalYield+1))
            {
                cout << "Error: counter is off at the location of total neutron yield data for isotope " << isoName << endl;
                cout << "Count= " << count << " should be equal to " << (-2*startNUBlock+1) << endl;
            }
            for(;count<(startNUBlock-startTotalYield+1); count++)
            {
                stream >> dummy;
            }

            totalYieldFlag=true;
            stream >> yieldDistType; count++;
            if(yieldDistType==1)
            {
                dNTotalYieldDist = new NYieldPolyFunc;
            }
            else
            {
                dNTotalYieldDist = new NYield1DTab;
            }
            dNTotalYieldDist->ExtractMCNPData(stream, count);
        }
    }

    //Extract the MTR block to determine reaction ordering
    if(numReactions!=0)
    {
        if((count!=startMTRBlock)&&(count+numCSEner!=startMTRBlock))
        {
            cout << "Error: counter is off at the start of the MTR Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be equal to " << (startMTRBlock) << " or " << (startMTRBlock-numCSEner) << endl;
        }
        for(;count<startMTRBlock; count++)
        {
            stream >> dummy;
        }

        intTemp=1;
        MTROrder[0]=0;
        for(int i=0; i<numReactions; i++, count++)
        {
            stream >> index;
            //make sure these MT numbers are right
            for(int j=1; j<numProcess; j++)
            {
                if(MTRList[j]==index)
                {
                    MTRListPos[j]=i;
                    lastMTRPos=i;
                    MTROrder[intTemp]=j;
                    intTemp++;
                    if(numProcess==intTemp)
                        break;
                }
            }
        }

        MakeReactionListFile(outDirName,isoName,MTRList,MTRListPos);

        //Extract the QVal block to get the reaction Q values
        if(count!=startQValBlock)
        {
            cout << "Error: counter is off at the start of the QVal Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be equal to " << startQValBlock << endl;
        }
        for(;count<startQValBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=1; i<numProcess; i++)
        {
            if(MTROrder[i]>0)
            {
                for(;count<MTRListPos[MTROrder[i]]+startQValBlock; count++)
                {
                    stream >> dummy;
                }
                stream >> temp; count++;
                reacQValue[MTROrder[i]]=temp;
            }
        }

        //Extract TYR block here to get the reaction yield data and reference frame
        if((count+numReactions-lastMTRPos-1)!=startTYRBlock)
        {
            cout << "Error: counter is off at the start of the TYR Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be equal to " << (startTYRBlock-(numReactions-lastMTRPos-1)) << endl;
        }
        for(;count<startTYRBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=1; i<numProcess; i++)
        {
            if(MTROrder[i]>0)
            {
                for(;count<MTRListPos[MTROrder[i]]+startTYRBlock; count++)
                {
                    stream >> dummy;
                }
                stream >> index;  count++;
                TYRList[MTROrder[i]]=index;
            }
        }

        //Extract the LSIG and the SIG blocks to get the cross-section data for reactions other than elastic scattering
        if((count+numReactions-lastMTRPos-1)!=startLSIGBlock)
        {
            cout << "Error: counter is off at the start of the LSIG Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be equal to " << (startLSIGBlock-(numReactions-lastMTRPos-1)) << endl;
        }
        for(;count<startLSIGBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=1; i<numProcess; i++)
        {
            if(MTROrder[i]>0)
            {
                for(;count<MTRListPos[MTROrder[i]]+startLSIGBlock; count++)
                {
                    stream >> dummy;
                }
                stream >> index; count++;
                LSIGList[MTROrder[i]]=index;
            }
        }

        for(int i=1; i<numProcess; i++)
        {
            for(int j=i+1; j<numProcess; j++)
            {
                if(LSIGList[extractOrder[j]]<LSIGList[extractOrder[i]])
                {
                    index = extractOrder[i];
                    extractOrder[i] = extractOrder[j];
                    extractOrder[j] = index;
                }
            }
        }

        for(int i=1; i<numProcess; i++)
        {
            if(LSIGList[extractOrder[i]]!=-1)
            {
                if(count>(LSIGList[extractOrder[i]]+startCSBlock-1))
                {
                    cout << "Error: counter is off at the location of cross-section data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << (LSIGList[extractOrder[i]]+startCSBlock-1) << endl;
                }
                for(;count<(LSIGList[extractOrder[i]]+startCSBlock-1); count++)
                {
                    stream >> dummy;
                }
                nCSVec[extractOrder[i]] = new CSDist1DTab(energyCSVec);
                nCSVec[extractOrder[i]]->ExtractMCNPData(stream, count);
            }
        }
        MakeCSDataFile(outDirName, isoName, MTRList, nCSVec, ascii);
    }

    if(onlyCS)
    {
        for(int i=0; i<numProcess; i++)
        {
            if(nCSVec[i])
                delete nCSVec[i];
        }
        if(dNTotalYieldDist)
            delete dNTotalYieldDist;

        if(dNPromptYieldDist)
            delete dNPromptYieldDist;

        if(energyCSVec)
            delete [] energyCSVec;

        return 0;
    }

    //### In this Section we extract the out-going neutron angular distribution data that is independant
    // of energy the rest of the angular data is extracted in the energy distribution

    //Extract the LAND and AND blocks to get the out-going neutron angular distribution data
    if(count>startLANDBlock)
    {
        cout << "Error: counter is off at the start of the LAND Block for isotope " << isoName << endl;
        cout << "Count= " << count << " should be less than equal to " << startLANDBlock << endl;
    }
    for(;count<startLANDBlock; count++)
    {
        stream >> dummy;
    }

    //elastic
    stream >> index; count++;
    LANDList[0]=index;

    for(int i=1; i<numProcess; i++)
    {
        if((MTROrder[i]>0)&&(TYRList[MTROrder[i]]!=0))
        {
            for(;count<MTRListPos[MTROrder[i]]+startLANDBlock+1; count++)
            {
                stream >> dummy;
            }
            stream >> index; count++;
            LANDList[MTROrder[i]]=index;
        }
    }

    for(int i=0; i<numProcess; i++)
    {
        for(int j=i+1; j<numProcess; j++)
        {
            if(LANDList[extractOrder[j]]<LANDList[extractOrder[i]])
            {
                index = extractOrder[i];
                extractOrder[i] = extractOrder[j];
                extractOrder[j] = index;
            }
        }
    }

    double *energyAngVec[numProcess]; // incoming neutron kinetic energy for Ang dist data
    for(int i=0; i<numProcess; i++)
    {
        energyAngVec[i]=NULL;
    }
    int *angTabPosVec[numProcess]; // the location of the angular probability table associated to the incoming neutron kinetic energy
    for(int i=0; i<numProcess; i++)
    {
        angTabPosVec[i]=NULL;
    }
    int *exOrderAng=NULL;
    int distType;
    bool angDistInEnDistFlag[numProcess];
    for(int i=0; i<numProcess; i++)
    {
        angDistInEnDistFlag[i]=false;
    }
    vector<AngularDist*> angDist[numProcess];

    for(int i=0; i<numProcess; i++)
    {
        // an isotropic distribution is assumed
        if(LANDList[extractOrder[i]]==0)
        {
            angDist[extractOrder[i]].push_back(new AngDistIso());
            angDist[extractOrder[i]].back()->SetTemperature(dataTemperature);
            numAngEner[extractOrder[i]]=2;
            continue;
        }

        //the angular distribution is contained in the energy distribution data
        if(LANDList[extractOrder[i]]==-1)
        {
            angDistInEnDistFlag[extractOrder[i]]=true;
            continue;
        }

        if(count>(LANDList[extractOrder[i]]+startANDBlock-1))
        {
            cout << "Error: counter is off at the location of primary neutron angular distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << (LANDList[extractOrder[i]]+startANDBlock-1) << endl;
        }
        for(;count<(LANDList[extractOrder[i]]+startANDBlock-1); count++)
        {
            stream >> dummy;
        }
        stream >> intTemp; count++;
        numAngEner[extractOrder[i]]=intTemp;

        energyAngVec[extractOrder[i]] = new double[numAngEner[extractOrder[i]]];
        for(int j=0; j<numAngEner[extractOrder[i]]; j++, count++)
        {
            stream >> temp;
            energyAngVec[extractOrder[i]][j] = temp;
        }

        angTabPosVec[extractOrder[i]] = new int[numAngEner[extractOrder[i]]];

        for(int j=0; j<numAngEner[extractOrder[i]]; j++, count++)
        {
            stream >> temp;
            angTabPosVec[extractOrder[i]][j] = temp;
        }

        exOrderAng = new int[numAngEner[extractOrder[i]]];
        for(int j=0; j<numAngEner[extractOrder[i]]; j++)
        {
            exOrderAng[j]=j;
        }
        //## we may be creating an error here by mixing up the order of the distributions with respect to energy
        for(int j=0; j<numAngEner[extractOrder[i]]; j++)
        {
            for(int k=j+1; k<numAngEner[extractOrder[i]]; k++)
            {
                if(abs(angTabPosVec[extractOrder[i]][exOrderAng[k]])<abs(angTabPosVec[extractOrder[i]][exOrderAng[j]]))
                {
                    index = exOrderAng[j];
                    exOrderAng[j] = exOrderAng[k];
                    exOrderAng[k] = index;
                }
            }
        }

        distType=0;
        for(int j=0; j<numAngEner[extractOrder[i]]; j++)
        {
            if(angTabPosVec[extractOrder[i]][exOrderAng[j]]==0)
            {
                //the distribution is isotropic and no further data is needed
                if(distType!=1)
                {
                    angDist[extractOrder[i]].push_back(new AngDistIso());
                    angDist[extractOrder[i]].back()->SetTemperature(dataTemperature);
                    distType=1;
                }
            }
            else if(angTabPosVec[extractOrder[i]][exOrderAng[j]]<0)
            {
                //the angular distribution is represented as a table of cosine and prob
                if(count>(startANDBlock-angTabPosVec[extractOrder[i]][exOrderAng[j]]-1))
                {
                    cout << "Error: counter is off at the location of prompt neutron angular distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << startANDBlock-angTabPosVec[extractOrder[i]][exOrderAng[j]]-1 << endl;
                }
                for(;count<(startANDBlock-angTabPosVec[extractOrder[i]][exOrderAng[j]]-1); count++)
                {
                    stream >> dummy;
                }

                if(distType!=2)
                {
                    angDist[extractOrder[i]].push_back(new AngDist2DTabular());
                    angDist[extractOrder[i]].back()->SetTemperature(dataTemperature);
                    distType=2;
                }
            }
            else
            {
                //the angular distribution is represented by a probability density function where each of the 32 points/bins are spaced along the entire cosine of the scattering angle axis
                //such that the integral of the probability density with the cosine of the scattering angle inbetween them is kept constant, meaning that the bins are in a equilprobable arrangement
                // the first bin probably corresponds to a cosine value of -1 and and the last bin being some where around 1 depending on the spacing
                if(count>(startANDBlock+angTabPosVec[extractOrder[i]][exOrderAng[j]]-1))
                {
                    cout << "Error: counter is off at the location of prompt neutron angular distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << startANDBlock+angTabPosVec[extractOrder[i]][exOrderAng[j]]-1 << endl;
                }
                for(;count<(startANDBlock+angTabPosVec[extractOrder[i]][exOrderAng[j]]-1); count++)
                {
                    stream >> dummy;
                }

                if(distType!=3)
                {
                    angDist[extractOrder[i]].push_back(new AngDist32EqualPBin());
                    angDist[extractOrder[i]].back()->SetTemperature(dataTemperature);
                    distType=3;
                }
            }
            angDist[extractOrder[i]].back()->SetPoint(stream, count, energyAngVec[extractOrder[i]][exOrderAng[j]]);
        }
        if(exOrderAng)
            delete [] exOrderAng;
    }

    for(int i=0; i<numProcess; i++)
    {
        for(int j=0; j<int(angDist[i].size()); j++)
        {
            if(!(angDist[i][j]->CheckData()))
            {
                cout << "Error in angular data ConvertMCNPtoG4NDL.cc:925" << endl;
            }
        }
    }

    //Create elastic FS files
    if((numReactions!=0)&&(LANDList[0]!=-1))
    {
        MakeElasFSFile(outDirName, isoName, isoMass, numAngEner, angDist, ascii);
    }

    //### In this section we extract the outgoing hadron/photon energy distribution data

    //Extract the LDLW, DLW blocks to get the out-going neutron energy distribtion data and the out-going neutron energy-angular distribtion data
    //have not iterated through this section yet
    if(numReacSecN!=0)
    {
        if(count>startLDLWBlock)
        {
            cout << "Error: counter is off at the start of the LDLW Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << startLDLWBlock << endl;
        }
        for(;count<startLDLWBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=1; i<numProcess; i++)
        {
            if((MTROrder[i]>0)&&(TYRList[MTROrder[i]]!=0))
            {
                for(;count<MTRListPos[MTROrder[i]]+startLDLWBlock; count++)
                {
                    stream >> dummy;
                }
                stream >> index; count++;
                LDLWList[MTROrder[i]]=index;
            }
        }

        for(int i=1; i<numProcess; i++)
        {
            for(int j=i+1; j<numProcess; j++)
            {
                if(LDLWList[extractOrder[j]]<LDLWList[extractOrder[i]])
                {
                    index = extractOrder[i];
                    extractOrder[i] = extractOrder[j];
                    extractOrder[j] = index;
                }
            }
        }
    }

    int nextLawPos, enDisLawDataPos;
    vector<int> enDisLaw[numProcess], enDisNumLawApplNReg[numProcess], enDisNumLawApplNEn[numProcess];
    vector<int*> enDisSchemeVec[numProcess], enDisRangeVec[numProcess];
    vector<double*> enDisEnApplVec[numProcess], enDisProbApplVec[numProcess];
    vector<EnergyDist*> enerDist[numProcess];
    vector<AngularEnergyDist*> angEnDist[numProcess];
    YieldDist *nYieldReac[numProcess];
    for(int i=0; i<numProcess; i++)
    {
        nYieldReac[i]=NULL;
    }

    //have not iterated through this section yet
    if(numReacSecN!=0)
    {
        for(int i=1; i<numProcess; i++)
        {
            if(abs(TYRList[extractOrder[i]])>100)
            {
                if(count>(abs(TYRList[extractOrder[i]])+startDLWBlock-101))
                {
                    cout << "Error: counter is off at the start of the energy dependant primary neutron yield Block for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << (abs(TYRList[extractOrder[i]])+startDLWBlock-101) << endl;
                }
                for(;count<(abs(TYRList[extractOrder[i]])+startDLWBlock-101); count++)
                {
                    stream >> dummy;
                }
                nYieldReac[extractOrder[i]] = new NYield1DTab;
                nYieldReac[extractOrder[i]]->ExtractMCNPData(stream, count);
            }
            if(LDLWList[extractOrder[i]]==-1)
                continue;
            if(count>(LDLWList[extractOrder[i]]+startDLWBlock-1))
            {
                cout << "Error: counter is off at the location of primary neutron energy distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                cout << "Count= " << count << " should be less than equal to " << (LDLWList[extractOrder[i]]+startDLWBlock-1) << endl;
            }
            for(;count<(LDLWList[extractOrder[i]]+startDLWBlock-1); count++)
            {
                stream >> dummy;
            }
            do
            {
                stream >> nextLawPos >> temp >> enDisLawDataPos; count=count+3;
                if(floor(temp)!=ceil(temp))
                {
                    cout << "Error: counter is off at the location of primary neutron energy distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                    cout << "NextLaw= " << temp << " should be an integer" << endl;
                }
                intTemp=int(temp);
                enDisLaw[extractOrder[i]].push_back(intTemp);
                stream >> intTemp; count++;
                enDisNumLawApplNReg[extractOrder[i]].push_back(intTemp);
                if(enDisNumLawApplNReg[extractOrder[i]].back()==0)
                {
                    enDisNumLawApplNReg[extractOrder[i]].back()=1;
                    enDisSchemeVec[extractOrder[i]].push_back(new int [enDisNumLawApplNReg[extractOrder[i]].back()]);
                    enDisRangeVec[extractOrder[i]].push_back(new int [enDisNumLawApplNReg[extractOrder[i]].back()]);

                    stream >> intTemp; count++;
                    enDisRangeVec[extractOrder[i]].back()[0]=intTemp;
                    enDisSchemeVec[extractOrder[i]].back()[0]=2;
                }
                else
                {
                    enDisSchemeVec[extractOrder[i]].push_back(new int [enDisNumLawApplNReg[extractOrder[i]].back()]);
                    enDisRangeVec[extractOrder[i]].push_back(new int [enDisNumLawApplNReg[extractOrder[i]].back()]);

                    for(int j=0; j<enDisNumLawApplNReg[extractOrder[i]].back(); j++, count++)
                    {
                        stream >> intTemp;
                        (enDisRangeVec[extractOrder[i]].back())[j]=intTemp;
                    }

                    for(int j=0; j<enDisNumLawApplNReg[extractOrder[i]].back(); j++, count++)
                    {
                        stream >> intTemp;
                        (enDisSchemeVec[extractOrder[i]].back())[j]=intTemp;
                    }

                    stream >> intTemp; count++;
                }

                enDisNumLawApplNEn[extractOrder[i]].push_back(intTemp);

                enDisEnApplVec[extractOrder[i]].push_back(new double [enDisNumLawApplNEn[extractOrder[i]].back()]);
                enDisProbApplVec[extractOrder[i]].push_back(new double [enDisNumLawApplNEn[extractOrder[i]].back()]);

                for(int j=0; j<enDisNumLawApplNEn[extractOrder[i]].back(); j++, count++)
                {
                    stream >> temp;
                    (enDisEnApplVec[extractOrder[i]].back())[j]=temp;
                }

                for(int j=0; j<enDisNumLawApplNEn[extractOrder[i]].back(); j++, count++)
                {
                    stream >> temp;
                    (enDisProbApplVec[extractOrder[i]].back())[j]=temp;
                }

                if(count>(enDisLawDataPos+startDLWBlock-1))
                {
                    cout << "Error: counter is off at the location of prompt neutron energy distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << (enDisLawDataPos+startDLWBlock-1) << endl;
                }
                for(;count<(enDisLawDataPos+startDLWBlock-1); count++)
                {
                    stream >> dummy;
                }

                // Now we extract the data depending on the format of the law using a different class for each distinct case
                // Don't have to worry about energy ordering of the laws, the law accuracy probability
                // ensures that the law will only be used in the right energy regime.

                if(enDisLaw[extractOrder[i]].back()==1)
                    enerDist[extractOrder[i]].push_back(new EnerDistEqPEnerBins());

                else if(enDisLaw[extractOrder[i]].back()==3)
                    enerDist[extractOrder[i]].push_back(new EnerDistLevScat((enDisEnApplVec[extractOrder[i]].back())[0], (enDisEnApplVec[extractOrder[i]].back())[enDisNumLawApplNEn[extractOrder[i]].back()-1]));

                else if(enDisLaw[extractOrder[i]].back()==4)
                    enerDist[extractOrder[i]].push_back(new EnerDistConTab());

                else if(enDisLaw[extractOrder[i]].back()==5)
                    enerDist[extractOrder[i]].push_back(new EnerDistGenEvapSpec());

                else if(enDisLaw[extractOrder[i]].back()==7)
                    enerDist[extractOrder[i]].push_back(new EnerDistMaxwellFisSpec());

                else if(enDisLaw[extractOrder[i]].back()==9)
                    enerDist[extractOrder[i]].push_back(new EnerDistEvapSpec());

                else if(enDisLaw[extractOrder[i]].back()==11)
                    enerDist[extractOrder[i]].push_back(new EnerDistWattSpec());

                else if(enDisLaw[extractOrder[i]].back()==22)
                    enerDist[extractOrder[i]].push_back(new EnerDistTabLinFunc());

                else if(enDisLaw[extractOrder[i]].back()==24)
                    enerDist[extractOrder[i]].push_back(new EnerDistTabMulti());

                else if(enDisLaw[extractOrder[i]].back()==44)
                    angEnDist[extractOrder[i]].push_back(new AngEnDistKallbach());

                else if(enDisLaw[extractOrder[i]].back()==61)
                    angEnDist[extractOrder[i]].push_back(new AngEnDist3DTab(startDLWBlock));

                else if(enDisLaw[extractOrder[i]].back()==66)
                    angEnDist[extractOrder[i]].push_back(new AngEnDistNBody());

                else if(enDisLaw[extractOrder[i]].back()==67)
                    angEnDist[extractOrder[i]].push_back(new AngEnDistLab3DTab());
                else
                {
                    cout << "\n### Error: Energy law not recognized! ###" << endl;
                    if((count>(nextLawPos+startDLWBlock-1))&&(nextLawPos!=0))
                    {
                        cout << "Error: counter is off at the location of prompt neutron energy distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                        cout << "Count= " << count << " should be less than equal to " << nextLawPos+startDLWBlock-1 << endl;
                    }
                    for(;count<(nextLawPos+startDLWBlock-1); count++)
                    {
                        stream >> dummy;
                    }
                    if(nextLawPos==0)
                        break;
                    else
                        continue;
                }

                if(enDisLaw[extractOrder[i]].back()<44)
                    enerDist[extractOrder[i]].back()->ExtractMCNPData(stream, count);
                else
                    angEnDist[extractOrder[i]].back()->ExtractMCNPData(stream, count);

                // go to next law position
                if((count>(nextLawPos+startDLWBlock-1))&&(nextLawPos!=0))
                {
                    cout << "Error: counter is off at the location of prompt neutron energy distribution data for reaction " << MTRList[extractOrder[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << nextLawPos+startDLWBlock-1 << endl;
                }
                for(;count<(nextLawPos+startDLWBlock-1); count++)
                {
                    stream >> dummy;
                }
            }
            while(nextLawPos>0);
        }
    }

    //### Now we extract the delayed info

    //Extract the DNU block to get the delayed neutron yield data
    NYield1DTab *dNYield=NULL;
    NDelayConstDist* nDelConst=NULL;
    if(startDelayedNUBlock!=0)
    {
        if(count>startDelayedNUBlock)
        {
            cout << "Error: counter is off at the start of the DNU Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << startDelayedNUBlock << endl;
        }
        for(;count<(startDelayedNUBlock); count++)
        {
            stream >> dummy;
        }
        stream >> dummy; count ++;

        dNYield = new NYield1DTab();
        dNYield->ExtractMCNPData(stream, count);

        nDelConst = new NDelayConstDist(numDNPrecursorFam);
        nDelConst->ExtractMCNPData(stream, count);
    }

    //### In this section we extract the delayed neutron energy distribution data

    //Extract the DNEDL and DNED blocks to get the delayed neutron energy distribution data
    int *dNEnDistPos=NULL;
    vector<int> enDisLawND, enDisNumLawApplNRegND, enDisNumLawApplNEnND;
    vector<int*> enDisSchemeVecND, enDisRangeVecND;
    vector<double*> enDisEnApplVecND, enDisProbApplVecND;
    vector<EnergyDist*> enerDistND;
    vector<AngularEnergyDist*> angEnDistND;
    if(startDNEDLBlock!=0)
    {
        if(count>startDNEDLBlock)
        {
            cout << "Error: counter is off at the start of the DNEDL Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << startDNEDLBlock << endl;
        }
        for(;count<(startDNEDLBlock); count++)
        {
            stream >> dummy;
        }

        dNEnDistPos = new int [numDNPrecursorFam];
        for(int i=0; i<numDNPrecursorFam; i++)
        {
            stream >> dNEnDistPos[i]; count++;
        }

        for(int i=0; i<numDNPrecursorFam; i++)
        {
            if(dNEnDistPos[i]!=-1)
            {
                if(count>(startDNEDBlock+dNEnDistPos[i]-1))
                {
                    cout << "Error: counter is off at the location of the delayed neutron energy distribution data for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << (startDNEDBlock+dNEnDistPos[i]-1) << endl;
                }
                for(;count<(startDNEDBlock+dNEnDistPos[i]-1); count++)
                {
                    stream >> dummy;
                }

                /*do
                {*/
                    stream >> nextLawPos >> intTemp >> enDisLawDataPos; count=count+3;
                    enDisLawND.push_back(intTemp);
                    stream >> intTemp; count++;
                    enDisNumLawApplNRegND.push_back(intTemp);

                    if(enDisNumLawApplNRegND.back()==0)
                    {
                        enDisNumLawApplNRegND.back()=1;
                        enDisSchemeVecND.push_back(new int [enDisNumLawApplNRegND.back()]);
                        enDisRangeVecND.push_back(new int [enDisNumLawApplNRegND.back()]);

                        stream >> intTemp; count++;
                        enDisRangeVecND.back()[0]=intTemp;
                        enDisSchemeVecND.back()[0]=2;
                    }
                    else
                    {
                        enDisSchemeVecND.push_back(new int [enDisNumLawApplNRegND.back()]);
                        enDisRangeVecND.push_back(new int [enDisNumLawApplNRegND.back()]);

                        for(int j=0; j<enDisNumLawApplNRegND.back(); j++, count++)
                        {
                            stream >> intTemp;
                            (enDisRangeVecND.back())[j]=intTemp;
                        }

                        for(int j=0; j<enDisNumLawApplNRegND.back(); j++, count++)
                        {
                            stream >> intTemp;
                            (enDisSchemeVecND.back())[j]=intTemp;
                        }
                        stream >> intTemp; count++;
                    }

                    enDisNumLawApplNEnND.push_back(intTemp);

                    enDisEnApplVecND.push_back(new double [enDisNumLawApplNEnND.back()]);
                    enDisProbApplVecND.push_back(new double [enDisNumLawApplNEnND.back()]);

                    for(int j=0; j<enDisNumLawApplNEnND.back(); j++, count++)
                    {
                        stream >> temp;
                        (enDisEnApplVecND.back())[j]=temp;
                    }

                    for(int j=0; j<enDisNumLawApplNEnND.back(); j++, count++)
                    {
                        stream >> temp;
                        (enDisProbApplVecND.back())[j]=temp;
                    }

                    if(count>(enDisLawDataPos+startDNEDBlock-1))
                    {
                        cout << "Error: counter is off at the location of delayed neutron energy distribution data for isotope " << isoName << endl;
                        cout << "Count= " << count << " should be less than equal to " << (enDisLawDataPos+startDNEDBlock-1) << endl;
                    }
                    for(;count<(enDisLawDataPos+startDNEDBlock-1); count++)
                    {
                        stream >> dummy;
                    }
                    // Now we extract the data depending on the format of the law using a different class for each distinct case
                    // Don't have to worry about energy ordering of the laws, the law accuracy probability
                    // ensures that the law will only be used in the right energy regime.

                    if(enDisLawND.back()==1)
                        enerDistND.push_back(new EnerDistEqPEnerBins());

                    else if(enDisLawND.back()==3)
                        enerDistND.push_back(new EnerDistLevScat((enDisEnApplVecND.back())[0], (enDisEnApplVecND.back())[enDisNumLawApplNEnND.back()-1]));

                    else if(enDisLawND.back()==4)
                        enerDistND.push_back(new EnerDistConTab());

                    else if(enDisLawND.back()==5)
                        enerDistND.push_back(new EnerDistGenEvapSpec());

                    else if(enDisLawND.back()==7)
                        enerDistND.push_back(new EnerDistMaxwellFisSpec());

                    else if(enDisLawND.back()==9)
                        enerDistND.push_back(new EnerDistEvapSpec());

                    else if(enDisLawND.back()==11)
                        enerDistND.push_back(new EnerDistWattSpec());

                    else if(enDisLawND.back()==22)
                        enerDistND.push_back(new EnerDistTabLinFunc());

                    else if(enDisLawND.back()==24)
                        enerDistND.push_back(new EnerDistTabMulti());

                    else if(enDisLawND.back()==44)
                        angEnDistND.push_back(new AngEnDistKallbach());

                    else if(enDisLawND.back()==61)
                        angEnDistND.push_back(new AngEnDist3DTab(startDNEDBlock));

                    else if(enDisLawND.back()==66)
                        angEnDistND.push_back(new AngEnDistNBody());

                    else if(enDisLawND.back()==67)
                        angEnDistND.push_back(new AngEnDistLab3DTab());
                    else
                    {
                        cout << "\n### Error: Energy law not recognized! ###" << endl;
                        if((count>(nextLawPos+startDNEDBlock-1))&&(nextLawPos!=0))
                        {
                            cout << "Error: counter is off at the location of delayed neutron energy distribution data for isotope " << isoName << endl;
                            cout << "Count= " << count << " should be less than equal to " << (nextLawPos+startDNEDBlock-1) << endl;
                        }
                        for(;count<(nextLawPos+startDNEDBlock-1); count++)
                        {
                            stream >> dummy;
                        }
                        continue;
                    }

                    if(enDisLawND.back()<44)
                        enerDistND.back()->ExtractMCNPData(stream, count);

                    else
                        angEnDistND.back()->ExtractMCNPData(stream, count);

                    // go to next law position
                    if((count>(nextLawPos+startDNEDBlock-1))&&(nextLawPos!=0))
                    {
                        cout << "Error: counter is off at the location of delayed neutron energy distribution data for isotope " << isoName << endl;
                        cout << "Count= " << count << " should be less than equal to " << (nextLawPos+startDNEDBlock-1) << endl;
                    }
                    for(;count<(nextLawPos+startDNEDBlock-1); count++)
                    {
                        stream >> dummy;
                    }
                /*}
                while(nextLawPos>0);*/
            }
        }
    }

    //### Extract photon data, follows same format as the neutron data above, except for a few minor differences
    // Probably only one energy distribution per photon production reaction, don't need vector containers

    vector<int> MTRPList;// this array sets the order of the processes
    vector<int> MTRPListPos; //position of the reaction MT in MT list

    if(numPReactions!=0)
    {
        if(count>startMTRPBlock)
        {
            cout << "Error: counter is off at the start of the MTRP Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << startMTRPBlock << endl;
        }
        for(;count<startMTRPBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=0; i<numPReactions; i++, count++)
        {
            stream >> index;
            //make sure these MT numbers are right
            for(int j=1; j<numProcess; j++)
            {
                if(MTROrder[j]>0)
                {
                    if(MTRList[MTROrder[j]]==int(index/1000))
                    {
                        MTRPList.push_back(index);
                        MTRPListPos.push_back(i);
                        break;
                    }
                }
            }
        }
    }

    int numPProcess = MTRPListPos.size();
    int *LSIGPList = new int [numPProcess]; //starting point of the reaction cross-section data in the data file
    for(int i=0; i<numPProcess; i++)
    {
        LSIGPList[i]=-1;
    }
    int *MFType = new int [numPProcess];
    for(int i=0; i<numPProcess; i++)
    {
        MFType[i]=-1;
    }
    CSDist **pCSVec = new CSDist *[numPProcess]; //cross section array for each reaction
    for(int i=0; i<numPProcess; i++)
    {
        pCSVec[i]=NULL;
    }
    int *LANDPList = new int [numPProcess]; //starting point of the reaction's out going neutron angular distribution data
    for(int i=0; i<numPProcess; i++)
    {
        LANDPList[i]=-1;
    }
    int *numAngEnerP = new int [numPProcess]; // the number of incoming neutron energy points for the angular distribution
    for(int i=0; i<numPProcess; i++)
    {
        numAngEnerP[i]=-1;
    }
    int *LDLWPList= new int [numPProcess]; //starting point of the reaction's out going neutron energy distribution data
    for(int i=0; i<numPProcess; i++)
    {
        LDLWPList[i]=-1;
    }
    vector <int> numYieldP;
    vector <int*> mtNums;
    int *extractOrderP = new int [numPProcess]; // the order in which the data for each reaction will be extracted, pos [0] will be the first to be extracted,
    //where the number it contains corresponds to the index in the above arrays that the reaction belongs to
    for(int i=0; i<numPProcess; i++)
    {
        extractOrderP[i]=i;
    }

    if(numPReactions!=0)
    {
        //### in this section we extract the CS data for the different photon production reactions
        //Extract the LSIGP and the SIGP blocks to get the cross-section data for photon reactions
        if(count!=startLSIGPBlock)
        {
            cout << "Error: counter is off at the start of the LSIGP Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be equal to " << (startLSIGPBlock) << endl;
        }
        for(;count<startLSIGPBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=0; i<numPProcess; i++, count++)
        {
            for(;count<MTRPListPos[i]+startLSIGPBlock; count++)
            {
                stream >> dummy;
            }
            stream >> index;
            LSIGPList[i]=index;
        }

        for(int i=0; i<numPProcess; i++)
        {
            for(int j=i+1; j<numPProcess; j++)
            {
                if(LSIGPList[extractOrderP[j]]<LSIGPList[extractOrderP[i]])
                {
                    index = extractOrderP[i];
                    extractOrderP[i] = extractOrderP[j];
                    extractOrderP[j] = index;
                }
            }
        }

        int MTRIndex=0;
        for(int i=0; i<numPProcess; i++)
        {
            if(LSIGPList[extractOrderP[i]]==-1)
                continue;
            if(count>(LSIGPList[extractOrderP[i]]+startCSPBlock-1))
            {
                cout << "Error: counter is off at the location of the cross-section data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                cout << "Count= " << count << " should be less than equal to " << (LSIGPList[extractOrderP[i]]+startCSPBlock-1) << endl;
            }
            for(;count<(LSIGPList[extractOrderP[i]]+startCSPBlock-1); count++)
            {
                stream >> dummy;
            }
            stream >> intTemp; count++;
            MFType[extractOrderP[i]] = intTemp;

            for(int j=1; j<numProcess; j++)
            {
                if(MTRList[j]==int(MTRPList[i]/1000))
                {
                    MTRIndex=j;
                    break;
                }
            }

            if(MFType[extractOrderP[i]]==13)
            {
                pCSVec[extractOrderP[i]] = new CSDist1DTabP();
            }
            else
            {
                pCSVec[extractOrderP[i]] = new CSDistYieldComp();
            }
            pCSVec[extractOrderP[i]]->ExtractMCNPData(stream, count);
            // ### here we use the MT value we have determined instead of that provided by MCNP for convience, check if they are equivalent
            pCSVec[extractOrderP[i]]->SetCSData(nCSVec[MTRIndex]);
        }

        //### In this Section we extract the out-going photon angular distribution data that is independant
        // of energy the rest of the angular data is extracted in the energy distribution

        //Extract the LANDP and ANDP blocks to get the out-going photon angular distribution data
        if(count>startLANDPBlock)
        {
            cout << "Error: counter is off at the start of the LANDP Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << (startLANDPBlock) << endl;
        }
        for(;count<startLANDPBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=0; i<numPProcess; i++, count++)
        {
            for(;count<MTRPListPos[i]+startLANDPBlock; count++)
            {
                stream >> dummy;
            }
            stream >> index;
            LANDPList[i]=index;
        }

        for(int i=0; i<numPProcess; i++)
        {
            for(int j=i+1; j<numPProcess; j++)
            {
                if(LANDPList[extractOrderP[j]]<LANDPList[extractOrderP[i]])
                {
                    index = extractOrderP[i];
                    extractOrderP[i] = extractOrderP[j];
                    extractOrderP[j] = index;
                }
            }
        }
    }

    double **energyAngVecP = new double *[numPProcess]; // incoming neutron kinetic energy for Ang dist data
    for(int i=0; i<numPProcess; i++)
    {
        energyAngVecP[i]=NULL;
    }
    int **angTabPosVecP= new int *[numPProcess]; // the location of the angular probability table associated to the incoming neutron kinetic energy
    for(int i=0; i<numPProcess; i++)
    {
        angTabPosVecP[i]=NULL;
    }
    int *exOrderAngP=NULL;
    vector<AngularDist*> *angDistP = new vector<AngularDist*> [numPProcess];

    if(numPReactions!=0)
    {
        for(int i=0; i<numPProcess; i++)
        {
            // an isotropic distribution is assumed
            if(LANDPList[extractOrderP[i]]==0)
            {
                angDistP[extractOrderP[i]].push_back(new AngDistIsoP());
                angDistP[extractOrderP[i]].back()->SetTemperature(dataTemperature);
                numAngEnerP[extractOrderP[i]]=2;
                continue;
            }

            if(LANDPList[extractOrderP[i]]==-1)
                continue;

            if(count>(LANDPList[extractOrderP[i]]+startANDPBlock-1))
            {
                cout << "Error: counter is off at the location of the out-going photon angular distribution data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                cout << "Count= " << count << " should be less than equal to " << (LANDPList[extractOrderP[i]]+startANDPBlock-1) << endl;
            }
            for(;count<(LANDPList[extractOrderP[i]]+startANDPBlock-1); count++)
            {
                stream >> dummy;
            }
            stream >> intTemp; count++;
            numAngEnerP[extractOrderP[i]]=intTemp;

            energyAngVecP[extractOrderP[i]] = new double[numAngEnerP[extractOrderP[i]]];
            for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++, count++)
            {
                stream >> temp;
                energyAngVecP[extractOrderP[i]][j] = temp;
            }

            angTabPosVecP[extractOrderP[i]] = new int[numAngEnerP[extractOrderP[i]]];
            for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++, count++)
            {
                stream >> temp;
                angTabPosVecP[extractOrderP[i]][j] = temp;
            }

            exOrderAngP = new int[numAngEnerP[extractOrderP[i]]];
            for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++)
            {
                exOrderAngP[j] = j;
            }
            //## we may be creating an error here by mixing up the order of the distributions with respect to energy
            for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++)
            {
                for(int k=j+1; k<numAngEnerP[extractOrderP[i]]; k++)
                {
                    if(abs(angTabPosVecP[extractOrderP[i]][exOrderAngP[k]])<abs(angTabPosVecP[extractOrderP[i]][exOrderAngP[j]]))
                    {
                        index = exOrderAngP[j];
                        exOrderAngP[j] = exOrderAngP[k];
                        exOrderAngP[k] = index;
                    }
                }
            }

            distType=0;
            for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++)
            {
            //Change this so that the data is inputted into the appropriate daughter of the mother class AngularDist
                //the angular distribution is represented by a probability density function where each of the 32 points/bins are spaced along the entire cosine of the scattering angle axis
                //such that the integral of the probability density with the cosine of the scattering angle inbetween them is kept constant, meaning that the bins are in a equilprobable arrangement
                // the first bin probably corresponds to a cosine value of -1 and and the last bin being some where around 1 depending on the spacing
                if(angTabPosVecP[extractOrderP[i]][exOrderAngP[j]]!=0)
                {
                    if(count>(startANDPBlock+angTabPosVecP[extractOrderP[i]][exOrderAngP[j]]-1))
                    {
                        cout << "Error: counter is off at the location of photon angular distribution data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                        cout << "Count= " << count << " should be less than equal to " << (startANDPBlock+angTabPosVecP[extractOrderP[i]][exOrderAngP[j]]-1) << endl;
                    }
                    for(;count<(startANDPBlock+angTabPosVecP[extractOrderP[i]][exOrderAngP[j]]-1); count++)
                    {
                        stream >> dummy;
                    }
                    if(distType!=1)
                    {
                        angDistP[extractOrderP[i]].push_back(new AngDist32EqPBinP());
                        distType=1;
                    }
                }
                else if(distType!=2)
                {
                    angDistP[extractOrderP[i]].push_back(new AngDistIsoP());
                    distType=2;
                }
                angDistP[extractOrderP[i]].back()->SetPoint(stream, count, energyAngVecP[extractOrderP[i]][exOrderAngP[j]]);

            }
            if(exOrderAngP)
                delete [] exOrderAngP;

        }

        for(int i=0; i<numPProcess; i++)
        {
            for(int j=0; j<int(angDistP[i].size()); j++)
            {
                if(!(angDistP[i][j]->CheckData()))
                {
                    cout << "Error in angular data ConvertMCNPtoG4NDL.cc:1674" << endl;
                }
            }
        }

        //### In this section we extract the outgoing photon energy distribution data

        //Extract the LDLWP, DLWP blocks to get the out-going photon energy distribtion data and the out-going photon energy-angular distribtion data
        if(count>startLDLWPBlock)
        {
            cout << "Error: counter is off at the start of the LDLWP Block for isotope " << isoName << endl;
            cout << "Count= " << count << " should be less than equal to " << (startLDLWPBlock) << endl;
        }
        for(;count<startLDLWPBlock; count++)
        {
            stream >> dummy;
        }

        for(int i=0; i<numPProcess; i++, count++)
        {
            for(;count<MTRPListPos[i]+startLDLWPBlock; count++)
            {
                stream >> dummy;
            }
            stream >> index;
            LDLWPList[i]=index;
        }

        for(int i=0; i<numPProcess; i++)
        {
            for(int j=i+1; j<numPProcess; j++)
            {
                if(LDLWPList[extractOrderP[j]]<LDLWPList[extractOrderP[i]])
                {
                    index = extractOrderP[i];
                    extractOrderP[i] = extractOrderP[j];
                    extractOrderP[j] = index;
                }
            }
        }
    }

    vector<int> *enDisLawP = new vector<int> [numPProcess], *enDisNumLawApplNRegP = new vector<int> [numPProcess], *enDisNumLawApplNEnP = new vector<int> [numPProcess];
    vector<int*> *enDisSchemeVecP = new vector<int*> [numPProcess], *enDisRangeVecP = new vector<int*> [numPProcess];
    vector<double*> *enDisEnApplVecP = new vector<double*> [numPProcess], *enDisProbApplVecP = new vector<double*> [numPProcess];
    vector<EnergyDist*> *enerDistP = new vector<EnergyDist*> [numPProcess];
    vector<AngularEnergyDist*> *angEnDistP = new vector<AngularEnergyDist*> [numPProcess];

    if(numPReactions!=0)
    {
        for(int i=0; i<numPProcess; i++)
        {
            if(LDLWPList[extractOrderP[i]]==-1)
                continue;
            if(count>(LDLWPList[extractOrderP[i]]+startDLWPBlock-1))
            {
                cout << "Error: counter is off at the location of the out-going photon energy distribution data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                cout << "Count= " << count << " should be less than equal to " << (LDLWPList[extractOrderP[i]]+startDLWPBlock-1) << endl;
            }
            for(;count<(LDLWPList[extractOrderP[i]]+startDLWPBlock-1); count++)
            {
                stream >> dummy;
            }
            do
            {
                stream >> nextLawPos >> intTemp >> enDisLawDataPos; count=count+3;
                enDisLawP[extractOrderP[i]].push_back(intTemp);
                stream >> intTemp; count++;
                enDisNumLawApplNRegP[extractOrderP[i]].push_back(intTemp);
                if(enDisNumLawApplNRegP[extractOrderP[i]].back()==0)
                {
                    enDisNumLawApplNRegP[extractOrderP[i]].back()=1;
                    enDisSchemeVecP[extractOrderP[i]].push_back(new int [enDisNumLawApplNRegP[extractOrderP[i]].back()]);
                    enDisRangeVecP[extractOrderP[i]].push_back(new int [enDisNumLawApplNRegP[extractOrderP[i]].back()]);

                    stream >> intTemp; count++;
                    enDisRangeVecP[extractOrderP[i]].back()[0]=intTemp;
                    enDisSchemeVecP[extractOrderP[i]].back()[0]=2;
                }
                else
                {
                    enDisSchemeVecP[extractOrderP[i]].push_back(new int [enDisNumLawApplNRegP[extractOrderP[i]].back()]);
                    enDisRangeVecP[extractOrderP[i]].push_back(new int [enDisNumLawApplNRegP[extractOrderP[i]].back()]);

                    for(int j=0; j<enDisNumLawApplNRegP[extractOrderP[i]].back(); j++, count++)
                    {
                        stream >> intTemp;
                        (enDisRangeVecP[extractOrderP[i]].back())[j]=intTemp;
                    }

                    for(int j=0; j<enDisNumLawApplNRegP[extractOrderP[i]].back(); j++, count++)
                    {
                        stream >> intTemp;
                        (enDisSchemeVecP[extractOrderP[i]].back())[j]=intTemp;
                    }

                    stream >> intTemp; count++;
                }

                enDisNumLawApplNEnP[extractOrderP[i]].push_back(intTemp);

                enDisEnApplVecP[extractOrderP[i]].push_back(new double [enDisNumLawApplNEnP[extractOrderP[i]].back()]);
                enDisProbApplVecP[extractOrderP[i]].push_back(new double [enDisNumLawApplNEnP[extractOrderP[i]].back()]);

                for(int j=0; j<enDisNumLawApplNEnP[extractOrderP[i]].back(); j++, count++)
                {
                    stream >> temp;
                    (enDisEnApplVecP[extractOrderP[i]].back())[j]=temp;
                }

                for(int j=0; j<enDisNumLawApplNEnP[extractOrderP[i]].back(); j++, count++)
                {
                    stream >> temp;
                    (enDisProbApplVecP[extractOrderP[i]].back())[j]=temp;
                }
                if(count>(enDisLawDataPos+startDLWPBlock-1))
                {
                    cout << "Error: counter is off at the location of photon energy distribution data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << (enDisLawDataPos+startDLWPBlock-1) << endl;
                }
                for(;count<(enDisLawDataPos+startDLWPBlock-1); count++)
                {
                    stream >> dummy;
                }

                // Now we extract the data depending on the format of the law using a different class for each distinct case
                // Don't have to worry about energy ordering of the laws, the law accuracy probability
                // ensures that the law will only be used in the right energy regime.
                /*
                if(enDisLawP[extractOrderP[i]].back()==1)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistEqPEnerBins());

                else */if(enDisLawP[extractOrderP[i]].back()==2)
                    enerDistP[extractOrderP[i]].push_back(new EnerDist1PhEner(isoMass, (enDisEnApplVecP[extractOrderP[i]].back())[0], (enDisEnApplVecP[extractOrderP[i]].back())[enDisNumLawApplNEnP[extractOrderP[i]].back()-1]));
                /*
                else if(enDisLawP[extractOrderP[i]].back()==3)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistLevScat());
                */
                else if(enDisLawP[extractOrderP[i]].back()==4)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistConTab());
                /*
                else if(enDisLawP[extractOrderP[i]].back()==5)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistGenEvapSpec());

                else if(enDisLawP[extractOrderP[i]].back()==7)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistMaxwellFisSpec());

                else if(enDisLawP[extractOrderP[i]].back()==9)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistEvapSpec());

                else if(enDisLawP[extractOrderP[i]].back()==11)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistWattSpec());

                else if(enDisLawP[extractOrderP[i]].back()==22)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistTabLinFunc(startDLWPBlock));

                else if(enDisLawP[extractOrderP[i]].back()==24)
                    enerDistP[extractOrderP[i]].push_back(new EnerDistTabMulti());
                */
                else if(enDisLawP[extractOrderP[i]].back()==44)
                    angEnDistP[extractOrderP[i]].push_back(new AngEnDistKallbach());

                else if(enDisLawP[extractOrderP[i]].back()==61)
                    angEnDistP[extractOrderP[i]].push_back(new AngEnDist3DTab(startDLWPBlock));

                else if(enDisLawP[extractOrderP[i]].back()==66)
                    angEnDistP[extractOrderP[i]].push_back(new AngEnDistNBody());

                else if(enDisLawP[extractOrderP[i]].back()==67)
                    angEnDistP[extractOrderP[i]].push_back(new AngEnDistLab3DTab());
                else
                {
                    cout << "\n### Error: Energy law not recognized! ###" << endl;
                    if((count>(nextLawPos+startDLWPBlock-1))&&(nextLawPos!=0))
                    {
                        cout << "Error: counter is off at the location of photon energy distribution data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                        cout << "Count= " << count << " should be less than equal to " << (nextLawPos+startDLWPBlock-1) << endl;
                    }
                    for(;count<(nextLawPos+startDLWPBlock-1); count++)
                    {
                        stream >> dummy;
                    }
                    continue;
                }

                if(enDisLawP[extractOrderP[i]].back()<44)
                    enerDistP[extractOrderP[i]].back()->ExtractMCNPData(stream, count);

                // these energy angular distributions should not occur for photon data, check
                else
                {
                    angEnDistP[extractOrderP[i]].back()->ExtractMCNPData(stream, count);
                    cout << "\n### Error: photon energy dependant angular data is being collected, should use it to create FSMF6/ data in the MakeCaptureFSFile function" << endl;
                }


                // go to next law position
                if((count>(nextLawPos+startDLWPBlock-1))&&(nextLawPos!=0))
                {
                    cout << "Error: counter is off at the location of photon energy distribution data for reaction " << MTRPList[extractOrderP[i]] << " for isotope " << isoName << endl;
                    cout << "Count= " << count << " should be less than equal to " << (nextLawPos+startDLWPBlock-1) << endl;
                }
                for(;count<(nextLawPos+startDLWPBlock-1); count++)
                {
                    stream >> dummy;
                }
            }
            while(nextLawPos>0);
        }
    }

    //Create capture FS file
    if(MTRListPos[1]!=-1)
    {
        MakeCaptureFSFile(outDirName, isoName, isoMass, dataTemperature, TYRList, pCSVec, angDistP, enDisLawP, enerDistP, MTRPList, extractOrderP, numAngEnerP, enDisNumLawApplNRegP,
                        enDisNumLawApplNEnP, enDisSchemeVecP, enDisRangeVecP, enDisEnApplVecP, enDisProbApplVecP, energyAngVecP, nCSVec[1], ascii);
    }
    // don't need this data
    /*
    for(;count<(startYPBlock); count++)
    {
        stream >> dummy;
    }

    while(count<startFisBlock)
    {
        stream >> intTemp; count++;
        numYieldP.push_back(intTemp);
        mtNums.push_back(new int [numYieldP.back()]);

        for(int i=0; i<numYieldP.back(); i++, count++)
        {
            stream >> intTemp;
            mtNums.back()[i] = intTemp;
        }
    }
    */

    //Create Fission FS files
    if(MTRListPos[2]==-1)
    {
        createMT2Flag=true;
        CreateMT2(MTRListPos, outDirName, isoName, isoNum, isoMass, dataTemperature, MTRList, nCSVec, enDisLaw, enDisNumLawApplNReg, enDisNumLawApplNEn,
                        enDisSchemeVec, enDisRangeVec, enDisEnApplVec, enDisProbApplVec, enerDist, nYieldReac, reacQValue, numAngEner, angDist,
                        angDistInEnDistFlag, angEnDist, &pCSVec, TYRList, &numAngEnerP, &angDistP, &enDisLawP, &enDisNumLawApplNRegP, &enDisNumLawApplNEnP, &enDisSchemeVecP,
                        &enDisRangeVecP, &enDisEnApplVecP, &enDisProbApplVecP, &enerDistP, MTRPList, &energyAngVecP, ascii);
    }
    if(MTRListPos[2]!=-1)
    {
        MakeFissionFSFile(MTRList, MTRListPos, outDirName, isoName, isoMass, dataTemperature, nCSVec, enDisLaw, enDisNumLawApplNReg, enDisNumLawApplNEn, enDisSchemeVec, enDisRangeVec,
                            enDisEnApplVec, enDisProbApplVec, enerDist, angEnDist, enDisLawND, enDisNumLawApplNRegND, enDisNumLawApplNEnND, enDisSchemeVecND, enDisRangeVecND,
                            enDisEnApplVecND, enDisProbApplVecND, enerDistND, reacQValue, numAngEner, angDist, angDistInEnDistFlag, pCSVec, TYRList, promptYieldFlag, totalYieldFlag,
                            dNPromptYieldDist, dNTotalYieldDist, dNYield, nDelConst, numAngEnerP, angDistP, enDisLawP, enDisNumLawApplNRegP, enDisNumLawApplNEnP, enDisSchemeVecP, enDisRangeVecP,
                            enDisEnApplVecP, enDisProbApplVecP, enerDistP, MTRPList, energyAngVecP, ascii);
    }

    if(MTRListPos[7]==-1)
    {
        createMT4Flag=true;
        CreateMT4(MTRListPos, outDirName, isoName, isoNum, isoMass, dataTemperature, MTRList, nCSVec, enDisLaw, enDisNumLawApplNReg, enDisNumLawApplNEn,
                        enDisSchemeVec, enDisRangeVec, enDisEnApplVec, enDisProbApplVec, enerDist, nYieldReac, reacQValue, numAngEner, angDist,
                        angDistInEnDistFlag, angEnDist, &pCSVec, TYRList, &numAngEnerP, &angDistP, &enDisLawP, &enDisNumLawApplNRegP, &enDisNumLawApplNEnP, &enDisSchemeVecP,
                        &enDisRangeVecP, &enDisEnApplVecP, &enDisProbApplVecP, &enerDistP, MTRPList, &energyAngVecP, ascii);
    }
    //Create Inelastic FS Files
    MakeInElasticFSFile(MTRListPos, outDirName, isoName, isoNum, isoMass, dataTemperature, MTRList, nCSVec, enDisLaw, enDisNumLawApplNReg, enDisNumLawApplNEn,
                        enDisSchemeVec, enDisRangeVec, enDisEnApplVec, enDisProbApplVec, enerDist, nYieldReac, reacQValue, numAngEner, angDist,
                        angDistInEnDistFlag, angEnDist, pCSVec, TYRList, numAngEnerP, angDistP, enDisLawP, enDisNumLawApplNRegP, enDisNumLawApplNEnP, enDisSchemeVecP,
                        enDisRangeVecP, enDisEnApplVecP, enDisProbApplVecP, enerDistP, MTRPList, energyAngVecP, ascii);

    //delete out going primary neutron data
    if(energyCSVec!=NULL)
        delete[] energyCSVec;

    for(int i=0; i<numProcess; i++)
    {
        if(angTabPosVec[i]!=NULL)
            delete[] angTabPosVec[i];
        if(nCSVec[i]!=NULL)
            delete nCSVec[i];
        if(energyAngVec[i]!=NULL)
            delete[] energyAngVec[i];
        if(nYieldReac[i]!=NULL)
            delete nYieldReac[i];

        for(int j=0; j<int(angDist[i].size()); j++)
        {
            if(angDist[i][j]!=NULL)
                delete angDist[i][j];
        }
        for(int j=0; j<int(enDisSchemeVec[i].size()); j++)
        {
            if(enDisSchemeVec[i][j]!=NULL)
                delete [] enDisSchemeVec[i][j];
        }
        for(int j=0; j<int(enDisRangeVec[i].size()); j++)
        {
            if(enDisRangeVec[i][j]!=NULL)
                delete [] enDisRangeVec[i][j];
        }
        for(int j=0; j<int(enDisEnApplVec[i].size()); j++)
        {
            if(enDisEnApplVec[i][j]!=NULL)
                delete [] enDisEnApplVec[i][j];
        }
        for(int j=0; j<int(enDisProbApplVec[i].size()); j++)
        {
            if(enDisProbApplVec[i][j]!=NULL)
                delete [] enDisProbApplVec[i][j];
        }
        for(int j=0; j<int(enerDist[i].size()); j++)
        {
            if((!(createMT4Flag&&(i>=numProcess2)))&&(!(createMT2Flag&&(i>2)&&(i<7))))
            {
                if(enerDist[i][j]!=NULL)
                    delete enerDist[i][j];
            }
        }
        for(int j=0; j<int(angEnDist[i].size()); j++)
        {
            if((!(createMT4Flag&&(i>=numProcess2)))&&(!(createMT2Flag&&(i>2)&&(i<7))))
            {
                if(angEnDist[i][j]!=NULL)
                    delete angEnDist[i][j];
            }
        }
    }

    if(dNPromptYieldDist!=NULL)
        delete dNPromptYieldDist;
    if(dNTotalYieldDist!=NULL)
        delete dNTotalYieldDist;

    //delete out going photon data
    if(numPProcess!=0)
    {
        if(LSIGPList!=NULL)
            delete [] LSIGPList;
        if(MFType!=NULL)
            delete [] MFType;
        if(LANDPList!=NULL)
            delete [] LANDPList;
        if(numAngEnerP!=NULL)
            delete [] numAngEnerP; // the number of incoming neutron energy points for the angular distribution
        if(LDLWPList!=NULL)
            delete [] LDLWPList;
        /*for(int i=0; i<int(mtNums.size());i++)
        {
            if(mtNums[i]!=NULL)
                delete [] mtNums[i];
        }*/
        if(extractOrderP!=NULL)
            delete [] extractOrderP;
        if(enDisLawP!=NULL)
            delete [] enDisLawP;
        if(enDisNumLawApplNRegP!=NULL)
            delete [] enDisNumLawApplNRegP;
        if(enDisNumLawApplNEnP!=NULL)
            delete [] enDisNumLawApplNEnP;

        for(int i=0; i<numPProcess; i++)
        {
            if(energyAngVecP[i]!=NULL)
                delete [] energyAngVecP[i];
            if(angTabPosVecP[i]!=NULL)
                delete [] angTabPosVecP[i];
            if(pCSVec[i]!=NULL)
                delete pCSVec[i];

            for(int j=0; j<int(angDistP[i].size()); j++)
            {
                if(angDistP[i][j]!=NULL)
                    delete angDistP[i][j];
            }
            for(int j=0; j<int(enDisSchemeVecP[i].size()); j++)
            {
                if(enDisSchemeVecP[i][j]!=NULL)
                    delete [] enDisSchemeVecP[i][j];
            }
            for(int j=0; j<int(enDisRangeVecP[i].size()); j++)
            {
                if(enDisRangeVecP[i][j]!=NULL)
                    delete [] enDisRangeVecP[i][j];
            }
            for(int j=0; j<int(enDisEnApplVecP[i].size()); j++)
            {
                if(enDisEnApplVecP[i][j]!=NULL)
                    delete [] enDisEnApplVecP[i][j];
            }
            for(int j=0; j<int(enDisProbApplVecP[i].size()); j++)
            {
                if(enDisProbApplVecP[i][j]!=NULL)
                    delete [] enDisProbApplVecP[i][j];
            }
            for(int j=0; j<int(enerDistP[i].size()); j++)
            {
                if(enerDistP[i][j]!=NULL)
                    delete enerDistP[i][j];
            }
            for(int j=0; j<int(angEnDistP[i].size()); j++)
            {
                if(angEnDistP[i][j]!=NULL)
                    delete angEnDistP[i][j];
            }
        }

        if(energyAngVecP!=NULL)
            delete [] energyAngVecP;
        if(angDistP!=NULL)
            delete [] angDistP;
        if(angTabPosVecP!=NULL)
            delete [] angTabPosVecP;
        if(pCSVec!=NULL)
            delete[] pCSVec;
        if(enDisSchemeVecP!=NULL)
            delete [] enDisSchemeVecP;
        if(enDisRangeVecP!=NULL)
            delete [] enDisRangeVecP;
        if(enDisProbApplVecP!=NULL)
            delete [] enDisProbApplVecP;
        if(enDisEnApplVecP!=NULL)
            delete [] enDisEnApplVecP;
        if(enerDistP!=NULL)
            delete [] enerDistP;
    }

    //delete out going delayed neutron data
    if(dNYield)
        delete dNYield;

    if(nDelConst)
        delete nDelConst;

    if(dNEnDistPos)
        delete [] dNEnDistPos;

    for(int i=0; i<int(enDisSchemeVecND.size());i++)
    {
        if(enDisSchemeVecND[i]!=NULL)
            delete [] enDisSchemeVecND[i];
    }
    for(int i=0; i<int(enDisRangeVecND.size());i++)
    {
        if(enDisRangeVecND[i]!=NULL)
            delete [] enDisRangeVecND[i];
    }
    for(int i=0; i<int(enDisEnApplVecND.size());i++)
    {
        if(enDisEnApplVecND[i]!=NULL)
            delete [] enDisEnApplVecND[i];
    }
    for(int i=0; i<int(enDisProbApplVecND.size());i++)
    {
        if(enDisProbApplVecND[i]!=NULL)
            delete [] enDisProbApplVecND[i];
    }
    for(int i=0; i<int(enerDistND.size());i++)
    {
        if(enerDistND[i]!=NULL)
            delete enerDistND[i];
    }
    for(int i=0; i<int(angEnDistND.size());i++)
    {
        if(angEnDistND[i]!=NULL)
            delete angEnDistND[i];
    }
    return 0;
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;
        closedir (pDir);
    }

    return bExists;
}

void MakeReactionListFile(string outDirName, string isoName, int MTRList[], int MTRListPos[])
{
    std::stringstream stream;
    string dirName = outDirName+"ReactionList/";
    if(!(DirectoryExists(dirName.c_str())))
    {
        system( ("mkdir -p -m=666 "+dirName).c_str());
        if(DirectoryExists(dirName.c_str()))
        {
            cout << "created directory " << dirName << "\n" << endl;
        }
        else
        {
            cout << "\nError: could not create directory " << dirName << "\n" << endl;
            return;
        }
    }
    stream << "Fission Reactions Present: ";
    for(int i=2; i<7; i++)
    {
        if(MTRListPos[i]!=-1)
        {
            stream << MTRList[i] << " ";
        }
    }
    stream << "\n\nFission Reactions Not Present: ";
    for(int i=2; i<7; i++)
    {
        if(MTRListPos[i]==-1)
        {
            stream << MTRList[i] << " ";
        }
    }

    stream << "\n\nInelastic Reactions Present: ";
    for(int i=7; i<numProcess; i++)
    {
        if(MTRListPos[i]!=-1)
        {
            stream << MTRList[i] << " ";
        }
    }
    stream << "\n\nInelastic Reactions Not Present: ";
    for(int i=7; i<numProcess; i++)
    {
        if(MTRListPos[i]==-1)
        {
            stream << MTRList[i] << " ";
        }
    }
    SetDataStream( dirName+isoName, stream, true);
}

void MakeCSDataFile(string outDirName, string isoName, int MTRNum[], CSDist **csDist, bool ascii)
{
    // creates the Cross-section data files
    string outDirNameCSProc[4];
    std::stringstream stream;
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    outDirNameCSProc[0] = outDirName +"Elastic/CrossSection/";
    outDirNameCSProc[1] = outDirName +"Capture/CrossSection/";
    outDirNameCSProc[2] = outDirName +"Fission/CrossSection/";
    outDirNameCSProc[3] = outDirName +"Inelastic/CrossSection/";

    // creates temperature directories for the converted files to go if they don't already exist
    for(int i=0; i<4; i++)
    {
        if(!(DirectoryExists((outDirNameCSProc[i]).c_str())))
        {
            system( ("mkdir -p -m=666 "+outDirNameCSProc[i]).c_str());
            if(DirectoryExists((outDirNameCSProc[i]).c_str()))
            {
                cout << "created directory " << outDirNameCSProc[i] << "\n" << endl;
            }
            else
            {
                cout << "\nError: could not create directory " << outDirNameCSProc[i] << "\n" << endl;
                return;
            }
        }
    }

    for(int i=0; i<2; i++)
    {
        if(csDist[i])
        {
            stream << std::setw(14) << std::right << MTRNum[i] << '\n';
            stream << std::setw(14) << std::right << '0' << '\n';

            csDist[i]->WriteG4NDLCSData(stream);
            SetDataStream( outDirNameCSProc[i]+isoName, stream, ascii);

            stream.clear();
            stream.str("");
        }
    }

    if(csDist[2])
    {
        stream << std::setw(14) << std::right << MTRNum[2] << '\n';
        stream << std::setw(14) << std::right << '0' << '\n';

        csDist[2]->WriteG4NDLCSData(stream);
        SetDataStream( outDirNameCSProc[2]+isoName, stream, ascii);

        stream.clear();
        stream.str("");
    }
    else
    {
        stream << std::setw(14) << std::right << MTRNum[2] << '\n';
        stream << std::setw(14) << std::right << '0' << '\n';

        bool first=true;
        CSDist1DTab* sumCSData = NULL;
        for(int i=3; i<7; i++)
        {
            if(csDist[i])
            {
                if(first)
                {
                    sumCSData = new CSDist1DTab(dynamic_cast<CSDist1DTab*>(csDist[i]));
                    first=false;
                }
                else
                    sumCSData->AddData(csDist[i]);
            }
        }

        if(sumCSData)
        {
            sumCSData->WriteG4NDLCSData(stream);
            SetDataStream( outDirNameCSProc[2]+isoName, stream, ascii);
            delete sumCSData;
        }

        stream.clear();
        stream.str("");
    }

    stream << std::setw(14) << std::right << 0 << '\n';
    stream << std::setw(14) << std::right << '0' << '\n';

    bool first=true;
    CSDist1DTab* sumCSData = NULL;
    if(csDist[7])
    {
        for(int i=7; i<numProcess2; i++)
        {
            if(csDist[i])
            {
                if(first)
                {
                    sumCSData = new CSDist1DTab(dynamic_cast<CSDist1DTab*>(csDist[i]));
                    first=false;
                }
                else
                    sumCSData->AddData(csDist[i]);
            }
        }
    }
    else
    {
        for(int i=8; i<numProcess; i++)
        {
            if(csDist[i])
            {
                if(first)
                {
                    sumCSData = new CSDist1DTab(dynamic_cast<CSDist1DTab*>(csDist[i]));
                    first=false;
                }
                else
                    sumCSData->AddData(csDist[i]);
            }
        }
    }

    if(sumCSData)
    {
        sumCSData->WriteG4NDLCSData(stream);
        SetDataStream( outDirNameCSProc[3]+isoName, stream, ascii);
        delete sumCSData;
    }

    stream.clear();
    stream.str("");

    if(isoName=="1_1_Hydrogen")
    {
        stream << std::setw(14) << std::right << 0 << '\n';
        stream << std::setw(14) << std::right << '0' << '\n';
        stream << std::setw(14) << std::right << '1' << '\n';
        stream << std::setw(14) << std::right << 2.000000e+07 << std::setw(14) << std::right << 0.0 << '\n';

        SetDataStream( outDirNameCSProc[3]+isoName, stream, ascii);

        stream.clear();
        stream.str("");
    }
}

void MakeElasFSFile(string outDirName, string isoName, double isoMass, int numAngEner[], vector<AngularDist*> *angDist, bool ascii)
{
    //Creates the elastic scattering final state files

    std::stringstream stream;
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    int frameFlag=2, repFlag=2; //assume that elastic data is collected in CM frame, may need to compare values to G4NDL to verify that this is correct

    //this section outputs the out-going neutron angular distribution data for the elastic scattering reaction
    stream.fill(' ');
    stream << std::setw(14) << std::right << repFlag;
    stream << std::setw(14) << std::right << isoMass;
    stream << std::setw(14) << std::right << frameFlag;
    stream << std::setw(14) << std::right << numAngEner[0] << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << numAngEner[0] << std::setw(14) << std::right << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

    for(int i=0; i<int(angDist[0].size()); i++)
    {
        angDist[0][i]->WriteG4NDLData(stream);
    }

    stream << '\n';

    outDirName = outDirName+"Elastic/FS/";

    if(!(DirectoryExists((outDirName).c_str())))
    {
        system( ("mkdir -p -m=666 "+outDirName).c_str());
        if(DirectoryExists((outDirName).c_str()))
        {
            cout << "created directory " << outDirName << "\n" << endl;
        }
        else
        {
            cout << "\nError: could not create directory " << outDirName << "\n" << endl;
            return;
        }
    }

    SetDataStream( outDirName+isoName, stream, ascii);

    stream.clear();
    stream.str("");
}

void MakeCaptureFSFile(string outDirName, string isoName, double isoMass, double temperature, int TYRList[], CSDist **pCSVec, vector<AngularDist*> *angDistP,
                        vector<int> *enDisLawP, vector<EnergyDist*> *enerDistP, vector<int> MTRPList, int *extractOrderP, int *numAngEnerP,
                        vector<int> *enDisNumLawApplNRegP, vector<int> *enDisNumLawApplNEnP, vector<int*> *enDisSchemeVecP, vector<int*> *enDisRangeVecP,
                        vector<double*> *enDisEnApplVecP, vector<double*> *enDisProbApplVecP, double **energyAngVecP, CSDist *nCSVec, bool ascii)
{
    //Creates the neutron capture final state files

    std::stringstream stream;
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    //int frameFlag=2;

    vector<int> capPReacIndex;

    //check the reference frame that the data has been gathered from
    //if(TYRList[1]>0)
        //frameFlag=1;

    stream.fill(' ');

    //create list of relevant photon production reactions
    for(int i=0; i<int(MTRPList.size()); i++)
    {
        if(int(MTRPList[i]/1000)==102)
        {
            capPReacIndex.push_back(i);
        }
    }

    // out going photon multiplicity
    // I use repFlag 1 no matter what, since we do not collect the partials
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << capPReacIndex.size() << '\n';
    for(int i=0; i<int(capPReacIndex.size()); i++)
    {
        //we always select a continous energy dist
        stream << std::setw(14) << std::right << 1;
        //This average energy is not used when the previous value is set to 1 (ie when there is an energy distribution)
        stream << std::setw(14) << std::right << enerDistP[capPReacIndex[i]].front()->GetAverageOutEnergy() << '\n';

        pCSVec[capPReacIndex[i]]->WriteG4NDLYieldData(stream);
    }
    stream << '\n';
    stream << '\n';

    // photon angular distribution
    // set repflag to 1 instead of 2 since we do not collect the partials
    stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << std::setw(14) << std::right << min(1,int(capPReacIndex.size())) << std::setw(14) << std::right << 0 << endl;

    if(capPReacIndex.size()>0)
    {
        vector<AngularDist*> angDistPVec;
        vector<CSDist*> pCSVecTemp;
        for(int j=0; j<int(capPReacIndex.size()); j++)
        {
            for(int k=0; k<int(angDistP[capPReacIndex[j]].size()); k++)
            {
                if(!(angDistP[capPReacIndex[j]][k]->CheckData()))
                {
                    cout << "Error in angular data ConvertMCNPtoG4NDL.cc:2464" << endl;
                }
                angDistPVec.push_back(angDistP[capPReacIndex[j]][k]);
                pCSVecTemp.push_back(pCSVec[capPReacIndex[j]]);
            }
        }

        int tempNumAngEner;
        AngularDist *tempAng = new AngDist2DTabularP();
        tempAng->SumAngularData(angDistPVec, pCSVecTemp, tempNumAngEner);
        if(!(tempAng->CheckData()))
        {
            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:2476" << endl;
        }
        stream << std::setw(14) << std::right << enerDistP[capPReacIndex[0]].front()->GetAverageOutEnergy() << std::setw(14) << std::right << 0 << endl;
        stream << std::setw(14) << std::right << tempNumAngEner << '\n';
        tempAng->WriteG4NDLData(stream);
        delete tempAng;

        stream << '\n';
    }

    stream << '\n';

    //photon energy distribution
    //Here we adjust the tables defining the probability of an energy dist applicablitity to be normalized to 1 and then weighted by the cross-section
    //this is an approximation that lets us put all the energy dist for every photon production reaction caused by this neutron production reaction
    //into one big list. the normalization ensures that the probability sum between photon production reactions is the same and the weighting of the cross-sections
    //ensures that the energy distributions from the more probable photon production reaction are more likely to be applied
    int numPartials=0, low, reg, enApplCount;
    double sum, energy;
    for(int i=0; i<int(capPReacIndex.size()); i++)
    {
        sum=0.;
        enApplCount=0;
        //put in while loop to guard against the case that there is an incoming energy where all of the distribution probabilities are zero causing the sum to be zero and the probabilities to be unnormalized
        while((sum==0.)&&(enDisNumLawApplNEnP[capPReacIndex[i]][0]>enApplCount))
        {
            energy=enDisEnApplVecP[capPReacIndex[i]][0][enApplCount];
            for(int j=0; j<int(enDisLawP[capPReacIndex[i]].size()); j++)
            {
                if(enDisLawP[capPReacIndex[i]][j]==2||enDisLawP[capPReacIndex[i]][j]==4)
                {
                    reg=0;
                    for(low=0; low<int(enDisNumLawApplNEnP[capPReacIndex[i]][j]-1); low++)
                    {
                        while((enDisRangeVecP[capPReacIndex[i]][j][reg]<=low)&&(enDisNumLawApplNRegP[capPReacIndex[i]][j]-1>reg))
                            reg++;
                        if(energy<enDisEnApplVecP[capPReacIndex[i]][j][low])
                        {
                            break;
                        }
                    }
                    low--;
                    if(low<0)
                        low=0;

                    if(enDisNumLawApplNEnP[capPReacIndex[i]][j]>1)
                        sum+=max(0.,Interpolate(enDisSchemeVecP[capPReacIndex[i]][j][reg], energy, enDisEnApplVecP[capPReacIndex[i]][j][low], enDisEnApplVecP[capPReacIndex[i]][j][low+1],
                                    enDisProbApplVecP[capPReacIndex[i]][j][low], enDisProbApplVecP[capPReacIndex[i]][j][low+1]));
                    else
                        sum+=enDisProbApplVecP[capPReacIndex[i]][j][0];
                }
            }
            enApplCount++;
        }
        for(int k=0; k<int(enDisLawP[capPReacIndex[i]].size()); k++)
        {
            if(enDisLawP[capPReacIndex[i]][k]==2||enDisLawP[capPReacIndex[i]][k]==4)
            {
                numPartials++;
            }
        }
        if(sum!=0.)
        {
            for(int j=0; j<int(enDisLawP[capPReacIndex[i]].size()); j++)
            {
                if(enDisLawP[capPReacIndex[i]][j]==2||enDisLawP[capPReacIndex[i]][j]==4)
                {
                    for(low=0; low<int(enDisNumLawApplNEnP[capPReacIndex[i]][j]); low++)
                    {
                        enDisProbApplVecP[capPReacIndex[i]][j][low]=enDisProbApplVecP[capPReacIndex[i]][j][low]*max(0.,pCSVec[capPReacIndex[i]]->GetAvgCS())/sum;
                    }
                }
            }
        }
    }

    // make sure that the probability of one of the reactions energy distributions is greater than zero for every possible incoming energy
    for(int i=0; i<int(capPReacIndex.size()); i++)
    {
        for(int j=0; j<int(enDisLawP[capPReacIndex[i]].size()); j++)
        {
            for(low=0; low<int(enDisNumLawApplNEnP[capPReacIndex[i]][j]); low++)
            {
                energy=enDisEnApplVecP[capPReacIndex[i]][j][low];
                sum=0.;
                for(int k=0; k<int(enDisLawP[capPReacIndex[i]].size()); k++)
                {
                    int m=0;
                    for(; m<int(enDisNumLawApplNEnP[capPReacIndex[i]][k]-1); m++)
                    {
                        if(enDisEnApplVecP[capPReacIndex[i]][k][m]>energy)
                        {
                            break;
                        }
                    }
                    if(m!=0)
                        m--;
                    if(enDisNumLawApplNEnP[capPReacIndex[i]][k]>1)
                        sum+=max(0.,Interpolate(2, energy, enDisEnApplVecP[capPReacIndex[i]][k][m], enDisEnApplVecP[capPReacIndex[i]][k][m+1],
                                    enDisProbApplVecP[capPReacIndex[i]][k][m], enDisProbApplVecP[capPReacIndex[i]][k][m+1]));
                    else
                        sum+=enDisProbApplVecP[capPReacIndex[i]][k][0];
                }
                if(sum==0.)
                {
                    enDisProbApplVecP[capPReacIndex[i]][j][low]=1;
                }
            }
        }
    }

    if(capPReacIndex.size()>0)
    {
        stream << std::setw(14) << std::right << numPartials << endl;
        for(int i=0; i<int(capPReacIndex.size()); i++)
        {
            for(int j=0, count=0; j<int(enDisLawP[capPReacIndex[i]].size()); j++, count++)
            {
                if(enDisLawP[capPReacIndex[i]][j]==2||enDisLawP[capPReacIndex[i]][j]==4)
                {
                    stream << std::setw(14) << std::right << 0 << endl;
                    stream << std::setw(14) << std::right << enDisNumLawApplNEnP[capPReacIndex[i]][j] <<'\n';
                    stream << std::setw(14) << std::right << enDisNumLawApplNRegP[capPReacIndex[i]][j] <<'\n';
                    for(int k=0; k<enDisNumLawApplNRegP[capPReacIndex[i]][j]; k++)
                    {
                        stream << std::setw(14) << std::right << enDisRangeVecP[capPReacIndex[i]][j][k] << std::setw(14) << std::right << enDisSchemeVecP[capPReacIndex[i]][j][k] << '\n';
                    }
                    for(int k=0; k<enDisNumLawApplNEnP[capPReacIndex[i]][j]; k++)
                    {
                        stream << std::setw(14) << std::right << enDisEnApplVecP[capPReacIndex[i]][j][k]*1000000 << std::setw(14) << std::right << enDisProbApplVecP[capPReacIndex[i]][j][k] << '\n';
                    }
                    stream << '\n';
                    enerDistP[capPReacIndex[i]][count]->WriteG4NDLData(stream);
                }
                else
                    count--;
            }
        }
        stream << '\n';
    }
    stream << '\n';

    string outDir = outDirName+"Capture/FS/";
    string fileName = outDirName+"Capture/FS/"+isoName;

    if(!(DirectoryExists((outDir).c_str())))
    {
        system( ("mkdir -p -m=666 "+outDir).c_str());
        if(DirectoryExists((outDir).c_str()))
        {
            cout << "created directory " << outDir << "\n" << endl;
        }
        else
        {
            cout << "\nError: could not create directory " << outDir << "\n" << endl;
            return;
        }
    }
    SetDataStream( fileName, stream, ascii);

    stream.clear();
    stream.str("");

    // we do not make the FSMF6 directory since the MCNP does not seem to contain 3D tables of energy, angle and probability for photon production
    // we must check to make sure that this is true
    /*
    stream << isoMass << frameFlag << nDiscrete << endl;

    for(int i=0; i<nDiscrete; i++)
    {
        stream << std::setw(14) << std::right << (isoMass*1.00867) << std::setw(14) << std::right << 0 << 0 << endl;
    }
    */

}

void CreateMT2(int* MTRListPos, string outDirName, string isoName, int isoNum, double isoMass, double temperature, int *MTRList, CSDist* nCSVec[],
                        vector<int>* enDisLaw, vector<int>* enDisNumLawApplNReg, vector<int>* enDisNumLawApplNEn, vector<int*>* enDisSchemeVec,
                        vector<int*>* enDisRangeVec, vector<double*>* enDisEnApplVec, vector<double*>* enDisProbApplVec, vector<EnergyDist*>* enerDist,
                        YieldDist* nYieldReac[], double* reacQValue, int* numAngEner, vector<AngularDist*>* angDist, bool* angDistInEnDistFlag,
                        vector<AngularEnergyDist*>* angEnDist, CSDist*** pCSVec, int *TYRList, int** numAngEnerP, vector<AngularDist*>** angDistP, vector<int>** enDisLawP,
                        vector<int>** enDisNumLawApplNRegP, vector<int>** enDisNumLawApplNEnP, vector<int*>** enDisSchemeVecP, vector<int*>** enDisRangeVecP,
                        vector<double*>** enDisEnApplVecP, vector<double*>** enDisProbApplVecP,vector<EnergyDist*>** enerDistP, vector<int> &MTRPList,
                        double*** energyAngVecP, bool ascii)
{
    bool first=true;
    CSDist1DTab* sumCSData = NULL;
    for(int i=3; i<7; i++)
    {
        if(nCSVec[i])
        {
            if(first)
            {
                sumCSData = new CSDist1DTab(dynamic_cast<CSDist1DTab*>(nCSVec[i]));
                first=false;
            }
            else
                sumCSData->AddData(nCSVec[i]);
        }
    }
    if(sumCSData)
    {
        TYRList[2] = 19;
        MTRListPos[2]=1;
        nCSVec[2] = sumCSData;
        numAngEner[2] = 0;

        for(int i=3; i<7; i++)
        {
            if(angDistInEnDistFlag[i])
            {
                for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                {
                    if(enDisLaw[i][j]>24)
                    {
                        enDisLaw[2].push_back(4);
                        enDisNumLawApplNEn[2].push_back(enDisNumLawApplNEn[i][j]);
                        enDisNumLawApplNReg[2].push_back(enDisNumLawApplNReg[i][j]);
                        enDisRangeVec[2].push_back(new int [enDisNumLawApplNReg[i][j]]);
                        enDisSchemeVec[2].push_back(new int [enDisNumLawApplNReg[i][j]]);
                        for(int k=0; k<enDisNumLawApplNReg[i][j]; k++)
                        {
                            enDisRangeVec[2].back()[k]=enDisRangeVec[i][j][k];
                            enDisSchemeVec[2].back()[k]=enDisSchemeVec[i][j][k];
                        }

                        enDisEnApplVec[2].push_back(new double [enDisNumLawApplNEn[i][j]]);
                        enDisProbApplVec[2].push_back(new double [enDisNumLawApplNEn[i][j]]);
                        for(int k=0; k<enDisNumLawApplNEn[i][j]; k++)
                        {
                            enDisEnApplVec[2].back()[k]=enDisEnApplVec[i][j][k];
                            enDisProbApplVec[2].back()[k]=enDisProbApplVec[i][j][k];
                        }
                        AngularDist *temAngDist = NULL;
                        EnergyDist *tempEnDist = NULL;
                        angDist[2].push_back(temAngDist);
                        enerDist[2].push_back(tempEnDist);
                        angEnDist[i][count]->ConvertToEnerAndAngDist(&(enerDist[2].back()), &(angDist[2].back()), numAngEner[2]);
                        if(!(angDist[2].back()->CheckData()))
                        {
                            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:2714" << endl;
                        }
                    }
                    else
                        count--;
                }
            }
        }

        angDistInEnDistFlag[2]=false;

        //Angular Distribution
        angDist[2].push_back(new AngDist2DTabular());
        angDist[2].back()->SumAngularData(angDist, nCSVec, 3, 7, numAngEner[2]);
        if(!(angDist[2].back()->CheckData()))
        {
            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:2730" << endl;
        }
        angDist[2].back()->SetTemperature(temperature);

        //Energy Distribution of prompt neutron
        int numPartials=0, low, reg, enApplCount;
        double sum, energy;
        for(int j=3; j<7; j++)
        {
            sum=0.;
            enApplCount=0;
            while((enDisNumLawApplNEn[j].size()>0)&&(sum==0.)&&(enDisNumLawApplNEn[j][0]>enApplCount))
            {
                energy=enDisEnApplVec[j][0][enApplCount];
                for(int k=0; k<int(enDisLaw[j].size()); k++)
                {
                    if(enDisLaw[j][k]<44)
                    {
                        reg=0;
                        for(low=0; low<int(enDisNumLawApplNEn[j][k]-1); low++)
                        {
                            while((enDisRangeVec[j][k][reg]<=low)&&(enDisNumLawApplNReg[j][k]-1>reg))
                                reg++;
                            if(energy<enDisEnApplVec[j][k][low])
                            {
                                break;
                            }
                        }
                        low--;
                        if(low<0)
                            low=0;

                        if(enDisNumLawApplNEn[j][k]>1)
                            sum+=max(0.,Interpolate(enDisSchemeVec[j][k][reg], energy, enDisEnApplVec[j][k][low], enDisEnApplVec[j][k][low+1],
                                        enDisProbApplVec[j][k][low], enDisProbApplVec[j][k][low+1]));
                        else
                            sum+=enDisProbApplVec[j][k][0];
                    }
                }
                enApplCount++;
            }
            for(int k=0; k<int(enDisLaw[j].size()); k++)
            {
                if(enDisLaw[j][k]<44)
                {
                    numPartials++;
                }
            }
            if(sum!=0.)
            {
                for(int k=0; k<int(enDisLaw[j].size()); k++)
                {
                    if(enDisLaw[j][k]<44)
                    {
                        for(low=0; low<int(enDisNumLawApplNEn[j][k]); low++)
                        {
                            enDisProbApplVec[j][k][low]=enDisProbApplVec[j][k][low]*max(0.,nCSVec[j]->GetAvgCS())/sum;
                        }
                    }
                }
            }
        }

        // make sure that the probability of one of the reactions energy distributions is greater than zero for every possible incoming energy
        for(int l=3; l<7; l++)
        {
            for(int j=0; j<int(enDisLaw[l].size()); j++)
            {
                for(low=0; low<int(enDisNumLawApplNEn[l][j]); low++)
                {
                    energy=enDisEnApplVec[l][j][low];
                    sum=0.;
                    for(int k=0; k<int(enDisLaw[l].size()); k++)
                    {
                        int m=0;
                        for(; m<int(enDisNumLawApplNEn[l][k]-1); m++)
                        {
                            if(enDisEnApplVec[l][k][m]>energy)
                            {
                                break;
                            }
                        }
                        if(m!=0)
                            m--;
                        if(enDisNumLawApplNEn[l][k]>1)
                            sum+=max(0.,Interpolate(2, energy, enDisEnApplVec[l][k][m], enDisEnApplVec[l][k][m+1],
                                        enDisProbApplVec[l][k][m], enDisProbApplVec[l][k][m+1]));
                        else
                            sum+=enDisProbApplVec[l][k][0];
                    }
                    if(sum==0.)
                    {
                        enDisProbApplVec[l][j][low]=1;
                    }
                }
            }
        }

        for(int j=3; j<7; j++)
        {
            for(int k=0; k<int(enDisLaw[j].size()); k++)
            {
                if(enDisLaw[j][k]<44)
                {
                    enDisNumLawApplNReg[2].push_back(enDisNumLawApplNReg[j][k]);
                    enDisNumLawApplNEn[2].push_back(enDisNumLawApplNEn[j][k]);
                    enDisRangeVec[2].push_back(new int[enDisNumLawApplNReg[j][k]]);
                    for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                    {
                        enDisRangeVec[2].back()[i]=enDisRangeVec[j][k][i];
                    }
                    enDisSchemeVec[2].push_back(new int[enDisNumLawApplNReg[j][k]]);
                    for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                    {
                        enDisSchemeVec[2].back()[i]=enDisSchemeVec[j][k][i];
                    }
                    enDisEnApplVec[2].push_back(new double[enDisNumLawApplNEn[j][k]]);
                    for(int i=0; i<enDisNumLawApplNEn[j][k]; i++)
                    {
                        enDisEnApplVec[2].back()[i]=enDisEnApplVec[j][k][i];
                    }
                    enDisProbApplVec[2].push_back(new double [enDisNumLawApplNEn[j][k]]);
                    for(low=0; low<int(enDisNumLawApplNEn[j][k]); low++)
                    {
                        (enDisProbApplVec[2].back())[low]=enDisProbApplVec[j][k][low];
                    }
                }
            }

        }

        int countProc=0;
        reacQValue[2]=0.;
        for(int i=3; i<7; i++)
        {
            // this is only used when a mean value from the energy dist cannot be extracted, we use the mean value here as they do in GEANT4
            reacQValue[2]+=reacQValue[i];
            if(reacQValue[i]!=0.)
                countProc++;
            for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
            {
                if(enDisLaw[i][j]<44)
                {
                    enDisLaw[2].push_back(enDisLaw[i][j]);

                    enerDist[2].push_back(enerDist[i][count]);
                }
                else
                    count--;
            }
        }
        reacQValue[2]/=countProc;

        // neutron energy-angular distribution
        // not necessary to convert since we don't use the  angular-energy dist in the fission FS files but leave in this code in case we ant to use it later
        /*
        if(angDistInEnDistFlag[2])
        {
            //this first part, adjusting the law applicability, is not needed or used so far but we keep it incase someone wants to use it later
            for(int j=3; j<7; j++)
            {
                sum=0.;
                energy=0.;
                first=true;
                for(int k=0; k<int(enDisLaw[j].size()); k++)
                {
                    if(enDisLaw[j][k]>24)
                    {
                        reg=0;
                        if(first)
                        {
                            energy=enDisEnApplVec[j][k][int((enDisNumLawApplNEn[j][k])/2)];
                            first=false;
                        }
                        numPartials++;
                        for(low=0; low<int(enDisNumLawApplNEn[j][k]-1); low++)
                        {
                            while((enDisRangeVec[j][k][reg]<=low)&&(enDisNumLawApplNReg[j][k]-1>reg))
                                reg++;
                            if(energy<enDisEnApplVec[j][k][low])
                            {
                                break;
                            }
                        }
                        low--;
                        if(low<0)
                            low=0;
                        sum+=max(0.,Interpolate(enDisSchemeVec[j][k][reg], energy, enDisEnApplVec[j][k][low], enDisEnApplVec[j][k][low+1],
                                        enDisProbApplVec[j][k][low], enDisProbApplVec[j][k][low+1]));
                    }
                }

                for(int l=0; l<int(enDisLaw[j].size()); l++)
                {
                    for(low=0; low<int(enDisNumLawApplNEn[j][l]); low++)
                    {
                        energy=enDisEnApplVec[j][l][low];
                        sum=0.;
                        for(int k=0; k<int(enDisLaw[j].size()); k++)
                        {
                            int m=0;
                            for(; m<int(enDisNumLawApplNEn[j][k]-1); m++)
                            {
                                if(enDisEnApplVec[j][k][m]>energy)
                                {
                                    break;
                                }
                            }
                            if(m!=0)
                                m--;
                            if(enDisNumLawApplNEn[j][k]>1)
                                sum+=max(0.,Interpolate(2, energy, enDisEnApplVec[j][k][m], enDisEnApplVec[j][k][m+1],
                                            enDisProbApplVec[j][k][m], enDisProbApplVec[j][k][m+1]));
                            else
                                sum+=enDisProbApplVec[j][k][0];
                        }
                        if(sum==0.)
                        {
                            enDisProbApplVec[j][l][low]=1;
                        }
                    }
                }

                for(int k=0; k<int(enDisLaw[j].size()); k++)
                {
                    if(enDisLaw[j][k]<24)
                    {
                        enDisNumLawApplNReg[2].push_back(enDisNumLawApplNReg[j][k]);
                        enDisNumLawApplNEn[2].push_back(enDisNumLawApplNEn[j][k]);
                        enDisRangeVec[2].push_back(new int[enDisNumLawApplNReg[j][k]]);
                        for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                        {
                            enDisRangeVec[2].back()[i]=enDisRangeVec[j][k][i];
                        }
                        enDisSchemeVec[2].push_back(new int[enDisNumLawApplNReg[j][k]]);
                        for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                        {
                            enDisSchemeVec[2].back()[i]=enDisSchemeVec[j][k][i];
                        }
                        enDisEnApplVec[2].push_back(new double[enDisNumLawApplNEn[j][k]]);
                        for(int i=0; i<enDisNumLawApplNEn[j][k]; i++)
                        {
                            enDisEnApplVec[2].back()[i]=enDisEnApplVec[j][k][i];
                        }
                        enDisProbApplVec[2].push_back(new double [enDisNumLawApplNEn[j][k]]);
                        for(low=0; low<int(enDisNumLawApplNEn[j][k]); low++)
                        {
                            if(sum!=0.)
                                (enDisProbApplVec[2].back())[low]=enDisProbApplVec[j][k][low]*max(0.,nCSVec[j]->GetAvgCS())/sum;
                            else
                                (enDisProbApplVec[2].back())[low]=0.;
                        }
                    }
                }

            }
        }
        */

        int pReacCount=0;
        for(int i=3; i<7; i++)
        {
            for(int j=0; j<int(MTRPList.size()); j++)
            {
                if(int(MTRPList[j]/1000)==MTRList[i])
                {
                    pReacCount++;
                }
            }
        }

        CSDist **pCSVecTemp = new CSDist* [MTRPList.size()+pReacCount];
        int *numAngEnerPTemp = new int [MTRPList.size()+pReacCount];
        vector<AngularDist*> *angDistPTemp = new vector<AngularDist*> [MTRPList.size()+pReacCount];
        vector<int> *enDisLawPTemp = new vector<int> [MTRPList.size()+pReacCount];
        vector<int> *enDisNumLawApplNRegPTemp = new vector<int> [MTRPList.size()+pReacCount];
        vector<int> *enDisNumLawApplNEnPTemp = new vector<int> [MTRPList.size()+pReacCount];
        vector<int*> *enDisSchemeVecPTemp = new vector<int*> [MTRPList.size()+pReacCount];
        vector<int*> *enDisRangeVecPTemp = new vector<int*> [MTRPList.size()+pReacCount];
        vector<double*> *enDisEnApplVecPTemp = new vector<double*> [MTRPList.size()+pReacCount];
        vector<double*> *enDisProbApplVecPTemp = new vector<double*> [MTRPList.size()+pReacCount];
        vector<EnergyDist*> *enerDistPTemp = new vector<EnergyDist*> [MTRPList.size()+pReacCount];
        double **energyAngVecPTemp = new double* [MTRPList.size()+pReacCount];

        for(int j=0; j<int(MTRPList.size()); j++)
        {
            pCSVecTemp[j] = pCSVec[0][j];
            numAngEnerPTemp[j] = numAngEnerP[0][j];
            angDistPTemp[j] = angDistP[0][j];
            enDisLawPTemp[j] = enDisLawP[0][j];
            enDisNumLawApplNRegPTemp[j] = enDisNumLawApplNRegP[0][j];
            enDisNumLawApplNEnPTemp[j] = enDisNumLawApplNEnP[0][j];
            enDisSchemeVecPTemp[j] = enDisSchemeVecP[0][j];
            enDisRangeVecPTemp[j] = enDisRangeVecP[0][j];
            enDisEnApplVecPTemp[j] = enDisEnApplVecP[0][j];
            enDisProbApplVecPTemp[j] = enDisProbApplVecP[0][j];
            enerDistPTemp[j] = enerDistP[0][j];
            energyAngVecPTemp[j] = energyAngVecP[0][j];
        }

        int pReacCount2=0;
        for(int i=3; i<7; i++)
        {
            for(int j=0; j<int(MTRPList.size()); j++)
            {
                if(int(MTRPList[j]/1000)==MTRList[i])
                {
                    pCSVecTemp[MTRPList.size()+pReacCount2] = pCSVec[0][j];
                    numAngEnerPTemp[MTRPList.size()+pReacCount2] = numAngEnerP[0][j];
                    angDistPTemp[MTRPList.size()+pReacCount2] = angDistP[0][j];
                    enDisLawPTemp[MTRPList.size()+pReacCount2] = enDisLawP[0][j];
                    enDisNumLawApplNRegPTemp[MTRPList.size()+pReacCount2] = enDisNumLawApplNRegP[0][j];
                    enDisNumLawApplNEnPTemp[MTRPList.size()+pReacCount2] = enDisNumLawApplNEnP[0][j];
                    enDisSchemeVecPTemp[MTRPList.size()+pReacCount2] = enDisSchemeVecP[0][j];
                    enDisRangeVecPTemp[MTRPList.size()+pReacCount2] = enDisRangeVecP[0][j];
                    enDisEnApplVecPTemp[MTRPList.size()+pReacCount2] = enDisEnApplVecP[0][j];
                    enDisProbApplVecPTemp[MTRPList.size()+pReacCount2] = enDisProbApplVecP[0][j];
                    enerDistPTemp[MTRPList.size()+pReacCount2] = enerDistP[0][j];
                    energyAngVecPTemp[MTRPList.size()+pReacCount2] = energyAngVecP[0][j];
                    pReacCount2++;
                }
            }
        }

        //delete out going photon data

        if(enDisLawP[0]!=NULL)
            delete [] enDisLawP[0];
        if(enDisNumLawApplNRegP[0]!=NULL)
            delete [] enDisNumLawApplNRegP[0];
        if(enDisNumLawApplNEnP[0]!=NULL)
            delete [] enDisNumLawApplNEnP[0];
        if(energyAngVecP[0]!=NULL)
            delete [] energyAngVecP[0];
        if(angDistP[0]!=NULL)
            delete [] angDistP[0];
        if(pCSVec[0]!=NULL)
            delete[] pCSVec[0];
        if(enDisSchemeVecP[0]!=NULL)
            delete [] enDisSchemeVecP[0];
        if(enDisRangeVecP[0]!=NULL)
            delete [] enDisRangeVecP[0];
        if(enDisProbApplVecP[0]!=NULL)
            delete [] enDisProbApplVecP[0];
        if(enDisEnApplVecP[0]!=NULL)
            delete [] enDisEnApplVecP[0];
        if(enerDistP[0]!=NULL)
            delete [] enerDistP[0];

        for(int j=0; j<pReacCount; j++)
        {
            MTRPList.push_back(18000+j);
        }

        pCSVec[0] = pCSVecTemp;
        numAngEnerP[0] = numAngEnerPTemp;
        angDistP[0] = angDistPTemp;
        enDisLawP[0] = enDisLawPTemp;
        enDisNumLawApplNRegP[0] = enDisNumLawApplNRegPTemp;
        enDisNumLawApplNEnP[0] = enDisNumLawApplNEnPTemp;
        enDisSchemeVecP[0] = enDisSchemeVecPTemp;
        enDisRangeVecP[0] = enDisRangeVecPTemp;
        enDisEnApplVecP[0] = enDisEnApplVecPTemp;
        enDisProbApplVecP[0] = enDisProbApplVecPTemp;
        enerDistP[0] = enerDistPTemp;
        energyAngVecP[0] = energyAngVecPTemp;

    }
}

void MakeFissionFSFile(int *MTRList, int *MTRListPos, string outDirName, string isoName, double isoMass, double temperature, CSDist **nCSVec, vector<int> *enDisLaw, vector<int> *enDisNumLawApplNReg, vector<int> *enDisNumLawApplNEn,
                        vector<int*> *enDisSchemeVec, vector<int*> *enDisRangeVec, vector<double*> *enDisEnApplVec, vector<double*> *enDisProbApplVec, vector<EnergyDist*> *enerDist,  vector<AngularEnergyDist*> *angEnDist,
                        vector<int> enDisLawND, vector<int> enDisNumLawApplNRegND, vector<int> enDisNumLawApplNEnND, vector<int*> enDisSchemeVecND, vector<int*> enDisRangeVecND,
                        vector<double*> enDisEnApplVecND, vector<double*> enDisProbApplVecND, vector<EnergyDist*> enerDistND, double *reacQValue, int *numAngEner, vector<AngularDist*> *angDist,
                        bool *angDistInEnDistFlag, CSDist **pCSVec, int *TYRList, bool promptYieldFlag, bool totalYieldFlag, YieldDist *dNPromptYieldDist, YieldDist *dNTotalYieldDist,
                        YieldDist *dNYield, NDelayConstDist *nDelConst, int *numAngEnerP, vector<AngularDist*> *angDistP, vector<int> *enDisLawP, vector<int> *enDisNumLawApplNRegP, vector<int> *enDisNumLawApplNEnP,
                        vector<int*> *enDisSchemeVecP, vector<int*> *enDisRangeVecP, vector<double*> *enDisEnApplVecP, vector<double*> *enDisProbApplVecP, vector<EnergyDist*> *enerDistP,
                        vector<int> MTRPList, double **energyAngVecP, bool ascii)
{
    //Creates the fission final state files

    std::stringstream stream;
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    int frameFlag=2, repFlag=2;

    vector<int> fisPReacIndex;

    //check the reference frame that the data has been gathered from
    if(TYRList[2]>0)
        frameFlag=1;

    stream.fill(' ');

    //create list of relevant photon production reactions
    for(int i=0; i<int(MTRPList.size()); i++)
    {
        if(int(MTRPList[i]/1000)==18)
        {
            fisPReacIndex.push_back(i);
        }
    }

    // FS/ directory

    //Convert Energy-Angular Dist to Energy and Angular Dist
    if(angDistInEnDistFlag[2])
    {
        if(numAngEner[2]<0)
            numAngEner[2]=0;
        for(int i=0, count=0; i<int(enDisLaw[2].size()); i++, count++)
        {
            if(enDisLaw[2][i]>24)
            {
                enDisLaw[2].push_back(4);
                enDisNumLawApplNEn[2].push_back(enDisNumLawApplNEn[2][i]);
                enDisNumLawApplNReg[2].push_back(enDisNumLawApplNReg[2][i]);
                enDisRangeVec[2].push_back(new int [enDisNumLawApplNReg[2][i]]);
                enDisSchemeVec[2].push_back(new int [enDisNumLawApplNReg[2][i]]);
                for(int j=0; j<enDisNumLawApplNReg[2][i]; j++)
                {
                    enDisRangeVec[2].back()[j]=enDisRangeVec[2][i][j];
                    enDisSchemeVec[2].back()[j]=enDisSchemeVec[2][i][j];
                }

                enDisEnApplVec[2].push_back(new double [enDisNumLawApplNEn[2][i]]);
                enDisProbApplVec[2].push_back(new double [enDisNumLawApplNEn[2][i]]);
                for(int j=0; j<enDisNumLawApplNEn[2][i]; j++)
                {
                    enDisEnApplVec[2].back()[j]=enDisEnApplVec[2][i][j];
                    enDisProbApplVec[2].back()[j]=enDisProbApplVec[2][i][j];
                }
                AngularDist *temAngDist = NULL;
                EnergyDist *tempEnDist = NULL;
                angDist[2].push_back(temAngDist);
                enerDist[2].push_back(tempEnDist);
                angEnDist[2][count]->ConvertToEnerAndAngDist(&(enerDist[2].back()), &(angDist[2].back()), numAngEner[2]);
                if(!(angDist[2].back()->CheckData()))
                {
                    cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3195" << endl;
                }
            }
            else
                count--;
        }
    }

    //prompt neutrons angular Distribution
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 4 << '\n';
    stream << std::setw(14) << std::right << repFlag;
    stream << std::setw(14) << std::right << isoMass;
    stream << std::setw(14) << std::right << frameFlag;
    stream << std::setw(14) << std::right << numAngEner[2] << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << numAngEner[2] << std::setw(14) << std::right << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

    //note angDist[0].size() does not necessarilly equal numAngEner[0] since if the distribution type is the same between adjacent
    //incoming energy points we don't add another object to angDist, instead we add the point to the previous object
    for(int i=0; i<int(angDist[2].size()); i++)
    {
        if(!(angDist[2][i]->CheckData()))
        {
            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3218" << endl;
        }
        angDist[2][i]->WriteG4NDLData(stream);
    }
    stream << '\n' << endl;

    //Energy Distribution of prompt neutron
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 5 << '\n';
    stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << enerDist[2].size() << endl;

    for(int i=0, count=0; i<int(enDisLaw[2].size()); i++, count++)
    {
        if(enDisLaw[2][i]<44)
        {
            if(enDisLaw[2][i]==1)
                stream << std::setw(14) << std::right << 1 << '\n';
            else if(enDisLaw[2][i]==3)
            {
                stream << std::setw(14) << std::right << 1 << '\n';
            }
            else if(enDisLaw[2][i]==4)
            {
                stream << std::setw(14) << std::right << 1 << '\n';
            }
            else if(enDisLaw[2][i]==5)
            {
                stream << std::setw(14) << std::right << 5 << '\n';
            }
            else if(enDisLaw[2][i]==7)
            {
                stream << std::setw(14) << std::right << 7 << '\n';
            }
            else if(enDisLaw[2][i]==9)
            {
                stream << std::setw(14) << std::right << 9 << '\n';
            }
            else if(enDisLaw[2][i]==11)
            {
                stream << std::setw(14) << std::right << 11 << '\n';
            }
            else if(enDisLaw[2][i]==22)
            {
                cout << "###: No direct translation for this law!" << endl;
                stream << std::setw(14) << std::right << 1 << '\n';
            }
            else if(enDisLaw[2][i]==24)
            {
                cout << "###: No direct translation for this law!" << endl;
                stream << std::setw(14) << std::right << 1 << '\n';
            }

            stream << std::setw(14) << std::right << enDisNumLawApplNEn[2][i] <<'\n';
            stream << std::setw(14) << std::right << enDisNumLawApplNReg[2][i] <<'\n';
            for(int j=0; j<enDisNumLawApplNReg[2][i]; j++)
            {
                stream << std::setw(14) << std::right << enDisRangeVec[2][i][j] << std::setw(14) << std::right << enDisSchemeVec[2][i][j] << '\n';
            }
            for(int j=0; j<enDisNumLawApplNEn[2][i]; j++)
            {
                stream << std::setw(14) << std::right << enDisEnApplVec[2][i][j]*1000000 << std::setw(14) << std::right << enDisProbApplVec[2][i][j] << '\n';
            }

            enerDist[2][count]->WriteG4NDLData(stream);
        }
        else
            count--;
    }
    stream << '\n';
    stream << '\n';

    // photon production distribution
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 12 << '\n';
    // I use repFlag 1 no matter what
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << fisPReacIndex.size() << '\n';
    for(int i=0; i<int(fisPReacIndex.size()); i++)
    {
        //we always select a continous energy dist
        stream << std::setw(14) << std::right << 1;
        //This average energy is not used when the previous value is set to 1 (ie when there is an energy distribution)
        stream << std::setw(14) << std::right << enerDistP[fisPReacIndex[i]].front()->GetAverageOutEnergy() << '\n';

        pCSVec[fisPReacIndex[i]]->WriteG4NDLYieldData(stream);
    }
    stream << '\n';
    stream << '\n';

    // photon angular distribution
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 14 << '\n';
    stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << std::setw(14) << std::right << min(1,int(fisPReacIndex.size())) << std::setw(14) << std::right << 0 << endl;

    if(fisPReacIndex.size()>0)
    {
        vector<AngularDist*> angDistPVec;
        vector<CSDist*> pCSVecTemp;
        for(int j=0; j<int(fisPReacIndex.size()); j++)
        {
            for(int k=0; k<int(angDistP[fisPReacIndex[j]].size()); k++)
            {
                if(!(angDistP[fisPReacIndex[j]][k]->CheckData()))
                {
                    cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3318" << endl;
                }
                angDistPVec.push_back(angDistP[fisPReacIndex[j]][k]);
                pCSVecTemp.push_back(pCSVec[fisPReacIndex[j]]);
            }
        }

        int tempNumAngEner;
        AngularDist *tempAng = new AngDist2DTabularP();
        tempAng->SumAngularData(angDistPVec, pCSVecTemp, tempNumAngEner);
        if(!(tempAng->CheckData()))
        {
            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3330" << endl;
        }
        stream << std::setw(14) << std::right << enerDistP[fisPReacIndex[0]].front()->GetAverageOutEnergy() << std::setw(14) << std::right << 0 << endl;
        stream << std::setw(14) << std::right << tempNumAngEner << '\n';
        tempAng->WriteG4NDLData(stream);
        delete tempAng;

        stream << '\n';
    }

    stream << '\n';

    // photon energy distribution
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 15 << '\n';
    //Here we adjust the tables defining the probability of an energy dist applicablitity to be normalized to 1 and then weighted by the cross-section
    //this is an approximation that lets us put all the energy dist for every photon production reaction caused by this neutron production reaction
    //into one big list. the normalization ensures that the probability sum between photon production reactions is the same and the weighting of the cross-sections
    //ensures that the energy distributions from the more probable photon production reaction are more likely to be applied
    int numPartials=0, low, reg, enApplCount;
    double sum, energy;
    for(int i=0; i<int(fisPReacIndex.size()); i++)
    {
        enApplCount=0;
        sum=0.;
        while((sum==0.)&&(enDisNumLawApplNEnP[fisPReacIndex[i]][0]>enApplCount))
        {
            energy=enDisEnApplVecP[fisPReacIndex[i]][0][enApplCount];
            for(int j=0; j<int(enDisLawP[fisPReacIndex[i]].size()); j++)
            {
                if(enDisLawP[fisPReacIndex[i]][j]==2||enDisLawP[fisPReacIndex[i]][j]==4)
                {
                    reg=0;
                    for(low=0; low<int(enDisNumLawApplNEnP[fisPReacIndex[i]][j]-1); low++)
                    {
                        while((enDisRangeVecP[fisPReacIndex[i]][j][reg]<=low)&&(enDisNumLawApplNRegP[fisPReacIndex[i]][j]-1>reg))
                            reg++;
                        if(energy<enDisEnApplVecP[fisPReacIndex[i]][j][low])
                        {
                            break;
                        }
                    }
                    low--;
                    if(low<0)
                        low=0;

                    if(enDisNumLawApplNEnP[fisPReacIndex[i]][j]>1)
                        sum+=max(0.,Interpolate(enDisSchemeVecP[fisPReacIndex[i]][j][reg], energy, enDisEnApplVecP[fisPReacIndex[i]][j][low], enDisEnApplVecP[fisPReacIndex[i]][j][low+1],
                                    enDisProbApplVecP[fisPReacIndex[i]][j][low], enDisProbApplVecP[fisPReacIndex[i]][j][low+1]));
                    else
                        sum+=enDisProbApplVecP[fisPReacIndex[i]][j][0];
                }
            }
            enApplCount++;
        }
        for(int k=0; k<int(enDisLawP[fisPReacIndex[i]].size()); k++)
        {
            if(enDisLawP[fisPReacIndex[i]][k]==2||enDisLawP[fisPReacIndex[i]][k]==4)
            {
                numPartials++;
            }
        }
        if(sum!=0.)
        {
            for(int j=0; j<int(enDisLawP[fisPReacIndex[i]].size()); j++)
            {
                if(enDisLawP[fisPReacIndex[i]][j]==2||enDisLawP[fisPReacIndex[i]][j]==4)
                {
                    for(low=0; low<int(enDisNumLawApplNEnP[fisPReacIndex[i]][j]); low++)
                    {
                        enDisProbApplVecP[fisPReacIndex[i]][j][low]=enDisProbApplVecP[fisPReacIndex[i]][j][low]*(max(0.,pCSVec[fisPReacIndex[i]]->GetAvgCS()))/sum;
                    }
                }
            }
        }
    }

    // make sure that the probability of one of the reactions energy distributions is greater than zero for every possible incoming energy
    for(int i=0; i<int(fisPReacIndex.size()); i++)
    {
        for(int j=0; j<int(enDisLawP[fisPReacIndex[i]].size()); j++)
        {
            for(low=0; low<int(enDisNumLawApplNEnP[fisPReacIndex[i]][j]); low++)
            {
                energy=enDisEnApplVecP[fisPReacIndex[i]][j][low];
                sum=0.;
                for(int k=0; k<int(enDisLawP[fisPReacIndex[i]].size()); k++)
                {
                    int m=0;
                    for(; m<int(enDisNumLawApplNEnP[fisPReacIndex[i]][k]-1); m++)
                    {
                        if(enDisEnApplVecP[fisPReacIndex[i]][k][m]>energy)
                        {
                            break;
                        }
                    }
                    if(m!=0)
                        m--;
                    if(enDisNumLawApplNEnP[fisPReacIndex[i]][k]>1)
                        sum+=max(0.,Interpolate(2, energy, enDisEnApplVecP[fisPReacIndex[i]][k][m], enDisEnApplVecP[fisPReacIndex[i]][k][m+1],
                                    enDisProbApplVecP[fisPReacIndex[i]][k][m], enDisProbApplVecP[fisPReacIndex[i]][k][m+1]));
                    else
                        sum+=enDisProbApplVecP[fisPReacIndex[i]][k][0];
                }
                if(sum==0.)
                {
                    enDisProbApplVecP[fisPReacIndex[i]][j][low]=1;
                }
            }
        }
    }

    if(fisPReacIndex.size()>0)
    {
        stream << std::setw(14) << std::right << numPartials << endl;
        for(int i=0; i<int(fisPReacIndex.size()); i++)
        {
            for(int j=0, count=0; j<int(enDisLawP[fisPReacIndex[i]].size()); j++, count++)
            {
                if(enDisLawP[fisPReacIndex[i]][j]==2||enDisLawP[fisPReacIndex[i]][j]==4)
                {
                    stream << std::setw(14) << std::right << 0 << endl;
                    stream << std::setw(14) << std::right << enDisNumLawApplNEnP[fisPReacIndex[i]][j] <<'\n';
                    stream << std::setw(14) << std::right << enDisNumLawApplNRegP[fisPReacIndex[i]][j] <<'\n';
                    for(int k=0; k<enDisNumLawApplNRegP[fisPReacIndex[i]][j]; k++)
                    {
                        stream << std::setw(14) << std::right << enDisRangeVecP[fisPReacIndex[i]][j][k] << std::setw(14) << std::right << enDisSchemeVecP[fisPReacIndex[i]][j][k] << '\n';
                    }
                    for(int k=0; k<enDisNumLawApplNEnP[fisPReacIndex[i]][j]; k++)
                    {
                        stream << std::setw(14) << std::right << enDisEnApplVecP[fisPReacIndex[i]][j][k]*1000000 << std::setw(14) << std::right << enDisProbApplVecP[fisPReacIndex[i]][j][k] << '\n';
                    }
                    stream << '\n';
                    enerDistP[fisPReacIndex[i]][count]->WriteG4NDLData(stream);
                }
                else
                    count--;
            }
        }
    }
    stream << '\n';
    stream << '\n';

    if(promptYieldFlag&&(dNYield!=NULL))
    {
        // delayed neutron production distribution
        stream << std::setw(14) << std::right << 3 << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << 2 << '\n';
        nDelConst->WriteG4NDLData(stream);
        dNYield->WriteG4NDLData(stream);
        // this approximation is not needed
        /*
        else
        {
            if((dNTotalYieldDist->IdentifyYourSelf()=="NYieldPolyFunc")&&(dNPromptYieldDist->IdentifyYourSelf()=="NYield1DTab"))
            {
                YieldDist *temp = new NYield1DTab(dNTotalYieldDist);
                delete dNTotalYieldDist;
                dNTotalYieldDist = temp;
            }
            dNTotalYieldDist->SubtractPrompt(dNPromptYieldDist);
            dNTotalYieldDist->WriteG4NDLData(stream);
        }*/
        stream << '\n';

        // prompt neutron production distribution
        stream << std::setw(14) << std::right << 4 << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << 2 << '\n';
        if(dNPromptYieldDist->IdentifyYourSelf()=="NYieldPolyFunc")
        {
            YieldDist *temp = new NYield1DTab(dNPromptYieldDist);
            delete dNPromptYieldDist;
            dNPromptYieldDist = temp;
        }
        dNPromptYieldDist->WriteG4NDLData(stream);
    }
    else if(totalYieldFlag)
    {
        // total neutron production distribution (creates all neutrons immediately as an approximation)
        stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << isoMass;
        if(dNTotalYieldDist->IdentifyYourSelf()=="NYieldPolyFunc")
            stream << std::setw(14) << std::right << 1 << '\n';
        else
            stream << std::setw(14) << std::right << 2 << '\n';
        dNTotalYieldDist->WriteG4NDLData(stream);
    }
    else
    {
        stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << 2 << '\n';
        stream << std::setw(14) << std::right << 2 << '\n';
        stream << std::setw(14) << std::right << 1 << '\n';
        stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2 << '\n';
        stream << std::setw(14) << std::right << 1.0e-5 << std::setw(14) << std::right << abs(TYRList[2]) << '\n';
        stream << std::setw(14) << std::right << 2.0e+7 << std::setw(14) << std::right << abs(TYRList[2]) << '\n';
        dNTotalYieldDist->WriteG4NDLData(stream);
    }
    stream << '\n';

    //Energy Distribution of delayed neutrons
    if(dNYield!=NULL)
    {
        stream << std::setw(14) << std::right << 3 << std::setw(14) << std::right << 5 << '\n';
        stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << enerDistND.size() << endl;

        for(int i=0, count=0; i<int(enDisLawND.size()); i++, count++)
        {
            if(enDisLawND[i]<44)
            {
                if(enDisLawND[i]==1)
                    stream << std::setw(14) << std::right << 1 << '\n';
                else if(enDisLawND[i]==3)
                {
                    stream << std::setw(14) << std::right << 1 << '\n';
                }
                else if(enDisLawND[i]==4)
                {
                    stream << std::setw(14) << std::right << 1 << '\n';
                }
                else if(enDisLawND[i]==5)
                {
                    stream << std::setw(14) << std::right << 5 << '\n';
                }
                else if(enDisLawND[i]==7)
                {
                    stream << std::setw(14) << std::right << 7 << '\n';
                }
                else if(enDisLawND[i]==9)
                {
                    stream << std::setw(14) << std::right << 9 << '\n';
                }
                else if(enDisLawND[i]==11)
                {
                    stream << std::setw(14) << std::right << 11 << '\n';
                }
                else if(enDisLawND[i]==22)
                {
                    cout << "###: No direct translation for this law!" << endl;
                    stream << std::setw(14) << std::right << 1 << '\n';
                }
                else if(enDisLawND[i]==24)
                {
                    cout << "###: No direct translation for this law!" << endl;
                    stream << std::setw(14) << std::right << 1 << '\n';
                }

                stream << std::setw(14) << std::right << enDisNumLawApplNEnND[i] <<'\n';
                stream << std::setw(14) << std::right << enDisNumLawApplNRegND[i] <<'\n';
                for(int j=0; j<enDisNumLawApplNRegND[i]; j++)
                {
                    stream << std::setw(14) << std::right << enDisRangeVecND[i][j] << std::setw(14) << std::right << enDisSchemeVecND[i][j] << '\n';
                }
                for(int j=0; j<enDisNumLawApplNEnND[i]; j++)
                {
                    stream << std::setw(14) << std::right << enDisEnApplVecND[i][j]*1000000 << std::setw(14) << std::right << enDisProbApplVecND[i][j] << '\n';
                }

                enerDistND[count]->WriteG4NDLData(stream);
            }
            else
                count--;
        }
        stream << '\n';
        stream << '\n';
    }

    //fission fragment data
    stream << std::setw(14) << std::right << 5 << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << reacQValue[2]*0.96; // this is a big approximation but is only relavant if fission fragment generation will be used in GEANT4,
    // could be fixed by subtracting Q by the sum of the average out-going neutron and photon energy
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0;
    stream << std::setw(14) << std::right << 0 << '\n';

    string fileName = outDirName+"Fission/FS/"+isoName;
    string outDir = outDirName+"Fission/FS/";
    if(!(DirectoryExists((outDir).c_str())))
    {
        system( ("mkdir -p -m=666 "+outDir).c_str());
        if(DirectoryExists((outDir).c_str()))
        {
            cout << "created directory " << outDir << "\n" << endl;
        }
        else
        {
            cout << "\nError: could not create directory " << outDir << "\n" << endl;
            return;
        }
    }
    SetDataStream( fileName, stream, ascii);

    stream.clear();
    stream.str("");

    // FC/, SC/, TC/, and LC/
    string fisDirName[4]={"FC/", "SC/","TC/","LC/"};
    for(int i=3; i<7; i++)
    {
        if(nCSVec[i])
        {
            //check the reference frame that the data has been gathered from
            if(TYRList[i]>0)
                frameFlag=1;
            else
                frameFlag=2;

            stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << endl;

            nCSVec[i]->WriteG4NDLCSData(stream);

            //Convert Energy-Angular Dist to Energy and Angular Dist
            if(angDistInEnDistFlag[i])
            {
                if(numAngEner[i]<0)
                    numAngEner[i]=0;
                for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                {
                    if(enDisLaw[i][j]>24)
                    {
                        enDisLaw[i].push_back(4);
                        enDisNumLawApplNEn[i].push_back(enDisNumLawApplNEn[i][j]);
                        enDisNumLawApplNReg[i].push_back(enDisNumLawApplNReg[i][j]);
                        enDisRangeVec[i].push_back(new int [enDisNumLawApplNReg[i][j]]);
                        enDisSchemeVec[i].push_back(new int [enDisNumLawApplNReg[i][j]]);
                        for(int k=0; k<enDisNumLawApplNReg[i][j]; k++)
                        {
                            enDisRangeVec[i].back()[k]=enDisRangeVec[i][j][k];
                            enDisSchemeVec[i].back()[k]=enDisSchemeVec[i][j][k];
                        }

                        enDisEnApplVec[i].push_back(new double [enDisNumLawApplNEn[i][j]]);
                        enDisProbApplVec[i].push_back(new double [enDisNumLawApplNEn[i][j]]);
                        for(int k=0; k<enDisNumLawApplNEn[i][j]; k++)
                        {
                            enDisEnApplVec[i].back()[k]=enDisEnApplVec[i][j][k];
                            enDisProbApplVec[i].back()[k]=enDisProbApplVec[i][j][k];
                        }

                        AngularDist *temAngDist = NULL;
                        EnergyDist *tempEnDist = NULL;
                        angDist[i].push_back(temAngDist);
                        enerDist[i].push_back(tempEnDist);
                        angEnDist[i][count]->ConvertToEnerAndAngDist(&(enerDist[i].back()), &(angDist[i].back()), numAngEner[i]);
                        if(!(angDist[i].back()->CheckData()))
                        {
                            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3670" << endl;
                        }
                    }
                    else
                        count--;
                }
            }
            stream << endl;

            //prompt neutrons angular Distribution
            stream << std::setw(14) << std::right << 6 << std::setw(14) << std::right << 4 << endl;
            stream << std::setw(14) << std::right << repFlag;
            stream << std::setw(14) << std::right << isoMass;
            stream << std::setw(14) << std::right << frameFlag;
            stream << std::setw(14) << std::right << numAngEner[i] << '\n';
            stream << std::setw(14) << std::right << 1 << '\n';
            stream << std::setw(14) << std::right << numAngEner[i]<< std::setw(14) << std::right << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

            //note angDist[0].size() does not necessarilly equal numAngEner[0] since if the distribution type is the same between adjacent
            //incoming energy points we don't add another object to angDist, instead we add the point to the previous object
            for(int j=0; j<int(angDist[i].size()); j++)
            {
                if(!(angDist[i][j]->CheckData()))
                {
                    cout << "Error in angular data ConvertMCNPtoG4NDL.cc:92" << endl;
                }
                angDist[i][j]->WriteG4NDLData(stream);
            }
            stream << '\n' << endl;

            //Energy Distribution of prompt neutron
            stream << std::setw(14) << std::right << 6 << std::setw(14) << std::right << 5 << endl;
            stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << enerDist[i].size() << endl;

            for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
            {
                if(enDisLaw[i][j]<44)
                {
                    if(enDisLaw[i][j]==1)
                        stream << std::setw(14) << std::right << 1 << '\n';
                    else if(enDisLaw[i][j]==3)
                    {
                        stream << std::setw(14) << std::right << 1 << '\n';
                    }
                    else if(enDisLaw[i][j]==4)
                    {
                        stream << std::setw(14) << std::right << 1 << '\n';
                    }
                    else if(enDisLaw[i][j]==5)
                    {
                        stream << std::setw(14) << std::right << 5 << '\n';
                    }
                    else if(enDisLaw[i][j]==7)
                    {
                        stream << std::setw(14) << std::right << 7 << '\n';
                    }
                    else if(enDisLaw[i][j]==9)
                    {
                        stream << std::setw(14) << std::right << 9 << '\n';
                    }
                    else if(enDisLaw[i][j]==11)
                    {
                        stream << std::setw(14) << std::right << 11 << '\n';
                    }
                    else if(enDisLaw[i][j]==22)
                    {
                        cout << "###: No direct translation for this law!" << endl;
                        stream << std::setw(14) << std::right << 1 << '\n';
                    }
                    else if(enDisLaw[i][j]==24)
                    {
                        cout << "###: No direct translation for this law!" << endl;
                        stream << std::setw(14) << std::right << 1 << '\n';
                    }

                    stream << std::setw(14) << std::right << enDisNumLawApplNEn[i][j] <<'\n';
                    stream << std::setw(14) << std::right << enDisNumLawApplNReg[i][j] <<'\n';
                    for(int k=0; k<enDisNumLawApplNReg[i][j]; k++)
                    {
                        stream << std::setw(14) << std::right << enDisRangeVec[i][j][k] << std::setw(14) << std::right << enDisSchemeVec[i][j][k] << '\n';
                    }
                    for(int k=0; k<enDisNumLawApplNEn[i][j]; k++)
                    {
                        stream << std::setw(14) << std::right << enDisEnApplVec[i][j][k]*1000000 << std::setw(14) << std::right << enDisProbApplVec[i][j][k] << '\n';
                    }

                    enerDist[i][count]->WriteG4NDLData(stream);
                }
                else
                    count--;
            }
            stream << '\n';
            stream << '\n';

            fileName = outDirName+"Fission/"+fisDirName[i-3]+isoName;
            outDir = outDirName+"Fission/"+fisDirName[i-3];
            if(!(DirectoryExists((outDir).c_str())))
            {
                system( ("mkdir -p -m=666 "+outDir).c_str());
                if(DirectoryExists((outDir).c_str()))
                {
                    cout << "created directory " << outDir << "\n" << endl;
                }
                else
                {
                    cout << "\nError: could not create directory " << outDir << "\n" << endl;
                    return;
                }
            }
            SetDataStream( fileName, stream, ascii);

            stream.clear();
            stream.str("");
        }
    }
}

void CreateMT4(int* MTRListPos, string outDirName, string isoName, int isoNum, double isoMass, double temperature, int *MTRList, CSDist* nCSVec[],
                        vector<int>* enDisLaw, vector<int>* enDisNumLawApplNReg, vector<int>* enDisNumLawApplNEn, vector<int*>* enDisSchemeVec,
                        vector<int*>* enDisRangeVec, vector<double*>* enDisEnApplVec, vector<double*>* enDisProbApplVec, vector<EnergyDist*>* enerDist,
                        YieldDist* nYieldReac[], double* reacQValue, int* numAngEner, vector<AngularDist*>* angDist, bool* angDistInEnDistFlag,
                        vector<AngularEnergyDist*>* angEnDist, CSDist*** pCSVec, int *TYRList, int** numAngEnerP, vector<AngularDist*>** angDistP, vector<int>** enDisLawP,
                        vector<int>** enDisNumLawApplNRegP, vector<int>** enDisNumLawApplNEnP, vector<int*>** enDisSchemeVecP, vector<int*>** enDisRangeVecP,
                        vector<double*>** enDisEnApplVecP, vector<double*>** enDisProbApplVecP,vector<EnergyDist*>** enerDistP, vector<int> &MTRPList,
                        double*** energyAngVecP, bool ascii)
{

    bool first=true;
    CSDist1DTab* sumCSData = NULL;
    for(int i=numProcess2; i<numProcess; i++)
    {
        if(nCSVec[i])
        {
            if(first)
            {
                sumCSData = new CSDist1DTab(dynamic_cast<CSDist1DTab*>(nCSVec[i]));
                first=false;
            }
            else
                sumCSData->AddData(nCSVec[i]);
        }
    }
    if(sumCSData)
    {
        MTRListPos[7]=1;
        nCSVec[7] = sumCSData;
        numAngEner[7] = 0;

        for(int i=numProcess2; i<numProcess; i++)
        {
            if(angDistInEnDistFlag[i])
            {
                for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                {
                    if(enDisLaw[i][j]>24)
                    {
                        enDisLaw[7].push_back(4);
                        enDisNumLawApplNEn[7].push_back(enDisNumLawApplNEn[i][j]);
                        enDisNumLawApplNReg[7].push_back(enDisNumLawApplNReg[i][j]);
                        enDisRangeVec[7].push_back(new int [enDisNumLawApplNReg[i][j]]);
                        enDisSchemeVec[7].push_back(new int [enDisNumLawApplNReg[i][j]]);
                        for(int k=0; k<enDisNumLawApplNReg[i][j]; k++)
                        {
                            enDisRangeVec[7].back()[k]=enDisRangeVec[i][j][k];
                            enDisSchemeVec[7].back()[k]=enDisSchemeVec[i][j][k];
                        }

                        enDisEnApplVec[7].push_back(new double [enDisNumLawApplNEn[i][j]]);
                        enDisProbApplVec[7].push_back(new double [enDisNumLawApplNEn[i][j]]);
                        for(int k=0; k<enDisNumLawApplNEn[i][j]; k++)
                        {
                            enDisEnApplVec[7].back()[k]=enDisEnApplVec[i][j][k];
                            enDisProbApplVec[7].back()[k]=enDisProbApplVec[i][j][k];
                        }
                        AngularDist *temAngDist = NULL;
                        EnergyDist *tempEnDist = NULL;
                        angDist[7].push_back(temAngDist);
                        enerDist[7].push_back(tempEnDist);
                        angEnDist[i][count]->ConvertToEnerAndAngDist(&(enerDist[7].back()), &(angDist[7].back()), numAngEner[7]);
                        if(!(angDist[7].back()->CheckData()))
                        {
                            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3848" << endl;
                        }
                    }
                    else
                        count--;
                }
            }
        }
        angDistInEnDistFlag[7]=false;

        //Angular Distribution
        angDist[7].push_back(new AngDist2DTabular());
        angDist[7].back()->SumAngularData(angDist, nCSVec, numProcess2, numProcess, numAngEner[7]);
        if(!(angDist[7].back()->CheckData()))
        {
            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:3863" << endl;
        }
        angDist[7].back()->SetTemperature(temperature);

        //Energy Distribution of prompt neutron
        int numPartials=0, low, reg, enApplCount;
        double sum, energy;
        for(int j=numProcess2; j<numProcess; j++)
        {
            sum=0.;
            enApplCount=0;
            while((enDisNumLawApplNEn[j].size()>0)&&(sum==0.)&&(enDisNumLawApplNEn[j][0]>enApplCount))
            {
                energy=enDisEnApplVec[j][0][enApplCount];
                for(int k=0; k<int(enDisLaw[j].size()); k++)
                {
                    if(enDisLaw[j][k]<44)
                    {
                        reg=0;
                        for(low=0; low<int(enDisNumLawApplNEn[j][k]-1); low++)
                        {
                            while((enDisRangeVec[j][k][reg]<=low)&&(enDisNumLawApplNReg[j][k]-1>reg))
                                reg++;
                            if(energy<enDisEnApplVec[j][k][low])
                            {
                                break;
                            }
                        }
                        low--;
                        if(low<0)
                            low=0;

                        if(enDisNumLawApplNEn[j][k]>1)
                            sum+=max(0.,Interpolate(enDisSchemeVec[j][k][reg], energy, enDisEnApplVec[j][k][low], enDisEnApplVec[j][k][low+1],
                                        enDisProbApplVec[j][k][low], enDisProbApplVec[j][k][low+1]));
                        else
                            sum+=enDisProbApplVec[j][k][0];
                    }
                }
                enApplCount++;
            }
            for(int k=0; k<int(enDisLaw[j].size()); k++)
            {
                if(enDisLaw[j][k]<44)
                {
                    numPartials++;
                }
            }
            if(sum!=0.)
            {
                for(int k=0; k<int(enDisLaw[j].size()); k++)
                {
                    if(enDisLaw[j][k]<44)
                    {
                        for(low=0; low<int(enDisNumLawApplNEn[j][k]); low++)
                        {
                            enDisProbApplVec[j][k][low]=enDisProbApplVec[j][k][low]*max(0.,nCSVec[j]->GetAvgCS())/sum;
                        }
                    }
                }
            }
        }

        // make sure that the probability of one of the reactions energy distributions is greater than zero for every possible incoming energy
        for(int l=numProcess2; l<numProcess; l++)
        {
            for(int j=0; j<int(enDisLaw[l].size()); j++)
            {
                for(low=0; low<int(enDisNumLawApplNEn[l][j]); low++)
                {
                    energy=enDisEnApplVec[l][j][low];
                    sum=0.;
                    for(int k=0; k<int(enDisLaw[l].size()); k++)
                    {
                        int m=0;
                        for(; m<int(enDisNumLawApplNEn[l][k]-1); m++)
                        {
                            if(enDisEnApplVec[l][k][m]>energy)
                            {
                                break;
                            }
                        }
                        if(m!=0)
                            m--;
                        if(enDisNumLawApplNEn[l][k]>1)
                            sum+=max(0.,Interpolate(2, energy, enDisEnApplVec[l][k][m], enDisEnApplVec[l][k][m+1],
                                        enDisProbApplVec[l][k][m], enDisProbApplVec[l][k][m+1]));
                        else
                            sum+=enDisProbApplVec[l][k][0];
                    }
                    if(sum==0.)
                    {
                        enDisProbApplVec[l][j][low]=1;
                    }
                }
            }
        }

        for(int j=numProcess2; j<numProcess; j++)
        {
            for(int k=0; k<int(enDisLaw[j].size()); k++)
            {
                if(enDisLaw[j][k]<44)
                {
                    enDisNumLawApplNReg[7].push_back(enDisNumLawApplNReg[j][k]);
                    enDisNumLawApplNEn[7].push_back(enDisNumLawApplNEn[j][k]);
                    enDisRangeVec[7].push_back(new int[enDisNumLawApplNReg[j][k]]);
                    for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                    {
                        enDisRangeVec[7].back()[i]=enDisRangeVec[j][k][i];
                    }
                    enDisSchemeVec[7].push_back(new int[enDisNumLawApplNReg[j][k]]);
                    for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                    {
                        enDisSchemeVec[7].back()[i]=enDisSchemeVec[j][k][i];
                    }
                    enDisEnApplVec[7].push_back(new double[enDisNumLawApplNEn[j][k]]);
                    for(int i=0; i<enDisNumLawApplNEn[j][k]; i++)
                    {
                        enDisEnApplVec[7].back()[i]=enDisEnApplVec[j][k][i];
                    }
                    enDisProbApplVec[7].push_back(new double [enDisNumLawApplNEn[j][k]]);
                    for(low=0; low<int(enDisNumLawApplNEn[j][k]); low++)
                    {
                        (enDisProbApplVec[7].back())[low]=enDisProbApplVec[j][k][low];
                    }
                }
            }

        }

        int countProc=0;
        reacQValue[7]=0.;
        for(int i=numProcess2; i<numProcess; i++)
        {
            // this is only used when a mean value from the energy dist cannot be extracted, we use the mean value here as they do in GEANT4
            reacQValue[7]+=reacQValue[i];
            if(reacQValue[i]!=0.)
                countProc++;
            for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
            {
                if(enDisLaw[i][j]<44)
                {
                    enDisLaw[7].push_back(enDisLaw[i][j]);

                    enerDist[7].push_back(enerDist[i][count]);
                }
                else
                    count--;
            }
        }
        reacQValue[7]/=countProc;

        // neutron energy-angular distribution

        //this first part, adjusting the law applicability, is not needed or used so far but we keep it incase someone wants to use it later
        for(int j=numProcess2; j<numProcess; j++)
        {
            sum=0.;
            energy=0.;
            first=true;
            for(int k=0; k<int(enDisLaw[j].size()); k++)
            {
                if(enDisLaw[j][k]>24)
                {
                    reg=0;
                    if(first)
                    {
                        energy=enDisEnApplVec[j][k][int((enDisNumLawApplNEn[j][k])/2)];
                        first=false;
                    }
                    numPartials++;
                    for(low=0; low<int(enDisNumLawApplNEn[j][k]-1); low++)
                    {
                        while((enDisRangeVec[j][k][reg]<=low)&&(enDisNumLawApplNReg[j][k]-1>reg))
                            reg++;
                        if(energy<enDisEnApplVec[j][k][low])
                        {
                            break;
                        }
                    }
                    low--;
                    if(low<0)
                        low=0;
                    sum+=max(0.,Interpolate(enDisSchemeVec[j][k][reg], energy, enDisEnApplVec[j][k][low], enDisEnApplVec[j][k][low+1],
                                    enDisProbApplVec[j][k][low], enDisProbApplVec[j][k][low+1]));
                }
            }

            for(int l=0; l<int(enDisLaw[j].size()); l++)
            {
                for(low=0; low<int(enDisNumLawApplNEn[j][l]); low++)
                {
                    energy=enDisEnApplVec[j][l][low];
                    sum=0.;
                    for(int k=0; k<int(enDisLaw[j].size()); k++)
                    {
                        int m=0;
                        for(; m<int(enDisNumLawApplNEn[j][k]-1); m++)
                        {
                            if(enDisEnApplVec[j][k][m]>energy)
                            {
                                break;
                            }
                        }
                        if(m!=0)
                            m--;
                        if(enDisNumLawApplNEn[j][k]>1)
                            sum+=max(0.,Interpolate(2, energy, enDisEnApplVec[j][k][m], enDisEnApplVec[j][k][m+1],
                                        enDisProbApplVec[j][k][m], enDisProbApplVec[j][k][m+1]));
                        else
                            sum+=enDisProbApplVec[j][k][0];
                    }
                    if(sum==0.)
                    {
                        enDisProbApplVec[j][l][low]=1;
                    }
                }
            }

            for(int k=0; k<int(enDisLaw[j].size()); k++)
            {
                if(enDisLaw[j][k]<24)
                {
                    enDisNumLawApplNReg[7].push_back(enDisNumLawApplNReg[j][k]);
                    enDisNumLawApplNEn[7].push_back(enDisNumLawApplNEn[j][k]);
                    enDisRangeVec[7].push_back(new int[enDisNumLawApplNReg[j][k]]);
                    for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                    {
                        enDisRangeVec[7].back()[i]=enDisRangeVec[j][k][i];
                    }
                    enDisSchemeVec[7].push_back(new int[enDisNumLawApplNReg[j][k]]);
                    for(int i=0; i<enDisNumLawApplNReg[j][k]; i++)
                    {
                        enDisSchemeVec[7].back()[i]=enDisSchemeVec[j][k][i];
                    }
                    enDisEnApplVec[7].push_back(new double[enDisNumLawApplNEn[j][k]]);
                    for(int i=0; i<enDisNumLawApplNEn[j][k]; i++)
                    {
                        enDisEnApplVec[7].back()[i]=enDisEnApplVec[j][k][i];
                    }
                    enDisProbApplVec[7].push_back(new double [enDisNumLawApplNEn[j][k]]);
                    for(low=0; low<int(enDisNumLawApplNEn[j][k]); low++)
                    {
                        (enDisProbApplVec[7].back())[low]=enDisProbApplVec[j][k][low];
                    }
                }
            }

        }

        first=true;
        for(int i=numProcess2; i<numProcess; i++)
        {
            if(angEnDist[i].size()>0)
            {
                TYRList[7]=101;
                // weight the yield by the reactions cross-section relative to the sum of cross-section at each energy
                if(first)
                {
                    if(abs(TYRList[i])>100)
                    {
                        if(nYieldReac[i]->IdentifyYourSelf()!="NYield1DTab")
                        {
                            YieldDist* temp = new NYield1DTab(nYieldReac[i]);
                            delete nYieldReac[i];
                            nYieldReac[i]=temp;
                        }
                        nYieldReac[7] = new NYield1DTab(nYieldReac[i], nCSVec, i, numProcess, numProcess2);
                        first=false;
                    }
                    else if(nCSVec[i])
                    {
                        nYieldReac[i] = new NYield1DTab(TYRList[i],0.,20., nCSVec, i, numProcess, numProcess2);
                        nYieldReac[7] = new NYield1DTab(nYieldReac[i], nCSVec, i, numProcess, numProcess2);
                        first=false;
                    }
                }
                else
                {
                    if(abs(TYRList[i])>100)
                        nYieldReac[7]->AddData(nYieldReac[i], nCSVec, i, numProcess, numProcess2);
                    else
                    {
                        nYieldReac[i] = new NYield1DTab(TYRList[i],0.,20., nCSVec, i, numProcess, numProcess2);
                        nYieldReac[7]->AddData(nYieldReac[i], nCSVec, i, numProcess, numProcess2);
                    }
                }
                for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                {
                    if(enDisLaw[i][j]>24)
                    {
                        enDisLaw[7].push_back(enDisLaw[i][j]);
                        angEnDist[7].push_back(angEnDist[i][count]);
                    }
                    else
                        count--;
                }
            }
        }

        // instead of merging all of the photon data we simply rename the photon reactions so that they will be called for MT=4

        int pReacCount=0;
        for(int i=numProcess2; i<numProcess; i++)
        {
            for(int j=0; j<int(MTRPList.size()); j++)
            {
                if(int(MTRPList[j]/1000)==MTRList[i])
                {
                    pReacCount++;
                }
            }
        }

        CSDist **pCSVecTemp = new CSDist* [MTRPList.size()+pReacCount];
        int *numAngEnerPTemp = new int [MTRPList.size()+pReacCount];
        vector<AngularDist*> *angDistPTemp = new vector<AngularDist*> [MTRPList.size()+pReacCount];
        vector<int> *enDisLawPTemp = new vector<int> [MTRPList.size()+pReacCount];
        vector<int> *enDisNumLawApplNRegPTemp = new vector<int> [MTRPList.size()+pReacCount];
        vector<int> *enDisNumLawApplNEnPTemp = new vector<int> [MTRPList.size()+pReacCount];
        vector<int*> *enDisSchemeVecPTemp = new vector<int*> [MTRPList.size()+pReacCount];
        vector<int*> *enDisRangeVecPTemp = new vector<int*> [MTRPList.size()+pReacCount];
        vector<double*> *enDisEnApplVecPTemp = new vector<double*> [MTRPList.size()+pReacCount];
        vector<double*> *enDisProbApplVecPTemp = new vector<double*> [MTRPList.size()+pReacCount];
        vector<EnergyDist*> *enerDistPTemp = new vector<EnergyDist*> [MTRPList.size()+pReacCount];
        double **energyAngVecPTemp = new double* [MTRPList.size()+pReacCount];

        for(int j=0; j<int(MTRPList.size()); j++)
        {
            pCSVecTemp[j] = pCSVec[0][j];
            numAngEnerPTemp[j] = numAngEnerP[0][j];
            angDistPTemp[j] = angDistP[0][j];
            enDisLawPTemp[j] = enDisLawP[0][j];
            enDisNumLawApplNRegPTemp[j] = enDisNumLawApplNRegP[0][j];
            enDisNumLawApplNEnPTemp[j] = enDisNumLawApplNEnP[0][j];
            enDisSchemeVecPTemp[j] = enDisSchemeVecP[0][j];
            enDisRangeVecPTemp[j] = enDisRangeVecP[0][j];
            enDisEnApplVecPTemp[j] = enDisEnApplVecP[0][j];
            enDisProbApplVecPTemp[j] = enDisProbApplVecP[0][j];
            enerDistPTemp[j] = enerDistP[0][j];
            energyAngVecPTemp[j] = energyAngVecP[0][j];
        }

        int pReacCount2=0;
        for(int i=numProcess2; i<numProcess; i++)
        {
            for(int j=0; j<int(MTRPList.size()); j++)
            {
                if(int(MTRPList[j]/1000)==MTRList[i])
                {
                    pCSVecTemp[MTRPList.size()+pReacCount2] = pCSVec[0][j];
                    numAngEnerPTemp[MTRPList.size()+pReacCount2] = numAngEnerP[0][j];
                    angDistPTemp[MTRPList.size()+pReacCount2] = angDistP[0][j];
                    enDisLawPTemp[MTRPList.size()+pReacCount2] = enDisLawP[0][j];
                    enDisNumLawApplNRegPTemp[MTRPList.size()+pReacCount2] = enDisNumLawApplNRegP[0][j];
                    enDisNumLawApplNEnPTemp[MTRPList.size()+pReacCount2] = enDisNumLawApplNEnP[0][j];
                    enDisSchemeVecPTemp[MTRPList.size()+pReacCount2] = enDisSchemeVecP[0][j];
                    enDisRangeVecPTemp[MTRPList.size()+pReacCount2] = enDisRangeVecP[0][j];
                    enDisEnApplVecPTemp[MTRPList.size()+pReacCount2] = enDisEnApplVecP[0][j];
                    enDisProbApplVecPTemp[MTRPList.size()+pReacCount2] = enDisProbApplVecP[0][j];
                    enerDistPTemp[MTRPList.size()+pReacCount2] = enerDistP[0][j];
                    energyAngVecPTemp[MTRPList.size()+pReacCount2] = energyAngVecP[0][j];
                    pReacCount2++;
                }
            }
        }

        //delete out going photon data

        if(enDisLawP[0]!=NULL)
            delete [] enDisLawP[0];
        if(enDisNumLawApplNRegP[0]!=NULL)
            delete [] enDisNumLawApplNRegP[0];
        if(enDisNumLawApplNEnP[0]!=NULL)
            delete [] enDisNumLawApplNEnP[0];
        if(energyAngVecP[0]!=NULL)
            delete [] energyAngVecP[0];
        if(angDistP[0]!=NULL)
            delete [] angDistP[0];
        if(pCSVec[0]!=NULL)
            delete[] pCSVec[0];
        if(enDisSchemeVecP[0]!=NULL)
            delete [] enDisSchemeVecP[0];
        if(enDisRangeVecP[0]!=NULL)
            delete [] enDisRangeVecP[0];
        if(enDisProbApplVecP[0]!=NULL)
            delete [] enDisProbApplVecP[0];
        if(enDisEnApplVecP[0]!=NULL)
            delete [] enDisEnApplVecP[0];
        if(enerDistP[0]!=NULL)
            delete [] enerDistP[0];

        for(int j=0; j<pReacCount; j++)
        {
            MTRPList.push_back(4000+j);
        }

        pCSVec[0] = pCSVecTemp;
        numAngEnerP[0] = numAngEnerPTemp;
        angDistP[0] = angDistPTemp;
        enDisLawP[0] = enDisLawPTemp;
        enDisNumLawApplNRegP[0] = enDisNumLawApplNRegPTemp;
        enDisNumLawApplNEnP[0] = enDisNumLawApplNEnPTemp;
        enDisSchemeVecP[0] = enDisSchemeVecPTemp;
        enDisRangeVecP[0] = enDisRangeVecPTemp;
        enDisEnApplVecP[0] = enDisEnApplVecPTemp;
        enDisProbApplVecP[0] = enDisProbApplVecPTemp;
        enerDistP[0] = enerDistPTemp;
        energyAngVecP[0] = energyAngVecPTemp;
    }
}

void MakeInElasticFSFile(int *MTRListPos, string outDirName, string isoName, int isoNum, double isoMass, double temperature, int *MTRList, CSDist *nCSVec[],
                        vector<int> *enDisLaw, vector<int> *enDisNumLawApplNReg, vector<int> *enDisNumLawApplNEn, vector<int*> *enDisSchemeVec,
                        vector<int*> *enDisRangeVec, vector<double*> *enDisEnApplVec, vector<double*> *enDisProbApplVec, vector<EnergyDist*> *enerDist,
                        YieldDist *nYieldReac[], double *reacQValue, int *numAngEner, vector<AngularDist*> *angDist, bool *angDistInEnDistFlag,
                        vector<AngularEnergyDist*> *angEnDist, CSDist **pCSVec, int *TYRList, int *numAngEnerP, vector<AngularDist*> *angDistP, vector<int> *enDisLawP,
                        vector<int> *enDisNumLawApplNRegP, vector<int> *enDisNumLawApplNEnP, vector<int*> *enDisSchemeVecP, vector<int*> *enDisRangeVecP,
                        vector<double*> *enDisEnApplVecP, vector<double*> *enDisProbApplVecP,vector<EnergyDist*> *enerDistP, vector<int> MTRPList,
                        double **energyAngVecP, bool ascii)
{
    //Creates the inelastic scattering final state files

    std::stringstream stream, numConv;
    stream.fill(' ');
    stream.precision(6);
    stream.setf(std::ios::scientific);

    int frameFlag, repFlag, dirNum=0, startProc=7, endProc=numProcess;
    vector<int> inEPReacIndex;
    vector<bool> firstPass(36,1);

    /*if(MTRListPos[7]<0)
    {
        startProc=8;
        endProc=numProcess;
    }*/
    for(int i=startProc; i<endProc; i++)
    {
        dirNum++;
        if((dirNum==12)||(dirNum==16)||(dirNum==29)||(dirNum==32))
            dirNum++;

        if(i>=numProcess2)
        {
            dirNum=1;
        }

        frameFlag=2;
        repFlag=2;
        stream.str("");
        stream.clear();
        numConv.str("");
        numConv.clear();
        inEPReacIndex.clear();

        if(MTRListPos[i]!=-1)
        {
            //check the reference frame that the data has been gathered from
            if(TYRList[MTRListPos[i]]>0)
                frameFlag=1;

            stream.fill(' ');

            //create list of relevant photon production reactions
            for(int j=0; j<int(MTRPList.size()); j++)
            {
                if(int(MTRPList[j]/1000)==MTRList[i])
                {
                    inEPReacIndex.push_back(j);
                }
            }

            // F02, F03, F04,  F05, F06, F07, F08, F09, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21, F22,  F28, F29, F30, F31, F32, F33, F34,
            // F35, F36 directories

            if((dirNum!=1)&&(!((dirNum>22)&&(dirNum<28))))
            {
                //probability of reaction occuring
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 3 << '\n';
                //dummy variables at the beginning of the file
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0;
                nCSVec[i]->WriteG4NDLCSData(stream);

                //prompt neutrons

                //Angular Distribution
                if(!angDistInEnDistFlag[i])
                {
                    stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 4 << '\n';
                    stream << std::setw(14) << std::right << repFlag;
                    stream << std::setw(14) << std::right << isoMass;
                    stream << std::setw(14) << std::right << frameFlag;
                    stream << std::setw(14) << std::right << numAngEner[i] << '\n';
                    stream << std::setw(14) << std::right << 1 << '\n';
                    stream << std::setw(14) << std::right << numAngEner[i] << std::setw(14) << std::right << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

                    //note angDist[i].size() does not necessarilly equal numAngEner[i] since if the distribution type is the same between adjacent
                    //incoming energy points we don't add another object to angDist, instead we add the point to the previous object
                    for(int j=0; j<int(angDist[i].size()); j++)
                    {
                        if(!(angDist[i][j]->CheckData()))
                        {
                            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:4344" << endl;
                        }
                        angDist[i][j]->WriteG4NDLData(stream);
                    }
                    stream << '\n' << endl;
                }
                stream << '\n' << endl;

                //Energy Distribution of prompt neutron
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 5 << '\n';
                stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << enerDist[i].size() << endl;

                for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                {
                    if(enDisLaw[i][j]<44)
                    {
                        if(enDisLaw[i][j]==1)
                            stream << std::setw(14) << std::right << 1 << '\n';
                        else if(enDisLaw[i][j]==3)
                        {
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }
                        else if(enDisLaw[i][j]==4)
                        {
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }
                        else if(enDisLaw[i][j]==5)
                        {
                            stream << std::setw(14) << std::right << 5 << '\n';
                        }
                        else if(enDisLaw[i][j]==7)
                        {
                            stream << std::setw(14) << std::right << 7 << '\n';
                        }
                        else if(enDisLaw[i][j]==9)
                        {
                            stream << std::setw(14) << std::right << 9 << '\n';
                        }
                        else if(enDisLaw[i][j]==11)
                        {
                            stream << std::setw(14) << std::right << 11 << '\n';
                        }
                        else if(enDisLaw[i][j]==22)
                        {
                            cout << "###: No direct translation for this law!" << endl;
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }
                        else if(enDisLaw[i][j]==24)
                        {
                            cout << "###: No direct translation for this law!" << endl;
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }

                        stream << std::setw(14) << std::right << enDisNumLawApplNEn[i][j] <<'\n';
                        stream << std::setw(14) << std::right << enDisNumLawApplNReg[i][j] <<'\n';
                        for(int k=0; k<enDisNumLawApplNReg[i][j]; k++)
                        {
                            stream << std::setw(14) << std::right << enDisRangeVec[i][j][k] << std::setw(14) << std::right << enDisSchemeVec[i][j][k] << '\n';
                        }
                        for(int k=0; k<enDisNumLawApplNEn[i][j]; k++)
                        {
                            stream << std::setw(14) << std::right << enDisEnApplVec[i][j][k]*1000000 << std::setw(14) << std::right << enDisProbApplVec[i][j][k] << '\n';
                        }

                        enerDist[i][count]->WriteG4NDLData(stream);
                    }
                    else
                        count--;
                }
                stream << '\n';

                // neutron energy-angular distribution
                if(angDistInEnDistFlag[i])
                {
                    stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 6 << '\n';
                    stream << std::setw(14) << std::right << isoMass;
                    stream << std::setw(14) << std::right << frameFlag;
                    stream << std::setw(14) << std::right << angEnDist[i].size() << '\n';

                    for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                    {
                        if(enDisLaw[i][j]>24)
                        {
                            //set the out-going particle type to always be a neutron
                            stream << std::setw(14) << std::right << 1;
                            stream << std::setw(14) << std::right << isoMass;
                            stream << std::setw(14) << std::right << 0;

                            if(enDisLaw[i][j]==44)
                            {
                                //cout << "###: minor translation used for this law!" << endl;
                                stream << std::setw(14) << std::right << 7 << '\n';
                            }
                            else if(enDisLaw[i][j]==61)
                            {
                                //cout << "###: minor translation used for this law!" << endl;
                                stream << std::setw(14) << std::right << 7 << '\n';
                            }
                            else if(enDisLaw[i][j]==66)
                            {
                                stream << std::setw(14) << std::right << 6 << '\n';
                            }
                            else if(enDisLaw[i][j]==67)
                            {
                                stream << std::setw(14) << std::right << 7 << '\n';
                            }

                            stream << std::setw(14) << std::right << reacQValue[i];
                            stream << std::setw(14) << std::right << reacQValue[i];

                            if(abs(TYRList[i])>100)
                                nYieldReac[i]->WriteG4NDLData(stream);
                            else
                            {
                                stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';
                                stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2 << '\n';

                                stream << std::setw(14) << std::right << 0.0;
                                stream << std::setw(14) << std::right << abs(TYRList[i]);

                                stream << std::setw(14) << std::right << 20000000.0;
                                stream << std::setw(14) << std::right <<  abs(TYRList[i]);
                            }

                            angEnDist[i][count]->WriteG4NDLData(stream);
                        }
                        else
                            count--;
                    }
                }
                stream << '\n' << endl;

                // photon production distribution
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 12 << '\n';
                // I use repFlag 1 no matter what
                stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << inEPReacIndex.size();
                for(int j=0; j<int(inEPReacIndex.size()); j++)
                {
                    //we always select a continous energy dist
                    stream << std::setw(14) << std::right << 1;
                    //This average energy is not used when the previous value is set to 1 (ie when there is an energy distribution)
                    stream << std::setw(14) << std::right << enerDistP[inEPReacIndex[j]].front()->GetAverageOutEnergy() << '\n';

                    pCSVec[inEPReacIndex[j]]->WriteG4NDLYieldData(stream);
                }
                stream << '\n';
                stream << '\n';

                // photon reaction cross-section
                // this format has not yet been fully implemented in GEANT4 and the data would not be used anyways since we set repFlag=1
                /*
                stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 13 << '\n';
                */

                // photon angular distribution
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 14 << '\n';
                stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << std::setw(14) << std::right << min(1,int(inEPReacIndex.size())) << std::setw(14) << std::right << 0 << endl;

                if(inEPReacIndex.size()>0)
                {
                    vector<AngularDist*> angDistPVec;
                    vector<CSDist*> pCSVecTemp;
                    for(int j=0; j<int(inEPReacIndex.size()); j++)
                    {
                        for(int k=0; k<int(angDistP[inEPReacIndex[j]].size()); k++)
                        {
                            if(!(angDistP[inEPReacIndex[j]][k]->CheckData()))
                            {
                                cout << "Error in angular data ConvertMCNPtoG4NDL.cc:4512" << endl;
                            }
                            angDistPVec.push_back(angDistP[inEPReacIndex[j]][k]);
                            pCSVecTemp.push_back(pCSVec[inEPReacIndex[j]]);
                        }
                    }

                    int tempNumAngEner;
                    AngularDist *tempAng = new AngDist2DTabularP();
                    tempAng->SumAngularData(angDistPVec, pCSVecTemp, tempNumAngEner);
                    if(!(tempAng->CheckData()))
                    {
                        cout << "Error in angular data ConvertMCNPtoG4NDL.cc:4524" << endl;
                    }
                    stream << std::setw(14) << std::right << enerDistP[inEPReacIndex[0]].front()->GetAverageOutEnergy() << std::setw(14) << std::right << 0 << endl;
                    stream << std::setw(14) << std::right << tempNumAngEner << '\n';
                    tempAng->WriteG4NDLData(stream);
                    delete tempAng;

                    stream << '\n';
                }
                stream << '\n';

                // photon energy distribution
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 15 << '\n';
                //Here we adjust the tables defining the probability of an energy dist applicablitity to be normalized to 1 and then weighted by the cross-section
                //this is an approximation that lets us put all the energy dist for every photon production reaction caused by this neutron production reaction
                //into one big list. the normalization ensures that the probability sum between photon production reactions is the same and the weighting of the cross-sections
                //ensures that the energy distributions from the more probable photon production reaction are more likely to be applied
                int numPartials=0, low, reg, enApplCount;
                double sum, energy;
                for(int j=0; j<int(inEPReacIndex.size()); j++)
                {
                    sum=0.;
                    enApplCount=0;
                    while((sum==0.)&&(enDisNumLawApplNEnP[inEPReacIndex[j]][0]>enApplCount))
                    {
                        energy=enDisEnApplVecP[inEPReacIndex[j]][0][enApplCount];
                        for(int k=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++)
                        {
                            if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                            {
                                reg=0;
                                for(low=0; low<int(enDisNumLawApplNEnP[inEPReacIndex[j]][k]-1); low++)
                                {
                                    while((enDisRangeVecP[inEPReacIndex[j]][k][reg]<=low)&&(enDisNumLawApplNRegP[inEPReacIndex[j]][k]-1>reg))
                                        reg++;
                                    if(energy<enDisEnApplVecP[inEPReacIndex[j]][k][low])
                                    {
                                        break;
                                    }
                                }
                                low--;
                                if(low<0)
                                    low=0;

                                if(enDisNumLawApplNEnP[inEPReacIndex[j]][k]>1)
                                    sum+=max(0.,Interpolate(enDisSchemeVecP[inEPReacIndex[j]][k][reg], energy, enDisEnApplVecP[inEPReacIndex[j]][k][low], enDisEnApplVecP[inEPReacIndex[j]][k][low+1],
                                                enDisProbApplVecP[inEPReacIndex[j]][k][low], enDisProbApplVecP[inEPReacIndex[j]][k][low+1]));
                                else
                                    sum+=enDisProbApplVecP[inEPReacIndex[j]][k][0];
                            }
                        }
                        enApplCount++;
                    }
                    for(int k=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++)
                    {
                        if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                        {
                            numPartials++;
                        }
                    }
                    if(sum!=0.)
                    {
                        for(int k=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++)
                        {
                            if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                            {
                                for(low=0; low<int(enDisNumLawApplNEnP[inEPReacIndex[j]][k]); low++)
                                {
                                    enDisProbApplVecP[inEPReacIndex[j]][k][low]=enDisProbApplVecP[inEPReacIndex[j]][k][low]*max(0.,pCSVec[inEPReacIndex[j]]->GetAvgCS())/sum;
                                }
                            }
                        }
                    }
                }

                // make sure that the probability of one of the reactions energy distributions is greater than zero for every possible incoming energy
                for(int l=0; l<int(inEPReacIndex.size()); l++)
                {
                    for(int j=0; j<int(enDisLawP[inEPReacIndex[l]].size()); j++)
                    {
                        for(low=0; low<int(enDisNumLawApplNEnP[inEPReacIndex[l]][j]); low++)
                        {
                            energy=enDisEnApplVecP[inEPReacIndex[l]][j][low];
                            sum=0.;
                            for(int k=0; k<int(enDisLawP[inEPReacIndex[l]].size()); k++)
                            {
                                int m=0;
                                for(; m<int(enDisNumLawApplNEnP[inEPReacIndex[l]][k]-1); m++)
                                {
                                    if(enDisEnApplVecP[inEPReacIndex[l]][k][m]>energy)
                                    {
                                        break;
                                    }
                                }
                                if(m!=0)
                                    m--;
                                if(enDisNumLawApplNEnP[inEPReacIndex[l]][k]>1)
                                    sum+=max(0.,Interpolate(2, energy, enDisEnApplVecP[inEPReacIndex[l]][k][m], enDisEnApplVecP[inEPReacIndex[l]][k][m+1],
                                                enDisProbApplVecP[inEPReacIndex[l]][k][m], enDisProbApplVecP[inEPReacIndex[l]][k][m+1]));
                                else
                                    sum+=enDisProbApplVecP[inEPReacIndex[l]][k][0];
                            }
                            if(sum==0.)
                            {
                                enDisProbApplVecP[inEPReacIndex[l]][j][low]=1;
                            }
                        }
                    }
                }

                if(inEPReacIndex.size()>0)
                {
                    stream << std::setw(14) << std::right << numPartials << endl;
                    for(int j=0; j<int(inEPReacIndex.size()); j++)
                    {
                        for(int k=0, count=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++, count++)
                        {
                            if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                            {
                                stream << std::setw(14) << std::right << 0 << endl;
                                stream << std::setw(14) << std::right << enDisNumLawApplNEnP[inEPReacIndex[j]][k] <<'\n';
                                stream << std::setw(14) << std::right << enDisNumLawApplNRegP[inEPReacIndex[j]][k] <<'\n';
                                for(int l=0; l<enDisNumLawApplNRegP[inEPReacIndex[j]][k]; l++)
                                {
                                    stream << std::setw(14) << std::right << enDisRangeVecP[inEPReacIndex[j]][k][l] << std::setw(14) << std::right << enDisSchemeVecP[inEPReacIndex[j]][k][l] << '\n';
                                }
                                for(int l=0; l<enDisNumLawApplNEnP[inEPReacIndex[j]][k]; l++)
                                {
                                    stream << std::setw(14) << std::right << enDisEnApplVecP[inEPReacIndex[j]][k][l]*1000000 << std::setw(14) << std::right << enDisProbApplVecP[inEPReacIndex[j]][k][l] << '\n';
                                }
                                stream << '\n';
                                enerDistP[inEPReacIndex[j]][count]->WriteG4NDLData(stream);
                            }
                            else
                                count--;
                        }
                    }
                }
                stream << '\n';
                stream << '\n';
            }

            // F01, F23, F24, F25, F26, F27, directories
            else
            {
                //probability of reaction occuring
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 3 << '\n';
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                //dummy variables at the beginning of the file
                stream << std::setw(14) << std::right << reacQValue[i] << std::setw(14) << std::right << 0;
                nCSVec[i]->WriteG4NDLCSData(stream);
                stream << '\n' << endl;

                //prompt neutrons

                //Angular Distribution
                if(!angDistInEnDistFlag[i])
                {
                    stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 4 << '\n';
                    stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                    stream << std::setw(14) << std::right << repFlag;
                    stream << std::setw(14) << std::right << isoMass;
                    stream << std::setw(14) << std::right << frameFlag;
                    stream << std::setw(14) << std::right << numAngEner[i] << '\n';
                    stream << std::setw(14) << std::right << 1 << '\n';
                    stream << std::setw(14) << std::right << numAngEner[i] << std::setw(14) << std::right << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

                    //note angDist[i].size() does not necessarilly equal numAngEner[i] since if the distribution type is the same between adjacent
                    //incoming energy points we don't add another object to angDist, instead we add the point to the previous object
                    for(int j=0; j<int(angDist[i].size()); j++)
                    {
                        if(!(angDist[i][j]->CheckData()))
                        {
                            cout << "Error in angular data ConvertMCNPtoG4NDL.cc:4696" << endl;
                        }
                        angDist[i][j]->WriteG4NDLData(stream);
                    }
                    stream << '\n' << endl;
                }

                //Energy Distribution of prompt neutron
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 5 << '\n';
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << enerDist[i].size() << endl;

                for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                {
                    if(enDisLaw[i][j]<44)
                    {
                        if(enDisLaw[i][j]==1)
                            stream << std::setw(14) << std::right << 1 << '\n';
                        else if(enDisLaw[i][j]==3)
                        {
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }
                        else if(enDisLaw[i][j]==4)
                        {
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }
                        else if(enDisLaw[i][j]==5)
                        {
                            stream << std::setw(14) << std::right << 5 << '\n';
                        }
                        else if(enDisLaw[i][j]==7)
                        {
                            stream << std::setw(14) << std::right << 7 << '\n';
                        }
                        else if(enDisLaw[i][j]==9)
                        {
                            stream << std::setw(14) << std::right << 9 << '\n';
                        }
                        else if(enDisLaw[i][j]==11)
                        {
                            stream << std::setw(14) << std::right << 11 << '\n';
                        }
                        else if(enDisLaw[i][j]==22)
                        {
                            cout << "###: No direct translation for this law!" << endl;
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }
                        else if(enDisLaw[i][j]==24)
                        {
                            cout << "###: No direct translation for this law!" << endl;
                            stream << std::setw(14) << std::right << 1 << '\n';
                        }

                        stream << std::setw(14) << std::right << enDisNumLawApplNEn[i][j] <<'\n';
                        stream << std::setw(14) << std::right << enDisNumLawApplNReg[i][j] <<'\n';
                        for(int k=0; k<enDisNumLawApplNReg[i][j]; k++)
                        {
                            stream << std::setw(14) << std::right << enDisRangeVec[i][j][k] << std::setw(14) << std::right << enDisSchemeVec[i][j][k] << '\n';
                        }
                        for(int k=0; k<enDisNumLawApplNEn[i][j]; k++)
                        {
                            stream << std::setw(14) << std::right << enDisEnApplVec[i][j][k]*1000000 << std::setw(14) << std::right << enDisProbApplVec[i][j][k] << '\n';
                        }

                        enerDist[i][count]->WriteG4NDLData(stream);
                        stream << '\n';
                    }
                    else
                        count--;
                }
                stream << '\n';

                // neutron energy-angular distribution
                if(angDistInEnDistFlag[i])
                {
                    stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 6 << '\n';
                    stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                    stream << std::setw(14) << std::right << isoMass;
                    stream << std::setw(14) << std::right << frameFlag;
                    stream << std::setw(14) << std::right << angEnDist[i].size() << '\n';

                    for(int j=0, count=0; j<int(enDisLaw[i].size()); j++, count++)
                    {
                        if(enDisLaw[i][j]>24)
                        {
                            //set the out-going particle type to always be a neutron
                            stream << std::setw(14) << std::right << 1;
                            stream << std::setw(14) << std::right << isoMass;
                            stream << std::setw(14) << std::right << 0;

                            if(enDisLaw[i][j]==44)
                            {
                                //cout << "###: minor translation used for this law!" << endl;
                                stream << std::setw(14) << std::right << 7 << '\n';
                            }
                            else if(enDisLaw[i][j]==61)
                            {
                                //cout << "###: minor translation used for this law!" << endl;
                                stream << std::setw(14) << std::right << 7 << '\n';
                            }
                            else if(enDisLaw[i][j]==66)
                            {
                                stream << std::setw(14) << std::right << 6 << '\n';
                            }
                            else if(enDisLaw[i][j]==67)
                            {
                                stream << std::setw(14) << std::right << 7 << '\n';
                            }

                            stream << std::setw(14) << std::right << reacQValue[i];
                            stream << std::setw(14) << std::right << reacQValue[i];

                            if(abs(TYRList[i])>100)
                                nYieldReac[i]->WriteG4NDLData(stream);
                            else
                            {
                                stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 1 << '\n';
                                stream << std::setw(14) << std::right << 2 << std::setw(14) << std::right << 2 << '\n';

                                stream << std::setw(14) << std::right << 0.0;
                                stream << std::setw(14) << std::right << abs(TYRList[i]);

                                stream << std::setw(14) << std::right << 20000000.0;
                                stream << std::setw(14) << std::right << abs(TYRList[i]);
                            }

                            angEnDist[i][count]->WriteG4NDLData(stream);
                        }
                        else
                            count--;
                    }
                    stream << '\n';
                }


                // photon production distribution
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 12 << '\n';
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                // I use repFlag 1 no matter what
                stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << inEPReacIndex.size();
                for(int j=0; j<int(inEPReacIndex.size()); j++)
                {
                    //we always select a continous energy dist
                    stream << std::setw(14) << std::right << 1;
                    //This average energy is not used when the previous value is set to 1 (ie when there is an energy distribution)
                    stream << std::setw(14) << std::right << enerDistP[inEPReacIndex[j]].front()->GetAverageOutEnergy() << '\n';

                    pCSVec[inEPReacIndex[j]]->WriteG4NDLYieldData(stream);
                }
                stream << '\n';
                stream << '\n';

                // photon reaction cross-section
                // this format has not yet been fully implemented in GEANT4 and the data would not be used anyways since we set repFlag=1
                /*
                stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 13 << '\n';
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                */

                // photon angular distribution
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 14 << '\n';
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << std::setw(14) << std::right << min(1,int(inEPReacIndex.size())) << std::setw(14) << std::right << 0 << endl;

                if(inEPReacIndex.size()>0)
                {
                    vector<AngularDist*> angDistPVec;
                    vector<CSDist*> pCSVecTemp;
                    for(int j=0; j<int(inEPReacIndex.size()); j++)
                    {
                        for(int k=0; k<int(angDistP[inEPReacIndex[j]].size()); k++)
                        {
                            if(!(angDistP[inEPReacIndex[j]][k]->CheckData()))
                            {
                                cout << "Error in angular data ConvertMCNPtoG4NDL.cc:4870" << endl;
                            }
                            angDistPVec.push_back(angDistP[inEPReacIndex[j]][k]);
                            pCSVecTemp.push_back(pCSVec[inEPReacIndex[j]]);
                        }
                    }

                    int tempNumAngEner;
                    AngularDist *tempAng = new AngDist2DTabularP();
                    tempAng->SumAngularData(angDistPVec, pCSVecTemp, tempNumAngEner);
                    if(!(tempAng->CheckData()))
                    {
                        cout << "Error in angular data ConvertMCNPtoG4NDL.cc:4882" << endl;
                    }
                    stream << std::setw(14) << std::right << enerDistP[inEPReacIndex[0]].front()->GetAverageOutEnergy() << std::setw(14) << std::right << 0 << endl;
                    stream << std::setw(14) << std::right << tempNumAngEner << '\n';
                    tempAng->WriteG4NDLData(stream);
                    delete tempAng;

                    stream << '\n';
                }
                stream << '\n';

                // photon energy distribution
                stream << std::setw(14) << std::right << dirNum-1 << std::setw(14) << std::right << 15 << '\n';
                stream << std::setw(14) << std::right << MTRList[i] << std::setw(14) << std::right << 0 << '\n';
                //Here we adjust the tables defining the probability of an energy dist applicablitity to be normalized to 1 and then weighted by the cross-section
                //this is an approximation that lets us put all the energy dist for every photon production reaction caused by this neutron production reaction
                //into one big list. the normalization ensures that the probability sum between photon production reactions is the same and the weighting of the cross-sections
                //ensures that the energy distributions from the more probable photon production reaction are more likely to be applied
                int numPartials=0, low, reg, enApplCount;
                double sum, energy;
                for(int j=0; j<int(inEPReacIndex.size()); j++)
                {
                    sum=0.;
                    enApplCount=0;
                    while((sum==0.)&&(enDisNumLawApplNEnP[inEPReacIndex[j]][0]>enApplCount))
                    {
                        energy=enDisEnApplVecP[inEPReacIndex[j]][0][enApplCount];
                        for(int k=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++)
                        {
                            if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                            {
                                reg=0;
                                for(low=0; low<int(enDisNumLawApplNEnP[inEPReacIndex[j]][k]-1); low++)
                                {
                                    while((enDisRangeVecP[inEPReacIndex[j]][k][reg]<=low)&&(enDisNumLawApplNRegP[inEPReacIndex[j]][k]-1>reg))
                                        reg++;
                                    if(energy<enDisEnApplVecP[inEPReacIndex[j]][k][low])
                                    {
                                        break;
                                    }
                                }
                                low--;
                                if(low<0)
                                    low=0;
                                if(enDisNumLawApplNEnP[inEPReacIndex[j]][k]>1)
                                    sum+=max(0.,Interpolate(enDisSchemeVecP[inEPReacIndex[j]][k][reg], energy, enDisEnApplVecP[inEPReacIndex[j]][k][low], enDisEnApplVecP[inEPReacIndex[j]][k][low+1],
                                                enDisProbApplVecP[inEPReacIndex[j]][k][low], enDisProbApplVecP[inEPReacIndex[j]][k][low+1]));
                                else
                                    sum+=enDisProbApplVecP[inEPReacIndex[j]][k][0];
                            }
                        }
                        enApplCount++;
                    }
                    for(int k=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++)
                    {
                        if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                        {
                            numPartials++;
                        }
                    }
                    if(sum!=0.)
                    {
                        for(int k=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++)
                        {
                            if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                            {
                                for(low=0; low<int(enDisNumLawApplNEnP[inEPReacIndex[j]][k]); low++)
                                {
                                    enDisProbApplVecP[inEPReacIndex[j]][k][low]=enDisProbApplVecP[inEPReacIndex[j]][k][low]*max(0.,pCSVec[inEPReacIndex[j]]->GetAvgCS())/sum;
                                }
                            }
                        }
                    }
                }

                // make sure that the probability of one of the reactions energy distributions is greater than zero for every possible incoming energy
                for(int l=0; l<int(inEPReacIndex.size()); l++)
                {
                    for(int j=0; j<int(enDisLawP[inEPReacIndex[l]].size()); j++)
                    {
                        for(low=0; low<int(enDisNumLawApplNEnP[inEPReacIndex[l]][j]); low++)
                        {
                            energy=enDisEnApplVecP[inEPReacIndex[l]][j][low];
                            sum=0.;
                            for(int k=0; k<int(enDisLawP[inEPReacIndex[l]].size()); k++)
                            {
                                int m=0;
                                for(; m<int(enDisNumLawApplNEnP[inEPReacIndex[l]][k]-1); m++)
                                {
                                    if(enDisEnApplVecP[inEPReacIndex[l]][k][m]>energy)
                                    {
                                        break;
                                    }
                                }
                                if(m!=0)
                                    m--;
                                if(enDisNumLawApplNEnP[inEPReacIndex[l]][k]>1)
                                    sum+=max(0.,Interpolate(2, energy, enDisEnApplVecP[inEPReacIndex[l]][k][m], enDisEnApplVecP[inEPReacIndex[l]][k][m+1],
                                                enDisProbApplVecP[inEPReacIndex[l]][k][m], enDisProbApplVecP[inEPReacIndex[l]][k][m+1]));
                                else
                                    sum+=enDisProbApplVecP[inEPReacIndex[l]][k][0];
                            }
                            if(sum==0.)
                            {
                                enDisProbApplVecP[inEPReacIndex[l]][j][low]=1;
                            }
                        }
                    }
                }

                if(inEPReacIndex.size())
                {
                    stream << std::setw(14) << std::right << numPartials << endl;
                    for(int j=0; j<int(inEPReacIndex.size()); j++)
                    {
                        for(int k=0, count=0; k<int(enDisLawP[inEPReacIndex[j]].size()); k++, count++)
                        {
                            if(enDisLawP[inEPReacIndex[j]][k]==2||enDisLawP[inEPReacIndex[j]][k]==4)
                            {
                                stream << std::setw(14) << std::right << 0 << endl;
                                stream << std::setw(14) << std::right << enDisNumLawApplNEnP[inEPReacIndex[j]][k] <<'\n';
                                stream << std::setw(14) << std::right << enDisNumLawApplNRegP[inEPReacIndex[j]][k] <<'\n';
                                for(int l=0; l<enDisNumLawApplNRegP[inEPReacIndex[j]][k]; l++)
                                {
                                    stream << std::setw(14) << std::right << enDisRangeVecP[inEPReacIndex[j]][k][l] << std::setw(14) << std::right << enDisSchemeVecP[inEPReacIndex[j]][k][l] << '\n';
                                }
                                for(int l=0; l<enDisNumLawApplNEnP[inEPReacIndex[j]][k]; l++)
                                {
                                    stream << std::setw(14) << std::right << enDisEnApplVecP[inEPReacIndex[j]][k][l]*1000000 << std::setw(14) << std::right << enDisProbApplVecP[inEPReacIndex[j]][k][l] << '\n';
                                }
                                stream << '\n';
                                enerDistP[inEPReacIndex[j]][count]->WriteG4NDLData(stream);
                            }
                            else
                                count--;
                        }
                    }
                    stream << '\n';
                    stream << '\n';
                }
            }

            if(dirNum<10)
                numConv << "Inelastic/F0" << dirNum << "/";
            else
                numConv << "Inelastic/F" << dirNum << "/";

            string fileName = outDirName+numConv.str()+isoName;
            string outDir = outDirName+numConv.str();
            if(!(DirectoryExists((outDir).c_str())))
            {
                system( ("mkdir -p -m=666 "+outDir).c_str());
                if(DirectoryExists((outDir).c_str()))
                {
                    cout << "created directory " << outDir << "\n" << endl;
                }
                else
                {
                    cout << "\nError: could not create directory " << outDir << "\n" << endl;
                    return;
                }
            }
            SetDataStream( fileName, stream, ascii, firstPass[dirNum-1]);
            firstPass[dirNum-1]=false;
        }
    }
    /*
        MCNP Data does not seem to contain the data for delayed photons so we don't create the /Inelastic/Gammas/ directory
    */

}

inline double Interpolate(int aScheme,
            double x, double x1, double x2, double y1, double y2)
{
  double result(0);
  int theScheme = aScheme;
  theScheme = theScheme%7;
  switch(theScheme)
  {
    case 1:
      //080809
      //result = Histogram(x, x1, x2, y1, y2);
      result = LinearLinear(x, x1, x2, y1, y2);
      break;
    case 2:
      result = LinearLinear(x, x1, x2, y1, y2);
      break;
    case 3:
      result = LinearLogarithmic(x, x1, x2, y1, y2);
      break;
    case 4:
      result = LogarithmicLinear(x, x1, x2, y1, y2);
      break;
    case 5:
      result = LogarithmicLogarithmic(x, x1, x2, y1, y2);
      break;
    default:
      cout << "Error: Unrecognized scheme = "<<theScheme<<endl;
      break;
  }
  return result;
}

inline double Histogram(double , double , double , double y1, double )
{
  double result;
  result = y1;
  return result;
}

inline double LinearLinear(double x, double x1, double x2, double y1, double y2)
{
  double slope=0, off=0;
  if(x2-x1==0) return (y2+y1)/2.;
  slope = (y2-y1)/(x2-x1);
  off = y2-x2*slope;
  double y = x*slope+off;
  return y;
}

inline double LinearLogarithmic(double x, double x1, double x2, double y1, double y2)
{
  double result;
  if(x==0) result = y1+y2/2.;
  else if(x1==0) result = y1;
  else if(x2==0) result = y2;
  else result = LinearLinear(std::log(x), std::log(x1), std::log(x2), y1, y2);
  return result;
}

inline double LogarithmicLinear(double x, double x1, double x2, double y1, double y2)
{
  double result;
  if(y1==0||y2==0) result = 0;
  else
  {
    result = LinearLinear(x, x1, x2, std::log(y1), std::log(y2));
    result = std::exp(result);
  }
  return result;
}

inline double LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2)
{
  double result;
  if(x==0) result = y1+y2/2.;
  else if(x1==0) result = y1;
  else if(x2==0) result = y2;
  if(y1==0||y2==0) result = 0;
  else
  {
    result = LinearLinear(std::log(x), std::log(x1), std::log(x2), std::log(y1), std::log(y2));
    result = std::exp(result);
  }
  return result;
}

void GetDataStream( string filename, std::stringstream& ss)
{
   string* data=NULL;
   std::ifstream* in=NULL;
   //string compfilename(filename);

   if(filename.substr((filename.length()-2),2)==".z")
   {
        in = new std::ifstream ( filename.c_str() , std::ios::binary | std::ios::ate );
   }

   if ( in!=NULL && in->good() )
   {
// Use the compressed file
      uLongf file_size = (uLongf)(in->tellg());
      in->seekg( 0 , std::ios::beg );
      Bytef* compdata = new Bytef[ file_size ];

      while ( *in )
      {
         in->read( (char*)compdata , file_size );
      }

      uLongf complen = (uLongf) ( file_size*4 );
      Bytef* uncompdata = new Bytef[complen];

      while ( Z_OK != uncompress ( uncompdata , &complen , compdata , file_size ) )
      {
         delete[] uncompdata;
         complen *= 2;
         uncompdata = new Bytef[complen];
      }
      delete [] compdata;
      //                                 Now "complen" has uncomplessed size
      data = new string ( (char*)uncompdata , (long)complen );
      delete [] uncompdata;
   }

   else
   {
// Use regular text file
      std::ifstream thefData( filename.c_str() , std::ios::in | std::ios::ate );
      if ( thefData.good() )
      {
         int file_size = thefData.tellg();
         thefData.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( thefData )
         {
            thefData.read( filedata , file_size );
         }
         thefData.close();
         data = new string ( filedata , file_size );
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );
         cout << endl << "### failed to open ascii file " << filename << " ###" << endl;
      }
   }
   if (data != NULL)
   {
        ss.str(*data);
        if(data->back()!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

   if(in!=NULL)
   {
        in->close();
        delete in;
   }

    if(data!=NULL)
        delete data;
}


void SetDataStream( string filename , std::stringstream& ss, bool ascii, bool overWrite )
{
    //bool cond=true;
   if (!ascii)
   {
        string compfilename(filename);

        if(compfilename.back()!='z')
            compfilename += ".z";

        std::ofstream* out;
        if(overWrite)
            out = new std::ofstream ( compfilename.c_str() , std::ios::binary | std::ios::trunc);
        else
            out = new std::ofstream ( compfilename.c_str() , std::ios::binary | std::ios::app );
       if ( ss.good() )
       {
       //
    // Create the compressed file
          ss.seekg( 0 , std::ios::end );
          uLongf file_size = (uLongf)(ss.tellg());
          ss.seekg( 0 , std::ios::beg );
          Bytef* uncompdata = new Bytef[ file_size ];

          while ( ss ) {
              ss.read( (char*)uncompdata , file_size );
          }

          uLongf complen = compressBound(file_size);

          Bytef* compdata = new Bytef[complen];

          if ( Z_OK == compress ( compdata , &complen , uncompdata , file_size ) )
          {
            out->write((char*)compdata, (long)complen);
            if (out->fail())
            {
                cout << endl << "writing the compressed data to the output file " << compfilename << " failed" << endl
                    << " may not have permission to delete an older version of the file" << endl;
            }
          }
          else
          {
            cout << endl << "compressing the data failed" << endl;
          }

          delete [] uncompdata;
          delete [] compdata;
       }
       else
       {
            cout << endl << "### failed to write to binary file ###" << endl;
       }

       out->close(); delete out;
   }
   else
   {
        // Use regular text file
        string compfilename(filename);

        if(compfilename.substr((compfilename.length()-2),2)==".z")
        {
            compfilename.pop_back();
            compfilename.pop_back();
        }

        std::ofstream *out;
        if(overWrite)
            out = new std::ofstream ( compfilename.c_str() , std::ios::out | std::ios::trunc );
        else
            out = new std::ofstream ( compfilename.c_str() , std::ios::out | std::ios::app );
        if ( ss.good() )
        {
             ss.seekg( 0 , std::ios::end );
             int file_size = ss.tellg();
             ss.seekg( 0 , std::ios::beg );
             char* filedata = new char[ file_size ];
             while ( ss ) {
                ss.read( filedata , file_size );
                if(!file_size)
                {
                    cout << "\n #### Error the size of the stringstream is invalid ###" << endl;
                    break;
                }
             }
             out->write(filedata, file_size);
            if (out->fail())
            {
                cout << endl << "writing the ascii data to the output file " << compfilename << " failed" << endl
                     << " may not have permission to delete an older version of the file" << endl;
            }

             delete [] filedata;
        }
        else
        {
        // found no data file
        //                 set error bit to the stream
         ss.setstate( std::ios::badbit );

         cout << endl << "### failed to write to ascii file " << compfilename << " ###" << endl;
        }
        out->close();
        delete out;
   }
   ss.str("");
}
