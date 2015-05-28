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

#define numProcess 4

int CreateIsoCSData(stringstream &stream, string outDirName, bool ascii, double setTemperature, bool limitTemp, bool onlyCS);
bool DirectoryExists( const char* pzPath );
void MakeCSDataFile(string fileName, int MTRNum, double energyCSVec[], CSDist *csDist, bool ascii);
void MakeElasFSFile(string fileName, double isoMass, double temperature, int TYRList[], int numAngEner[], double *energyAngVec[],
                    int *angDistType[], int *intSchemeAng[], int *numAngProb[], double **angVec[], double **angProbVec[], bool ascii);
void GetDataStream( string, std::stringstream&);
void SetDataStream( string, std::stringstream&, bool);

// we will have to check to make sure that the distributions are being used in the right energy regime, they may have to be listed in the G$NDL file in order of increasing energy
// the chronological ordering of the process extraction is probably unnessary since they seem to already be in order
// Maybe an error in the extraction of inelastic, only a few of the isotopes seem to have the identifier it is looking for
// X Note: range values (NBT) seem to be absolute positions in the array, test to make sure this is true
// X create code interpreter using the functions from ExtractMaterialComposition to make sure that every dynamic variable is deleted
// in AngEnDit3D classes check to make sure that the locators are being used properly, we assumed a -1 add on based off other examples
// X give the option to only extract CS data
// go throught the write functions and check that they match the format of G4NDL
// optimize data usage by only passing the sections of data containers and arrays needed within a particular function

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
                    // Gets data from the file and stores it into a data stream
                    GetDataStream(fileName, stream);
                    while(stream)
                    {
                        stream >> word;
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
        GetDataStream(inFileName, stream);
        lib = inFileName[int(inFileName.length()-3)];

        while(stream)
        {
            stream >> word;
            check1=word[int(word.find_last_of('.')+1)];
            check2=word[int(word.find_last_of('.')+2)];
            check3=word[int(word.find_last_of('.')+3)];

            //checks whether the word matches the beggining of an isotope data set identifier
            if((check1==lib)&&((check2>='0')&&(check2<='9'))/*&&((check3=='c')||(check3=='d'))*/)
            {
                if((check3=='c')||(check3=='d'))
                {
                    // gets the elastic, inelastic, fission and capture CS data for the isotope
                    result += CreateIsoCSData(stream, outDirName, ascii);
                }
            }
        }
    }

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
    string outDirNameCSProc[4];
    stringstream numConv;
    int Znum, isoNum;
    double isoMass;
    int yieldDistType;
    bool promptYieldFlag=false, totalYieldFlag=false;

    //fixed location data
    int numCSEner=0, numReactions=0, numPReactions=0, numDNPrecursorFam=0;
    int startEnerTable, startElasticBlock, startMTRBlock;
    int startNUBlock, startQValBlock, startTYRBlock;
    int startLSIGBlock, startCSBlock;
    int startLANDBlock, startANDBlock;
    int startLDLWBlock, startDLWBlock;
    int startMTRPBlock, startLSIGPBlock, startCSPBlock;
    int startLANDPBlock, startANDPBlock;
    int startLDLWPBlock, startDLWBPlock;
    int startYPBlock;
    int startDelayedNUBlock, startDNEDLBlock, startDNEDBlock;

    //temperary data
    double temp;
    int intTemp;
    int index, index2, count=1;
    char line[256];

    //fixed data arrays
    double *energyCSVec=NULL; // incoming neutron kinetic energy for CS data

    //for the following arrays, the indicies 0,1,2,3 correspond to elastic, capture, fission and inelastic respectively
    //location arrays are set to -1 to show that the data has not been set
    int MTRList[numProcess]={2,102,18,4}; // this array sets the order of the processes, must be manually enetered
    int MTRListPos[numProcess]; //position of the reaction MT in MT list
    for(int i=0; i<numProcess; i++)
    {
        MTRListPos[i]=-1;
    }
    double reacQValue[numProcess]; // Q-Value of the reaction
    int TYRList[numProcess]; // the average number of neutrons released from the reaction,
    // the sign determines the the reference frame that the angular and energy dist where measured from, - means CM and + means Lab
    TYRList[0]=-1;
    int LSIGList[numProcess]; //starting point of the reaction cross-section data in the data file
    for(int i=0; i<numProcess; i++)
    {
        LSIGList[i]=-1;
    }
    CSDist nCSVec[numProcess]; //cross section array for each reaction
    int LANDList[numProcess]; //starting point of the reaction's out going neutron angular distribution data
    for(int i=0; i<numProcess; i++)
    {
        LANDList[i]=-1;
    }
    int numAngEner[numProcess]; // the number of incoming neutron energy points for the angular distribution
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


    // first 12 lines contain miscellaneous data and the NXS and JXS arrays page 118 of the MCNP manual
    //### in this section we extract the isotope naming, temperature and mass data

    stream >> isoMass >> temp;
    isoMass=isoMass; //the isotope mass / neutron mass
    temp=temp*1000000/(8.6173324*(pow(10, -5)));
    numConv << temp;
    numConv >> temperature;

    if(limitTemp&&(temperature!=setTemperature))
    {
        return;
    }

    numConv.clear();
    numConv.str("");

    if(outDirName.back!='/')
        outDirName+='/';

    outDirName = outDirName+temperature+'/';

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
                cout << "created temperature directory " << outDirNameCSProc[i] << "\n" << endl;
            }
            else
            {
                cout << "\nError: could not create temperature Directory " << outDirNameCSProc[i] << "\n" << endl;
                return 1;
            }
        }
    }

    for(int i=0; i<6; i++)
    {
        stream.getline(line,256);
    }
    // Extracts the isotope name
    //start of NXS block
    stream >> dummy >> isoNum;
    Znum = floor(isoNum/1000);
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

    stream >> numCSEner >> numReactions >> dummy >> numPReactions >> dummy >> numDNPrecursorFam;

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

    //### In this section we extract the cross-section data

    for(;count<startEnerTable; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numCSEner; i++, count++)
    {
        stream >> temp;
        energyCSVec[i]=temp*1000000.; //Converting from MeV to eV
    }

    for(;count<startElasticBlock; count++)
    {
        stream >> dummy;
    }

    nCSVec[0] = new CSDist1DTab(energyCSVec, numCSEner, startEnerTable);
    nCSVec[0]->ExtractMCNPData(stream, count);

    if(numCSEner!=0)
            MakeCSDataFile(outDirNameCSProc[0]+isoName, MTRList[0], energyCSVec, nCSVec[0], ascii);

    //Extract NU block here for fission yield
    YieldDist *dNPromptYieldDist;
    YieldDist *dNTotalYieldDist;
    if(startNUBlock>0)
    {
        for(;count<startNUBlock; count++)
        {
            stream >> dummy;
        }
        //###Check this: if the delayed neutron yield exists, assume the NU block contains the prompt yield
        (startDelayedNUBlock!=0) ? promptYieldFlag=true : totalYieldFlag=true;

        if(promptYieldFlag)
        {
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
        }
        else
        {
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
    else if(startNUBlock<0)
    {
        for(;count<(-1*startNUBlock+1); count++)
        {
            stream >> dummy;
        }
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

        for(;count<(-2*startNUBlock+1); count++)
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

    for(;count<startMTRBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        //make sure these MT numbers are right
        for(int j=1; j<numProcess; j++)
        {
            if(MTRList[j]==index)
            {
                MTRListPos[j]=i;
            }
        }
    }

    for(;count<startQValBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> temp;
        for(int j=1; j<numProcess; j++)
        {
            if(MTRListPos[j]==i)
            {
                reacQValue[j]=temp;
            }
        }
    }


    for(;count<startTYRBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numProcess; j++)
        {
            if(MTRListPos[j]==i)
            {
                TYRList[j]=index;
            }
        }
    }

    for(;count<startLSIGBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numProcess; j++)
        {
            if(MTRListPos[j]==i)
            {
                LSIGList[j]=index;
            }
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
        for(;count<(LSIGList[extractOrder[i]]+startCSBlock-1); count++)
        {
            stream >> dummy;
        }
        nCSVec[[extractOrder[i]]] = new CSDist1DTab(energyCSVec);
        nCSVec[[extractOrder[i]]]->ExtractMCNPData(stream, count);
        if((CSVecSize[extractOrder[i]]!=0)&&(MTRList[extractOrder[i]]==102||MTRList[extractOrder[i]]==18||MTRList[extractOrder[i]]==4))
            MakeCSDataFile(outDirNameCSProc[extractOrder[i]]+isoName, MTRList[extractOrder[i]], energyCSVec, nCSVec[[extractOrder[i]]], ascii);
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

        return;
    }

    //### In this Section we extract the angular distribution data that is independant
    // of energy the rest of the angular data is extracted in the energy distribution


    for(;count<startLANDBlock; count++)
    {
        stream >> dummy;
    }

    //elastic
    stream >> index; count++;
    LANDList[0]=index;

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numProcess; j++)
        {
            if(MTRListPos[j]==i)
            {
                LANDList[j]=index;
            }
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
    double *angTabPosVec[numProcess]; // the location of the angular probability table associated to the incoming neutron kinetic energy
    int *exOrderAng;
    int distType=0;
    vector<AngularDist*> angDist[numProcess];

    for(int i=0; i<numProcess; i++)
    {
        // an isotropic distribution is assumed
        if(LANDList[extractOrder[i]]==0)
        {
            angDist[extractOrder[i]].push_back(new AngDistIso());
            angDist[extractOrder[i]].back()->SetTemperature(temperature);
            continue;
        }

        //the angular distribution is contained in the energy distribution data
        if(LANDList[extractOrder[i]]==-1)
        {
            angDistInEnDistFlag[extractOrder[i]]=true;
            continue;
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

        angTabPosVec[extractOrder[i]] = new double[numAngEner[extractOrder[i]]];

        for(int j=0; j<numAngEner[extractOrder[i]]; j++, count++)
        {
            stream >> temp;
            angTabPosVec[extractOrder[i]][j] = temp;
        }

        exOrderAng = new int[numAngEner[extractOrder[i]]]
        //## we may be creating an error here by mixing up the order of the distributions with respect to energy
        for(int j=0; j<numAngEner[extractOrder[i]]; j++)
        {
            for(int k=j+1; k<numAngEner[extractOrder[i]]; k++)
            {
                if(angTabPosVec[extractOrder[i]][exOrderAng[k]]<angTabPosVec[extractOrder[i]][exOrderAng[j]])
                {
                    index = exOrderAng[j];
                    exOrderAng[j] = exOrderAng[k];
                    exOrderAng[k] = index;
                }
            }
        }

        for(int j=0; j<numAngEner[extractOrder[i]]; j++)
        {
        //Change this so that the data is inputted into the appropriate daughter of the mother class AngularDist
            if(angTabPosVec[extractOrder[i]][exOrderAng[j]]=0)
            {
                //the distribution is isotropic and no further data is needed
                if(distType!=1)
                {
                    angDist[extractOrder[i]].push_back(new AngDistIso());
                    angDist[extractOrder[i]].back()->SetTemperature(temperature);
                    distType=1;
                }
            }
            else if(angTabPosVec[extractOrder[i]][exOrderAng[j]]<0)
            {
                //the angular distribution is represented as a table of cosine and prob
                for(;count<(startANDBlock-angTabPosVec[extractOrder[i]][exOrderAng[j]]-1); count++)
                {
                    stream >> dummy;
                }

                if(distType!=2)
                {
                    angDist[extractOrder[i]].push_back(new AngDist2DTabular());
                    angDist[extractOrder[i]].back()->SetTemperature(temperature);
                    distType=2;
                }
            }
            else
            {
                //the angular distribution is represented by a probability density function where each of the 32 points/bins are spaced along the entire cosine of the scattering angle axis
                //such that the integral of the probability density with the cosine of the scattering angle inbetween them is kept constant, meaning that the bins are in a equilprobable arrangement
                // the first bin probably corresponds to a cosine value of -1 and and the last bin being some where around 1 depending on the spacing
                for(;count<(startANDBlock+angTabPosVec[extractOrder[i]][exOrderAng[j]]-1); count++)
                {
                    stream >> dummy;
                }

                if(distType!=3)
                {
                    angDist[extractOrder[i]].push_back(new AngDist32EqualPBin());
                    angDist[extractOrder[i]].back()->SetTemperature(temperature);
                    distType=3;
                }
            }
            angDist[extractOrder[i]].back()->SetPoint(stream, count, energyAngVec[extractOrder[i]][exOrderAng[j]]);
        }
        if(exOrderAng)
            delete [] exOrderAng;

    }

    //### In this section we extract the outgoing hadron/photon energy distribution data

    for(;count<startLDLWBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numProcess; j++)
        {
            if(MTRListPos[j]==i)
            {
                LDLWList[j]=index;
            }
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

    int nextLawPos, enDisLawDataPos;
    vector<int> enDisLaw[numProcess], enDisNumLawApplNReg[numProcess], enDisNumLawApplNEn[numProcess];
    vector<int*> enDisSchemeVec[numProcess], enDisRangeVec[numProcess];
    vector<double*> enDisEnApplVec[numProcess], enDisProbApplVec[numProcess];
    vector<EnergyDist*> enerDist[numProcess];
    NYield1DTab *nYieldReac[numProcess];

    for(int i=1; i<numProcess; i++)
    {
        for(;count<(LDLWList[extractOrder[i]]+startDLWBlock-1); count++)
        {
            stream >> dummy;
        }
        do
        {
            stream >> nextLawPos >> intTemp >> enDisLawDataPos; count=count+3;
            enDisLaw[extractOrder[i]].push_back(intTemp);
            stream >> intTemp; count++;
            enDisNumLawApplNReg[extractOrder[i]].push_back(intTemp);

            enDisSchemeVec[extractOrder[i]].push_back(new int [enDisNumLawApplNReg[extractOrder[i]].back()]);
            enDisRangeVec[extractOrder[i]].push_back(new int [enDisNumLawApplNReg[extractOrder[i]].back()]);

            for(int j=0; j<enDisNumLawApplNReg[extractOrder[i]].back(); j++, count++)
            {
                stream >> intTemp;
                (enDisSchemeVec[extractOrder[i]].back())[j]=intTemp;
            }

            for(int j=0; j<enDisNumLawApplNReg[extractOrder[i]].back(); j++, count++)
            {
                stream >> intTemp;
                (enDisRangeVec[extractOrder[i]].back())[j]=intTemp;
            }

            stream >> intTemp; count++;
            enDisNumLawApplNEn[extractOrder[i]].push_back(intTemp);

            enDisEnApplVec[extractOrder[i]].push_back(new double [enDisNumLawApplNEn[extractOrder[i]].back()]);
            enDisProbApplVec[extractOrder[i]].push_back(new double [enDisNumLawApplNEn[extractOrder[i]].back()]);

            for(int j=0; j<enDisNumLawApplNEn[extractOrder[i]].back(); j++, count++)
            {
                stream >> temp;
                (enDisEnApplVec[extractOrder[i]].back())[j]=temp;
            }

            for(int j=0; j<enDisNumLawApplNEn.back(); j++, count++)
            {
                stream >> temp;
                (enDisProbApplVec[extractOrder[i]].back())[j]=temp;
            }

            for(;count<(enDisLawDataPos+startDLWBlock-1); count++)
            {
                stream >> dummy;
            }

            // Now we extract the data depending on the format of the law using a different class for each distinct case
            // Don't have to worry about energy ordering of the laws, the law accuracy probability
            // ensures that the law will only be used in the right energy regime.

            if(enDisLaw[extractOrder[i]]==1)
                enerDist[extractOrder[i]].push_back(new EnerDistEqPEnerBins());

            else if(enDisLaw[extractOrder[i]]==2)
                enerDist[extractOrder[i]].push_back(new EnerDist1PhEner(isoMass));

            else if(enDisLaw[extractOrder[i]]==3)
                enerDist[extractOrder[i]].push_back(new EnerDistLevScat());

            else if(enDisLaw[extractOrder[i]]==4)
                enerDist[extractOrder[i]].push_back(new EnerDistConTab());

            else if(enDisLaw[extractOrder[i]]==5)
                enerDist[extractOrder[i]].push_back(new EnerDistGenEvapSpec());

            else if(enDisLaw[extractOrder[i]]==7)
                enerDist[extractOrder[i]].push_back(new EnerDistMaxwellFisSpec());

            else if(enDisLaw[extractOrder[i]]==9)
                enerDist[extractOrder[i]].push_back(new EnerDistEvapSpec());

            else if(enDisLaw[extractOrder[i]]==11)
                enerDist[extractOrder[i]].push_back(new EnerDistWattSpec());

            else if(enDisLaw[extractOrder[i]]==22)
                enerDist[extractOrder[i]].push_back(new EnerDistTabLinFunc(startDLWBlock));

            else if(enDisLaw[extractOrder[i]]==24)
                enerDist[extractOrder[i]].push_back(new EnerDistTabMulti());

            else if(enDisLaw[extractOrder[i]]==44)
                angDist[extractOrder[i]].push_back(new AngEnDistKallbach(startDLWBlock));

            else if(enDisLaw[extractOrder[i]]==61)
                angDist[extractOrder[i]].push_back(new AngEnDist3DTab(startDLWBlock));

            else if(enDisLaw[extractOrder[i]]==66)
                angDist[extractOrder[i]].push_back(new AngEnDistNBody());

            else if(enDisLaw[extractOrder[i]]==67)
                angDist[extractOrder[i]].push_back(new AngEnDistLab3DTab());
            else
            {
                cout << "\n### Error: Energy law not recognized! ###" << endl;
                continue;
            }

            if(enDisLaw[extractOrder[i]]<44)
                enerDist[extractOrder[i]].back()->ExtractMCNPData(stream, count);

            else
                angDist[extractOrder[i]].back()->ExtractMCNPData(stream, count);

            // go to next law position
            for(;count<(nextLawPos+startDLWBlock-1); count++)
            {
                stream >> dummy;
            }
        }
        //this assumes that the energy distribution data collection is terminated when nextLawPos=0, may not be true check
        while(nextLawPos!=0)

        if(abs(TYRList[extractOrder[i]])>100);
        {
            for(;count<(abs(TYRList[extractOrder[i]])+startDLWBlock-101); count++)
            {
                stream >> dummy;
            }
            nYieldReac[extractOrder[i]] = new NYield1DTab;
            nYieldReac[extractOrder[i]]->ExtractMCNPData(stream, count);
        }

    }

    //Create elastic FS files
    // may get some of the angular dist data in the energy dist data so the elastic final state files can only be produced after the energy dist has been extracted
    outDirNameElasFS = outDirName +"Elastic/FS/";
    MakeElasFSFile(outDirNameElasFS+isoName, isoMass, temperature, TYRList, numAngEner, energyAngVec, angDistType, intSchemeAng, numAngProb, angVec, angProbVec, ascii);



    //### Extract photon data, follows same format as the neutron data above, except for a few minor differences
    // Probably only one energy distribution per photon production reaction, don't need vector containers

    vector<int> MTRPList;// this array sets the order of the processes
    vector<int> MTRPListPos; //position of the reaction MT in MT list

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
            if(MTRList[j]==int(index/1000))
            {
                MTRPList.push_back(index);
                MTRPListPos.push_back(i);
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
    CSDist **pCSVec = new CSDist *[numPProcess]; //cross section array for each reaction
    int *LANDPList = new int [numPProcess]; //starting point of the reaction's out going neutron angular distribution data
    for(int i=0; i<numPProcess; i++)
    {
        LANDPList[i]=-1;
    }
    int *numAngEnerP = new int [numPProcess]; // the number of incoming neutron energy points for the angular distribution
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

    //### in this section we extract the CS data for the different photon production reactions

    for(;count<startLSIGPBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numPReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numPProcess; j++)
        {
            if(MTRPListPos[j]==i)
            {
                LSIGPList[j]=index;
            }
        }
    }

    for(int i=1; i<numPProcess; i++)
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

    for(int i=1; i<numPProcess; i++)
    {
        for(;count<(LSIGPList[extractOrderP[i]]+startCSPBlock-1); count++)
        {
            stream >> dummy;
        }
        stream >> intTemp; count++;
        MFType[extractOrderP[i]] = intTemp;

        if(MFType[extractOrderP[i]]==13)
        {
            pCSVec[extractOrderP[i]] = new CSDist1DTabP(energyCSVec);
            pCSVec[extractOrderP[i]]->ExtractMCNPData(stream, count);
        }
        else
        {
            pCSVec[extractOrderP[i]] = new CSDistYieldComp();
            pCSVec[extractOrderP[i]]->ExtractMCNPData(stream, count);
        }
    }

    //### In this Section we extract the angular distribution data that is independant
    // of energy the rest of the angular data is extracted in the energy distribution


    for(;count<startLANDPBlock; count++)
    {
        stream >> dummy;
    }

    //elastic
    stream >> index; count++;
    LANDPList[0]=index;

    for(int i=0; i<numPReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numPProcess; j++)
        {
            if(MTRPListPos[j]==i)
            {
                LANDPList[j]=index;
            }
        }
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

    double **energyAngVecP = new double *[numPProcess]; // incoming neutron kinetic energy for Ang dist data
    double **angTabPosVecP= new double *[numPProcess]; // the location of the angular probability table associated to the incoming neutron kinetic energy
    int *exOrderAngP;
    int distType=0;
    vector<AngularDist*> *angDistP = new vector<AngularDist*> [numPProcess];

    for(int i=0; i<numPProcess; i++)
    {
        // an isotropic distribution is assumed
        if(LANDPList[extractOrderP[i]]==0)
        {
            angDistP[extractOrderP[i]].push_back(new AngDistIsoP());
            angDistP[extractOrderP[i]].back()->SetTemperature(temperature);
            continue;
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

        angTabPosVecP[extractOrderP[i]] = new double[numAngEnerP[extractOrderP[i]]];

        for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++, count++)
        {
            stream >> temp;
            angTabPosVecP[extractOrderP[i]][j] = temp;
        }

        exOrderAngP = new int[numAngEnerP[extractOrderP[i]]]
        //## we may be creating an error here by mixing up the order of the distributions with respect to energy
        for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++)
        {
            for(int k=j+1; k<numAngEnerP[extractOrderP[i]]; k++)
            {
                if(angTabPosVecP[extractOrderP[i]][exOrderAngP[k]]<angTabPosVecP[extractOrderP[i]][exOrderAngP[j]])
                {
                    index = exOrderAngP[j];
                    exOrderAngP[j] = exOrderAngP[k];
                    exOrderAngP[k] = index;
                }
            }
        }

        for(int j=0; j<numAngEnerP[extractOrderP[i]]; j++)
        {
        //Change this so that the data is inputted into the appropriate daughter of the mother class AngularDist
            //the angular distribution is represented by a probability density function where each of the 32 points/bins are spaced along the entire cosine of the scattering angle axis
            //such that the integral of the probability density with the cosine of the scattering angle inbetween them is kept constant, meaning that the bins are in a equilprobable arrangement
            // the first bin probably corresponds to a cosine value of -1 and and the last bin being some where around 1 depending on the spacing
            if(angTabPosVecP[extractOrderP[i]]!=0)
            {
                for(;count<(startANDPBlock+angTabPosVecP[extractOrderP[i]][exOrderAngP[j]]-1); count++)
                {
                    stream >> dummy;
                }
                if(distType!=1)
                {
                    angDistP[extractOrderP[i]].push_back(new AngDist32EqualPBin());
                    distType==1;
                }
            }
            else if(distType!=2)
            {
                angDistP[extractOrderP[i]].push_back(new AngDistIsoP());
                distType==2;
            }
            angDistP[extractOrderP[i]].back()->SetPoint(stream, count, energyAngVecP[extractOrderP[i]][exOrderAngP[j]]);

        }
        if(exOrderAngP)
            delete [] exOrderAngP;

    }

    //### In this section we extract the outgoing hadron/photon energy distribution data

    for(;count<startLDLWPBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numPReactions; i++, count++)
    {
        stream >> index;
        for(int j=1; j<numPProcess; j++)
        {
            if(MTRPListPos[j]==i)
            {
                LDLWPList[j]=index;
            }
        }
    }

    for(int i=1; i<numPProcess; i++)
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

    int nextLawPos, enDisLawDataPos;
    vector<int> *enDisLawP = new vector<int> [numPProcess], *enDisNumLawApplNRegP = new vector<int> [numPProcess], *enDisNumLawApplNEnP = new vector<int> [numPProcess];
    vector<int*> *enDisSchemeVecP = new vector<int*> [numPProcess], *enDisRangeVecP = new vector<int*> [numPProcess];
    vector<double*> *enDisEnApplVecP = new vector<double*> [numPProcess], *enDisProbApplVecP = new vector<double*> [numPProcess];
    vector<EnergyDist*> *enerDistP = new vector<EnergyDist*> [numPProcess];

    for(int i=1; i<numPProcess; i++)
    {
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

            enDisSchemeVecP[extractOrderP[i]].push_back(new int [enDisNumLawApplNRegP[extractOrderP[i]].back()]);
            enDisRangeVecP[extractOrderP[i]].push_back(new int [enDisNumLawApplNRegP[extractOrderP[i]].back()]);

            for(int j=0; j<enDisNumLawApplNRegP[extractOrderP[i]].back(); j++, count++)
            {
                stream >> intTemp;
                (enDisSchemeVecP[extractOrderP[i]].back())[j]=intTemp;
            }

            for(int j=0; j<enDisNumLawApplNRegP[extractOrderP[i]].back(); j++, count++)
            {
                stream >> intTemp;
                (enDisRangeVecP[extractOrderP[i]].back())[j]=intTemp;
            }

            stream >> intTemp; count++;
            enDisNumLawApplNEnP[extractOrderP[i]].push_back(intTemp);

            enDisEnApplVecP[extractOrderP[i]].push_back(new double [enDisNumLawApplNEnP[extractOrderP[i]].back()]);
            enDisProbApplVecP[extractOrderP[i]].push_back(new double [enDisNumLawApplNEnP[extractOrderP[i]].back()]);

            for(int j=0; j<enDisNumLawApplNEnP[extractOrderP[i]].back(); j++, count++)
            {
                stream >> temp;
                (enDisEnApplVecP[extractOrderP[i]].back())[j]=temp;
            }

            for(int j=0; j<enDisNumLawApplNEnP.back(); j++, count++)
            {
                stream >> temp;
                (enDisProbApplVecP[extractOrderP[i]].back())[j]=temp;
            }

            for(;count<(enDisLawDataPos+startDLWPBlock-1); count++)
            {
                stream >> dummy;
            }

            // Now we extract the data depending on the format of the law using a different class for each distinct case
            // Don't have to worry about energy ordering of the laws, the law accuracy probability
            // ensures that the law will only be used in the right energy regime.

            if(enDisLawP[extractOrderP[i]]==1)
                enerDistP[extractOrderP[i]].push_back(new EnerDistEqPEnerBins());

            else if(enDisLawP[extractOrderP[i]]==2)
                enerDistP[extractOrderP[i]].push_back(new EnerDist1PhEner(isoMass));

            else if(enDisLawP[extractOrderP[i]]==3)
                enerDistP[extractOrderP[i]].push_back(new EnerDistLevScat());

            else if(enDisLawP[extractOrderP[i]]==4)
                enerDistP[extractOrderP[i]].push_back(new EnerDistConTab());

            else if(enDisLawP[extractOrderP[i]]==5)
                enerDistP[extractOrderP[i]].push_back(new EnerDistGenEvapSpec());

            else if(enDisLawP[extractOrderP[i]]==7)
                enerDistP[extractOrderP[i]].push_back(new EnerDistMaxwellFisSpec());

            else if(enDisLawP[extractOrderP[i]]==9)
                enerDistP[extractOrderP[i]].push_back(new EnerDistEvapSpec());

            else if(enDisLawP[extractOrderP[i]]==11)
                enerDistP[extractOrderP[i]].push_back(new EnerDistWattSpec());

            else if(enDisLawP[extractOrderP[i]]==22)
                enerDistP[extractOrderP[i]].push_back(new EnerDistTabLinFunc(startDLWPBlock));

            else if(enDisLawP[extractOrderP[i]]==24)
                enerDistP[extractOrderP[i]].push_back(new EnerDistTabMulti());

            else if(enDisLawP[extractOrderP[i]]==44)
                angDistP[extractOrderP[i]].push_back(new AngEnDistKallbach(startDLWPBlock));

            else if(enDisLawP[extractOrderP[i]]==61)
                angDistP[extractOrderP[i]].push_back(new AngEnDist3DTab(startDLWPBlock));

            else if(enDisLawP[extractOrderP[i]]==66)
                angDistP[extractOrderP[i]].push_back(new AngEnDistNBody());

            else if(enDisLawP[extractOrderP[i]]==67)
                angDistP[extractOrderP[i]].push_back(new AngEnDistLab3DTab());
            else
            {
                cout << "\n### Error: Energy law not recognized! ###" << endl;
                continue;
            }

            if(enDisLawP[extractOrderP[i]]<44)
                enerDistP[extractOrderP[i]].back()->ExtractMCNPData(stream, count);

            // these energy angular distributions should not occur for photon data, check
            else
            {
                angDistP[extractOrderP[i]].back()->ExtractMCNPData(stream, count);
                cout << "\n### Error: photon energy dependant angular data is being collected" << endl;
            }


            // go to next law position
            for(;count<(nextLawPos+startDLWPBlock-1); count++)
            {
                stream >> dummy;
            }
        }
        //this assumes that the energy distribution data collection is terminated when nextLawPos=0, may not be true check
        while(nextLawPos!=0)
    }

    //Create capture FS files
    MakeCaptureFSFile(outDirName, isoName, isoMass, temperature, pCSVec, angDistP, enerDistP, MTRPList, extractOrderP, energyAngVecP, ascii);
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

    // now extract delayed info
    NYield1DTab *dNYield;

    if(startDelayedNUBlock!=0)
    {
        for(;count<(startDelayedNUBlock); count++)
        {
            stream >> dummy;
        }
        stream >> dummy; count ++;

        dNYield = new NYield1DTab();
        dNYield->ExtractMCNPData(stream, count);
    }

    int dNEnDistPos;
    int nextLawPos, enDisLawDataPos;
    vector<int> enDisLawND, enDisNumLawApplNRegND, enDisNumLawApplNEnND;
    vector<int*> enDisSchemeVecND, enDisRangeVecND;
    vector<double*> enDisEnApplVecND, enDisProbApplVecND;
    vector<EnergyDist*> enerDistND;
    vector<AngularDist*> angDistND;
    if(startDNEDLBlock!=0)
    {
        for(;count<(startDNEDLBlock); count++)
        {
            stream >> dummy;
        }
        stream >> dNEnDistPos; count ++;

        for(;count<(startDNEDBlock+dNEnDistPos-1); count++)
        {
            stream >> dummy;
        }

        do
        {
            stream >> nextLawPos >> intTemp >> enDisLawDataPos; count=count+3;
            enDisLawND[extractOrderP[i]].push_back(intTemp);
            stream >> intTemp; count++;
            enDisNumLawApplNRegND.push_back(intTemp);

            enDisSchemeVecND.push_back(new int [enDisNumLawApplNRegND.back()]);
            enDisRangeVecND.push_back(new int [enDisNumLawApplNRegND.back()]);

            for(int j=0; j<enDisNumLawApplNRegND.back(); j++, count++)
            {
                stream >> intTemp;
                (enDisSchemeVecND.back())[j]=intTemp;
            }

            for(int j=0; j<enDisNumLawApplNRegND.back(); j++, count++)
            {
                stream >> intTemp;
                (enDisRangeVecND.back())[j]=intTemp;
            }

            stream >> intTemp; count++;
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

            for(;count<(enDisLawDataPos+startDNEDBlock-1); count++)
            {
                stream >> dummy;
            }
            // Now we extract the data depending on the format of the law using a different class for each distinct case
            // Don't have to worry about energy ordering of the laws, the law accuracy probability
            // ensures that the law will only be used in the right energy regime.

            if(enDisLawND==1)
                enerDistND.push_back(new EnerDistEqPEnerBins());

            else if(enDisLawND==2)
                enerDistND.push_back(new EnerDist1PhEner(isoMass));

            else if(enDisLawND==3)
                enerDistND.push_back(new EnerDistLevScat());

            else if(enDisLawND==4)
                enerDistND.push_back(new EnerDistConTab());

            else if(enDisLawND==5)
                enerDistND.push_back(new EnerDistGenEvapSpec());

            else if(enDisLawND==7)
                enerDistND.push_back(new EnerDistMaxwellFisSpec());

            else if(enDisLawND==9)
                enerDistND.push_back(new EnerDistEvapSpec());

            else if(enDisLawND==11)
                enerDistND.push_back(new EnerDistWattSpec());

            else if(enDisLawND==22)
                enerDistND.push_back(new EnerDistTabLinFunc(startDNEDBlock));

            else if(enDisLawND==24)
                enerDistND.push_back(new EnerDistTabMulti());

            else if(enDisLawND==44)
                angDistND.push_back(new AngEnDistKallbach(startDNEDBlock));

            else if(enDisLawND==61)
                angDistND.push_back(new AngEnDist3DTab(startDNEDBlock));

            else if(enDisLawND==66)
                angDistND.push_back(new AngEnDistNBody());

            else if(enDisLawND==67)
                angDistND.push_back(new AngEnDistLab3DTab());
            else
            {
                cout << "\n### Error: Energy law not recognized! ###" << endl;
                continue;
            }

            if(enDisLawND<44)
                enerDistND.back()->ExtractMCNPData(stream, count);

            else
                angDistND.back()->ExtractMCNPData(stream, count);

            // go to next law position
            for(;count<(nextLawPos+startDNEDBlock-1); count++)
            {
                stream >> dummy;
            }
        }
        //this assumes that the energy distribution data collection is terminated when nextLawPos=0, may not be true check
        while(nextLawPos!=0)
    }

    //Create Fission FS files
    MakeFissionFSFile(outDirName, isoName, isoMass, temperature, enDisLaw, enDisNumLawApplNReg, enDisNumLawApplNEn, enDisSchemeVec, enDisRangeVecpCSVec, enDisEnApplVec, enDisProbApplVec, enerDist,
                        angDistP, enerDistP, MTRPList, extractOrderP, energyAngVecP, ascii);

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
            if(enerDist[i][j]!=NULL)
                delete enerDist[i][j];
        }
    }

    if(dNPromptYieldDist!=NULL)
        delete dNPromptYieldDist;
    if(dNTotalYieldDist!=NULL)
        delete dNTotalYieldDist;
    if(exOrderAng!=NULL)
        delete[] exOrderAng;

    //delete out going photon data
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
    if(exOrderAngP!=NULL)
        delete [] exOrderAngP;
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
        for(int j=0; j<int(enDisEnApplVecP[i].size()); j++)
        {
            if(enDisEnApplVecP[i][j]!=NULL)
                delete [] enDisEnApplVecP[i][j];
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

    //delete out going delayed neutron data
    if(dNYield)
        delete dNYield;

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
    for(int i=0; i<int(angDistND.size());i++)
    {
        if(angDistND[i]!=NULL)
            delete angDistND[i];
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

void MakeCSDataFile(string fileName, int MTRNum, CSDist *csDist, bool ascii)
{
    std::stringstream stream;

    stream.fill(' ');
    stream << std::setw(14) << std::right << MTRNum << '\n';
    WriteG4NDLData(stream);
    SetDataStream( fileName, stream, ascii);

    stream.clear();
    stream.str("");
}

void MakeElasFSFile(string fileName, double isoMass, double temperature, int TYRList[], int numAngEner[], double *energyAngVec[], int *angDistType[], int *intSchemeAng[], int *numAngProb[], double **angVec[], double **angProbVec[], double energyCSVec[], bool ascii)
{
    std::stringstream stream;
    int frameFlag=1, repFlag=2; //assume that elastic data is collected in LAB frame, may need to compare values to G4NDL to verify that this is correct
    double cosine=-1;

    stream.fill(' ');
    stream << std::setw(14) << std::right << repFlag;
    stream << std::setw(14) << std::right << isoMass;
    stream << std::setw(14) << std::right << frameFlag;
    stream << std::setw(14) << std::right << numAngEner[0] << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << numAngEner[0] << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

    for(int i=0; i<int(angDist[0].size()); i++)
    {
        angDist[0][i].WriteG4NDLData(stream)
    }

    stream << '\n';

    SetDataStream( fileName, stream, ascii);

    stream.clear();
    stream.str("");
}

void MakeCaptureFSFile(string outDirName, string isoName, double isoMass, double temperature, int TYRList[], CSDist* pCSVec, vector<AngularDist*> *angDistP, vector<EnergyDist*> *enerDistP, vector<int> MTRPList, int *extractOrderP, int *numAngEnerP, CSDist nCSVec, bool ascii)
{
    std::stringstream stream;
    int frameFlag=2, repFlag=2, nDiscrete=0;
    double cosine=-1;
    vector<int> capPReacIndex;

    if(TYRList[1]>0)
        frameFlag=1;

    stream.fill(' ');

    for(int i=0; i<int(MTRPList.size()); i++)
    {
        if(int(MTRPList[i]/1000)==102)
        {
            capPReacIndex.push_back(i);
            nDiscrete++;
        }
    }
    // I use repFlag 1 no matter what
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << isoMass << std::setw(14) << std::right << nDiscrete;
    for(int i=0; i<nDiscrete; i++)
    {
        if(enerDistP[extractOrderP[capPReacIndex[i]]].front()->IdentifyYourself()=="EnerDist1PhEner")
            stream << std::setw(14) << std::right << 0;
        else
            stream << std::setw(14) << std::right << 1;

        stream << std::setw(14) << std::right << enerDistP[extractOrderP[capPReacIndex[i]]].front()->GetAverageOutEnergy();

        pCSVec[extractOrderP[capPReacIndex[i]]].WriteG4NDLData(stream);
    }

    stream << '\n';
    stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << 2 << std::setw(14) << std::right << nDiscrete << std::setw(14) << std::right << 0 << endl;

    for(int i=0; i<nDiscrete; i++)
    {
        stream << std::setw(14) << std::right << enerDistP[extractOrderP[capPReacIndex[i]]].front()->GetAverageOutEnergy() << std::setw(14) << std::right << 0. << endl;
        stream << std::setw(14) << std::right << numAngEnerP[extractOrderP[capPReacIndex[i]]] << '\n';
        for(int j=0; j<int(angDistP[extractOrderP[capPReacIndex[i]]].size()); j++)
        {
            angDistP[extractOrderP[capPReacIndex[i]]][j]->WriteG4NDLData(stream);
        }
    }

    stream << '\n';

    if(enerDistP[extractOrderP[capPReacIndex[i]]].front()->IdentifyYourself()!="EnerDist1PhEner")
    {
        stream << std::setw(14) << std::right << nDiscrete << std::setw(14) << std::right << 0 << endl;
        pCSVec[extractOrderP[capPReacIndex[i]]].SetCSData(nCSVec.GetEnergyVec(), nCSVec.GetCSVec(), nCSVec.GetCSSize());
        pCSVec[extractOrderP[capPReacIndex[i]]].WriteG4NDLCSData(stream);
        stream << '\n';
        enerDistP[extractOrderP[capPReacIndex[i]]].front()->WriteG4NDLData(stream);
    }
    string fileName = outDirName+"Capture/FS/"+isoName;
    SetDataStream( fileName, stream, ascii);

    stream.clear();
    stream.str("");

    // we do not make the FSMF6 directory since the MCNP does not seem to contain 3D tables of energy, angle and probability for photon production
    // we must check to make sure that this is true
    /*
    stream << isoMass << frameFlag << nDiscrete << endl;

    for(int i=0; i<nDiscrete; i++)
    {
        stream << std::setw(14) << std::right << (isoMass*1.00867) << std::setw(14) << std::right << 0. << 0 << endl;
    }
    */

}

void MakeFissionFSFile(outDirName, isoName, isoMass, temperature, enDisLaw, enDisNumLawApplNReg, enDisNumLawApplNEn, enDisSchemeVec, enDisRangeVec, enDisEnApplVec, enDisProbApplVec, enerDist,
                    numAngEner, angDist, pCSVec, angDistP, enerDistP, MTRPList, extractOrderP, energyAngVecP, ascii);
{
    std::stringstream stream;
    int frameFlag=2, repFlag=2, nDiscrete=0;
    double cosine=-1;
    vector<int> fisPReacIndex;

    if(TYRList[2]>0)
        frameFlag=1;

    stream.fill(' ');

    for(int i=0; i<int(MTRPList.size()); i++)
    {
        if(int(MTRPList[i]/1000)==18)
        {
            fisPReacIndex.push_back(i);
            nDiscrete++;
        }
    }

    // FS/ directory

    //prompt neutrons

    //Angular Distribution
    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 4;
    stream << std::setw(14) << std::right << repFlag;
    stream << std::setw(14) << std::right << isoMass;
    stream << std::setw(14) << std::right << frameFlag;
    stream << std::setw(14) << std::right << numAngEner[2] << '\n';
    stream << std::setw(14) << std::right << 1 << '\n';
    stream << std::setw(14) << std::right << numAngEner[2] << 2 << '\n'; // assuming linear scheme here based off G4NDL examples

    for(int i=0; i<int(angDist[2].size()); i++)
    {
        angDist[2][i].WriteG4NDLData(stream)
    }
    stream << '\n' << endl;
    //Energy Distribution

    stream << std::setw(14) << std::right << 1 << std::setw(14) << std::right << 5;
    stream << std::setw(14) << std::right << 0 << std::setw(14) << std::right << enDisLaw[2].size() << endl;

    for(int i=0; i<enDisLaw[2].size(); i++)
    {
        if(enDisLaw[2]==1)
            stream << std::setw(14) << std::right << 1 << '\n';
        else if(enDisLaw[2]==3)
        {
            stream << std::setw(14) << std::right << 1 << '\n';
        }
        else if(enDisLaw[2]==4)
        {
            stream << std::setw(14) << std::right << 1 << '\n';
            cout << "###: not sure if this distribution is meant for out-going neutron energies" << endl;
        }
        else if(enDisLaw[2]==5)
        {
            stream << std::setw(14) << std::right << 5 << '\n';
        }
        else if(enDisLaw[2]==7)
        {
            stream << std::setw(14) << std::right << 7 << '\n';
        }
        else if(enDisLaw[2]==9)
        {
            stream << std::setw(14) << std::right << 9 << '\n';
        }
        else if(enDisLaw[2]==11)
        {
            stream << std::setw(14) << std::right << 11 << '\n';
        }
        else if(enDisLaw[2]==22)
        {
            cout << "###: No direct translation for this law!" << endl;
            stream << std::setw(14) << std::right << 1 << '\n';
        }
        else if(enDisLaw[2]==24)
        {
            cout << "###: No direct translation for this law!" << endl;
            stream << std::setw(14) << std::right << 1 << '\n';
        }

        stream << std::setw(14) << std::right << enDisNumLawApplNEn[2][i] <<'\n';
        stream << std::setw(14) << std::right << enDisNumLawApplNReg[2][i] <<'\n';
        for(int j=0; j<enDisNumLawApplNReg[2][i]; j++)
        {
            stream << std::setw(14) << std::right << enDisRangeVec[2][i][j] << std::right << enDisSchemeVec[2][i][j] << '\n';
        }
        for(int j=0; j<enDisNumLawApplNEn[2][i]; j++)
        {
            stream << std::setw(14) << std::right << enDisEnApplVec[2][i][j] << std::right << enDisProbApplVec[2][i][j] << '\n';
        }

        enerDist[2][i].WriteG4NDLData(stream);
    }
    stream << '\n';


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


void SetDataStream( string filename , std::stringstream& ss, bool ascii )
{
    //bool cond=true;
   if (!ascii)
   {
        string compfilename(filename);

        if(compfilename.back()!='z')
            compfilename += ".z";

       std::ofstream* out = new std::ofstream ( compfilename.c_str() , std::ios::binary | std::ios::trunc);
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

      std::ofstream out( compfilename.c_str() , std::ios::out | std::ios::trunc );
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
         out.write(filedata, file_size);
         if (out.fail())
        {
            cout << endl << "writing the ascii data to the output file " << compfilename << " failed" << endl
                 << " may not have permission to delete an older version of the file" << endl;
        }
         out.close();
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );

         cout << endl << "### failed to write to ascii file " << compfilename << " ###" << endl;
      }
   }
   ss.str("");
}
