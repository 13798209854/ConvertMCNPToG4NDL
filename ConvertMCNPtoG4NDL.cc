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

int CreateIsoCSData(stringstream &stream, string outDirName, bool ascii);
bool DirectoryExists( const char* pzPath );
void MakeCSDataFile(string fileName, std::stringstream& stream, int MTRNum, double energyVec[], double CSVec[], int vecSize, int startEnergy, bool ascii);
void GetDataStream( string, std::stringstream&);
void SetDataStream( string, std::stringstream&, bool);

//Maybe an error in the extraction of inelastic, only a few of the isotopes seem to have the identifier it is looking for

//Takes in a directory of MCNP cross-section libraries, converts the data into the G4NDL format and then outputs the information in a given directory
int main(int argc, char **argv)
{
    ElementNames elementNames;
    elementNames.SetElementNames();

    string outputType="ascii";
    bool ascii=true;
    string word, libName="endf";
    char lib, version='7', check1, check2, check3;
    string inFileName, outDirName, fileName;
    int result=0;

    stringstream stream;

    //Extracts user Inputs
    if(argc==5)
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
                                result += CreateIsoCSData(stream, outDirName, ascii);
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
int CreateIsoCSData(stringstream &stream, string outDirName, bool ascii)
{
    ElementNames *elementNames;
    char line[256];
    string temperature, isoName;
    string dummy, Z, A, outDirNameFis, outDirNameElas, outDirNameInE, outDirNameCap;
    string outDirNameProc[4];
    stringstream numConv;
    int numEner=0, numReactions=0;
    int startEnerTable, startMTRBlock;
    int startLSIGBlock, startCSBlock;
    int startElasticBlock;
    int startInElasEner, startFisEner, startCapEner;
    int Znum, isoNum;
    int index, index2, count=1;
    int fisCSVecSize=0, inECSVecSize=0, capCSVecSize=0;
    double temp;
    double *energyVec=NULL, *fisCSVec=NULL, *elasCSVec=NULL, *inECSVec=NULL, *capCSVec=NULL;
    std::vector<int> MTRList, LSIGList;

    //extracts the temperature
    stream >> dummy >> temp;
    temp=temp*1000000/(8.6173324*(pow(10, -5)));
    numConv << temp;
    numConv >> temperature;
    numConv.clear();
    numConv.str("");
    outDirName = outDirName+temperature+'/';

    outDirNameProc[2] = outDirName +"Fission/CrossSection/";
    outDirNameProc[1] = outDirName +"Elastic/CrossSection/";
    outDirNameProc[3] = outDirName +"Inelastic/CrossSection/";
    outDirNameProc[0] = outDirName +"Capture/CrossSection/";

    // creates temperature directories for the converted files to go if they don't already exist
    for(int i=0; i<4; i++)
    {
        if(!(DirectoryExists((outDirNameProc[i]).c_str())))
        {
            system( ("mkdir -p -m=666 "+outDirNameProc[i]).c_str());
            if(DirectoryExists((outDirNameProc[i]).c_str()))
            {
                cout << "created temperature directory " << outDirNameProc[i] << "\n" << endl;
            }
            else
            {
                cout << "\nError: could not create temperature Directory " << outDirNameProc[i] << "\n" << endl;
                if(energyVec!=NULL)
                    delete[] energyVec;
                if(elasCSVec!=NULL)
                    delete[] elasCSVec;
                if(capCSVec!=NULL)
                    delete[] capCSVec;
                if(inECSVec!=NULL)
                    delete[] inECSVec;
                if(fisCSVec!=NULL)
                    delete[] fisCSVec;
                return 1;
            }
        }
    }

    for(int i=0; i<6; i++)
    {
        stream.getline(line,256);
    }
    // Extracts the isotope name
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

    stream >> numEner >> numReactions;

    energyVec = new double[numEner];
    elasCSVec = new double[numEner];

    for(int i=0; i<2; i++)
    {
        stream.getline(line,256);
    }

    stream >> startEnerTable >> dummy >> startMTRBlock >> dummy >> dummy >> startLSIGBlock >> startCSBlock;

    startElasticBlock = startEnerTable+3*numEner;

    for(int i=0; i<4; i++)
    {
        stream.getline(line,256);
    }

    for(;count<startEnerTable; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numEner; i++, count++)
    {
        stream >> temp;
        energyVec[i]=temp;
    }

    for(;count<startElasticBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numEner; i++, count++)
    {
        stream >> temp;
        elasCSVec[i] = temp;
    }

    for(;count<startMTRBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        MTRList.push_back(index);
    }

    for(;count<startLSIGBlock; count++)
    {
        stream >> dummy;
    }

    for(int i=0; i<numReactions; i++, count++)
    {
        stream >> index;
        LSIGList.push_back(index);
    }

    for(int i=0; i<int(MTRList.size()); i++)
    {
        if((MTRList[i]==18)||(MTRList[i]==102)||(MTRList[i]==4))
        {

        }
        else
        {
            MTRList.erase(MTRList.begin()+i);
            LSIGList.erase(LSIGList.begin()+i);
            i--;
        }
    }

    for(int i=0; i<int(MTRList.size()); i++)
    {
        for(int j=i+1; j<int(MTRList.size()); j++)
        {
            if(LSIGList[j]<LSIGList[i])
            {
                index = MTRList[i];
                index2 = LSIGList[i];

                MTRList[i] = MTRList[j];
                LSIGList[i] = LSIGList[j];

                MTRList[j] = index;
                LSIGList[j] = index2;
            }
        }
    }

    for(int i=0; i<int(MTRList.size()); i++)
    {
        for(;count<(LSIGList[i]+startCSBlock-1); count++)
        {
            stream >> dummy;
        }
        if(MTRList[i]==18)
        {
            stream >> startFisEner;
            count++;
            stream >> fisCSVecSize;
            count++;
            fisCSVec = new double[fisCSVecSize];

            for(int j=0; j<fisCSVecSize; j++, count++)
            {
                stream >> temp;
                fisCSVec[j] = temp;
            }

        }
        else if(MTRList[i]==4)
        {
            stream >> startInElasEner;
            count++;
            stream >> inECSVecSize;
            count++;
            inECSVec = new double[inECSVecSize];

            for(int j=0; j<inECSVecSize; j++, count++)
            {
                stream >> temp;
                inECSVec[j] = temp;
            }
        }
        else if(MTRList[i]==102)
        {
            stream >> startCapEner;
            count++;
            stream >> capCSVecSize;
            count++;
            capCSVec = new double[capCSVecSize];

            for(int j=0; j<capCSVecSize; j++, count++)
            {
                stream >> temp;
                capCSVec[j] = temp;
            }
        }
    }

    stringstream ssCSFile;

    if(numEner!=0)
        MakeCSDataFile(outDirNameProc[1]+isoName, ssCSFile, 2, energyVec, elasCSVec, numEner, startEnerTable, ascii);
    if(fisCSVecSize!=0)
        MakeCSDataFile(outDirNameProc[2]+isoName, ssCSFile, 18, energyVec, fisCSVec, fisCSVecSize, startFisEner, ascii);
    if(inECSVecSize!=0)
        MakeCSDataFile(outDirNameProc[3]+isoName, ssCSFile, 4, energyVec, inECSVec, inECSVecSize, startInElasEner, ascii);
    if(capCSVecSize!=0)
        MakeCSDataFile(outDirNameProc[0]+isoName, ssCSFile, 102, energyVec, capCSVec, capCSVecSize, startCapEner, ascii);

    if(energyVec!=NULL)
        delete[] energyVec;
    if(elasCSVec!=NULL)
        delete[] elasCSVec;
    if(capCSVec!=NULL)
        delete[] capCSVec;
    if(inECSVec!=NULL)
        delete[] inECSVec;
    if(fisCSVec!=NULL)
        delete[] fisCSVec;

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

void MakeCSDataFile(string fileName, std::stringstream& stream, int MTRNum, double energyVec[], double CSVec[], int vecSize, int startEnergy, bool ascii)
{
    stream.clear();
    stream.str("");

    stream << MTRNum << '\n';
    stream << '0' << '\n';
    stream << vecSize << '\n';

    stream.fill(' ');

    for(int i=0; i<vecSize; i++)
    {
        if(floor(i/2)==ceil(double(i)/2))
        {
            stream << '\n';
        }
        stream << std::setw(14) << std::left << energyVec[i+startEnergy-1];
        stream << std::setw(14) << std::left << CSVec[i];
    }
    stream << '\n';

    SetDataStream( fileName, stream, ascii);

    stream.clear();
    stream.str("");
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
