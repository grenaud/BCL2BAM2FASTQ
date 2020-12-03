/*
 * bcl2bam
 * Date: Aug-16-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

//TODO


#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <errno.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "libgab.h"

#define MAXCYCLES 10000

using namespace std;
using namespace BamTools;


//constants taken from http://www.umanitoba.ca/afs/Plant_Science/psgendb/local/pkg/CASAVA_v1.8.2/src/c++/include/alignment/BclReader.hh
static const int blockSize      = 25;
static const int imageWidth     = 2048;
static const int blocksPerlLine = (imageWidth + blockSize - 1) / blockSize;
static const int imageHeight    = 20000;

const uint32_t flagSingleReads =  4; // 00000100
const uint32_t flagFirstPair   = 77; // 01001101
const uint32_t flagSecondPair  =141; // 10001101

static const int offseqQualScores=33;

typedef struct {
    unsigned int x;
    unsigned int y;    
} tileCoords;
	





inline void reverseCstr(char * inputString,int sizeS){
    char  tempbuffer [sizeS-1];

    if(sizeS >1 ){	

	for(int i=(sizeS-1);i>=1;i--){
	    
	    tempbuffer[sizeS-1-i] =( inputString[i-1] );
	    // cout<<"tempbuffer["<<(sizeS-1-i)<<"]  inputString["<<i-1<<"] )\t"<< inputString[i-1]<<endl;
	}

	for(int i=0;i<(sizeS-1);i++){	    
	    inputString[i]  = tempbuffer[i];
	    // cout<<"tempbuffer["<<(i)<<"] \t( inputString["<<i<<"] )"<<tempbuffer[i]<<endl;
	}
	inputString[sizeS-1]='\0';
	// cout<<string(inputString)<<endl;
    }
    //    return toReturn;
}


inline char complementBase(const char c){
    if(c ==    'A')
	return 'T';

    if(c ==    'C')
	return 'G';

    if(c ==    'G')
	return 'C';

    if(c ==    'T')
	return 'A';

    if(c ==    'N')
	return 'N';


    if(c ==    'a')
	return 't';

    if(c ==    'c')
	return 'g';

    if(c ==    'g')
	return 'c';

    if(c ==    't')
	return 'a';




    if(c ==    'n')
	return 'n';

    cerr<<"Complement: Invalid base pair="<<c<<endl;
    exit(1);
}



inline void reverseCstrRC(char * inputString,int sizeS){
    char  tempbuffer [sizeS-1];

    if(sizeS >1 ){	

	for(int i=(sizeS-1);i>=1;i--){
	    
	    tempbuffer[sizeS-1-i] =complementBase( inputString[i-1] );
	    // cout<<"tempbuffer["<<(sizeS-1-i)<<"]  inputString["<<i-1<<"] )\t"<< inputString[i-1]<<endl;
	}

	for(int i=0;i<(sizeS-1);i++){	    
	    inputString[i]  = tempbuffer[i];
	    // cout<<"tempbuffer["<<(i)<<"] \t( inputString["<<i<<"] )"<<tempbuffer[i]<<endl;
	}
	inputString[sizeS-1]='\0';
	// cout<<string(inputString)<<endl;
    }
    //    return toReturn;
}


typedef struct { 
    char baseC;
    int  qualC;
} base;

bool str2intcmp (string i,string j) { 
    int iInt=destringify<int>(i);
    int jInt=destringify<int>(j);

    return (iInt<jInt); 
}

inline base bin2base(unsigned char c ){
    base toreturn;
    toreturn.baseC="ACGT"[c % 4];
    c /= 4;
    toreturn.qualC=(c % 64);
    return toreturn;
}
// char *chartobin ( unsigned char c ){
//     static char bin[8 + 1] = {0};
//     int         i;
//     for ( i = 8 - 1; i >= 0; i-- ){
//         bin[i] = (c % 2) + '0';
// 	c /= 2;
//     }
//     return bin;
// }


inline double round(double n){
    if(n < 0.0)
	return ceil(n - 0.5);
    else
	return floor(n + 0.5);
}

//we should not support format formats anymore
//inline int do_round( int coordtype_, float x ){
inline int do_round( float x ){
    // switch( coordtype_ )
    // {
    //     case 1: return int( round( x ) ) ;
    //     case 2: return int( floor( fabs( x ) ) ) ;
    //     case 3:
    return int( round( 10 * x + 1000 ) ) ;
    //     default: throw "unknown coordinate format" ;
    // }
}

int main (int argc, char *argv[]) {
    vector<int> lanesToUse;
    int forwardCycles=0;
    int reverseCycles=0;
    int index1Cycles=0;
    int index2Cycles=0;
    string nameOfExperiment="SOLEXA";
    string pathToPositionFiles="";
    string bclDirectory   = "";
    string bamfiletowrite = "";
    bool onlyIndex=false;
    bool fiir=false; //if the order was forward, i1, i2 and reverse

    bool noFinishedFlag=false;
    string suffixFinishedFlag=".finished";

    string indexToExtract="";
    vector<int> tilesToUse;
    bool force=false;
    bool fakePos=false;
    
    string usage=string(""+string(argv[0])+" options  "+
			"\nThis program reads a BCL directory and produces a BAM file\n"+
			"\nOptions:\n"+
			"\tMandatory options:\n"+
			"\t\t"+"-f"+" " "--forward"+"\t\t\t"+"Number of cycles for the forward read ex:76 (default: none)\n"+
			"\t\t"+"-r"+" " "--reverse"+"\t\t\t"+"Number of cycles for the reverse read ex:76, use 0 if none (default: none)\n"+

			"\t\t"+"-i"+" " "--index1" +"\t\t\t"+"Number of cycles for the first index  ex:7 use 0 if none   (default: none)\n"+
			"\t\t"+"-j"+" " "--index2" +"\t\t\t"+"Number of cycles for the second index ex:7 use 0 if none   (default: none)\n"+

			"\t\t"+"-p"+" " "--path"+"\t\t\t"+"Path to position files (default: none)\n"+
			"\t\t"+"-b"+" " "--bcl"+"\t\t\t"+"Path to bcl files (default: none)\n"+
			"\t\t"+"-o"+" " "--outbam"+"\t\t\t"+"output file (BAM format)  (default: none)\n"+
			
			"\tOptional:\n"+
			"\t\t"+""+" "+"--fiir"+"\t\t\t"+"If the run was forward, i1, i2 (reverse) we also rev. comp. i2  \n"+
			"\t\t"+""+" "+""+"\t\t\t"+"The default is : forward, i1, reverse, i2 \n"+

			"\t\t"+""+" " "--fakepos"+"\t\t\t"+"Output fake positions, not recommended (default: "+booleanAsString(fakePos)+")\n"+
			"\t\t"+"-e"+" " "--exp"+"\t\t\t"+"Name of experiment (default: "+nameOfExperiment+")\n"+
			"\t\t"+"-l"+" " "--lanes"+"\t\t\t"+"Colon separated list of lanes to use ex:1,3 (default: all)\n"+
			"\t\t"+"-t"+" " "--tiles"+"\t\t\t"+"Colon separated list of tiles to use ex:1112,3212 (default: all)\n"+
			"\t\t"+"-s"+" " "--indexseq"+"\t\t\t"+"Only extract the subset of sequences with this string as first index (default: return everything)\n"+
			"\t\t"+"--noflag"+"\t\t\t"+"Do not print a "+suffixFinishedFlag+" file when the program finished gracefully\n"+
			"\t\t"+"-z"+"\t\t\t"+"(force) Do not exit if the number of cycles found is not the one that the user entered, just warn\n"+

			"\n");
    
    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //parse options
    //use destringify
    for(int i=1;i<(argc);i++){ 

	if( (strcmp(argv[i],"--noflag") == 0)  ){
	    suffixFinishedFlag=true;
            continue;
	}

	if( (strcmp(argv[i],"--fiir") == 0)  ){
	    fiir=true;
            continue;
	}

	if( (strcmp(argv[i],"--fakepos") == 0)  ){
	    fakePos=true;
            continue;
	}
	
	if( (strcmp(argv[i],"-s") == 0) || (strcmp(argv[i],"--indexseq") == 0)  ){
	    onlyIndex=true;
	    indexToExtract=string(argv[i+1]);
	    i++;
            continue;
	}


	if( (strcmp(argv[i],"-l") == 0) || (strcmp(argv[i],"--lanes") == 0)  ){
	    vector<string> temp=allTokens(string(argv[i+1]),',');
	    for(int k=0;k<int(temp.size());k++)
		lanesToUse.push_back( destringify<int>( temp[k]) );
	    i++;
            continue;
        }

	if( (strcmp(argv[i],"-t") == 0) || (strcmp(argv[i],"--tiles") == 0)  ){
	    vector<string> temp=allTokens(string(argv[i+1]),',');
	    for(int k=0;k<int(temp.size());k++)
		tilesToUse.push_back( destringify<int>( temp[k]) );
	    i++;
            continue;
        }


	if( (strcmp(argv[i],"-o") == 0) || (strcmp(argv[i],"--outbam") == 0)  ){
	    bamfiletowrite=string(argv[i+1]);
	    i++;
            continue;
        }

	if( (strcmp(argv[i],"-e") == 0) || (strcmp(argv[i],"--exp") == 0)  ){
	    nameOfExperiment=string(argv[i+1]);
	    i++;
            continue;
        }

	//to force is zwingen in german, hidden option to transform the error of number 
	//of cycles into a warning, USE AT YOUR OWN RISK
	if( strcmp(argv[i],"-z") == 0  ){
	    force=true;
            continue;
        }



	if( (strcmp(argv[i],"-b") == 0) || (strcmp(argv[i],"--bcl") == 0)  ){
	    bclDirectory=string(argv[i+1]);
	    i++;
            continue;
        }

        if( (strcmp(argv[i],"-p") == 0) || (strcmp(argv[i],"--path") == 0)  ){
	    pathToPositionFiles=string(argv[i+1]);
	    i++;
            continue;
        }

        if( (strcmp(argv[i],"-f") == 0) || (strcmp(argv[i],"--forward") == 0)  ){
	    forwardCycles=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

        if( (strcmp(argv[i],"-r") == 0) || (strcmp(argv[i],"--reverse") == 0)  ){
	    reverseCycles=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

        if( (strcmp(argv[i],"-i") == 0) || (strcmp(argv[i],"--index1") == 0)  ){
	    index1Cycles=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

        if( (strcmp(argv[i],"-j") == 0) || (strcmp(argv[i],"--index2") == 0)  ){
	    index2Cycles=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }


        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           
    }


    if(bamfiletowrite == "" ){ cerr<<"The output bam file must be specified"<<endl;    return 1;  }
    if(bclDirectory == "" ){   cerr<<"The folder of the BCL files must be specified"<<endl;    return 1;  }
    if(pathToPositionFiles == "" ){ if(!fakePos){ cerr<<"The folder of the position files must be specified"<<endl;    return 1; } }
    if(forwardCycles == 0 ){ cerr<<"The number of cycles for the forward read must be specified"<<endl;    return 1;  }
    // if(reverseCycles == -1 ){ cerr<<"The number of cycles for the reverse read must be specified"<<endl;    return 1;  }
    // if(index1Cycles  == -1 ){ cerr<<"The number of cycles for the first index must be specified"<<endl;    return 1;  }
    // if(index2Cycles  == 0 ){ cerr<<"The number of cycles for the second index must be specified"<<endl;    return 1;  }
    if(fiir){
	 if(index1Cycles  == 0 ){ cerr<<"The number of cycles for the first index must be specified when --fiir is used"<<endl;    return 1;  }
	 if(index2Cycles  == 0 ){ cerr<<"The number of cycles for the second index must be specified when --fiir is used"<<endl;    return 1;  }
    }

    if(onlyIndex){
	if(int(indexToExtract.length()) != index1Cycles){ cerr<<"The number of cycles for the first index is not equal to the size of the index you entered"<<endl;    return 1;  }
    }
    //parse cmd line arguments

    bool isPairedEnd=false;
    if(reverseCycles > 0)
	isPairedEnd=true;
    

    if(!isDirectory(bclDirectory)){
	cerr<<"Directory  "<<bclDirectory <<" is not a directory, exiting"<<endl;
        return 1;    
    }

    if(!fakePos)
    if(!isDirectory(pathToPositionFiles)){
	cerr<<"The folder of the position files "<<pathToPositionFiles <<" is not a directory, exiting"<<endl;
        return 1;    
    }

    //try to detect RTA version
    string versionRTA="unknown";
    string configXML=pathToPositionFiles+"/config.xml";
    if(isFile(configXML)){
	string lineXML;
	ifstream myRTAconfFile;
	myRTAconfFile.open(configXML.c_str(), ios::in);

	if (myRTAconfFile.is_open()){
	    while ( getline (myRTAconfFile,lineXML)){
		vector<string> testxml= allTokens(lineXML,' '); 
		for(unsigned int indexxml=0;indexxml<testxml.size();indexxml++){

		    //looking for     <Software Name="RTA" Version="1.16.18.0" />
		    if( (testxml[indexxml]  == "<Software") &&
		        (        indexxml   == (testxml.size()-4)) ){
			string versionstring=testxml[testxml.size()-2];

			if(strBeginsWith(versionstring,"Version=")){

			    vector<string> testxml2= allTokens(versionstring,'='); 
			    if(testxml2.size() == 2){

				versionRTA=testxml2[1].substr(1, testxml2[1].size() -2);
			    }
			}
		    }
		}
	    }
	    myRTAconfFile.close();
	}else{
	    cerr << "Unable to open file "<<configXML<<" will leave the RTA version as "<<versionRTA<<endl;
	}
    }

    //initialize BAM writer
    BamWriter writer;
    SamProgram bcl2bamProg ("BCL2BAM");
    bcl2bamProg.Version=versionRTA;
    SamProgramChain spcToAdd;
    spcToAdd.Add(bcl2bamProg);
    SamHeader header;
    header.Programs=spcToAdd;

    RefVector references;
    if( !writer.Open(bamfiletowrite,header,references) ) {
	cerr << "Could not open output BAM file "<<bamfiletowrite << endl;
	return 1;
    }



    //DETECT NUMBER OF CYCLES
    int numberOfCycles=-1;
    for(int cycle=1;cycle<=MAXCYCLES;cycle++){
	string lane="1";
	if(lanesToUse.size()!=0){
	    lane=stringify( lanesToUse[0] );
	}
	string bclFile=bclDirectory+"/L00"+lane+"/C"+ stringify(cycle)+".1/";

	if(isDirectory(bclFile)){
	    numberOfCycles=cycle;
	}else{
	    break;
	}
    }

    cerr<<"Found "<<numberOfCycles<<" cycles"<<endl;

    if(numberOfCycles != (forwardCycles +  reverseCycles + index1Cycles + index2Cycles)){
	
	cerr<<"The total number of cycles specified "<<(forwardCycles +  reverseCycles + index1Cycles + index2Cycles)<<" does not correspond to the total number of cycles detected "<<numberOfCycles<<endl;
	if(!force){
	    return 1;    
	}else{
	    numberOfCycles = (forwardCycles +  reverseCycles + index1Cycles + index2Cycles);
	    cerr<<"Using hidden -z, setting the number of cycles to "<<numberOfCycles<<endl;
	}
    }


    //mark every cycle as being forward, reverse, index1, index2
    //not the sexiest C++ code but should work
    char cycleCode [numberOfCycles];
    if(fiir){
	for(int i=0;i<forwardCycles;i++){
	    cycleCode[i]='f';
	}

	for(int i=forwardCycles;i<(forwardCycles+index1Cycles);i++){
	    cycleCode[i]='i';
	}

	for(int i=(forwardCycles+index1Cycles);i<(forwardCycles+index1Cycles+index2Cycles);i++){
	    cycleCode[i]='j';
	}
    
	for(int i=(forwardCycles+index1Cycles+index2Cycles);i<(forwardCycles+index1Cycles+index2Cycles+reverseCycles);i++){
	    cycleCode[i]='r';
	}

    }else{
	for(int i=0;i<forwardCycles;i++){
	    cycleCode[i]='f';
	}

	for(int i=forwardCycles;i<(forwardCycles+index1Cycles);i++){
	    cycleCode[i]='i';
	}

	for(int i=(forwardCycles+index1Cycles);i<(forwardCycles+index1Cycles+reverseCycles);i++){
	    cycleCode[i]='r';
	}
    
	for(int i=(forwardCycles+index1Cycles+reverseCycles);i<(forwardCycles+index1Cycles+reverseCycles+index2Cycles);i++){
	    cycleCode[i]='j';
	}
    }
    // cout<<cycle<<endl;
    // return 1;
    // for(int i=(0);i<(numberOfCycles);i++){
    // 	cout<<i<<"\t"<<cycleCode[i]<<endl;
    // }
    // return 1;


    //DETECTING TILE "ARRANGEMENT"
    string laneToUseForDetection ="1";  //using first lane , should always exists
    string cycleToUseForDetection="1"; //using first cycle, should always exists
    if(lanesToUse.size()!=0){
	laneToUseForDetection=stringify( lanesToUse[0] );	
    }
    string bclDir=bclDirectory+"/L00"+laneToUseForDetection+"/C"+cycleToUseForDetection+".1/"; 

    vector<string> testTiles=getdir(bclDir);

    if(tilesToUse.empty()){ //if the user has not defined which tiles to use, autodetect, otherwise trust the user
	for(int i=0;i<int(testTiles.size());i++){
	    string requiredSuffix=".bcl";
	    if(strEndsWith(testTiles[i],requiredSuffix)){
		string requiredPrefix="s_"+laneToUseForDetection+"_";
		if(strBeginsWith(testTiles[i],requiredPrefix)){
		    string tileTocheck=testTiles[i].substr(requiredPrefix.length(), testTiles[i].length()-requiredSuffix.length()-requiredPrefix.length());
		    int tileNumber=destringify<int>(tileTocheck);
		    tilesToUse.push_back(tileNumber);
		}else{
		    cerr<<"File  "<< testTiles[i] <<" in "<<bclDir<<" should begin with "<<requiredPrefix<<endl;
		    return 1;
		}
	    }else{

		requiredSuffix=".bcl.gz";
		if(strEndsWith(testTiles[i],requiredSuffix)){
		    string requiredPrefix="s_"+laneToUseForDetection+"_";
		    if(strBeginsWith(testTiles[i],requiredPrefix)){
			string tileTocheck=testTiles[i].substr(requiredPrefix.length(), testTiles[i].length()-requiredSuffix.length()-requiredPrefix.length());
			int tileNumber=destringify<int>(tileTocheck);
			tilesToUse.push_back(tileNumber);
		    }else{
			cerr<<"File  "<< testTiles[i] <<" in "<<bclDir<<" should begin with "<<requiredPrefix<<endl;
			return 1;
		    }
		}


	    }
	}

	sort (tilesToUse.begin(), tilesToUse.end());
	cerr<<"Found the following tiles: "<<vectorToString(tilesToUse)<<endl;

    }else{ //user defined tiles
	sort (tilesToUse.begin(), tilesToUse.end());
	cerr<<"Will use the following tiles: "<<vectorToString(tilesToUse)<<endl;
    }





    //DETECTING LANES TO USE
    if( lanesToUse.empty() ){ //The user has not specified which lanes to use
	for(int putativeLane=1;putativeLane<=8;putativeLane++){

	    bclDir=bclDirectory+"/L00"+stringify(putativeLane)+"/C"+cycleToUseForDetection+".1/"; 
	    if(isDirectory(bclDir)){
		lanesToUse.push_back(putativeLane);
	    }

	}
    }else{

	for(int i=0;i<int(lanesToUse.size());i++){
	    if(lanesToUse[i] < 1 || lanesToUse[i] > 8 ){
		cerr<<"Lanes must be numerical and between 1 and 8"<<endl;
		return 1;
	    }
	}
	     
    }
    

    for(int laneIndex=0;laneIndex<int(lanesToUse.size());laneIndex++){ //for each lane
	//cerr<<"Processing lane #"<<lanesToUse[laneIndex]<<endl;
	for(int tileIndex=0;tileIndex<int(tilesToUse.size());tileIndex++){ //for each tile
	    // cerr<<"Processing tile #"<<tilesToUse[tileIndex]<<endl;

	    igzstream * mybclfile= new igzstream[numberOfCycles];
	    unsigned int numberOfClustersFirst=0 ;

	    //Opening the BCL files and opening file pointers
	    for(int cycle=1;cycle<=numberOfCycles;cycle++){ //for each cycle
		string bclFile=bclDirectory+"/L00"+ stringify(lanesToUse[laneIndex]) +"/C"+stringify(cycle)+".1/s_"+stringify(lanesToUse[laneIndex])+"_"+stringify(tilesToUse[tileIndex])+".bcl"; 

		if(isFile(bclFile)){
		
		    mybclfile[cycle-1].open(bclFile.c_str(),ios::in|ios::binary);
		    if (!mybclfile[cycle-1]) {
			cerr<<"Unable to read BCL file "<<bclFile<<" "<<strerror(errno)<<endl;
			return 1;
		    }
		    
		    unsigned int numberOfClusters ;
		    mybclfile[cycle-1].read((char*)&numberOfClusters, sizeof (int));

		    if(numberOfClustersFirst == 0 ){
			numberOfClustersFirst=numberOfClusters;
			//cerr<"Found "<<numberOfClustersFirst<<" clusters "<<endl;
			cerr<<"Processing lane #"<<lanesToUse[laneIndex]<<" tile #"<<tilesToUse[tileIndex]<<" with "<<numberOfClustersFirst<<" clusters"<<endl;
		    }else{
			if(numberOfClustersFirst!=numberOfClusters){
			    cerr<<"Number of clusters in "<<bclFile<<" is different from the other cycles, exiting\n";
			    return 1;
			}
		    }

	
		}else{
		    bclFile=bclDirectory+"/L00"+ stringify(lanesToUse[laneIndex]) +"/C"+stringify(cycle)+".1/s_"+stringify(lanesToUse[laneIndex])+"_"+stringify(tilesToUse[tileIndex])+".bcl.gz"; 

		    if(isFile(bclFile)){
			
			mybclfile[cycle-1].open(bclFile.c_str(),ios::in|ios::binary);
			if (!mybclfile[cycle-1]) {
			    cerr<<"Unable to read BCL file "<<bclFile<<" "<<strerror(errno)<<endl;
			    return 1;
			}
			
			unsigned int numberOfClusters ;
			mybclfile[cycle-1].read((char*)&numberOfClusters, sizeof (int));
			
			if(numberOfClustersFirst == 0 ){
			    numberOfClustersFirst=numberOfClusters;
			    //cerr<"Found "<<numberOfClustersFirst<<" clusters "<<endl;
			    cerr<<"Processing lane #"<<lanesToUse[laneIndex]<<" tile #"<<tilesToUse[tileIndex]<<" with "<<numberOfClustersFirst<<" clusters"<<endl;
			}else{
			    if(numberOfClustersFirst!=numberOfClusters){
				cerr<<"Number of clusters in "<<bclFile<<" is different from the other cycles, exiting\n";
				return 1;
			    }
			}
		    }else{			
			cerr<<"BCL File "<< bclFile<<" does not exist"<<endl;
			return 1;
		    }
		}
	    }//each cycle

	    //READ POSITION FILES
	    vector<tileCoords> coordinatesFound;

	    if( fakePos ){

		//if(coordinatesFound.size() != numberOfClustersFirst){
		for(unsigned int i = 0;i<numberOfClustersFirst;i++){
		    tileCoords toadd;
		    toadd.x=i;
		    toadd.y=i;
		    coordinatesFound.push_back(toadd);
		}
		
	    }else{
		string posFile;
		//POS FILE
		posFile=pathToPositionFiles+"/s_"+ stringify(lanesToUse[laneIndex]) +"_"+zeroPad(tilesToUse[tileIndex],4)+"_pos.txt";

		if(isFile(posFile)){
		    //cerr<<"exists "<<posFile<<endl;
		    ifstream myFile;
		    string line;

		    myFile.open(posFile.c_str(), ios::in);
		    if (myFile.is_open()){
			while ( getline (myFile,line)){
			    istringstream ss( line ) ;
			    tileCoords toadd;
			    float x;
			    float y;
			    if( ss >> x >> y ) {
				toadd.x=do_round(x);
				toadd.y=do_round(y);
				coordinatesFound.push_back(toadd);
			    }else{
				cerr<<"Invalid line "<<line<<" in file:  "<<posFile<<endl;
				return 1;
			    }

			}
		    }else{
			cerr<<"Cannot open file  "<<posFile<<endl;
			return 1;
		    }
		    myFile.close();

		}else{
		    //checking locs
		    posFile=pathToPositionFiles+"/L00"+ stringify(lanesToUse[laneIndex]) +"/s_"+stringify(lanesToUse[laneIndex])+"_"+stringify(tilesToUse[tileIndex])+".locs";

		    if(isFile(posFile)){
			ifstream myFile;
			string line;

			myFile.open(posFile.c_str(),ios::in|ios::binary);
			if (myFile.is_open()){
			    //header
			    uint32_t headerData[3] ;
			    myFile.read( (char*)headerData, 12 ) ;

			    if(headerData[2] != numberOfClustersFirst){
				cerr<<"Found "<<numberOfClustersFirst <<" clusters in the bcl files but found "<<headerData[2]<<" coordinates in "<<posFile<<endl;
				return 1;
			    }

			
			    tileCoords toadd;
			
			    float cc[2] ;
			    while(myFile.read( (char*)cc, 8 )){
				toadd.x=do_round( cc[0] ) ;
				toadd.y=do_round( cc[1] ) ;
				// cout<<toadd.x<<"\t"<<toadd.y<<endl;
				coordinatesFound.push_back(toadd);
			    }
			
			

			    //			nclusters_ = gunk[2] ;
			}else{
			    cerr<<"Cannot open file  "<<posFile<<endl;
			    return 1;
			}
			myFile.close();


		    }else{
			//checking clocs
			posFile=pathToPositionFiles+"/L00"+ stringify(lanesToUse[laneIndex]) +"/s_"+stringify(lanesToUse[laneIndex])+"_"+stringify(tilesToUse[tileIndex])+".clocs";
			//cerr<<posFile<<endl;
			if(isFile(posFile)){
			

			    unsigned long totalClusters=0;
			    fstream myclocsfile (posFile.c_str(),ios::in|ios::binary);
			    if (!myclocsfile) {
				cout<<"Unable to read position file "<<posFile<<" "<<strerror(errno)<<endl;
			    }
			    char toread;
			    myclocsfile.read(&toread,sizeof(char)); //version
			    unsigned int numberOfTotalBins ;
			    myclocsfile.read((char*)&numberOfTotalBins, sizeof (int)); //number of bins
			    unsigned int binIndex=0 ;
			    double xOffset = 0;
			    double yOffset = 0;


			    while(binIndex<numberOfTotalBins){ //for each bin

				unsigned char numberOfrecordsInbin ;
				myclocsfile.read((char*)&numberOfrecordsInbin, sizeof (char));

				totalClusters+=int(numberOfrecordsInbin);
				if(numberOfrecordsInbin == 0){
				    binIndex++;
				    continue;
				}

				xOffset = float(blockSize) * (binIndex % blocksPerlLine);
				yOffset = float(blockSize) * (binIndex / blocksPerlLine);

				for(int i=0;i<int(numberOfrecordsInbin);i++){ //for each block in the bin
				    unsigned char x ;
				    myclocsfile.read((char*)&x, sizeof (char));
				    unsigned char y ;
				    myclocsfile.read((char*)&y, sizeof (char));
				    tileCoords toadd;
			
				    toadd.x=do_round(xOffset + float(int(x))/10.0);
				    toadd.y=do_round(yOffset + float(int(y))/10.0);
				    coordinatesFound.push_back(toadd);

				}

				binIndex++;
	
			    }
			    //cout<<totalClusters<<endl;
			    myclocsfile.close();



			}else{
			    cerr<<"Cannot seem to find the position files, exiting"<<endl;
			    return 1;
			}
		    
		    }

		}
	    }	    

	    if(coordinatesFound.size() != numberOfClustersFirst){
		cerr<<"Found "<<numberOfClustersFirst <<" clusters in the bcl files but "<<coordinatesFound.size()<<" coordinates exiting"<<endl;
		return 1;
	    }
	    // return 1;
	    

		//bclDirectory+"/L00"+ stringify(lanesToUse[laneIndex]) +"/C"+stringify(cycle)+".1/s_"+stringify(lanesToUse[laneIndex])+"_"+stringify(tilesToUse[tileIndex])+".bcl"; 


	    //For sequences
	    vector<char *> forwardS;
	    vector<char *> reverseS;
	    vector<char *> index1S;
	    vector<char *> index2S;
	    //for qualities
	    vector<char *> forwardQ;
	    vector<char *> reverseQ;
	    vector<char *> index1Q;
	    vector<char *> index2Q;

	    
	    for(unsigned int i=0;i<numberOfClustersFirst;i++){

		char * forwardST;
		char * reverseST;
		char * index1ST ;
		char * index2ST ;
	
		char * forwardQT;
		char * reverseQT;
		char * index1QT ;
		char * index2QT ;

		try{
		    forwardST = new char[forwardCycles+1];
		    reverseST = new char[reverseCycles+1];
		    index1ST  = new char [index1Cycles+1];
		    index2ST  = new char[index2Cycles+1];
	
		    forwardQT = new char[forwardCycles+1];
		    reverseQT = new char[reverseCycles+1];
		    index1QT  = new char[index1Cycles+1];
		    index2QT  = new char[index2Cycles+1];
		}catch(bad_alloc& ba){
		    cerr << "Cannot allocate memory: " << ba.what() <<" try a machine with more RAM"<<endl;
		    return 1;
		}
		forwardS.push_back(forwardST);
		reverseS.push_back(reverseST);
		index1S.push_back( index1ST  );
		index2S.push_back( index2ST  );

		forwardQ.push_back(forwardQT);
		reverseQ.push_back(reverseQT);
		index1Q.push_back( index1QT  );
		index2Q.push_back( index2QT  );

	    }


	    
	    for(int cycle=0;cycle<numberOfCycles;cycle++){ //for each cycle
		vector<char *> * ptr2vectorS;
		vector<char *> * ptr2vectorQ;

		int cycleToUse;
		switch(cycleCode[cycle]){
		case 'f':
		    ptr2vectorS=&forwardS;
		    ptr2vectorQ=&forwardQ;
		    cycleToUse=cycle;
		    break;
		case 'i':
		    ptr2vectorS=&index1S;
		    ptr2vectorQ=&index1Q;
		    cycleToUse=cycle-forwardCycles;
		    break;
		case 'r': //reverse
		    ptr2vectorS=&reverseS;
		    ptr2vectorQ=&reverseQ;
		    if(fiir)
			cycleToUse=cycle-forwardCycles-index1Cycles-index2Cycles;
		    else
			cycleToUse=cycle-forwardCycles-index1Cycles;
		    break;
		case 'j': //i2
		    ptr2vectorS=&index2S;
		    ptr2vectorQ=&index2Q;
		    if(fiir)
			cycleToUse=cycle-forwardCycles-index1Cycles;
		    else
			cycleToUse=cycle-forwardCycles-index1Cycles-reverseCycles;
		    break;
		default:
		    cerr<<"Internal error, unable to determine cycle type"<<endl;
		    return 1;
		}

		for(unsigned int i=0;i<numberOfClustersFirst;i++){
		    char toread;
		    mybclfile[cycle].read(&toread, sizeof (char));
		    base returned=bin2base(toread);
		    if(returned.baseC  == 'A' && returned.qualC==0){
			(*ptr2vectorS)[i][cycleToUse]='N';
			(*ptr2vectorQ)[i][cycleToUse]=char(offseqQualScores);
		    }else{
			(*ptr2vectorS)[i][cycleToUse]=returned.baseC;
			(*ptr2vectorQ)[i][cycleToUse]=char(offseqQualScores+returned.qualC);
		    }
		}
      
	    }


	    //capping the char *
	    for(unsigned int i=0;i<numberOfClustersFirst;i++){
		forwardS[i][forwardCycles]='\0';
		reverseS[i][reverseCycles]='\0';
		index1S[i][index1Cycles]='\0';
		index2S[i][index2Cycles]='\0';

		forwardQ[i][forwardCycles]='\0';
		reverseQ[i][reverseCycles]='\0';
		index1Q[i][index1Cycles]='\0';
		index2Q[i][index2Cycles]='\0';

		if(fiir){
		    // cout<<"init"<<endl;
		    // cout<<string(index2S[i])<<endl;		    
		    // cout<<string(index2Q[i])<<endl;

		    reverseCstrRC( index2S[i], index2Cycles+1);
		    reverseCstr( index2Q[i], index2Cycles+1);

		    // cout<<string(index2S[i])<<endl;		    
		    // cout<<string(index2Q[i])<<endl;

		    
		}

	    }

	   
	    for(unsigned int i=0;i<numberOfClustersFirst;i++){

		BamAlignment  al;
		BamAlignment al2; //reverse read

		if(isPairedEnd){
		    string readName=nameOfExperiment+":"+stringify(lanesToUse[laneIndex])+":"+stringify(tilesToUse[tileIndex])+":"+stringify(coordinatesFound[i].x)+":"+stringify(coordinatesFound[i].y);

		    al.Name  =readName;
		    al2.Name =readName;
		    al.SetIsMapped (false);
		    al.SetIsMateMapped (false);
		    al2.SetIsMapped (false);
		    al2.SetIsMateMapped (false);


		    al.SetIsPaired(true);
		    al2.SetIsPaired(true);
		    al.SetIsFirstMate(true);
		    al2.SetIsSecondMate(true);	  		  

		    al.AlignmentFlag =  flagFirstPair;
		    al2.AlignmentFlag =  flagSecondPair;
		    
		    if(onlyIndex ){
			if(index1S[i] != indexToExtract)
			    continue;
		    }

		    if(index1Cycles != 0){
			if( !al.AddTag("XI","Z",string(index1S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if( !al.AddTag("YI","Z",string(index1Q[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(!al2.AddTag("XI","Z",string(index1S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(!al2.AddTag("YI","Z",string(index1Q[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		    }

		    if(index2Cycles != 0){
			if(! al.AddTag("XJ","Z",string(index2S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(! al.AddTag("YJ","Z",string(index2Q[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(!al2.AddTag("XJ","Z",string(index2S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(!al2.AddTag("YJ","Z",string(index2Q[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		    }


		    al.QueryBases  = forwardS[i];
		    al2.QueryBases = reverseS[i];
		    al.Qualities   = forwardQ[i];
		    al2.Qualities  = reverseQ[i];
		    writer.SaveAlignment(al);
		    writer.SaveAlignment(al2);

		}else{
		    //toadd
		    
		    string readName=nameOfExperiment+":"+stringify(lanesToUse[laneIndex])+stringify(tilesToUse[tileIndex])+":"+stringify(coordinatesFound[i].x)+":"+stringify(coordinatesFound[i].y);

		    al.Name  =readName;

		    al.SetIsMapped (false);
		    al.SetIsPaired(false);
		    al.AlignmentFlag =  flagSingleReads;

		    if(onlyIndex ){			
			if(index1S[i] != indexToExtract)
			    continue;
		    }

		    if(index1Cycles != 0){
			if(!al.AddTag("XI", "Z",string(index1S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(!al.AddTag("YI", "Z",string(index1Q[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		    }

		    if(index2Cycles != 0){
			if(!al.AddTag("XJ", "Z",string(index2S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
			if(!al.AddTag("YJ", "Z",string(index2Q[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }		    }


		    al.QueryBases  = forwardS[i];
		    al.Qualities   = forwardQ[i];

		    writer.SaveAlignment(al);


		}



		// cout<<<<endl;
		// cout<<forwardQ[i]<<endl;
	    }




	    //closing the BCL files
	    for(int cycle=1;cycle<=numberOfCycles;cycle++){ //for each cycle
		mybclfile[cycle-1].close();
	    }


	    //deallocating memory
	    for(unsigned int i=0;i<numberOfClustersFirst;i++){
		delete  forwardS[i];
		delete  reverseS[i];
		delete  index1S[i];
		delete  index2S[i];
		delete  forwardQ[i];
		delete  reverseQ[i];
		delete  index1Q[i];
		delete  index2Q[i];
	    }

	    


	} //end for each tile
    } //end for each lane

    writer.Close();
    
    //graciously finished flag
    if(!noFinishedFlag){
	ofstream fileoutFlag;
	fileoutFlag.open( string(bamfiletowrite+suffixFinishedFlag).c_str() );
	fileoutFlag<<""<<endl;
	fileoutFlag.close();
    }
    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

