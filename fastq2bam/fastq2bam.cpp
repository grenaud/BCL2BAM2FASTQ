/*
 * fastq2bam
 * Date: Feb-07-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>

#include "PutProgramInHeader.h"
#include "utils.h"
#include "FastQObj.h"
#include "FastQParser.h"

const uint32_t flagSingleReads =  4; // 00000100
const uint32_t flagFirstPair   = 77; // 01001101
const uint32_t flagSecondPair  =141; // 10001101

using namespace std;

int main (int argc, char *argv[]) {

    string bamoutfile;
    string index1    ="";
    string index2    ="";
    bool hasSpeciedIndex=false;

    string readgroup ="";
    bool isFasta=false;
    bool singleIndexFromDefline = false;
    bool doubleIndexFromDefline = false;
    bool indexFromDefline       = false;

    bool singleEndMode=false;

    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<"This program converts fasta/q files into unaligned bam\nUsage: "<<string(argv[0])<<" [options] [fastq in r1] (fastq in r2)"<<endl;
	cout<<"The second fastq represents the reverse reads, do not specify it for single-end reads"<<endl;
        cout<<"Options:"<<endl;
        cout<<"\tMandatory:"<<endl;
        cout<<"\t\t-o [output bam]"<<endl;
        cout<<"\tOptional:"<<endl;
        cout<<"\t\t-i1 [index1]"<<"\t\t"<<"Use this index sequence as the first  index field (XI)"<<endl;
        cout<<"\t\t-i2 [index2]"<<"\t\t"<<"Use this index sequence as the second index field (XJ)"<<endl;
        cout<<"\t\t-r  [read group]"<<endl;
        cout<<"\t\t-a  If input is fasta"<<endl;
        cout<<"\t\t-si Use the single index specified in the defline (as usually provided by Illumina)"<<endl;
        cout<<"\t\t-di As above but for double index"<<endl;

        return 1;
    }

    int lastIndex=1;
    for(int i=1;i<(argc);i++){
	if(strcmp(argv[i],"-o") == 0  ){
            bamoutfile=string(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-i1") == 0  ){
            index1=string(argv[i+1]);
            i++;
	    hasSpeciedIndex=true;
            continue;
        }

	if(strcmp(argv[i],"-i2") == 0  ){
            index2=string(argv[i+1]);
	    hasSpeciedIndex=true;
            i++;
            continue;
        }

	if(strcmp(argv[i],"-r") == 0  ){
            readgroup=string(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-a") == 0  ){
	    isFasta=true;
            continue;
        }

	if(strcmp(argv[i],"-si") == 0  ){
	    singleIndexFromDefline=true;
            continue;
        }

	if(strcmp(argv[i],"-di") == 0  ){
	    doubleIndexFromDefline=true;
            continue;
        }
	
	lastIndex=i;
	break;
    }
    
    if(lastIndex == (argc-1)){
	singleEndMode=true;
    }else{
	if(lastIndex == (argc-2)){
	    singleEndMode=false;
	}else{
	    cerr << "ERROR the option "<<argv[lastIndex] <<" is unknown" << endl;
	    return 1;	    
	}
    }
    
    if(bamoutfile.empty()){
	cerr << "ERROR -o option is mandatory " << endl;
        return 1;
    }

    if( (singleIndexFromDefline || doubleIndexFromDefline) &&
	hasSpeciedIndex ){
	cerr << "ERROR cannot specify options (-si or -di) and (-i1 or -i2)" << endl;
        return 1;
    }
    indexFromDefline = (singleIndexFromDefline || doubleIndexFromDefline);

    string fastqin1;
    string fastqin2;

    if(singleEndMode){
	fastqin1=argv[argc-1];
    }else{
	fastqin1=argv[argc-2];
	fastqin2=argv[argc-1];
    }

    BamWriter writer;

    SamHeader header ;
    string pID          = string(argv[0]);
    string pName        = string(argv[0]);
    string pCommandLine = "";
    
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));

    //adding @RG in read group
    if(!readgroup.empty()){	
	SamReadGroupDictionary  srgd;
	SamReadGroup srg ( readgroup );	
	srg.Sample        = readgroup;
	srgd.Add( srg );       
	header.ReadGroups=srgd;
    }

    RefVector references;

    if ( !writer.Open(bamoutfile,header,references) ) {
        cerr << "Could not open output BAM file "<<bamoutfile << endl;
        return 1;
    }

    FastQParser * fqp1;
    FastQParser * fqp2;

    if(singleEndMode){
	fqp1 = new FastQParser (fastqin1,isFasta);
    }else{
	fqp1 = new FastQParser (fastqin1,isFasta);
	fqp2 = new FastQParser (fastqin2,isFasta);
    }

    unsigned int totalSeqs=0;
    while(fqp1->hasData()){

	FastQObj * fo1=fqp1->getData();
	vector<string> def1=allTokens( *(fo1->getID()), ' '  );
	string def1s=def1[0];
	
	string ext1s;
	if(indexFromDefline){
	    if(def1.size() != 2){
		cerr << "ERROR: The following record does not have 2 fields for index assignment " <<  *(fo1->getID()) <<endl;
		return 1;
	    }
	    ext1s=def1[1];
	}
	// cout<<def1s<<endl;
	// cout<<ext1s<<endl;

	FastQObj * fo2;
	string def2s;
	string ext2s;

	if(!singleEndMode){
	    if(!fqp2->hasData()){
		cerr << "ERROR: Discrepency between fastq files at record " <<  *(fo1->getID()) <<endl;
		return 1;
	    }

	    fo2=fqp2->getData();
	    vector<string> def2=allTokens( *(fo2->getID()), ' ' );
	    def2s=def2[0];

	    if(indexFromDefline){
		if(def2.size() != 2){
		    cerr << "ERROR: The following record does not have 2 fields for index assignment " <<  *(fo2->getID()) <<endl;
		    return 1;
		}
		ext2s=def2[1];
	    }



	    if(strEndsWith(def1s,"/1")){
		def1s=def1s.substr(0,def1s.size()-2);
	    }
	    if(strEndsWith(def2s,"/2")){
		def2s=def2s.substr(0,def2s.size()-2);
	    }

	    if(strBeginsWith(def1s,"@")){
		def1s=def1s.substr(1,def1s.size()-1);
	    }
	    if(strBeginsWith(def2s,"@")){
		def2s=def2s.substr(1,def2s.size()-1);
	    }


	    if(def1s != def2s){
		cerr << "ERROR: Discrepency between fastq files, different names " << *(fo1->getID()) <<" and "<< *(fo2->getID()) <<endl;
		return 1;
	    }
	}
	// cout<<def1s<<endl;
	// cout<<ext1s<<endl;

	if(indexFromDefline){
	    vector<string> tokensSemiCol=allTokens( ext1s, ':' );
	    
	    if(singleIndexFromDefline){
		index1 = tokensSemiCol[ tokensSemiCol.size() -1 ]; //last element should be the only index

		for(unsigned int indexi=0;indexi<index1.size();indexi++){
		    if( !isValidDNA( index1[ indexi ] ) ){
			cerr << "ERROR: index found for " << index1 <<" in sequence "<<def1s<<" is not a valid DNA sequence" <<endl;
			return 1;
		    }
		}
	    }

	    if(doubleIndexFromDefline){
		index1 = tokensSemiCol[ tokensSemiCol.size() -2 ]; //second to last element should be the second index
		index2 = tokensSemiCol[ tokensSemiCol.size() -2 ]; //          last element should be the first index

		for(unsigned int indexi=0;indexi<index1.size();indexi++){
		    if( !isValidDNA( index1[ indexi ] ) ){
			cerr << "ERROR: index found for " << index1 <<" in sequence "<<def1s<<" is not a valid DNA sequence" <<endl;
			return 1;
		    }
		}

		for(unsigned int indexi=0;indexi<index2.size();indexi++){
		    if( !isValidDNA( index2[ indexi ] ) ){
			cerr << "ERROR: index found for " << index2 <<" in sequence "<<def1s<<" is not a valid DNA sequence" <<endl;
			return 1;
		    }
		}


	    }
	    
	}
	// cout<<def1s<<endl;
	// cout<<ext1s<<endl;

	BamAlignment toWrite1;
	BamAlignment toWrite2;

	toWrite1.Name=def1s;
	toWrite1.MapQuality=0;
	toWrite1.QueryBases =  *(fo1->getSeq());

	if(isFasta){
	    toWrite1.Qualities  =  string(toWrite1.QueryBases.length(),'!');	
	}else{
	    toWrite1.Qualities  =  *(fo1->getQual());
	}


	if(singleEndMode){
	    toWrite1.AlignmentFlag=flagSingleReads;	    
	}else{
	    toWrite1.AlignmentFlag=flagFirstPair;
	    toWrite2.Name=def2s;

	    toWrite2.AlignmentFlag=flagSecondPair;
	    toWrite2.MapQuality=0;
	    toWrite2.QueryBases =  *(fo2->getSeq());

	    if(isFasta)
		toWrite2.Qualities  =  string(toWrite2.QueryBases.length(),'!');
	    else
		toWrite2.Qualities  =  *(fo2->getQual());
	}
	
	

	//add tags for indices and fake qualities for the indices
	if(!index1.empty()){
	    if(!toWrite1.AddTag("XI", "Z",index1) )                      {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite1.AddTag("YI", "Z",string(index1.length(),'!')) ) {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!singleEndMode){
		if(!toWrite2.AddTag("XI", "Z",index1) )                      {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		if(!toWrite2.AddTag("YI", "Z",string(index1.length(),'!')) ) {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    }
	}

	if(!index2.empty()){
	    if(!toWrite1.AddTag("XJ", "Z",index2) )                      {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite1.AddTag("YJ", "Z",string(index2.length(),'!')) ) {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!singleEndMode){
		if(!toWrite2.AddTag("XJ", "Z",index2) )                      {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		if(!toWrite2.AddTag("YJ", "Z",string(index2.length(),'!')) ) {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    }
	}

	if(!readgroup.empty()){
	    if(!toWrite1.AddTag("RG", "Z",readgroup) )                      { cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!singleEndMode){
		if(!toWrite2.AddTag("RG", "Z",readgroup) )                      { cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    }
	}

	

	writer.SaveAlignment(toWrite1);
	if(!singleEndMode){
	    writer.SaveAlignment(toWrite2);
	}

	totalSeqs++;
    }
    
    delete fqp1;
    if(!singleEndMode){
	delete fqp2;
    }
    
    writer.Close();

    cerr<<"Wrote "<<totalSeqs<<" sequences successfully"<<endl;


    return 0;
}

