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
    string index1     ="";
    string index2     ="";
    string index1q    ="";
    string index2q    ="";

    bool hasSpeciedIndex=false;
    bool hasSpeciedIndexFile =false;
    bool hasSpeciedIndexFile2=false;

    string index1filename;
    string index2filename;

    string readgroup ="";
    bool isFasta=false;
    bool singleIndexFromDefline = false;
    bool doubleIndexFromDefline = false;
    bool indexFromDefline       = false;

    bool singleEndMode=false;
    int qualForFasta=0;
    int qualForIndices=0;

    bool qualScoresCapBool=false;
    bool baseQual64       =false;
    
    int qualScoresCap=60;

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
        cout<<"\t\t-a  If input is fasta"<<endl;
	cout<<"\t\t-q  [qual]"<<"\t\t"<<"If input is fasta, use this qual as the quality (Default : "+stringify(qualForFasta)+")"<<endl;
	cout<<"\t\t-qi  [qual]"<<"\t\t"<<"Value for indices quality scores if the qual scores are missing (Default : "+stringify(qualForIndices)+")"<<endl;
	
	cout<<"\t\t-m  [qual]"<<"\t\t"<<"Cap quality scores at this value (Default : "+stringify(qualScoresCapBool)+")"<<endl;

	cout<<"\t\t-b64  [qual]"<<"\t\t"<<"Quality scores are on the 64 offset instead of 33 (Default : "+boolStringify(baseQual64)+")"<<endl;


        cout<<"\n\t\tIndex options:"<<endl;
	
        cout<<"\t\t\t-if1 [file index1]"<<"\t\t"<<"Use this fastq file as first index sequence  (XI)"<<endl;
        cout<<"\t\t\t-if2 [file index2]"<<"\t\t"<<"Use this fastq file as second index sequence (XJ)"<<endl;

        cout<<"\t\t\t-i1 [index1 string]"<<"\t\t"<<"Use this index sequence as the first  index field (XI)"<<endl;
        cout<<"\t\t\t-i2 [index2 string]"<<"\t\t"<<"Use this index sequence as the second index field (XJ)"<<endl;


        cout<<"\t\t\t-si Use the single index specified in the defline (as usually provided by Illumina)"<<endl;
        cout<<"\t\t\t-di As above but for double index"<<endl;

        return 1;
    }

    int lastIndex=1;
    bool specifiedQual=false;
    bool specifiedQualIndices=false;

    for(int i=1;i<(argc);i++){
	if(strcmp(argv[i],"-o") == 0  ){
            bamoutfile=string(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-m") == 0  ){
            qualScoresCapBool=true;
	    qualScoresCap=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-b64") == 0  ){
            baseQual64=true;
            continue;
        }

	
	if(strcmp(argv[i],"-q") == 0  ){
            qualForFasta=destringify<int>(argv[i+1]);
            i++;
	    specifiedQual=true;
            continue;
        }

	if(strcmp(argv[i],"-qi") == 0  ){
            qualForIndices=destringify<int>(argv[i+1]);
            i++;
	    specifiedQualIndices=true;
            continue;
        }

	if(strcmp(argv[i],"-i1") == 0  ){
            index1  = string(argv[i+1]);
	    index1q = string(index1.length(),char(qualForFasta+33));

            i++;
	    hasSpeciedIndex=true;
            continue;
        }

	if(strcmp(argv[i],"-i2") == 0  ){
            index2  = string(argv[i+1]);
	    index2q = string(index2.length(),char(qualForFasta+33));

	    hasSpeciedIndex=true;
            i++;
            continue;
        }

	if(strcmp(argv[i],"-if1") == 0  ){
            index1filename=string(argv[i+1]);
            i++;
	    hasSpeciedIndexFile=true;
            continue;
        }

	if(strcmp(argv[i],"-if2") == 0  ){
	    index2filename=string(argv[i+1]);
	    hasSpeciedIndexFile2=true;
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

    if( (specifiedQual && !isFasta) ){
	cerr << "Cannot specify quality if you do not use fasta" << endl;
        return 1;
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

    if( (hasSpeciedIndexFile || hasSpeciedIndexFile2) &&
	indexFromDefline){
	cerr << "ERROR cannot specify index from the defline and specify index files at once" << endl;
        return 1;
    }

    if( (!hasSpeciedIndexFile && hasSpeciedIndexFile2) ){
	cerr << "ERROR cannot specify second index and not first" << endl;
        return 1;
    }

    if( (hasSpeciedIndex && hasSpeciedIndexFile) ){
	cerr << "ERROR cannot specify both index file and index from the command line" << endl;
        return 1;
    }
    

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
    FastQParser * fqpindx1;
    FastQParser * fqpindx2;

    if(singleEndMode){
	fqp1 = new FastQParser (fastqin1,isFasta);
    }else{
	fqp1 = new FastQParser (fastqin1,isFasta);
	fqp2 = new FastQParser (fastqin2,isFasta);
    }

    if(hasSpeciedIndexFile){
	fqpindx1      = new FastQParser (index1filename,isFasta);
	if(hasSpeciedIndexFile2){
	    fqpindx2  = new FastQParser (index2filename,isFasta);
	}
    }
    


    unsigned int totalSeqs=0;
    while(fqp1->hasData()){

	FastQObj * fo1=fqp1->getData();
	vector<string> def1=allTokens( *(fo1->getID()), ' '  );
	string def1s=def1[0];

	string ext1s;
	if(indexFromDefline){
	    if(def1.size() != 2){
		vector<string> def1WithPound=allTokens( *(fo1->getID()), '#'  );
		if(def1WithPound.size() != 2){
		    cerr << "ERROR: The following record does not have 2 fields for index assignment " <<  *(fo1->getID()) <<endl;
		    return 1;
		}

		def1s=def1WithPound[0];
		ext1s=def1WithPound[1];
		
		if(strEndsWith(ext1s,"/1")){
		    ext1s=ext1s.substr(0,ext1s.size()-2);
		}

		//cout<<ext1s<<endl;
		// return 1;

		
	    }else{
		ext1s=def1[1];
	    }
	}
	// cout<<def1s<<endl;
	// cout<<ext1s<<endl;

	FastQObj * fo2;
	FastQObj * foi1;
	FastQObj * foi2;


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
		    vector<string> def2WithPound=allTokens( *(fo2->getID()), '#'  );
		    if(def2WithPound.size() != 2){
			cerr << "ERROR: The following record does not have 2 fields for index assignment " <<  *(fo2->getID()) <<endl;
			return 1;
		    }

		    def2s=def2WithPound[0];
		    ext2s=def2WithPound[1];
		    // cerr << "ERROR: The following record does not have 2 fields for index assignment " <<  *(fo2->getID()) <<endl;
		    // return 1;
		}else{
		    ext2s=def2[1];
		}
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
		cerr << " we tried names " << def1s <<" and "<< def2s <<endl;
		return 1;
	    }
	}


	if( hasSpeciedIndexFile ){
	    if(!fqpindx1->hasData()){
		cerr << "ERROR: Discrepency between fastq files at record " <<  *(fo1->getID()) <<endl;
		return 1;
	    }

	    foi1=fqpindx1->getData();
		
	    if( hasSpeciedIndexFile2 ){
		if(!fqpindx2->hasData()){
		    cerr << "ERROR: Discrepency between fastq files at record " <<  *(fo1->getID()) <<endl;
		    return 1;
		}
		    
		foi2=fqpindx2->getData();
	    }

	}
	// cout<<def1s<<endl;
	// cout<<ext1s<<endl;

	if(indexFromDefline){
	    vector<string> tokensSemiCol=allTokens( ext1s, ':' );
	    
	    if(singleIndexFromDefline){
		index1  = tokensSemiCol[ tokensSemiCol.size() -1 ]; //last element should be the only index
		index1q = string(index1.length(),char(qualForIndices+33));

		for(unsigned int indexi=0;indexi<index1.size();indexi++){
		    if( !isValidDNA( index1[ indexi ] ) ){
			cerr << "ERROR: index found for " << index1 <<" in sequence "<<def1s<<" is not a valid DNA sequence" <<endl;
			return 1;
		    }
		}
	    }

	    if(doubleIndexFromDefline){
		index1  = tokensSemiCol[ tokensSemiCol.size() -2 ];             // second to last element should be the second index
		index1q = string(index1.length(),char(qualForIndices+33));

		index2  = tokensSemiCol[ tokensSemiCol.size() -2 ];             // last element should be the first index
		index2q = string(index2.length(),char(qualForIndices+33));
		
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
	    
	}//end index from defline


	if( hasSpeciedIndexFile ){	    

	    index1      = *(foi1->getSeq());
	    index1q     = *(foi1->getQual());

	    if( hasSpeciedIndexFile2 ){		
		index2  = *(foi2->getSeq());
		index2q = *(foi2->getQual());		
	    }
	}else{
	    
	}

	// cout<<def1s<<endl;
	// cout<<ext1s<<endl;

	BamAlignment toWrite1;
	BamAlignment toWrite2;

	toWrite1.Name=def1s;
	toWrite1.MapQuality=0;
	toWrite1.QueryBases =  *(fo1->getSeq());

	if(isFasta){
	    toWrite1.Qualities  =  string(toWrite1.QueryBases.length(),char(qualForFasta+33));	
	}else{
	    toWrite1.Qualities  =  *(fo1->getQual());
	    if(baseQual64){
		for(unsigned int indexS=0;indexS<toWrite1.Qualities.size();indexS++){
		    // cerr<<toWrite1.Qualities[indexS]<<endl;
		    // cerr<<int(toWrite1.Qualities[indexS])<<endl;
		    // cerr<< ( (int(toWrite1.Qualities[indexS])-31) ) <<endl;
		    // cerr<< char( (int(toWrite1.Qualities[indexS])-31) ) <<endl;
		    
		    toWrite1.Qualities[indexS] = char ( (int(toWrite1.Qualities[indexS])-31) ) ; //-64+33 = 31
		    //return 1;
		}
	    }
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
		toWrite2.Qualities  =  string(toWrite2.QueryBases.length(),char(qualForFasta+33));
	    else{
		toWrite2.Qualities  =  *(fo2->getQual());
		if(baseQual64){
		    for(unsigned int indexS=0;indexS<toWrite2.Qualities.size();indexS++){
			toWrite2.Qualities[indexS] = char ( (int(toWrite2.Qualities[indexS])-31) ) ; //-64+33 = 31
		    }
		}
	    }
	}

	
	if(qualScoresCapBool){

	    for(unsigned int indexS=0;indexS<toWrite1.Qualities.size();indexS++){
		toWrite1.Qualities[indexS] = char ( min(qualScoresCap, int(toWrite1.Qualities[indexS])-33) +33) ;
	    }

	    for(unsigned int indexS=0;indexS<toWrite2.Qualities.size();indexS++){
		toWrite2.Qualities[indexS] = char ( min(qualScoresCap, int(toWrite2.Qualities[indexS])-33) +33) ;
	    }

	    for(unsigned int indexS=0;indexS<index1q.size();indexS++){
		index1q[indexS]            = char ( min(qualScoresCap, int(index1q[indexS])-33) +33) ;
	    }	

	    for(unsigned int indexS=0;indexS<index2q.size();indexS++){
		index2q[indexS]            = char ( min(qualScoresCap, int(index2q[indexS])-33) +33) ;
	    }	

	}

	//add tags for indices and fake qualities for the indices
	if(!index1.empty()){
	    if(!toWrite1.AddTag("XI", "Z",index1) ) {        cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite1.AddTag("YI", "Z",index1q)) {        cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!singleEndMode){
		if(!toWrite2.AddTag("XI", "Z",index1)  ) {   cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		if(!toWrite2.AddTag("YI", "Z",index1q) ) {   cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    }
	}

	if(!index2.empty()){
	    if(!toWrite1.AddTag("XJ", "Z",index2)  ) {       cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite1.AddTag("YJ", "Z",index2q) ) {       cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!singleEndMode){
		if(!toWrite2.AddTag("XJ", "Z",index2) ) {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
		if(!toWrite2.AddTag("YJ", "Z",index2q)) {    cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    }
	}

	if(!readgroup.empty()){
	    if(!toWrite1.AddTag("RG", "Z",readgroup) )     {  cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!singleEndMode){
		if(!toWrite2.AddTag("RG", "Z",readgroup) ) {  cerr<<"Internal error, cannot add tag"<<endl; return 1; }
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

    if(hasSpeciedIndexFile){
	delete fqpindx1;
	if(hasSpeciedIndexFile2){
	    delete fqpindx2;
	}
    }


    
    writer.Close();

    cerr<<"Wrote "<<totalSeqs<<" sequences successfully"<<endl;


    return 0;
}

