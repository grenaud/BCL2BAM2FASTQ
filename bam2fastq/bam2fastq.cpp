/*
 * bam2fastq
 * Date: May-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>


#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>

#include "utils.h"

using namespace std;
using namespace BamTools;

typedef struct { 

    ogzstream single;
    ogzstream pairr1;
    ogzstream pairr2;
    //index1
    ogzstream singlei1;
    ogzstream pairi1;
    //index2
    ogzstream singlei2;
    ogzstream pairi2;

    //
    // FAILED
    //
    ogzstream singlef;
    ogzstream pairr1f;
    ogzstream pairr2f;
    
    //index 1
    ogzstream singlei1f;
    ogzstream pairi1f;

    //index 2
    ogzstream singlei2f;
    ogzstream pairi2f;

 } fqwriters;


int main (int argc, char *argv[]) {
    bool splitRG=false;

    string usage=string(""+string(argv[0])+" <options> [bam file] [output prefix]"+		     
			"\n\nThis program will read a bam file sorted by name and produce 6 files:\n"+
			"\n\t[output prefix].fq.gz\t\tFor single reads"+
			"\n\t[output prefix]_r1.fq.gz\tFor first mates"+
			"\n\t[output prefix]_r2.fq.gz\tFor second mates\n"
			"\n\t[output prefix].fq.fail.gz\t\tFor single reads that failed"+
			"\n\t[output prefix]_r1.fq.fail.gz\tFor first mates for clusters that failed"+
			"\n\t[output prefix]_r2.fq.fail.gz\tFor second mates for clusters that failed\n"

			"Options:\n"+
			"\t\t"+"--rg"+"\t"+"Split according to read group (Default: "+booleanAsString(splitRG)+" )\n"
			);		      


    if(argc  < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    for(int i=1;i<(argc-2);i++){ 

	if(string(argv[i]) == "--rg" ) {
	    splitRG=true;
	    continue;
	}

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    string bamfile = string(argv[argc-2]);
    string outdir  = string(argv[argc-1]);



    fqwriters onereadgroup;
    map<string,fqwriters *> rg2fqwriters;


    BamReader reader;
    if ( !reader.Open(bamfile) ) {
    	cerr << "Could not open input BAM file  "<<bamfile << endl;
    	return 1;
    }
    unsigned int totalReads = 0;
    unsigned int singlereads = 0;
    unsigned int pairedreads = 0;

    BamAlignment al;
    BamAlignment al2;
    bool al2Null=true;
    bool firstRead=true;
    bool indexUnknown=true;
    bool i1present=false;
    bool i2present=false;

    while ( reader.GetNextAlignment(al) ) {
	totalReads++;
	
	string rgTag;

	if(indexUnknown){
	    i1present=al.HasTag("XI");
	    i2present=al.HasTag("XJ");
	    cerr<<i1present<<i2present<<endl;
	    indexUnknown=false;
	}

	if(splitRG){
 
	    if(al.HasTag("RG")){
		al.GetTag("RG",rgTag);
	    }else{
		rgTag="unknown";
	    }


	    if(rg2fqwriters.find(rgTag) == rg2fqwriters.end()){ //new
		cerr<<"Found new RG "<<rgTag<<endl;

		rg2fqwriters[rgTag] = new  fqwriters();
		string outdirs   = outdir+"rg_"+rgTag+""+".fq.gz";
		string outdir1   = outdir+"rg_"+rgTag+""+"_r1.fq.gz";
		string outdir2   = outdir+"rg_"+rgTag+""+"_r2.fq.gz";
		string outdirsf  = outdir+"rg_"+rgTag+""+".fq.fail.gz";
		string outdir1f  = outdir+"rg_"+rgTag+""+"_r1.fail.fq.gz";
		string outdir2f  = outdir+"rg_"+rgTag+""+"_r2.fail.fq.gz";

		rg2fqwriters[rgTag]->single.open(outdirs.c_str(), ios::out);
		rg2fqwriters[rgTag]->pairr1.open(outdir1.c_str(), ios::out);
		rg2fqwriters[rgTag]->pairr2.open(outdir2.c_str(), ios::out);

		rg2fqwriters[rgTag]->singlef.open(outdirsf.c_str(), ios::out);
		rg2fqwriters[rgTag]->pairr1f.open(outdir1f.c_str(), ios::out);
		rg2fqwriters[rgTag]->pairr2f.open(outdir2f.c_str(), ios::out);

		
		if(i1present){

		    string outdirsi1    = outdir+"rg_"+rgTag+".i1.fq.gz";
		    string outdiri1     = outdir+"rg_"+rgTag+"_i1.fq.gz";

		    rg2fqwriters[rgTag]->singlei1.open(outdirsi1.c_str(), ios::out);
		    rg2fqwriters[rgTag]->pairi1.open(outdiri1.c_str(), ios::out);

		    string outdirsi1f    = outdir+"rg_"+rgTag+".i1.fail.fq.gz";
		    string outdiri1f     = outdir+"rg_"+rgTag+"_i1.fail.fq.gz";

		    rg2fqwriters[rgTag]->singlei1f.open(outdirsi1f.c_str(), ios::out);
		    rg2fqwriters[rgTag]->pairi1f.open(outdiri1f.c_str(), ios::out);
		}

		if(i2present){

		    string outdirsi2    = outdir+"rg_"+rgTag+".i2.fq.gz";
		    string outdiri2     = outdir+"rg_"+rgTag+"_i2.fq.gz";

		    rg2fqwriters[rgTag]->singlei2.open(outdirsi2.c_str(), ios::out);
		    rg2fqwriters[rgTag]->pairi2.open(outdiri2.c_str(), ios::out);

		    string outdirsi2f    = outdir+"rg_"+rgTag+".i2.fail.fq.gz";
		    string outdiri2f     = outdir+"rg_"+rgTag+"_i2.fail.fq.gz";

		    rg2fqwriters[rgTag]->singlei2f.open(outdirsi2f.c_str(), ios::out);
		    rg2fqwriters[rgTag]->pairi2f.open(outdiri2f.c_str(), ios::out);
		}


	    }

	}else{

	    if(firstRead){
		string outdirs   = outdir+".fq.gz";
		string outdir1   = outdir+"_r1.fq.gz";
		string outdir2   = outdir+"_r2.fq.gz";
		string outdirsf  = outdir+".fq.fail.gz";
		string outdir1f  = outdir+"_r1.fail.fq.gz";
		string outdir2f  = outdir+"_r2.fail.fq.gz";


		onereadgroup.single.open(outdirs.c_str(), ios::out);
		onereadgroup.pairr1.open(outdir1.c_str(), ios::out);
		onereadgroup.pairr2.open(outdir2.c_str(), ios::out);

		onereadgroup.singlef.open(outdirsf.c_str(), ios::out);
		onereadgroup.pairr1f.open(outdir1f.c_str(), ios::out);
		onereadgroup.pairr2f.open(outdir2f.c_str(), ios::out);

		if(i1present){

		    string outdirsi1   = outdir+".i1.fq.gz";
		    string outdiri1     = outdir+"_i1.fq.gz";

		    onereadgroup.singlei1.open(outdirsi1.c_str(), ios::out);
		    onereadgroup.pairi1.open(outdiri1.c_str(), ios::out);

		    string outdirsi1f   = outdir+".i1.fail.fq.gz";
		    string outdiri1f     = outdir+"_i1.fail.fq.gz";

		    onereadgroup.singlei1f.open(outdirsi1f.c_str(), ios::out);
		    onereadgroup.pairi1f.open(outdiri1f.c_str(), ios::out);
		}

		if(i2present){

		    string outdirsi2   = outdir+".i2.fq.gz";
		    string outdiri2     = outdir+"_i2.fq.gz";

		    onereadgroup.singlei2.open(outdirsi2.c_str(), ios::out);
		    onereadgroup.pairi2.open(outdiri2.c_str(), ios::out);

		    string outdirsi2f   = outdir+".i2.fail.fq.gz";
		    string outdiri2f     = outdir+"_i2.fail.fq.gz";

		    onereadgroup.singlei2f.open(outdirsi2f.c_str(), ios::out);
		    onereadgroup.pairi2f.open(outdiri2f.c_str(), ios::out);
		}

		firstRead=false;
	    }

	}


	if(al.IsPaired() ){
	    if( al2Null ){
		al2=al;
		al2Null=false;
		continue;
	    }else{
		if(al.Name != al2.Name ){
		    cerr << "Seq#1 " <<al.Name <<" has a different id than seq #2 " <<al2.Name <<" exiting " << endl;
		    return 1;
		}


		if(splitRG){

		    
		    if(!al.IsFailedQC() && !al2.IsFailedQC() ){//read passed
			rg2fqwriters[rgTag]->pairr2 <<"@"<< al.Name<<"/2" <<endl <<al.QueryBases<<endl<<"+"<<endl <<al.Qualities<<endl;
			rg2fqwriters[rgTag]->pairr1 <<"@"<<al2.Name<<"/1"<<endl <<al2.QueryBases<<endl<<"+"<<endl<<al2.Qualities<<endl;

			
			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->pairi1 <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;

			}

			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->pairi2 <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;

			}

		    }else{
			rg2fqwriters[rgTag]->pairr2f<<"@"<< al.Name<<"/2" <<endl <<al.QueryBases<<endl<<"+"<<endl <<al.Qualities<<endl;
			rg2fqwriters[rgTag]->pairr1f<<"@"<<al2.Name<<"/1"<<endl <<al2.QueryBases<<endl<<"+"<<endl<<al2.Qualities<<endl;


			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->pairi1f <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;
			}


			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->pairi2f <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;
			}
		    }


		}else{
		
		    if(!al.IsFailedQC() && !al2.IsFailedQC() ){//read passed
			onereadgroup.pairr2 <<"@"<< al.Name<<"/2" <<endl <<al.QueryBases<<endl<<"+"<<endl <<al.Qualities<<endl;
			onereadgroup.pairr1 <<"@"<<al2.Name<<"/1"<<endl <<al2.QueryBases<<endl<<"+"<<endl<<al2.Qualities<<endl;

			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.pairi1 <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;

			}

			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.pairi2 <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;

			}



		    }else{
			onereadgroup.pairr2f<<"@"<< al.Name<<"/2" <<endl <<al.QueryBases<<endl<<"+"<<endl <<al.Qualities<<endl;
			onereadgroup.pairr1f<<"@"<<al2.Name<<"/1"<<endl <<al2.QueryBases<<endl<<"+"<<endl<<al2.Qualities<<endl;

			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.pairi1f <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;
			}


			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.pairi2f <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;
			}

		    }
		}
		pairedreads++;
	    }

	}else{
	    if( al2Null ){





		if(splitRG){

			
		    if(!al.IsFailedQC() ){//read passed
			rg2fqwriters[rgTag]->single <<"@"<<al.Name<<endl<<al.QueryBases<<endl<<"+"<<endl<<al.Qualities<<endl;
			    			    
			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->singlei1 <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;
			}

			    
			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->singlei2 <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;
			}

		    }else{//read failed
			rg2fqwriters[rgTag]->singlef<<"@"<<al.Name<<endl<<al.QueryBases<<endl<<"+"<<endl<<al.Qualities<<endl;

			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->singlei1f <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;
			}

			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    rg2fqwriters[rgTag]->singlei2f <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;
			}

		    }

		}else{ //split rg

		    if(!al.IsFailedQC() ){//read passed
			onereadgroup.single <<"@"<<al.Name<<endl<<al.QueryBases<<endl<<"+"<<endl<<al.Qualities<<endl;
			    
			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.singlei1 <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;
			}

			    
			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.singlei2 <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;
			}

			    

		    }else{ //read failed
			onereadgroup.singlef<<"@"<<al.Name<<endl<<al.QueryBases<<endl<<"+"<<endl<<al.Qualities<<endl;
			    
			if(i1present){
			    string i1stw;
			    string i1qtw;

			    if(!al.GetTag("XI",i1stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YI",i1qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.singlei1f <<"@"<< al.Name<<"" <<endl <<i1stw<<endl<<"+"<<endl <<i1qtw<<endl;
			}

			if(i2present){
			    string i2stw;
			    string i2qtw;

			    if(!al.GetTag("XJ",i2stw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    if(!al.GetTag("YJ",i2qtw)){
				cerr<<"Cannot get index tag from read "<<al.Name<<endl;
				return 1;
			    }

			    onereadgroup.singlei2f <<"@"<< al.Name<<"" <<endl <<i2stw<<endl<<"+"<<endl <<i2qtw<<endl;
			}

		    }

		}//not split rg

		singlereads++;
	    }else{ //if al2 is not null
		cerr << "Seq#1 " <<al.Name <<" was found next to a lone paired reads " <<al2.Name <<" exiting " << endl;
		return 1;
	    }

	}
	al2Null=true;
    }

    if(!splitRG){
	
	onereadgroup.single.close();
	onereadgroup.pairr1.close();
	onereadgroup.pairr2.close();
	onereadgroup.singlef.close();
	onereadgroup.pairr1f.close();
	onereadgroup.pairr2f.close();

	if(i1present){
	    onereadgroup.singlei1.close();
	    onereadgroup.pairi1.close();
   
	    onereadgroup.singlei1f.close();
	    onereadgroup.pairi1f.close();
	}


	if(i2present){
	    onereadgroup.singlei2.close();
	    onereadgroup.pairi2.close();
   
	    onereadgroup.singlei2f.close();
	    onereadgroup.pairi2f.close();
	}

    }else{

	vector<string> vsKey = allKeysMap(rg2fqwriters);
	
	for(unsigned int i=0;i<vsKey.size();i++){

	    rg2fqwriters[ vsKey[i] ]->single.close();
	    rg2fqwriters[ vsKey[i] ]->pairr1.close();
	    rg2fqwriters[ vsKey[i] ]->pairr2.close();
	    rg2fqwriters[ vsKey[i] ]->singlef.close();
	    rg2fqwriters[ vsKey[i] ]->pairr1f.close();
	    rg2fqwriters[ vsKey[i] ]->pairr2f.close();

	    if(i1present){
		rg2fqwriters[ vsKey[i] ]->singlei1.close();
		rg2fqwriters[ vsKey[i] ]->pairi1.close();
   
		rg2fqwriters[ vsKey[i] ]->singlei1f.close();
		rg2fqwriters[ vsKey[i] ]->pairi1f.close();
	    }


	    if(i2present){
		rg2fqwriters[ vsKey[i] ]->singlei2.close();
		rg2fqwriters[ vsKey[i] ]->pairi2.close();
   
		rg2fqwriters[ vsKey[i] ]->singlei2f.close();
		rg2fqwriters[ vsKey[i] ]->pairi2f.close();
	    }


	}
	
    }
	
    cerr<<"Program finished successfully, wrote "<<totalReads<<" reads, "<<singlereads<<" single reads and "<<pairedreads<<" pairs"<<endl;
    return 0;
}

