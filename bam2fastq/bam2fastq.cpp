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

int main (int argc, char *argv[]) {
    // bool splitQC=false;

    string usage=string(""+string(argv[0])+" <options> [bam file] [output prefix]"+		     
			"\n\nThis program will read a bam file sorted by name and produce 6 files:\n"+
			"\n\t[output prefix].fq.gz\t\tFor single reads"+
			"\n\t[output prefix]_r1.fq.gz\tFor first mates"+
			"\n\t[output prefix]_r2.fq.gz\tFor second mates\n"
			"\n\t[output prefix].fq.fail.gz\t\tFor single reads that failed"+
			"\n\t[output prefix]_r1.fq.fail.gz\tFor first mates for clusters that failed"+
			"\n\t[output prefix]_r2.fq.fail.gz\tFor second mates for clusters that failed\n"

			// "Options:\n"+
			// "\t\t"+"--qc"+"\t"+"Split quality failed and non-quality failed (Default: "+booleanAsString(splitQC)+" )\n"
			);		      


    if(argc != 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    for(int i=1;i<(argc-2);i++){ 

	// if(string(argv[i]) == "--qc" ) {
	//     splitQC=true;
	//     continue;
	// }

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    string bamfile = string(argv[argc-2]);
    string outdir  = string(argv[argc-1]);

    string outdirs   = outdir+".fq.gz";
    string outdir1   = outdir+"_r1.fq.gz";
    string outdir2   = outdir+"_r2.fq.gz";
    string outdirsf  = outdir+".fq.fail.gz";
    string outdir1f  = outdir+"_r1.fail.fq.gz";
    string outdir2f  = outdir+"_r2.fail.fq.gz";


    ogzstream single;
    single.open(outdirs.c_str(), ios::out);
    ogzstream pairr1;
    pairr1.open(outdir1.c_str(), ios::out);
    ogzstream pairr2;
    pairr2.open(outdir2.c_str(), ios::out);

    ogzstream singlef;
    singlef.open(outdirsf.c_str(), ios::out);
    ogzstream pairr1f;
    pairr1f.open(outdir1f.c_str(), ios::out);
    ogzstream pairr2f;
    pairr2f.open(outdir2f.c_str(), ios::out);

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
    while ( reader.GetNextAlignment(al) ) {
	totalReads++;

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

		
		if(!al.IsFailedQC() && !al2.IsFailedQC() ){//read passed
		    pairr2 <<"@"<< al.Name<<"/2" <<endl <<al.QueryBases<<endl<<"+"<<endl <<al.Qualities<<endl;
		    pairr1 <<"@"<<al2.Name<<"/1"<<endl <<al2.QueryBases<<endl<<"+"<<endl<<al2.Qualities<<endl;
		}else{
		    pairr2f<<"@"<< al.Name<<"/2" <<endl <<al.QueryBases<<endl<<"+"<<endl <<al.Qualities<<endl;
		    pairr1f<<"@"<<al2.Name<<"/1"<<endl <<al2.QueryBases<<endl<<"+"<<endl<<al2.Qualities<<endl;
		}
		pairedreads++;
	    }

	}else{
	    if( al2Null ){


		if(!al.IsFailedQC() ){//read passed
		    single <<"@"<<al.Name<<endl<<al.QueryBases<<endl<<"+"<<endl<<al.Qualities<<endl;
		}else{
		    singlef<<"@"<<al.Name<<endl<<al.QueryBases<<endl<<"+"<<endl<<al.Qualities<<endl;
		}

		singlereads++;
	    }else{
		cerr << "Seq#1 " <<al.Name <<" was found next to a lone paired reads " <<al2.Name <<" exiting " << endl;
		return 1;
	    }
	}
	al2Null=true;
    }

    single.close();
    pairr1.close();
    pairr2.close();
    singlef.close();
    pairr1f.close();
    pairr2f.close();
	
    cerr<<"Program finished successfully, wrote "<<totalReads<<" reads, "<<singlereads<<" single reads and "<<pairedreads<<" pairs"<<endl;
    return 0;
}

