/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>

#include "uvtDetectDistCorr.h"
#include<glog/logging.h>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    time_t st ,et;
    st = time(NULL);//Computes execution start time(in seconds)
    int status=0;//Flag to store return status of functions
     google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Settings path to create log files. It reads the environment variable 'GLOG_log_dir'
      checkParFile(argv[0]);//Checks for existence of parameter file and 'PFILES' environment variable
    uvtDetectDistCorr obj;//Creating object for DetectDIstCorr class
    status=obj.read(argc,argv);//Read parameters from uvtDetectDistCorr.par file and return status
      if(status){
        cerr<<"***Error in Reading the Parameter ***"<<endl;
        return(EXIT_FAILURE);
    }
    obj.display();//Display detectDistCorr input parameters
    status=obj.uvtDetectDistCorrProcess ();//Performing Detect Dist Corr computation
    if(status){
        cerr<<"***Error in uvtDetectDistortion Module***"<<endl;
        return(EXIT_FAILURE);
    }
    et = time(NULL);//Computes execution end time(in seconds)
    cerr<<"Execution Time "<<et-st<<" Seconds"<<endl;
    google::ShutdownGoogleLogging ();//Close logging library
    
    return 0;
}

