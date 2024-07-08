/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>

#include<iostream>
#include<ctime>

//#include "uvtCentroidBias.h"
#include "uvtCentroidCorr.h"
#include<glog/logging.h>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    time_t et,st;
    st=time(NULL);//Computes execution start time(in seconds)
    
    int status=0;//FLag to store return status of functions
      google::InitGoogleLogging(argv[0]);//Initializing of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files.It reads the environment variable 'GLOG_log_dir'
   // uvtFilterBadpix obj;
    //status=obj.read(argc,argv);
    //
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtCentroidCorr obj;//Creating object for centroidCorr class
    status=obj.read(argc,argv);//Read parameters from uvtCentroidCorr.par file and return status
     if(status){
        cerr<<endl<<"***Error in reading parameters for Centroid Correction module*** ";
        return (EXIT_FAILURE);
     }
    obj.display();//Display centroidCorr input parameters
     status=obj.uvtCentroidCorrProcess();//Performing Centroid COrr computation
     if(status){
         cerr<<endl<<"Error in process for uvtCentroidCorr"<<endl;
         return(EXIT_FAILURE);
     }
   
    cerr<<endl<<"::::::::::::::::::::::::::Centroid Correction module completed successfully:::::::::::::::::::::::::::";
    et=time(NULL);//Computes execution end time(in seconds)
    LOG(INFO)<<endl<<"Execution Time : "<<et-st<<" seconds"<<endl;
      google::ShutdownGoogleLogging ();//Close logging library
    return(EXIT_SUCCESS);

}
