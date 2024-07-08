/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on October 24, 2013, 4:27 PM
 */

 

#include <cstdlib>

#include<iostream>
#include<ctime>

#include "uvtCentroidBias.h"
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
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtCentroidBias obj;//Creating object for centroidBias class
    status=obj.read(argc,argv);//Read parameters from uvtCentroidBias.par file and return status
    status =obj.uvtCentroidBiasProcess();//Performing Centroid Bias computation
     if(status){
        LOG(INFO)<<endl<<"***Error in  Centroid bias module*** ";
        return (EXIT_FAILURE);
     }
    
    LOG(INFO)<<endl<<"::::::::::::::::::::::::::Centroid Bias module completed successfully:::::::::::::::::::::::::::";
    et=time(NULL);//COmputes execution end time(in seconds)
    LOG(INFO)<<endl<<"Execution Time : "<<et-st<<" seconds"<<endl;
     google::ShutdownGoogleLogging ();//Close logging library
    return(EXIT_SUCCESS);

}
