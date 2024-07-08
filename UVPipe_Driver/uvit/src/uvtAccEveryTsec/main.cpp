/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 * Created on July 25, 2013, 1:29 PM
 */

#include <cstdlib>
#include<iostream>
#include<ctime>
#include"uvtAccEveryTsec.h"
#include<glog/logging.h>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

     time_t et,st;
    st=time(NULL);//Computes execution start time(in seconds)
    
    int status=0;//Flag to store return status of functions
    google::InitGoogleLogging(argv[0]);//Initializing of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtAccEveryTsec obj;//Creating object for AccEveryTsec class
    status=obj.read(argc,argv);//Read parameters from uvtAccEveryTsec.par file and return status
     if(status)
        cerr<<endl<<"***Error in Accumulation module*** ";
    obj.display();//Display accEveryTsec input parameters
    status=obj.uvtAccEveryTsecProcess();//Performing  Acc Every Tsc Computation
    if(status){
        cerr<<endl<<"***Error in Accumulation module*** ";
        return(EXIT_FAILURE);
    }
    else{
        cerr<<endl<<"::::::::::::::::::::::::::Accumulation module completed successfully:::::::::::::::::::::::::::";
           return(EXIT_FAILURE);
    }
    et=time(NULL);//Computes execution time end time(in seconds)
    cerr<<endl<<"Execution Time : "<<et-st<<" seconds"<<endl;
      google::ShutdownGoogleLogging ();//Close logging library
    return(EXIT_SUCCESS);
    
    return 0;
}

