
/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>
#include <uvtRefFrameCal.h>
#include<iostream>
#include<glog/logging.h>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    time_t st,et;                                 
    st =time(NULL);//Computes execution start time (in seconds)
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'      
     
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtRefFrameCal obj;//Creating object for refFrameCal class
  
    int status =obj.read(argc,argv);;//Read parameters from uvtRefFrameCal.par file and return status
    if(status){                          
        cerr<<"***Error in  reading parameters for star detection***"<<endl;
        return(EXIT_FAILURE);
    }
   
   obj.display();//Display Ref Frame Cal input parameters
  
   status=obj.uvtRefFrameCalProcess();//Performing Ref Frame Cal computation
   if(status){
        cerr<<"***Error in star detection***"<<endl;
        return(EXIT_FAILURE);
   }
  
   et=time(NULL);//Computes execution end time(in seconds)
   cerr<<"Execution Time::"<<et-st<<" Seconds"<<endl;
   
   google::ShutdownGoogleLogging ();//Close logging library
   return 0;
}

