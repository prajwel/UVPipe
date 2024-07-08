/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on June 27, 2013, 5:31 PM
 */

#include <cstdlib>
#include<iostream>
#include<ctime>
#include "uvtComputeThermal.h"
#include<glog/logging.h>
using namespace std;

int main(int argc, char** argv) {
    time_t et,st;
   st=time(NULL);//Computes execution start time (in seconds)
    int status=0;//FLag to store return status of functions
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir' 
     checkParFile(argv[0]); //Check for existence of parameter file and 'PFILES' environment variable   
    uvtComputeThermal obj;//Creates object for computeThermal class
    /*Reading PIL file */
    status=obj.read(argc,argv);//Reads parameter from uvtComputeThermal.par file and return status
    if(status)
    {
        cerr<<"***Error reading input parameters***"<<endl;;
        return(EXIT_FAILURE);
    }
  obj.display();//Display computeThermal input parameters
  if(obj.uvtThermalCalcProcess())//Performing COmpute thermal computation
  {
      cerr<<"***Error in uvtComputeThermal Module***"<<endl;
      return(EXIT_FAILURE);
  }
  else{
      cerr<<"\nuvtComputeThermal Module completed Successfully"<<endl;
  }
    et=time(NULL);//Computes execution end time(in seconds) 
    cerr<<"\nExecution Time:::"<< et-st <<" Seconds"<<endl;
       google::ShutdownGoogleLogging ();//Close logging library
    return (EXIT_SUCCESS);
}




