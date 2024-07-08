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
#include "uvtPixPadding.h"
#include<glog/logging.h>
using namespace std;



/*
 * 
 */
int main(int argc, char** argv) {
    time_t et,st;                                                        
    st=time(NULL); //Computes execution start time(in seconds)           
    int status=0; //Flag to store return status of functions
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtPixPadding obj;//Creating object for pixPadding class
    status=obj.read(argc,argv);//Read parameters from uvtPixPadding.par file and return status                          
    if(status)                      
    {
        LOG(INFO)<<"***Error reading input parameters***"<<endl;;
        return(EXIT_FAILURE);
    }
  obj.display();//Display Pixel  padding input parameters
  
  if(obj.uvtPixPaddingProcess())//Performing Pix Padding computation 
  {
      LOG(INFO)<<"***Error in PixPadding Module***"<<endl;
      return(EXIT_FAILURE);
  }
  else{
      LOG(INFO)<<"\nPixPadding Module completed Successfully"<<endl;
  }
    
    et=time(NULL);//Computes execution time(in seconds)
    LOG(INFO)<<"\nExecution Time:::"<< et-st <<" Seconds"<<endl;
   
    google::ShutdownGoogleLogging ();//Close logging library
    return (EXIT_SUCCESS);
}



