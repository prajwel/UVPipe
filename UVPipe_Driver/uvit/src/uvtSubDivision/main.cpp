/* 
 * File:   main.cpp
 * Authors:: Dhruv,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */


#include <cstdlib>

#include<iostream>
#include<ctime>
#include<glog/logging.h>

#include "uvtSubDivision.h"


using namespace std;



/*
 * 
 */
int main(int argc, char** argv) {
    time_t et,st;                                 
   st=time(NULL);//Computes execution start time in seconds 
    int status=0;//Flag to store return status of functions
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
      checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtSubDivision obj;//Creating object for subDivision class 
    status=obj.read(argc,argv);//Read parameters from uvtSubDivision.par file and return status
    if(status)                        
    {
        cerr<<endl<<"Error in  Reading The Content of the Parameter Files";
        return(EXIT_FAILURE);
    }
    
    obj.display();//Display subDivision input parameters
    
   if(obj.uvtSubDivisionProcess())//Performing Sub Division calculation
   {
       cerr<<"Error in SubDivision Module"<<endl;
       return(EXIT_FAILURE);
   }
   else
   {
       cerr<<"SubDivision Module Completed Successfully"<<endl;
   }
    
   et =time(NULL);//Computes execution end time (in seconds)
   cout<<"\nExecution Time "<<et-st<< "seconds"<<endl;
   google::ShutdownGoogleLogging ();//Close logging library
    
    return (EXIT_SUCCESS);
}



