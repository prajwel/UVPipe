/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */


#include<uvtQEMCPCorr.h>
#include<cstdlib>
#include<iostream>
#include<ctime>
#include<glog/logging.h>

using namespace std;

int main(int argc, char **argv){
    time_t et,st;
    st=time(NULL);//Computes execution start time (in seconds)
     int status=0;//Flag to store return status of functions
      google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'

    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
   
      uvtQEMCPCorr obj;//Creating object for QEMCPCorr class
    
      status=obj.read(argc,argv);//Read parameters from uvtQEMCPCorr.par file and return status
      if(status)
        LOG(ERROR)<<endl<<"***Error in QE MCP correction module*** ";
    obj.display();//Display QEMPCorr input parameters 
    status=obj.uvtQEMCPCorrectionProcess();//Performing QEMCP Corr computation
    if(status){
        LOG(ERROR)<<endl<<"***Error in QE MCP correction module*** ";
        return(EXIT_FAILURE);
    }
    else
    {
        LOG(INFO)<<endl<<"::::::::::::::::::::::::::Filter QE MCP correction completed successfully:::::::::::::::::::::::::::";
       
    }
    et=time(NULL);//Computes execution end time(in seconds)
    LOG(INFO)<<endl<<"Execution Time : "<<et-st<<" seconds"<<endl;
    google::ShutdownGoogleLogging ();//CLose logging library
    return(EXIT_SUCCESS);
}
