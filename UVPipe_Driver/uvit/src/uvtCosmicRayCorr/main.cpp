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
#include<uvtCosmicRayCorr.h>
#include<glog/logging.h>

using namespace std;
/*
 * 
 */
int main(int argc, char** argv) {
    time_t et,st;
  st=time(NULL);//Computes execution start time (in seconds)
    int status=0;//Flag to store return status of functions
     google::InitGoogleLogging(argv[0]);//Initializing of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files.It reads the environment variable 'GLOG_log_dir'
      checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtCosmicRayCorr obj;//Creating object for uvtCosmicRayCorr class
    /*Reading PIL file */
    status=obj.read(argc,argv);//Read parameters from uvtCosmicRayCorr.par file and return status
      if(status)
    {
        cerr<<"Error in  Reading The Content of the Parameter Files";
        return(EXIT_FAILURE);
    }
    obj.display();//Display cosmicRayCorr input parameters
    status=obj.uvtCosmicRayCorrProcess();//Performing Cosmic Ray Corr Computation
    if(status)
    {
         cerr<<"\nError in CosmicRayCorrection Module"<<endl;
         return(EXIT_FAILURE);
    }
    else{
        cerr<<"\nCosmicRayCorrection Module Completed Successfully"<<endl;
    }
    et =time(NULL);//Computes execution end time(in seconds)
    LOG(INFO)<<"Execution Time ::"<<et-st<<" seconds"<<endl;
    google::ShutdownGoogleLogging ();//close logging library
    return (EXIT_SUCCESS);
}

