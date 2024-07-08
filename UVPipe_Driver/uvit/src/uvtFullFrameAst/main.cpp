
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
#include<glog/logging.h>


#include <uvtFullFrameAst.h> 

using namespace std;


int main(int argc, char** argv) {
    time_t et,st;
    st=time(NULL);//Computes execution start time(in seconds)
    int status=0;//FLag to store return status of functions
    
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    
    uvtFullFrameAst obj;//Creating object for fullFrameAst_new class
  
    status=obj.read(argc,argv);//Read parameters from uvtFullFrameAst_new.par file and return status
    if(status)
    {
        LOG(ERROR)<<"***Error in reading parameters***";
        return(EXIT_FAILURE);
    }
    obj.display();//Display fullFrameAst_new input parameters
  
    status = obj.uvtFullFrmAstProcess();//Performing fullFrameAst_new computation
    if(status)
   {
      LOG(ERROR)<<"\033[1;31m***Error in full frame astrometry***\033[0m";
      LOG(ERROR)<<"CRASH ASTROMETRY FAILED (main.cpp in uvtFullFrameAst)";
       return(EXIT_FAILURE);
   }
   else
   {
       LOG(ERROR)<<"\033[1;34mFull Frame Astrometry completed successfully\033[0m";
   }
   et =time(NULL);//Computes execution end time(in seconds)
   LOG(INFO)<<"..................Execution Time "<<et-st<< "  seconds...................";;
    
   google::ShutdownGoogleLogging();//CLose logging library
   
    return (EXIT_SUCCESS);
}



