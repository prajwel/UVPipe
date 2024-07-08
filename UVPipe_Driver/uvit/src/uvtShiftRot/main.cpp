/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>

#include "uvtShiftRot.h"
#include<glog/logging.h>
//#include "uvtShiftRot.h"
using namespace std;
/*
 * 
 */
int main(int argc, char** argv) {

    time_t st ,et;
    st = time(NULL);//Computes execution start time (in seconds)
    int status=0;//Flag to store return status of functions
     google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'  
     checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable    
    uvtShiftRot obj;//Creating object for shiftRot class
    status=obj.read(argc,argv);//Read parameters from uvtShiftRot.par and return status
      if(status){
        cerr<<"***Error in Reading the Parameter ***"<<endl;
        return(EXIT_FAILURE);
    }
    obj.display();//Display shiftRot input parameters
    status=obj.uvtShiftRotProcess();//Performing shift rot computation
    if(status){
        cerr<<"***Error in uvtShiftRot Module***"<<endl;
        return(EXIT_FAILURE);
    }
    et = time(NULL);//COmputes executin end time(in seconds)
    cerr<<"Execution Time "<<et-st<<" Seconds"<<endl;
     google::ShutdownGoogleLogging ();//Close logging library 
    return 0;
}


