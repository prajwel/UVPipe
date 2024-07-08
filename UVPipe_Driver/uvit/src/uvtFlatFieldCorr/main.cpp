
/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>
//#include<cstdlib>
#include<iostream>
#include<ctime>
#include<uvtFlatFieldCorr.h>
#include<glog/logging.h>

using namespace std;
/*
 * 
 */
int main(int argc, char** argv) {
    time_t et,st;
    st=time(NULL);//Computes execution start time(in seconds)
    int status=0;//Flag to store return status of functions
     google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to craete log files. It reads the environment variable 'GLOG_log_dir'
      checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    uvtFlatFieldCorr obj;//Creating object for flatFieldCorr class
          
        status=obj.read(argc,argv);//Read parameters fom uvtFlatFieldCorr.par and return status
        if(status)
        {
            cerr<<endl<<"***Error in Parameter Reading*** ";
            return (EXIT_FAILURE);
        }
        //displaying the PIL content which has been read.
        obj.display();//Display uvtFlatFieldCorr input parameters 
        
        /*the process of FlatField Correction */
        status= obj.uvtFlatFieldCorrProcess ();//Performing Flat Field Corr calculation
        if(status){
            cerr<<endl<<"***Error in FlatFieldCorrection module *** ";
        return(EXIT_FAILURE);
        }
        else
        {
            cerr<<endl<<"---------- FlatFieldCorrection module completed successfully----------"<<endl;
            
        }
        et=time(NULL);//Computes execution end time(in seconds)
        cout<<"\nExecution Time : "<<et-st<<" seconds "<<endl;
        google::ShutdownGoogleLogging ();//Close logging library
        return (EXIT_SUCCESS);
}

