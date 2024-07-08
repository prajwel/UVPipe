/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include<uvtMaskBadPix.h>
#include<cstdlib>
#include<iostream>
#include<ctime>
#include<glog/logging.h>
using namespace std;

int main(int argc, char **argv){
    time_t et,st;                                                              
    st=time(NULL);//Computes execution start time(in seconds)                                                     
    
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir
    
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    //checkParFile(argv[0]);             
    int status=0;//Flag to store return status of functions
    uvtMaskBadPix obj;//Creating object for maskBadPix class
     LOG(INFO)<<"uvtMaskBadPix  STARTED....................";
    status=obj.read(argc,argv);//read parameters from uvtMaskBadPix.par file and return status
     if(status){                               
        cerr<<endl<<"Error in reading parameters for Filter Bad Pixel module ";
        return (EXIT_FAILURE);
     }
    obj.display();//Display Maskbadpix input parameters
    status=obj.uvtMaskBadPixProcess(); //Performing MaskBadPix computation
    if(status){
        cerr<<endl<<"***Error in Filter Bad Pixel module*** ";
        return (EXIT_FAILURE);
    }
    cerr<<endl<<"::::::::::::::::::::::::::Filter Bad Pixel module completed successfully:::::::::::::::::::::::::::";
    et=time(NULL);//Computes execution endtime(in seconds)
    cerr<<endl<<"Execution Time : "<<et-st<<" seconds"<<endl;     
       google::ShutdownGoogleLogging ();//Close logging library
    return(EXIT_SUCCESS);
}
