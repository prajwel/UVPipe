/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on July 4, 2013, 4:16 PM
 */

#include<cstdlib>
#include<iostream>
#include<ctime>
#include<DataIngest.h>
#include<glog/logging.h>
#include<uvtUtils.h>

using namespace std;

int main(int argc, char **argv)
{
    time_t et,st;                                       
    st=time(NULL);//Computes execution start time(in seconds)
    
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    
    checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
    
    int status=0;//Flag to store return status of functions
    DataIngest obj;//Creating object for DataIngest class
    
    LOG(INFO)<<"DATAINGEST STARTED....................";              
    
    status=obj.read(argc,argv);//Read parameters from DataIngest.par file and return status
    if(status){                                              
         LOG(ERROR)<<"***Error in DataIngest module*** ";
        return(EXIT_FAILURE);
     }
    
    obj.display();//Display DataIngest input parameters
    
    status=obj.DataIngestProcess();//Performing Data Ingest Computation
    if(status){                                             
        //cerr<<endl<<"***Error in DataIngest module*** ";
        LOG(ERROR)<<"***Error in DataIngest module*** ";
LOG(ERROR)<<"CRASH   DATAINGEST FAILED (main.cpp in DataIngest)";
        return (EXIT_FAILURE);
    }else{
        LOG(INFO)<<"::::::::::::::::::::::::::DataIngest completed successfully:::::::::::::::::::::::::::";
    }
    
    et=time(NULL);//Computes execution end time(in seconds)
    //cerr<<endl<<"Execution Time : "<<et-st<<" seconds"<<endl;
    LOG(INFO)<<"Execution Time : "<<et-st<<" seconds"<<endl;  
    
    google::ShutdownGoogleLogging ();//Close logging library
    
    return(EXIT_SUCCESS);
}
