
/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on August 20, 2013, 10:22 AM
 */
#include <cstdlib>
#include <uvtFindWtdMean.h>
#include<iostream>
#include<glog/logging.h>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    time_t st,et;
    st =time(NULL);//Comptes execution start time(inseconds)
    google::InitGoogleLogging(argv[0]);//Initialization of google logging library
    google::SetStderrLogging(google::INFO);//Setting path to create log files.It reads the environment variable 'GLOG_log_dir'
      checkParFile(argv[0]);//Check for existence of parameter file and 'PFILES' environment variable
   uvtFindWtdMean obj;//Creating object for findWtdMean class
   int status =obj.read(argc,argv);;//Read parameters from uvtFindWtdMean.par file and return status
   if(status){
        LOG(ERROR)<<"***Error in  reading parameters for star detection***"<<endl;
        return(EXIT_FAILURE);
    }
   
   obj.display();//Display findWtdMean input paramters
   
   status=obj.uvtFindWeightedMeanProcess();//Performing findWtdMean computation
   if(status){
        LOG(ERROR)<<"***Error in star detection***"<<endl;
        return(EXIT_FAILURE);
   }
   et=time(NULL);//Computes execution end time(in seconds)
   LOG(ERROR)<<"Execution Time::"<<et-st<<" Seconds"<<endl;
    google::ShutdownGoogleLogging ();//Close logging library
   return 0;
}


