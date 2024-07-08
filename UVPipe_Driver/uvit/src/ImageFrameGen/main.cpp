/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on December 9, 2013, 4:04 PM
 */

#include <cstdlib>

#include "ImageFrameGen.h"
#include<glog/logging.h>
#include<csignal>
using namespace std;

/*
 * 
 */
void signalHandler( int signum )
{
    char temp;
    cout << "\033[1;31mDo you Want to close  the Execution??[y/n].\n\033[0m";//Do your handling here
    cin>>temp;
    if (temp=='y'){
       exit(signum);  
    }   
}
int main(int argc, char** argv) {

     signal(SIGINT, signalHandler);  //Interrupt will be captured here.
     while(1){
     time_t st,et;
      google::InitGoogleLogging(argv[0]);                //Initialization of google logging library
    google::SetStderrLogging(google::INFO);      //Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    st=time(NULL);
    uvtImRa_commonArray obj;
int     status=obj.read(argc,argv);
    if(status){
        cerr<<endl<<"Error in reading parameters"<<endl;
        return (EXIT_FAILURE);
    }
    obj.uvtImRacomArrProcess();
    
   et=time(NULL);
   google::ShutdownGoogleLogging ();        //Close logging library
   cout<<endl<<"Execution time :"<<et-st<<"  seconds"<<endl;
    return (EXIT_SUCCESS);
     }
}

