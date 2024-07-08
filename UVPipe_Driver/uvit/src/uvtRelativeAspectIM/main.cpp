/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on December 9, 2013, 4:04 PM
 */

#include <cstdlib>

#include "uvtRelativeAspectIM.h"
#include<glog/logging.h>
#include<csignal>
#include<unistd.h>

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
         
           pid_t parentPid=getppid();
//    const char* cmdline=get_process_name_by_pid(parentPid);
    string cmdline=get_process_name_by_pid(parentPid);
    //cout<<cmdline;
//      if(strcmp(cmdline.c_str(),"UVIT_DriverModule")!=0){
//    //LOG(INFO)<<cmdline;
//          
//    
//    if(strstr(cmdline.c_str(),"-")==NULL){
//        LOG(ERROR)<<"\033[1;31m***SHELL used is NOT C-SHELL!!,Please change the shell to C-shell (just type 'csh' on terminal ) and try again!!***\033[0m";
//                return(EXIT_FAILURE);
//    }
//    }
     time_t st,et;
      google::InitGoogleLogging(argv[0]);                //Initialization of google logging library
    google::SetStderrLogging(google::INFO);      //Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
    st=time(NULL);
    uvtRelativeAspectIM  obj;
int     status=obj.readPILParameters(argc,argv);
    if(status){
        cerr<<endl<<"Error in reading parameters"<<endl;
        return (EXIT_FAILURE);
    }
    status=obj.uvtRelativeAspectIMProcess();
    if (status)
        {
            LOG (ERROR) << "Error in relative aspect chain for IM mode" ;
	    LOG(ERROR)<<"CRASH RA_IM CHAIN FAILED (main.cpp in uvtRelativeAspectIM)";
            return (EXIT_FAILURE) ;
        }
   et=time(NULL);
   google::ShutdownGoogleLogging ();        //Close logging library
   cout<<endl<<"Execution time :"<<et-st<<"  seconds"<<endl;
    return (EXIT_SUCCESS);
     }
}

