/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on December 9, 2013, 4:04 PM
 */

#include <cstdlib>

#include "uvtLevel2PC.h"
#include<glog/logging.h>
#include<csignal>
#include <fstream>
using namespace std;
void signalHandler( int signum )
{
    char temp;
    cout << "\033[1;31mDo you Want to close  the Execution??[y/n].\n\033[0m";//Do your handling here
    cin>>temp;
    if (temp=='y'){
       exit(signum);  
    }   
}
/*
 * 
 */
int main(int argc, char** argv) {
    
    
    
    
    signal(SIGINT, signalHandler);  //Interrupt will be captured here.
    while(1){
         pid_t parentPid=getppid();
//    const char* cmdline=get_process_name_by_pid(parentPid);
    string cmdline=get_process_name_by_pid(parentPid);
//    if(strcmp(cmdline.c_str(),"UVIT_DriverModule")!=0){
//    //LOG(INFO)<<cmdline;
//    
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
    
    uvtLevel2PC obj;
int status=obj.readPILParameters(argc,argv);
st=time(NULL);
    if(status)
    {
        cerr<<endl<<"Error in reading parameters"<<endl;
        return (EXIT_FAILURE);
    }
    status=obj.uvtLevel2PCprocess ();
    if(status){
        LOG(ERROR)<<"**Error in process***";
	LOG(ERROR)<<" CRASH L2PC FAILED (main.cpp in uvtLevel2PC)";
        return(EXIT_FAILURE);
    }
   et=time(NULL);
   google::ShutdownGoogleLogging ();        //Close logging library
   
 
   vector<string> filename_scienceDataFail=obj.getFailedScienceDataFilename();
   //write into information file.
//   fitsfile *fout;
   LOG(INFO)<<filename_scienceDataFail.size()<<" "<<obj.orbnum;
  // char outfile[FLEN_FILENAME];
  // int tfields = 1 ;
    //char *ttype[] = {"SCIENCEDATA_FILENAME " } ;
    //char *tform[] = {"A256"} ;
    ofstream of1;
    string name_File=TIME_MISMATCH_FILENAME;
    of1.open((char*)name_File.c_str(),ios::app);
    for (int i=0;i<filename_scienceDataFail.size();i++)
    {
        of1<<filename_scienceDataFail[i]<<endl;        
    }
    
    of1.close();
    
    
    
//    sprintf (outfile , "%s%s_%s.fits" , "Driver_",obj.orbnum , "_NO_VIS_FOUND" ) ;
//        fits_create_file (&fout , outfile , &status) ;
//        printError (status , "Error creating the output File " , outfile) ;
//        fits_create_tbl (fout , ASCII_TBL , 0 , tfields , ttype , tform , NULL , "FILENAMES" , &status) ;
//        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
//        fits_write_col (fout , TSTRING , colnum , firstrow , firstelem , nrows , (void *) filename_scienceDataFail , &status) ;
//        printError (status , "***Error writing filename***" , outfile) ;
//        fits_close_file (fout , &status) ;
//        printError (status , "Error closing the file " , outfile) ;

      
      
   cout<<endl<<"Execution time :"<<et-st<<"  seconds"<<endl;
   return (EXIT_SUCCESS);
    }
    
}

