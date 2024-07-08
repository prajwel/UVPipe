/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>

#include "uvtFrameIntegration.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    time_t st ,et;
    st = time(NULL);//Computes execution start time(in seconds)
    int status=0;//Flag to store return status of functions
    uvtFrameIntegration obj;//Creating object for frameIntegration class
    status=obj.read(argc,argv);//Read parameters from uvtFrameIntegration.par file and return status
      if(status){
        cerr<<"***Error in Reading the Parameter ***"<<endl;
        return(EXIT_FAILURE);
    }
    obj.display();//Display frameIntegration input parameters
   status=obj.uvtFrameIntProcess();//Performing Frame Integration Computation
    if(status){
        cerr<<"***Error in uvtFrameIntegration Module***"<<endl;
        return(EXIT_FAILURE);
    }
    et = time(NULL);//Computes execution end time(in seonds)
    cerr<<"Execution Time "<<et-st<<" Seconds"<<endl;
    
    return 0;
}

