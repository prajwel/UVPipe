/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <cstdlib>

#include "uvtComputeJitter.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    time_t st ,et;
    st = time(NULL);//Computes execution start time(i9n seconds)
    int status=0;//FLag to store return status of functions
    uvtComputeJitter obj;//Creating object for computeJitter class
    status=obj.read(argc,argv);//Read parameters from uvtComputeJitter.par file and return status
      if(status){
        cerr<<"***Error in Reading the Parameter ***"<<endl;
        return(EXIT_FAILURE);
    }
    obj.display();//Display computeJitter input parameters
    status=obj.uvtComputeJitterProcess();//Performing Jitter computation
    if(status){
        cerr<<"***Error in uvtDetectDestortion Module***"<<endl;
        return(EXIT_FAILURE);
    }
    et = time(NULL);//Computes execution end time(in seconds)
    cerr<<"Execution Time "<<et-st<<" Seconds"<<endl;
    
    return 0;
}


