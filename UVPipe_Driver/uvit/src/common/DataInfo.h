/* 
 * File:   DataInfo.h
 * Authors: Preeti Tahlani,  Dhruv Shelat
 * Technical Guidance: Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *Description: In UVIT software, an information file acts as a source of auxiliary information about the 
 *                      data at every processing step. Every module reads information file and also generates 
 *                      an information file with updated data.
 *                      This file contains a class having data structures to hold information contained in information 
 *                      file and also interfaces to read and write data from info file. 
 * Created on September 2, 2013, 3:26 PM
 */


#ifndef DATAINFO_H
#define	DATAINFO_H

#include<fitsio.h>
#include<cstring>
#include<glog/logging.h>
#define PC 0
#define IM 1

#define FF_COLNAME "EFFECTIVE_NUM_PHOTONS"

class DataInfo
{
    double tstarti,tstartf,tstopi,tstopf;            //Integer and fractional parts of start and stop time for observation data
    char *obsmode;                                        //Observation mode
    bool modeflag;                                          //modeflag= 1 for IM /  0 for PC
    char *filter;                                              //Filter used for observation 
    char *source;                                        //data source        
    double integrationTime;                      //Integration time for data   
    char *detector;                                  //channels - FUV,NUV,VIS
    int xsize;                                             //Width of frame
    int ysize;                                             //Height of frame
 
 public:
    int x , y , xoff , yoff ;                      // for window - x->xsize, y->ysize, xoff ->X offset, yoff -> Y offset    
    DataInfo();
    ~DataInfo();
    
    void setTime(double tstarti,double tstartf,double tstopi,double tstopf)  { this->tstarti=tstarti;  this->tstopi=tstopi; 
                                                                                                                                this->tstartf=tstartf;   this->tstopf=tstopf;    }
    void setObsMode(char *obsmode);
    void setObsMode(bool modeflag);
    void setFilter(char *filter){   strcpy(this->filter,filter);  }
    void setSource(char *source)  {  strcpy(this->source,source); }
    void setIntegrationTime(double t)   { integrationTime=t;  }
    void setDetector(char *d)  { strcpy(detector,d);  }
    void setXsize(int xsize)   { this->xsize=xsize;  }
    void setYsize(int ysize)   { this->ysize=ysize;  }
    
    void getTime(double *tstarti, double *tstartf, double *tstopi, double *tstopf) {  *tstarti=this->tstarti;   *tstartf=this->tstartf;   
                                                                                                                                        *tstopi=this->tstopi;  *tstopf=this->tstopf;  }
    char *getObsMode() const  { return obsmode; }
    bool getModeFlag()   const  { return modeflag; }
    char *getFilter()  const { return filter;}
    char *getsource()  const { return source; }
    double getIntegrationTime() const  { return integrationTime;  }
    char *getDetector()   const { return detector;  }
    double getTimeDuration()  const {   return ((tstopi+tstopf) - (tstarti+tstartf));   }
    int getXsize() const { return xsize; }
    int getYsize() const { return ysize; }
    
    double getTstart()  { return (tstarti+tstartf); }
    double getTstop()  { return (tstopi+tstopf); }
    
    void write(fitsfile *fptr);             //function to write information to a fits file (already created)
    void getInfo(fitsfile *fptr);            //function to read information from a fits file (already created)
};

#endif	/* DATAINFO_H */

