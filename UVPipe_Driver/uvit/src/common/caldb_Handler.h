/* 
 * File:   caldb_Handler.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#ifndef CALDB_HANDLER_H
#define CALDB_HANDLER_H

#define CALDB_BADPIXEL_DIR  "BAD_PIXELS" 
#define CALDB_TEMPLATE_EXP_DIR  "EXPOSURE_TEMPLATE" 
#define CALDB_FLATFIELD_DIR  "FLAT_FIELDS_FILTER"
#define  CALDB_DETECTORDISTO_DIR "DISTORTION/DETECTOR/"
#define  CALDB_OPTICALDISTO_DIR "DISTORTION/OPTICS/"
#define  CALDB_CENTROIDBIAS_DIR "BIAS_CORR"
#define  CALDB_CENTROIDCORR_DIR "AV_PH_ENERGY"
#define  CALDB_QE_DIR "QE_TEMP"
#define  CALDB_THERMAL_DIR "THERMAL_ASPECT"
#define  CALDB_PLATESCALE_DIR "PLATE_SCALE"
#define CALDB_DIST_SIZE 512
#define CALDB_TELDEF_DIR "TELDEF"
#define CALDB_CATALOGUEFILE_PATH "LOOKUP_FOR_CATALOGUE"

class caldb_Handler {
    // const char * cal_badpixeldir=""
    
    char caldbpath[NAMESIZE];
   // char detector[5];
public:
   string  getFlatFieldFile(char *detector, char *mod, char *filter, char *path);
    string getBadPixelFile(char *detector, char *mode, int xsize, int ysize, char *path);
    string  getDetectorFile(char *detector, char *path);
    string  getOpticalDetectorFile(char *detector, char *path, char *filter);
    string getCentroidBiasFile(char *detector, char *path, int winsize);
    string  getCentroidEAFile(char *detector, char *path);
    string getQEFile(char *detector, char *mod, char *path);
    string  getThermalFile(char *detector, char *path);
    int readCaldbDistFile(float *Xdist, float *Ydist, char *distortionCorrfile);
    int readCaldbOpticDistFile(float *Xdist, float *Ydist, char *distortionCorrfile) ;
    char* getTelDefFile(char *path ,char *detector,char *filter);
     char *getPlatScaleFile(char *path ,char *detector);
     string getCatalogue_LookupFile(char *path);
      string  getTemplateFileForExposure(char *detector,char *mode,int xsize, int ysize, char *path);
    // void setXdist( float *xdist);

//     caldb_Handler(char *path, char *detector){
//         strcpy(caldbpath,path);
//         //strcpy(this->detector,detector);
//     }
     
//    float * getXdist() {
//        return x_dist;
//    }
//    //   void setYdist(float *ydist);
//
//    float * getYdist() {
//        return y_dist;
//    }
    
};


struct TelDef{
    float coef_x0_a, coef_x0_b, coef_x0_c;
    float coef_y0_a, coef_y0_b,coef_y0_c;
    float det_xpix1,det_ypix1,det_xsiz,det_ysiz;
    float det_xcen,det_ycen;
    float int_xcen,int_ycen,det_xoff,det_yoff;
    int detxflip, detyflip;
    float det_scal, focal_len;
    float m11,m12,m13,m21,m22,m23,m31,m32,m33;
    float det_xscl,det_yscl;
      
    int read(char *teldeffile);
    
};

#endif	/* CALDB_HANDLER_H */
