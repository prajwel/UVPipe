/* 
 * File:   transform.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna K Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef TRANSFORM_H
#define	TRANSFORM_H

#include<caldb_Handler.h>
/*
#define XSCALE_ARCSEC_VIS 3.357
#define YSCALE_ARCSEC_VIS  3.311
#define XSCALE_ARCSEC_NUV  3.3307
#define YSCALE_ARCSEC_NUV  3.3307
#define XSCALE_ARCSEC_FUV 3.3373
#define YSCALE_ARCSEC_FUV 3.3373
*/

// NEW VALUES for xscale /yscale(6 lines below) : 26-July-2017

#define XSCALE_ARCSEC_VIS 3.357
#define YSCALE_ARCSEC_VIS  3.311
#define XSCALE_ARCSEC_NUV  3.33
#define YSCALE_ARCSEC_NUV  3.33
#define XSCALE_ARCSEC_FUV 3.33
#define YSCALE_ARCSEC_FUV 3.33

/*
#define ANGLE_NUV 31.26
#define ANGLE_FUV 0.733 //original value
*/


 //NEW values (+0.25 degree for NUV & possibilities for FUV) 19-June-2017 :

#define ANGLE_NUV 31.51
// #define ANGLE_FUV 0.733 //original value
// #define ANGLE_FUV 0.983 // if +0.25 deg
// #define ANGLE_FUV 0.483 // if -0.25 deg
#define ANGLE_FUV 0.483


//#define ANGLE_FUV 60.733
#define ANGLE_VIS 34.13
#define VIS_TO_NUV_COEFF_1 1.03458
#define VIS_TO_NUV_COEFF_2 -0.04366
#define VIS_TO_NUV_COEFF_3 -0.04366
#define VIS_TO_NUV_COEFF_4 -1.03458

#define  FUV_TO_NUV_COEFF_1 0.854078
#define  FUV_TO_NUV_COEFF_2 0.52394
#define  FUV_TO_NUV_COEFF_3 0.52394
#define  FUV_TO_NUV_COEFF_4 -0.854078

#define VIS_TO_FUV_COEFF_1  0.85734 
#define VIS_TO_FUV_COEFF_2 -0.57706
#define VIS_TO_FUV_COEFF_3 0.57706
#define VIS_TO_FUV_COEFF_4 0.85734

#define RPYTOXYTHETA_COEFF_1_VIS  0.1627
#define RPYTOXYTHETA_COEFF_2_VIS 0.2400
#define RPYTOXYTHETA_COEFF_3_VIS  0.2400
#define RPYTOXYTHETA_COEFF_4_VIS -0.1627

/*
     4 lines below have coefficients from V1.5 release 17-June-2016
#define RPYTOXYTHETA_COEFF_1_NUV  0.15785
#define RPYTOXYTHETA_COEFF_2_NUV 0.2554
#define RPYTOXYTHETA_COEFF_3_NUV  -0.2554
#define RPYTOXYTHETA_COEFF_4_NUV 0.15785

*/
/*
// 4 lines below are coefficients deduced on 01-Jul-2016 from SNT's note
#define RPYTOXYTHETA_COEFF_1_NUV  0.15313
#define RPYTOXYTHETA_COEFF_2_NUV 0.2522
#define RPYTOXYTHETA_COEFF_3_NUV  -0.2522
#define RPYTOXYTHETA_COEFF_4_NUV 0.15313
*/

//++++++
// NEW values (+0.25 degree) 13-June-2017 :


#define RPYTOXYTHETA_COEFF_1_NUV 0.15423
#define RPYTOXYTHETA_COEFF_2_NUV 0.25153
#define RPYTOXYTHETA_COEFF_3_NUV -0.25153
#define RPYTOXYTHETA_COEFF_4_NUV 0.15423


/*
     4 lines below have coefficients from V1.55 release 13-July-2016
#define RPYTOXYTHETA_COEFF_1_FUV  0.001
#define RPYTOXYTHETA_COEFF_2_FUV 0.29965
#define RPYTOXYTHETA_COEFF_3_FUV  0.29965
#define RPYTOXYTHETA_COEFF_4_FUV -0.001
*/

/* four lines below : being COMMENTED on 19-Jun-2017 (being replaced by 
   NEW values from Shyam+Prajwel tests ....

// 4 lines below are coefficients from SNT's note (based on
// Joe's findings of 21-Jan-2016 [inserted on 19-July-2016 by SKG]
#define RPYTOXYTHETA_COEFF_1_FUV  -0.0037911
#define RPYTOXYTHETA_COEFF_2_FUV  0.29635469
#define RPYTOXYTHETA_COEFF_3_FUV  0.29635469
#define RPYTOXYTHETA_COEFF_4_FUV  0.0037911

*/

// Added NEW values (FUV) 4 lines below (19-Jun-2017)

#define RPYTOXYTHETA_COEFF_1_FUV  -0.00249798
#define RPYTOXYTHETA_COEFF_2_FUV  0.29636841
#define RPYTOXYTHETA_COEFF_3_FUV  0.29636841
#define RPYTOXYTHETA_COEFF_4_FUV  0.00249798


//Added on 23august testing 

//#define RPYTOXYTHETA_COEFF_1_FUV  0.00383
//#define RPYTOXYTHETA_COEFF_2_FUV  0.2996
//#define RPYTOXYTHETA_COEFF_3_FUV  -0.2996
//#define RPYTOXYTHETA_COEFF_4_FUV  0.00383

//
//#define SHIFTS_TO_PY_NUV_COEFF_ONE  1.03458
//#define SHIFTS_TO_PY_FUV_COEFF_ONE  0.85734
//#define SHIFTS_TO_PY_VIS_COEFF_1  1.9353
//#define SHIFTS_TO_PY_NUV_COEFF_TWO 0.04366
//#define SHIFTS_TO_PY_FUV_COEFF_TWO 0.57706
//#define SHIFTS_TO_PY_VIS_COEFF_2 2.8547
#define ROLL_TERM 0


 static bool  Flag_NUVtoVIS=FALSE,Flag_FUVtoVIS=FALSE,Flag_dtermtoRPY_vis=FALSE,Flag_dtermtoRPY_nuv=FALSE,Flag_dtermtoRPY_fuv=FALSE,Flag_NUVtoFUV=FALSE;
 static double nuvTovis_inverse1,nuvTovis_inverse2,nuvTovis_inverse3,nuvTovis_inverse4;
 static double fuvTovis_inverse1,fuvTovis_inverse2,fuvTovis_inverse3,fuvTovis_inverse4;
 static double nuvTofuv_inverse1, nuvTofuv_inverse2, nuvTofuv_inverse3, nuvTofuv_inverse4;
 static double dtermToRPY_VIS_inverse1,dtermToRPY_VIS_inverse2,dtermToRPY_VIS_inverse3,dtermToRPY_VIS_inverse4;
 static double dtermToRPY_NUV_inverse1,dtermToRPY_NUV_inverse2,dtermToRPY_NUV_inverse3,dtermToRPY_NUV_inverse4;
 static double dtermToRPY_FUV_inverse1,dtermToRPY_FUV_inverse2,dtermToRPY_FUV_inverse3,dtermToRPY_FUV_inverse4;
 
using namespace std;

/**
 * Function to transform drift series,jitter series,thermal series from CH1 to CH2 frame
 * @param scaleNUV-scale for NUV 
 * @param scaleVIS-scale for VIS
 * @param angle_bet_NUVnVIS -Angle between NUV and VIS 
 * @param tx - Translation between NUV and VIS in X direction
 * @param ty - Translation between NUV and VIS in y direction
 * @param deltaX- input delta X array which will be updated accordingly
 * @param deltaY-input delta Y array which will be updated accordingly
 * @param nsamples- number of samples at input.basically it is the size of DeltaX array
 * @return 
 */
int transformCH1toCH2 ( float scalexNUV,float scalexVIS,float scaleyNUV,float scaleyVIS,float angle_bet_NUVnVIS,float tx,float ty,double  *deltaX, double *deltaY,int nsamples);


/**
 * Function to transform Roll, Pitch, Yaw to dx, dy and dtheta
 * @param perdegreePix- number of pixels in 1 degree
 * @param yaw-yaw array
 * @param pitch-pitch array
 * @param nsamples- number of samples at input ,basically it is the size of yaw/pitch array
 * @return 
 */
int RPY2dxdy(float  perdegreePix_x,float perdegreePix_y,double *yaw,double *pitch ,int nsamples);
/**
 * Function to transform spacecraft to VIS
 * @param perdegpix-number of pixels in 1 degree
 * @param scale-
 * @param eulerangle_One-Euler angle one of rotation from spacecraft frame to VIS frame
 * @param eulerangle_Two-Euler angle two  of rotation from spacecraft frame to VIS frame
 * @param eulerangle_Three-Euler angle three  of rotation from spacecraft frame to VIS frame
 * @param deltaX-input delta X array which will be updated accordingly
 * @param deltaY- input delta Y array which will be updated accordingly
 * @param nsamples- number of samples at input ,basically it is the size of yaw/pitch array
 * @return 
 */
 int transformSpacecraft2VIS(float perdegpix_x,float perdegpix_y, float scale, float eulerangle_One,float eulerangle_Two,float eulerangle_Three,double *deltaX, double *deltaY, int nsamples);
// int  transformNUVtoVIS(double  *Xnuv,double *Ynuv,double  *thetanuv,double *Xvis,double *Yvis,double *thetavis);

 int transformNUVtoVIS(double  &Xnuv,double &Ynuv,double  &thetanuv,double &Xvis,double &Yvis,double &thetavis);
 int transformVIStoNUV(float *Xvis,float *Yvis,float *thetavis,float *Xnuv,float *Ynuv,float *thetanuv);
int transformFUVtoNUV (double &Xfuv , double  &Yfuv , double  &thetafuv ,double  &Xnuv , double &Ynuv , double &thetanuv );
 int transformVIStoFUV(float *Xvis,float *Yvis,float *thetavis,float *Xfuv,float *Yfuv,float *thetafuv);
int transformNUVtoFUV (double  &Xnuv , double &Ynuv , double  &thetanuv , double  &Xfuv , double  &Yfuv , double  &thetafuv );

int transformRPYtoDXDYDTHETA_VIS(double &roll,double  &pitch ,double  &yaw,double &dx,double &dy,double  &dtheta);
int transformRPYtoDXDYDTHETA_NUV(double &roll,double  &pitch ,double  &yaw,double &dx,double &dy,double  &dtheta);
int transformRPYtoDXDYDTHETA_FUV(double &roll,double  &pitch ,double  &yaw,double &dx,double &dy,double  &dtheta);

/**
 * Functions to transform Xshift,yshift and rotation to Roll,pitch and yaw.in case of  VIS/NUV/FUV channel.
 * @param dx-x shift 
 * @param dy-yshift
 * @param dtheta -rotation angle in degree
 * @param roll- theta to be converted in to roll
 * @param pitch-Xshift to be converted in  to pitch
 * @param yaw -yshift to be converted in to yaw
 * @return 
 */
int transformDXDYDTHETAtoRPY_VIS(double &dx,double &dy ,double  &dtheta,double  &roll,double  &pitch ,double  &yaw);
int transformDXDYDTHETAtoRPY_NUV(double &dx,double &dy ,double  &dtheta,double  &roll,double  &pitch ,double  &yaw);
int transformDXDYDTHETAtoRPY_FUV(double &dx,double &dy ,double  &dtheta,double  &roll,double  &pitch ,double  &yaw);

int readParamsFrmteldef(char *caldbfile,float *shiftx,float *shifty);
int readPlateScaleFile ( char* filename, float *scaleX,float *scaleY);
//
//struct TelDef{
//    float coef_x0_a, coef_x0_b, coef_x0_c;
//    float coef_y0_a, coef_y0_b,coef_y0_c;
//    float det_xpix1,det_ypix1,det_xsiz,det_ysiz;
//    float det_xcen,det_ycen;
//    float int_xcen,int_ycen,det_xoff,det_yoff;
//    int detxflip, detyflip;
//    float det_scal, focal_len;
//    float m11,m12,m13,m21,m22,m23,m31,m32,m33;
//      
//    int read(char *teldeffile);
//    
//};

//transformation from raw to satellite coordinates
int RawToSat (double *deltaX, double *deltaY, int nsamples, char *teldef);

int transform_uvitRPY_to_spacecraftRPY();

#endif	/* TRANSFORM_H */
