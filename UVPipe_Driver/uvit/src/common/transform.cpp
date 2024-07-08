/* 
 * File:   transform.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include<iostream>
#include<stdlib.h>
#include<cstdlib>
#include<vector>
#include<math.h>
#include<fitsio.h>
#include<uvtUtils.h>
#include<fitsio.h>
#include<spMatrix.h>
#include "transform.h"

using namespace std ;
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function
#define lj 89

int transformCH1toCH2 (float scalexNUV , float scalexVIS , float scaleyNUV , float scaleyVIS , float angle_bet_NUVnVIS , float tx , float ty , double *deltaX , double *deltaY , int nsamples)
{
    float scaleFactorx = scalexVIS / scalexNUV ; //scalefactor
    float scaleFactory = scaleyVIS / scaleyNUV ; //scalefactor
    float angle_rad = angle_bet_NUVnVIS * M_PI / 180 ;

    for (int i = 0 ; i < nsamples ; i ++)
    {
        deltaX[i] = scaleFactorx * (deltaX[i] * cos (angle_rad) + deltaY[i] * sin (angle_rad)) + tx ;
        deltaY[i] = scaleFactory * (deltaY[i] * cos (angle_rad) - deltaX[i] * sin (angle_rad)) + ty ;
    }
    return (EXIT_SUCCESS) ;
}


int RPY2dxdy (float perdegreePix_x , float perdegreePix_y , double *yaw , double *pitch , int nsamples)
{
    for (int i = 0 ; i < nsamples ; i ++)
    {

        yaw[i] = yaw[i] / perdegreePix_x ;
        pitch[i] = pitch[i] / perdegreePix_y ;

    }

    return (EXIT_SUCCESS) ;
}
int transform_uvitRPY_to_spacecraftRPY(){
    
}

int transformSpacecraft2VIS (float perdegpix_x , float perdegpix_y , float scale , float eulerangle_One , float eulerangle_Two , float eulerangle_Three , double *deltaX , double *deltaY , int nsamples)
{
    float Eul_One_rad , Eul_Two_rad , Eul_Three_rad ;

    Eul_One_rad = eulerangle_One * M_PI / 180 ;
    Eul_Two_rad = eulerangle_Two * M_PI / 180 ;
    Eul_Three_rad = eulerangle_Three * M_PI / 180 ;

    for (int i = 0 ; i < nsamples ; i ++)
    {
        deltaX[i] = scale * (deltaX[i] * cos (Eul_Three_rad) + deltaY[i] * sin (Eul_Three_rad)) + Eul_One_rad / perdegpix_x ;
        deltaY[i] = scale * (deltaY[i] * cos (Eul_Three_rad) - deltaX[i] * sin (Eul_Three_rad)) + Eul_Two_rad / perdegpix_y ;
    }

    return (EXIT_SUCCESS) ;
}


int readParamsFrmteldef (char *caldbfile , float *shiftx , float *shifty)
{
    fitsfile *fptr ;
    int status = 0 ;
    fits_open_file (&fptr , caldbfile , READONLY , &status) ;
    printError (status , "***input File reading Fails***") ;
    fits_read_key (fptr , TFLOAT , "TX_NV" , shiftx , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , caldbfile) ;
    fits_read_key (fptr , TFLOAT , "TY_NV" , shifty , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , caldbfile) ;
    return (EXIT_SUCCESS) ;
}


int readPlateScaleFile (char* filename , float *scaleX , float *scaleY)
{
    fitsfile *fptr ;
    int status = 0 ;

    fits_open_file (&fptr , filename , READONLY , &status) ;
    printError (status , " Error in opening the file " , filename) ;
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "Error in Moving the 2nd HDU" , filename) ;
    fits_read_col (fptr , TFLOAT , 1 , 1 , 1 , 1 , NULL , scaleX , NULL , &status) ;
    printError (status , "Error in reading the column " , filename) ;
    fits_read_col (fptr , TFLOAT , 2 , 1 , 1 , 1 , NULL , scaleY , NULL , &status) ;
    printError (status , "Error in reading the column " , filename) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file" , filename) ;

    return (EXIT_SUCCESS) ;
}

// int TelDef::read(char* teldeffile){
//       
//   int flag = readKeywords(teldeffile,1,27,TFLOAT,"COE_X0_A",&coef_x0_a,
//                                                TFLOAT,"COE_X0_B",&coef_x0_b,
//                                                TFLOAT,"COE_X0_C",&coef_x0_c,
//                                                TFLOAT,"COE_Y0_A",&coef_y0_a,
//                                                TFLOAT,"COE_Y0_B",&coef_y0_b,
//                                                TFLOAT,"COE_Y0_C",&coef_y0_c,
//                                                TFLOAT,"DETXPIX1",&det_xpix1,
//                                                TFLOAT,"DETYPIX1",&det_ypix1,
//                                                TFLOAT,"DET_XSIZ",&det_xsiz,
//                                                TFLOAT,"DET_YSIZ",&det_ysiz,
//                                                TFLOAT,"INT_XCEN",&int_xcen,
//                                                TFLOAT,"INT_YCEN",&int_ycen,
//                                                TINT,"DETXFLIP",&detxflip,
//                                                TINT,"DETYFLIP",&detyflip,
//                                                TFLOAT,"DET_XOFF",&det_xoff,
//                                                TFLOAT,"DET_YOFF",&det_yoff,
//                                                TFLOAT,"DET_SCAL",&det_scal,
//                                                TFLOAT,"FOCALLEN",&focal_len,
//                                                TFLOAT,"ALIGNM11",&m11,
//                                                TFLOAT,"ALIGNM12",&m12,
//                                                TFLOAT,"ALIGNM13",&m13,
//                                                TFLOAT,"ALIGNM21",&m21,
//                                                TFLOAT,"ALIGNM22",&m22,
//                                                TFLOAT,"ALIGNM23",&m23,
//                                                TFLOAT,"ALIGNM31",&m31,
//                                                TFLOAT,"ALIGNM32",&m32,
//                                                TFLOAT,"ALIGNM33",&m33);
//   
//   if(flag)  return (EXIT_FAILURE);
//   
//   det_xcen = det_xpix1+ (det_xsiz-1)/2.0;
//   det_ycen = det_ypix1 + (det_ysiz-1)/2.0;
//              
//   return (EXIT_SUCCESS);
//   
// }


int RawToSat (double *deltaX , double *deltaY , int nsamples , char *teldeffile)
{
    TelDef teldef ;
    teldef.read (teldeffile) ;

    for (int i = 0 ; i < nsamples ; i ++)
    {
        float xint = teldef.coef_x0_a + teldef.coef_x0_b * deltaX[i] + teldef.coef_x0_c * deltaY[i] ;
        float yint = teldef.coef_y0_a + teldef.coef_y0_b * deltaX[i] + teldef.coef_y0_c * deltaY[i] ;


        float detx = teldef.det_xcen + teldef.detxflip + (xint - teldef.int_xcen - teldef.det_xoff) / teldef.det_scal ;
        float dety = teldef.focal_len ;
        float detz = teldef.det_ycen + teldef.detyflip + (yint - teldef.int_ycen - teldef.det_yoff) / teldef.det_scal ;
        //         cout<<teldef.det_xcen<<" "<<teldef.detxflip<<" "<<teldef.int_xcen<<" "<<teldef.det_xoff<<" "<<teldef.det_scal<<endl;
        //         cout<<teldef.det_ycen<<" "<<teldef.detyflip<<" "<<teldef.int_ycen<<" "<<teldef.det_yoff<<" "<<teldef.focal_len<<endl;
        //         exit(1);

        //assuming detector plane in X-Z plane of satellite
        deltaX[i] = teldef.m11 * detx + teldef.m12 * dety + teldef.m13*detz ;
        deltaY[i] = teldef.m31 * detx + teldef.m32 * dety + teldef.m33*detz ;

    }

    return (EXIT_SUCCESS) ;
}


int transformNUVtoVIS (double &Xnuv , double &Ynuv , double &thetanuv , double &Xvis , double &Yvis , double &thetavis)
{

    spMatrix NUV_mat (2 , 1) ;
    spMatrix VIS_mat (2 , 1) ;

    spMatrix NUV_inverse (2 , 2) ;
    NUV_mat (0 , 0) = Xnuv ;
    NUV_mat (1 , 0) = Ynuv ;

    if (! Flag_NUVtoVIS)
    {
        spMatrix A (2 , 2) ;
        spMatrix B (2 , 2) ;
        //cout<<"inside"<<endl;
        A (0 , 0) = VIS_TO_NUV_COEFF_1 ;
        A (0 , 1) = VIS_TO_NUV_COEFF_2 ;
        A (1 , 0) = VIS_TO_NUV_COEFF_3 ;
        A (1 , 1) = VIS_TO_NUV_COEFF_4 ;
        Flag_NUVtoVIS = TRUE ;
        A = A.Inverse () ;

        nuvTovis_inverse1 = A (0 , 0) ;
        nuvTovis_inverse2 = A (0 , 1) ;
        nuvTovis_inverse3 = A (1 , 0) ;
        nuvTovis_inverse4 = A (1 , 1) ;
    }
    NUV_inverse (0 , 0) = nuvTovis_inverse1 ;
    NUV_inverse (0 , 1) = nuvTovis_inverse2 ;
    NUV_inverse (1 , 0) = nuvTovis_inverse3 ;
    NUV_inverse (1 , 1) = nuvTovis_inverse4 ;

    VIS_mat = NUV_inverse*NUV_mat ;
    Xvis = VIS_mat (0 , 0) ;
    Yvis = VIS_mat (1 , 0) ;
    thetavis = - thetanuv ;

    return (EXIT_SUCCESS) ;
}
//to be changed


int transformVIStoNUV (double &Xvis , double &Yvis , double &thetavis , double &Xnuv , double &Ynuv , float &thetanuv)
{
    Xnuv = VIS_TO_NUV_COEFF_1 * Xvis + VIS_TO_NUV_COEFF_2* Yvis ;
    Ynuv = VIS_TO_NUV_COEFF_3 * Xvis + VIS_TO_NUV_COEFF_4 *Yvis ;
    thetanuv = - thetavis ;
    return (EXIT_SUCCESS) ;
}


int transformFUVtoVIS (double &Xfuv , double &Yfuv , double &thetafuv , double &Xvis , double &Yvis , double &thetavis)
{
    spMatrix FUV_mat (2 , 1) ;
    spMatrix VIS_mat (2 , 1) ;

    spMatrix FUV_inverse (2 , 2) ;
    FUV_mat (0 , 0) = Xfuv ;
    FUV_mat (1 , 0) = Yfuv ;

    if (! Flag_FUVtoVIS)
    {
        spMatrix A (2 , 2) ;
        A (0 , 0) = VIS_TO_FUV_COEFF_1 ;
        A (0 , 1) = VIS_TO_FUV_COEFF_2 ;
        A (1 , 0) = VIS_TO_FUV_COEFF_3 ;
        A (1 , 1) = VIS_TO_FUV_COEFF_4 ;
        Flag_FUVtoVIS = TRUE ;
        A = A.Inverse () ;
        fuvTovis_inverse1 = A (0 , 0) ;
        fuvTovis_inverse2 = A (0 , 1) ;
        fuvTovis_inverse3 = A (1 , 0) ;
        fuvTovis_inverse4 = A (1 , 1) ;
    }

    FUV_inverse (0 , 0) = fuvTovis_inverse1 ;
    FUV_inverse (0 , 1) = fuvTovis_inverse2 ;
    FUV_inverse (1 , 0) = fuvTovis_inverse3 ;
    FUV_inverse (1 , 1) = fuvTovis_inverse4 ;


    VIS_mat = FUV_inverse*FUV_mat ;
    Xvis = VIS_mat (0 , 0) ;
    Yvis = VIS_mat (1 , 0) ;
    thetavis = thetafuv ;

    return (EXIT_SUCCESS) ;
}


int transformVIStoFUV (double &Xvis , double &Yvis , double &thetavis , double &Xfuv , double &Yfuv , double &thetafuv)
{

    Xfuv = VIS_TO_FUV_COEFF_1 * Xvis + VIS_TO_FUV_COEFF_2 * Yvis ;
    Yfuv = VIS_TO_FUV_COEFF_3 * Xvis + VIS_TO_FUV_COEFF_4*Yvis ;
    thetafuv = thetavis ;

    return (EXIT_SUCCESS) ;
}


int transformFUVtoNUV (double &Xfuv , double &Yfuv , double &thetafuv , double &Xnuv , double &Ynuv , double &thetanuv)
{

    Xnuv = FUV_TO_NUV_COEFF_1 * Xfuv + FUV_TO_NUV_COEFF_2 * Yfuv ;
    Ynuv = FUV_TO_NUV_COEFF_3 * Xfuv + FUV_TO_NUV_COEFF_4*Yfuv ;
    thetanuv = thetafuv ;
    return (EXIT_SUCCESS) ;
}


int transformNUVtoFUV (double &Xnuv , double &Ynuv , double &thetanuv , double &Xfuv , double &Yfuv , double &thetafuv)
{
    spMatrix NUV_mat (2 , 1) ;
    spMatrix FUV_mat (2 , 1) ;

    spMatrix inverse_mat (2 , 2) ;
    NUV_mat (0 , 0) = Xnuv ;
    NUV_mat (1 , 0) = Ynuv ;
    if (! Flag_NUVtoFUV)
    {
        spMatrix A (2 , 2) ;
        A (0 , 0) = FUV_TO_NUV_COEFF_1 ;
        A (0 , 1) = FUV_TO_NUV_COEFF_2 ;
        A (1 , 0) = FUV_TO_NUV_COEFF_3 ;
        A (1 , 1) = FUV_TO_NUV_COEFF_4 ;
        Flag_NUVtoFUV = TRUE ;
        A = A.Inverse () ;
        nuvTofuv_inverse1 = A (0 , 0) ;
        nuvTofuv_inverse2 = A (0 , 1) ;
        nuvTofuv_inverse3 = A (1 , 0) ;
        nuvTofuv_inverse4 = A (1 , 1) ;
    }
    
    inverse_mat (0 , 0) = nuvTofuv_inverse1 ;
    inverse_mat (0 , 1) = nuvTofuv_inverse2 ;
    inverse_mat (1 , 0) = nuvTofuv_inverse3 ;
    inverse_mat (1 , 1) = nuvTofuv_inverse4 ;


    FUV_mat = inverse_mat*NUV_mat ;
    Xfuv = FUV_mat (0 , 0) ;
    Yfuv = FUV_mat (1 , 0) ;
    thetafuv = thetanuv ;

    return (EXIT_SUCCESS) ;
}


int transformRPYtoDXDYDTHETA_VIS (double &roll , double &pitch , double &yaw , double &dx , double &dy , double &dtheta)
{
    dx = RPYTOXYTHETA_COEFF_1_VIS * yaw + RPYTOXYTHETA_COEFF_2_VIS*pitch ;
    dy = RPYTOXYTHETA_COEFF_3_VIS * yaw + RPYTOXYTHETA_COEFF_4_VIS*pitch ;
    //dtheta = roll * (ROLL_TERM - roll) ;
    dtheta = roll ;
    cout<<dx<<" "<<dy<<" "<<dtheta<<endl;
    return (EXIT_SUCCESS) ;
}


int transformDXDYDTHETAtoRPY_VIS (double &dx , double &dy , double &dtheta , double &roll , double &pitch , double &yaw)
{

    spMatrix Dterm_mat (2 , 1) ;
    spMatrix RPY_mat (2 , 1) ;

    spMatrix inverse_vis (2 , 2) ;
    Dterm_mat (0 , 0) = dx ;
    Dterm_mat (1 , 0) = dy ;

    if (! Flag_dtermtoRPY_vis)
    {
        spMatrix A (2 , 2) ;
        A (0 , 0) = RPYTOXYTHETA_COEFF_1_VIS ;
        A (0 , 1) = RPYTOXYTHETA_COEFF_2_VIS ;
        A (1 , 0) = RPYTOXYTHETA_COEFF_3_VIS ;
        A (1 , 1) = RPYTOXYTHETA_COEFF_4_VIS ;
        Flag_dtermtoRPY_vis = TRUE ;
        A = A.Inverse () ;
        dtermToRPY_VIS_inverse1 = A (0 , 0) ;
        dtermToRPY_VIS_inverse2 = A (0 , 1) ;
        dtermToRPY_VIS_inverse3 = A (1 , 0) ;
        dtermToRPY_VIS_inverse4 = A (1 , 1) ;
       //cout<<A;exit(1);
    }
    inverse_vis (0 , 0) = dtermToRPY_VIS_inverse1 ;
    inverse_vis (0 , 1) = dtermToRPY_VIS_inverse2 ;
    inverse_vis (1 , 0) = dtermToRPY_VIS_inverse3 ;
    inverse_vis (1 , 1) = dtermToRPY_VIS_inverse4 ;

    RPY_mat = inverse_vis*Dterm_mat ;
    yaw = RPY_mat (0 , 0) ;
    pitch = RPY_mat (1 , 0) ;
    roll = dtheta ;

   return (EXIT_SUCCESS) ;
}


int transformRPYtoDXDYDTHETA_NUV (double &roll , double &pitch , double &yaw , double &dx , double &dy , double &dtheta)
{
    dx = RPYTOXYTHETA_COEFF_1_NUV * yaw + RPYTOXYTHETA_COEFF_2_NUV*pitch ;
    dy = RPYTOXYTHETA_COEFF_3_NUV * yaw + RPYTOXYTHETA_COEFF_4_NUV*pitch ;
  //  dtheta[i]=roll *(ROLL_TERM -roll);
    dtheta = roll ;
    return (EXIT_SUCCESS) ;
}


int transformDXDYDTHETAtoRPY_NUV (double &dx , double &dy , double &dtheta , double &roll , double &pitch , double &yaw)
{

    spMatrix Dterm_mat (2 , 1) ;
    spMatrix RPY_mat (2 , 1) ;

    spMatrix inverse_nuv (2 , 2) ;
    Dterm_mat (0 , 0) = dx ;
    Dterm_mat (1 , 0) = dy ;

    if (! Flag_dtermtoRPY_nuv)
    {
        spMatrix A (2 , 2) ;
        A (0 , 0) = RPYTOXYTHETA_COEFF_1_NUV ;
        A (0 , 1) = RPYTOXYTHETA_COEFF_2_NUV ;
        A (1 , 0) = RPYTOXYTHETA_COEFF_3_NUV ;
        A (1 , 1) = RPYTOXYTHETA_COEFF_4_NUV ;
        Flag_dtermtoRPY_nuv = TRUE ;
        A = A.Inverse () ;
        dtermToRPY_NUV_inverse1 = A (0 , 0) ;
        dtermToRPY_NUV_inverse2 = A (0 , 1) ;
        dtermToRPY_NUV_inverse3 = A (1 , 0) ;
        dtermToRPY_NUV_inverse4 = A (1 , 1) ;
        
    }
    inverse_nuv (0 , 0) = dtermToRPY_NUV_inverse1 ;
    inverse_nuv (0 , 1) = dtermToRPY_NUV_inverse2 ;
    inverse_nuv (1 , 0) = dtermToRPY_NUV_inverse3 ;
    inverse_nuv (1 , 1) = dtermToRPY_NUV_inverse4 ;

    RPY_mat = inverse_nuv*Dterm_mat ;
    yaw = RPY_mat (0 , 0) ;
    pitch = RPY_mat (1 , 0) ;
    roll = dtheta ;
    return (EXIT_SUCCESS) ;
}


int transformRPYtoDXDYDTHETA_FUV (double &roll , double &pitch , double &yaw , double &dx , double &dy , double &dtheta)
{
    dx = RPYTOXYTHETA_COEFF_1_FUV * yaw + RPYTOXYTHETA_COEFF_2_FUV*pitch ;
    dy = RPYTOXYTHETA_COEFF_3_FUV * yaw + RPYTOXYTHETA_COEFF_4_FUV*pitch ;
//     dx = RPYTOXYTHETA_COEFF_1_FUV *pitch + RPYTOXYTHETA_COEFF_2_FUV*yaw ;
//    dy = RPYTOXYTHETA_COEFF_3_FUV * pitch + RPYTOXYTHETA_COEFF_4_FUV*yaw ;
    //dtheta[i]=roll[i] *(ROLL_TERM -roll[i]);
//--------------------------------------------------------------------
//  dtheta = roll ;   sign flipped in the line below (01-Oct-2016)
    dtheta = -roll ;
//-------------------------------------------------------------------------
    return (EXIT_SUCCESS) ;
}


int transformDXDYDTHETAtoRPY_FUV (double &dx , double &dy , double &dtheta , double &roll , double &pitch , double &yaw)
{
    spMatrix Dterm_mat (2 , 1) ;
    spMatrix RPY_mat (2 , 1) ;

    spMatrix inverse_fuv (2 , 2) ;
    Dterm_mat (0 , 0) = dx ;
    Dterm_mat (1 , 0) = dy ;
    if (! Flag_dtermtoRPY_fuv)
    {
        spMatrix A (2 , 2) ;
        A (0 , 0) = RPYTOXYTHETA_COEFF_1_FUV ;
        A (0 , 1) = RPYTOXYTHETA_COEFF_2_FUV ;
        A (1 , 0) = RPYTOXYTHETA_COEFF_3_FUV ;
        A (1 , 1) = RPYTOXYTHETA_COEFF_4_FUV ;
        Flag_dtermtoRPY_fuv = TRUE ;
        A = A.Inverse () ;
        dtermToRPY_FUV_inverse1 = A (0 , 0) ;
        dtermToRPY_FUV_inverse2 = A (0 , 1) ;
        dtermToRPY_FUV_inverse3 = A (1 , 0) ;
        dtermToRPY_FUV_inverse4 = A (1 , 1) ;
    }

    inverse_fuv (0 , 0) = dtermToRPY_FUV_inverse1 ;
    inverse_fuv (0 , 1) = dtermToRPY_FUV_inverse2 ;
    inverse_fuv (1 , 0) = dtermToRPY_FUV_inverse3 ;
    inverse_fuv (1 , 1) = dtermToRPY_FUV_inverse4 ;

    RPY_mat = inverse_fuv*Dterm_mat ;
    yaw = RPY_mat (0 , 0) ;
    pitch = RPY_mat (1 , 0) ;
    roll = dtheta ;
    //yaw = dtermToRPY_FUV_inverse1 * dx + dtermToRPY_FUV_inverse2*dy ;
    //pitch = dtermToRPY_FUV_inverse2 * dx + dtermToRPY_FUV_inverse1*dy ;
    //roll =dtheta;
    return (EXIT_SUCCESS) ;
}
//  int transformDXDYDTHETAtoPY_VIS(double &xvis,double  &yvis ,double  &thetavis,double &pitch ,double &yaw)
//  {
//          yaw =  SHIFTS_TO_PY_VIS_COEFF_1*xvis + SHIFTS_TO_PY_VIS_COEFF_2*yvis;
//          pitch =   SHIFTS_TO_PY_VIS_COEFF_2*xvis- SHIFTS_TO_PY_VIS_COEFF_1*yvis; 
//     
//      return(EXIT_SUCCESS);
//  }
//  int transformdxdydthetatoPYNUV(float *xnuv,float *ynuv ,float *thetanuv,float* pitch ,float *yaw, int nsample)
//  {      
//    for(int i=0;i<nsample;i++)
//      {
//        yaw[i] =  SHIFTS_TO_PY_NUV_COEFF_ONE*xnuv[i] -SHIFTS_TO_PY_NUV_COEFF_TWO*ynuv[i];
//        pitch[i] =   SHIFTS_TO_PY_NUV_COEFF_TWO*xnuv[i]- SHIFTS_TO_PY_NUV_COEFF_ONE*ynuv[i]; 
//      }    
//          
//      
//      return(EXIT_SUCCESS);
//  }
//  int transformdxdydthetatoPYFUV(float *xfuv,float *yfuv ,float *thetafuv,float* pitch ,float *yaw, int nsample)
//  {      
//    for(int i=0;i<nsample;i++)
//      {
//        yaw[i] =  SHIFTS_TO_PY_FUV_COEFF_ONE*xfuv[i] -SHIFTS_TO_PY_FUV_COEFF_TWO*yfuv[i];
//        pitch[i] =   SHIFTS_TO_PY_FUV_COEFF_TWO*xfuv[i]- SHIFTS_TO_PY_FUV_COEFF_ONE*yfuv[i]; 
//      }       
//      
//      return(EXIT_SUCCESS);
//  }
