/* 
 * File:   uvtFullFrameAst.h
 * Author: uvit
 *
 * Created on December 4, 2013, 9:35 AM
 */

#ifndef UVTFULLFRAMEAST_H
#define	UVTFULLFRAMEAST_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<Attitude.h>
#include<caldb_Handler.h>
#include<glog/logging.h>
#include<macro_def.h>
#include<Database.h>
#define MODULENAME "uvtFullFrameAst"
#define caldbfinalsize 600 
typedef struct Star1
{
    float intensity ;
    float x , y ;
};
class uvtFullFrameAst{
    private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char RADECfileinput[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
        char attitudefile[PIL_LINESIZE];
        char shiftRotfile[PIL_LINESIZE];
        char expfile[PIL_LINESIZE];
        char noiceMapFile[PIL_LINESIZE];
        char catalogpath[PIL_LINESIZE];
        char att_timecol[NAMESIZE];
        char att_qcol[NAMESIZE];
        char databasename[NAMESIZE];
          int algo_flag,minimum_No_of_Stars;  
        int clobber;
        int history;
        char mode[10];            //pil mode
        float center_ra,center_dec,twist,center_ra_prev,center_dec_prev;
          float star_ra,star_dec;
        char imagefile_in[PIL_LINESIZE];
        char imagefile_out[PIL_LINESIZE];
        char imagefile_out_expName[PIL_LINESIZE];
        char imagefile_out_noicemapName[PIL_LINESIZE];
        char infofile_in[PIL_LINESIZE];   //full path
        char infofile_out[PIL_LINESIZE];
        caldb_Handler caldb_handler;
        DataInfo datainfo;
        float dayRefMJD;
        int xsize,ysize;
        float roll_angle;
        char nameprefix[FLEN_FILENAME]; 
        char moduleoutdir[PIL_LINESIZE];
          int centroid_Winsize,refine_Winsize;
         vector<int> Fx,Fy,Rx,Ry;
         vector<float> Fval,Rval;
         vector<double> Cx,Cy, Ci;
        Q uvitAlign;                //for UVIT alignment quaternion
         vector<Star1> star_track;
        TelDef teldef;
            float sd_mul_factor,sd_multi_factor_default;
       bool flag_NOT_FOUND_CATA_MATCH;
       bool flag_Optic_Catalogue;
       bool flag_UV_Catalogue;
       float center_RA_UV,center_DEC_UV,center_ROLL_UV;
       int numStars_frmCatamatch_optic,numStars_frmCatamatch_UV;
       float *bmag_rmag ;
       float *umag_rmag;
        float diff_ra_add,diff_dec_add;
        vector<float> track_raback,track_decback;
        float DX_UVIT,DYshift_UVIT,DTheta_UVIT;
        
         /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
       int numrows_lookupTable;
    int search_algo_ctlg;     
    char len_a[PIL_LINESIZE] ,len_b[PIL_LINESIZE],rad_search[PIL_LINESIZE];
    int getHistory(vector<string> &vhistory);
    int readcatalogueFile();             //to read catalog file
    int getAttitudeQ();
    int getRaDecTwist();
    int addWCS(fitsfile *fptr);
    int findStar_algo1 (float *inputArray, float *expdata ,float peakOFexp);
    void doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w);
     vector<int> valid_Star_index,validStar_index_filtered;
    
public:
    uvtFullFrameAst();
    ~uvtFullFrameAst();
    bool flag_Seconditeration;
    bool flag_out_From_First_Iteration;
    int read(int argc, char **argv);
    float getRAVAL(){
        return this->center_ra;
    }
    float getDECVAL(){
        return this->center_dec;
    }
    int flag_inputImage;
    double rapnt_att,decpnt_att;
    double rapnt_obs,decpnt_obs;
  //  int read(char *inputdatadir, char * caldbDir, char *outdir, char *attitudefile,
     //                                           char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int clobber, int history) ;
    
  //  int read(char *inputdatadir, char * caldbDir, char *outdir, char *attitudefile,
    //                                            char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char*len,char*wid,int clobber, int history);
//    int read(char *inputdatadir, char * caldbDir, char *shiftNRotfile,char *outdir, char *attitudefile,
//                                                char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char *len,char *wid, char *rad,int clobber, int history) ;
                                              //  char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char*len,char*wid, char* rad,int clobber, int history);
    
    int read(char *inputdatadir, char * caldbDir, char *shiftNRotfile, char *NoiceMapfile,char *exposureFile,char *outdir, char *attitudefile,
                                                char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char *len,char *wid, char *rad,int clobber, int history,int flag_Imageinput);
    int read(char *inputdatadir, char * caldbDir, char *shiftNRotfile, char *outdir, char *attitudefile,
                                                char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char *len,char *wid, char *rad,int clobber, int history,int flag_Imageinput);
    int read(char *inputdatadir, char* RADECfilename,char * caldbDir, char *shiftNRotfile, char *NoiceMapfile,char *exposureFile,char *outdir, char *attitudefile,
                                                char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char *len,char *wid, char *rad,int clobber, int history,int flag_Imageinput);
    void display();
    int uvtFullFrmAstProcess();
    int convertToUVITAxes (float &x ,float  &y,double &new_x,double &new_y,double &new_z);
    int toNormalize (double &x, double &y, double &z,double &nor_x,double &nor_y,double &nor_z);
    const char *getModuleOutdir() const {
        return moduleoutdir;
    } 
    //string getRaDECmatch(vector<float> &RA_img_Stars,vector<float> &DEC_img_Stars,int search_algo_ctlg,string len_a,string len_b,string rad_search,int no_of_stars,double &Max_Ra_value, double &Max_Dec_value);
//string getRaDECmatch(vector<float> &RA_img_Stars,vector<float> &DEC_img_Stars,int search_algo_ctlg,string len_a,string len_b,string rad_search,int no_of_stars,double &Max_Ra_value, double &Max_Dec_value ,string nameFile);    

//string getRaDECmatch(vector<float> &RA_img_Stars,vector<float> &DEC_img_Stars,int search_algo_ctlg,string len_a,string len_b,string rad_search,int no_of_stars,
  //      double &Max_Ra_value, double &Max_Dec_value ,string nameFile,int *numStarsmatch);
string getRaDECmatch(vector<float> &RA_img_Stars,vector<float> &DEC_img_Stars,int search_algo_ctlg,string len_a,string len_b,string rad_search,int no_of_stars,
        double &Abs_Val_DiffRA, double &Abs_Val_DiffDEC ,string nameFile,int *numStarsmatch,double decangle,bool flag_channel);
float diff_ra_add_opt,diff_dec_add_opt;
float getDiffRAval(){
        return this->diff_ra_add_opt;
    }
    float getDiffDECval(){
        return this->diff_dec_add_opt;
    }
    float getXshift(){
        return this->DX_UVIT;
    }
    float getYshift(){
        return this->DYshift_UVIT;
    }
    float getTheta(){
        return this->DTheta_UVIT;
    }
    
    
int calculateShiftsAndRoll( int totalelements,vector<double> &Xone,vector<double> &Yone,vector<double> &Xtwo,vector<double> &Ytwo,int xsize,int ysize,double &Xdx,double  &Ydy,double &Theta);
};


#endif	/* UVTFULLFRAMEAST_H */


