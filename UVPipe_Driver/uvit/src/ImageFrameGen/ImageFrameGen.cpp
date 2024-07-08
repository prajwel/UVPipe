/* 
 * File:   uvtImRa_commonArray.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <stdlib.h>

#include "ImageFrameGen.h"
#include<DataIngest.h>
//#include "uvtUnitConversion.h"
//#include "caldb_Handler.h"
//#include "uvtFilterBadpix.h"
//#include "uvtPixPadding.h"
//#include "uvtSubDivision.h"
//#include<uvtComputeDrift.h>
//#include "uvtDetectStar.h"
#include<Directory.h>
#include<DataInfo.h>
#include<DataIngest.h>
#include <algorithm>
//#include<uvtQEMCPCorr.h>
#include<set>
//#include<uvtComputeJitter.h>
//#include<uvtComputeThermal.h>
//#include<uvtRelAspCal.h>
//#include<uvtDetectDistCorr.h>
#include <vector>
#include<memory.h>
#include<spMatrix.h>
#include<glog/logging.h>
#include<uvtUtils.h>
#include<macro_def.h>
//#include<uvtComputeDrift.h>
#define  NBHD_RADIUS 5
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function
#define moduleoutdir_unit "uvtUnitConvertion"
#define moduleoutdir_badpix "uvtMaskBadPix"
#define moduleoutdir_fltfield "uvtFlatFieldCorrection"
#define moduleoutdir_qe "uvtQEMCPCorr"
#define moduleoutdir_pixpad "uvtPixPadding"
#define moduleoutdir_subdiv "uvtSubDivision"
#define moduleoutdir_cosmicray "uvtCosmicRayCorrection"
#define moduleoutdir_AccEverytsec "uvtAccEveryTsec"
#define moduleoutdir_findstarcentroid "uvtFindStarCentroid"
#define moduleoutdir_detectordistortion "uvtDetectDistortion"
#define moduleoutdir_opticaldistortion "uvtOpticalDistortion"
#define moduleoutdir_driftExercise "uvtComputeDrift"
#define moduleoutdir_refFrameCal "uvtRefFrameCal"
#define moduleoutdir_FrameIntegration "uvtFrameIntegration"
//#define IMAGE_ARRAYSIZE  4800
//#define subDivision_size  2400




uvtImRa_commonArray::uvtImRa_commonArray ()
{
//    sprintf (moduleoutdir_bp , "%s_%s" , moduleoutdir_badpix , VERSION) ;
//    sprintf (moduleoutdir_uc , "%s_%s" , moduleoutdir_unit , VERSION) ;
//    sprintf (moduleoutdir_ff , "%s_%s" , moduleoutdir_fltfield , VERSION) ;
//    sprintf (moduleoutdir_pp , "%s_%s" , moduleoutdir_pixpad , VERSION) ;
//    sprintf (moduleoutdir_sd , "%s_%s" , moduleoutdir_subdiv , VERSION) ;
//    sprintf (moduleoutdir_cr , "%s_%s" , moduleoutdir_cosmicray , VERSION) ;
//    sprintf (moduleoutdir_ac , "%s_%s" , moduleoutdir_AccEverytsec , VERSION) ;
//    sprintf (moduleoutdir_sc , "%s_%s" , moduleoutdir_findstarcentroid , VERSION) ;
//    sprintf (moduleoutdir_dd , "%s_%s" , moduleoutdir_detectordistortion , VERSION) ;
//    sprintf (moduleoutdir_od , "%s_%s" , moduleoutdir_opticaldistortion , VERSION) ;
//    sprintf (moduleoutdir_de , "%s_%s" , moduleoutdir_driftExercise , VERSION) ;
//    sprintf (moduleoutdir_rfc , "%s_%s" , moduleoutdir_refFrameCal , VERSION) ;
}

uvtImRa_commonArray::~uvtImRa_commonArray ()
{
    
//    delete[] frame_Data_subdivided ;
//    delete[] frame_ExpData_subdivided ;
//    delete[] frame_Data_Padded ;
//    delete[] frame_ExpoData_padded ;
//    delete[] frame_fc_data ;
    //delete[] frame_Data;
   // delete[] flatfielddata;
}

int uvtImRa_commonArray::read (int argc , char** argv)
{
    int status = 0 ;
    char temp[PIL_LINESIZE] ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG (INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("level1indir" , temp)))
    {
        LOG (INFO) << endl << "***Error reading input level1 directory name***" ;
        return status ;
    }
    level1indir.assign (temp) ;
    if (PIL_OK != (status = PILGetFname ("caldbdir" , temp)))
    {
        LOG (INFO) << endl << "***Error reading caldb directory name***" ;
        return status ;
    }
    caldbindir.assign (temp) ;
    if (PIL_OK != (status = PILGetFname ("level2outdir" , temp)))
    {
        LOG (INFO) << endl << "***Error reading output level2 directory name***" ;
        return status ;
    }
    level2outdir.assign (temp) ;

    if (PIL_OK != (status = PILGetString ("channel" , temp)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    channel.assign (temp);

    //dataingest param
    if (PIL_OK != (status = PILGetBool ("dropframe" , &dropframe)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }

  

    if (PIL_OK != (status = PILGetBool ("GTI_FLAG" , &gti_flag)))
    {
        LOG (INFO) << endl << "***Error Reading darkframeFlag:" <<gti_flag<< "***" ;
        return status ;
    }
    if (gti_flag)
    {
        all_or_custom = 0 ;
        valid_bit = 1 ;
    }
    

    if (PIL_OK != (status = PILGetInt ("Nacc" , &Nacc)))
    {
        LOG (INFO) << endl << "***Error Reading Nacc :" << Nacc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesDiscard" , &nFrameDiscard_fi)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesCompute" , &nFrameIntegrate_fi)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("frameintegrateDim" , &subDivision_size)))
    {
        LOG (INFO) << endl << "***Error Reading padding_dim :" << padding_dim << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("history" , &history)))
    {
        LOG (INFO) << endl << "***Error Reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    {
        LOG (INFO) << "***Error Reading mode parameter:" << history << "***" ;
        return status ;
    }  
   
    if (PIL_OK != (status = PILGetInt ("Write_todiskac" , &wtd_ac)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskfi" , &wtd_fi)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
        return status ;
    }
  
 
    PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtImRa_commonArray::uvtImRacomArrProcess ()
{
    int status = 0 ; //status flag to check return status of functions
    string nameOf_Tar,orbno;
    nameOf_Tar.assign (level1indir);
    cout<<"temp_dir "<<nameOf_Tar<<endl;
    level1indir="";
//     status= extractTars (nameOf_Tar,level1indir);
//   if(status){
//       LOG(INFO)<<"Error in extracting tar";
//       return(EXIT_FAILURE);
//   }
   status= extractTars (nameOf_Tar,level1indir,orbno);
   if(status){
       LOG(INFO)<<"Error in extracting tar";
       return(EXIT_FAILURE);
   }
   //ut<<"level1  tar ->"<<level1indir<<endl;
    //check Whether  level -1 directory Exist or not.
    if (!DirExists ((char *) level1indir.c_str ()))
    {
        LOG(ERROR) << endl << "Input level 1 directory not found " << endl ;
        return (EXIT_FAILURE) ;
    }

    //creating output directory if it does not exists
    if (!DirExists ((char *) level2outdir.c_str ()))
    {
        string cmd = "mkdir -p " + level2outdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists ((char*)level2outdir.c_str ()) && clobber==YES){
         string cmd = "rm  -r " + level2outdir ;
        system (cmd.c_str ()) ;
         cmd = "mkdir -p " + level2outdir ;
        system (cmd.c_str ()) ;
    }
   int temp_channel_id=-1;
   if(strcmp ((char*)channel.c_str (),"NUV")==0){
       temp_channel_id=NUV;
       
               
   }
   else if(strcmp ((char*)channel.c_str (),"FUV")==0){
       temp_channel_id=FUV;
   }
   else if(strcmp ((char*)channel.c_str (),"VIS")==0){
       temp_channel_id=VIS;
   }
   else{
       LOG(INFO)<<"***Invalid channel***";
       return(EXIT_FAILURE);
   }
    
   //LOG(INFO)<<"channel id"<<temp_channel_id;
    //function to get level 1 files
    Directory dirobj ;
    dirobj.setOrbitNo (orbno);
    /*---NOTE : Change 'NUV' to 'VIS'  for valid chain operation, presently using NUV due to non availability of data*/
    if (dirobj.setup (level1indir , level2outdir , temp_channel_id))
    {
        //IM RA  needs VIS data only //setup done only if uvitV directory found
        LOG(ERROR) << endl << "***Error in directory set up***" << endl ;
        LOG(ERROR) << endl << "***This chain runs for VIS data only...check for uvitV directory in input*** " << endl ;
        return (EXIT_FAILURE) ;
    }
    int numofsciencedatafiles = dirobj.sciencedatafile.size () ;//total number of science data file
    if (numofsciencedatafiles <= 0)
    {
        LOG(ERROR) << endl << "No science data file found" << endl ;
        return (EXIT_FAILURE) ;
    }

    char obsmode[FLEN_VALUE] ;
    char moduleIndir[FLEN_FILENAME] ; //hold path for input directory for every module, will be updated after every module run
    char outputdir[FLEN_FILENAME] ;
  //  int division_fact;

//    uvtDetectStar obj_sc ;//creating the object for finding peaks for the image
//incase of subdivision  not to be done
//   if (subdivisionFlag == 0)       subDivision_size = padding_dim ;
//   //Array for storing  pixel padding 
//    frame_Data_Padded = new float[padding_dim * padding_dim] ;
//    frame_ExpoData_padded = new float[padding_dim * padding_dim] ;
//    //Array for storing subdivided pixels incase of subdivision to be done
//    if(subdivisionFlag)
//    {
//    frame_Data_subdivided = new float[subDivision_size * subDivision_size] ;
//    frame_ExpData_subdivided = new float[subDivision_size * subDivision_size] ;
//       division_fact = subDivision_size / padding_dim ;
//        division_fact = division_fact*division_fact ;
//    }
     
    char errstr[500] ;
    unsigned short frameno ;
    double frametime , integrationtime ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
   
      //loop for number of science data file
   for (int dataindex = 0 ; dataindex < numofsciencedatafiles ; dataindex++)
    {
        name_track.clear ();
        LOG(ERROR) << endl << "------------------Data Set " << dataindex + 1 << " : " << dirobj.sciencedatafile[dataindex] << "-----------------------" << endl ;
        /*---finding mode from science data file---*/
        getKeywordVal ((char *) dirobj.sciencedatafile[dataindex].c_str () , "OBS_MODE" , 1 , obsmode) ;
        if (strcasecmp (obsmode , "IM") != 0 &&  strcasecmp (obsmode , "PC") != 0)
        {
            LOG(ERROR) << endl << "Observation mode is " << obsmode << "  in file  " << dirobj.sciencedatafile[dataindex] ;
            LOG(ERROR) << endl << "Checking next file...." << endl ;
            continue ; //go to next file if obs mode is not IM
        }

        LOG(INFO) << endl << "Data Mode is " << obsmode << endl ;

        strcpy (outputdir , dirobj.level2path[dataindex].c_str ()) ;
    //     if(wtd_bp==1)    sprintf (moduleoutdir_bp , "%s/%s_%s", outputdir, moduleoutdir_badpix,VERSION) ;
  //   if(wtd_uc==1)   sprintf (moduleoutdir_uc , "%s/%s_%s" , outputdir ,moduleoutdir_unit,VERSION) ;
   // if(flatfieldFlag==1 && wtd_ff==1)  sprintf (moduleoutdir_ff , "%s/%s_%s" ,outputdir, moduleoutdir_fltfield,VERSION ) ;
   // if(qemcpFlag==1 && wtd_qemcp==1)  sprintf (moduleoutdir_qemcp , "%s/%s_%s" ,outputdir, moduleoutdir_qe,VERSION ) ;
   // if(wtd_pp==1)  sprintf (moduleoutdir_pp , "%s/%s_%s" , outputdir ,moduleoutdir_pixpad,VERSION) ;
   // if(subdivisionFlag==1 && wtd_sd)  sprintf (moduleoutdir_sd , "%s/%s_%s" , outputdir , moduleoutdir_subdiv,VERSION) ;
  //  if(wtd_cr==1)sprintf (moduleoutdir_cr , "%s/%s_%s" , outputdir , moduleoutdir_cosmicray,VERSION) ;
    if (wtd_ac==1)sprintf (moduleoutdir_ac , "%s/%s_%s" , outputdir , moduleoutdir_AccEverytsec,VERSION) ;
      
 //   if(wtd_fsc==1)sprintf (moduleoutdir_sc , "%s/%s_%s" , outputdir , moduleoutdir_findstarcentroid,VERSION) ;
//    if(wtd_dd==1)sprintf (moduleoutdir_dd , "%s/%s_%s" ,outputdir ,moduleoutdir_detectordistortion,VERSION) ;
  //  if(wtd_od==1)sprintf (moduleoutdir_od , "%s/%s_%s" ,outputdir , moduleoutdir_opticaldistortion,VERSION) ;
   // sprintf (moduleoutdir_de , "%s/%s_%s" , outputdir , moduleoutdir_driftExercise,VERSION) ;
   // sprintf (moduleoutdir_rfc , "%s/%s_%s" , outputdir , moduleoutdir_refFrameCal,VERSION) ;
    //LOG(INFO)<<"DHDHH";
  //  strcpy();
        //----------DATAINGEST----------//
        LOG(ERROR) << endl << "===================================DATAINGEST===================================================" << endl ;
        //strcpy (lbtfile,(char *) dirobj.lbtfile.c_str ());
        DataIngest di_obj ;
        di_obj.read ((char *) dirobj.sciencedatafile[dataindex].c_str () ,(char*)caldbindir.c_str (),"" , "" ,"", "" , gti_flag , valid_bit , all_or_custom , outputdir , dropframe , parity_flag , clobber , history) ;
        di_obj.display () ;
        status = di_obj.DataIngestProcess () ;
        if (status)
        {
            LOG(ERROR) <<  "***Error in Data Ingest Process***"  ;
            return (EXIT_FAILURE) ;
        }
      
        strcpy (moduleIndir , di_obj.getModuleOutdir ()) ;
         strcpy (Indir_dataIngest , di_obj.getModuleOutdir ()) ;
       
        
     char infofile_in[PIL_LINESIZE] ;
    string  tempfilepath = searchFile (moduleIndir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << "***Information file not found in " << moduleIndir << "***" ;
        continue;
    }
    sprintf (infofile_in , "%s/%s" , moduleIndir , tempfilepath.c_str()) ;
    LOG(INFO) << "Information file :" << infofile_in;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "***Input FileList not Found at Specified PATH,Check Input Direcrory***" ;
        continue;
    }
    /*
 open the .info FITS file  and read the header information from the second HDU.
     */
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "") ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "") ;
    datainfo.getInfo(finfo_in) ; //reading basic information for data from information file
    //xsize = datainfo.getXsize () ;
    //ysize = datainfo.getYsize () ;
    xsize = IMG_DIM_DI ;
    ysize = IMG_DIM_DI ;
    int nframes ;
    char nameprefix[PIL_LINESIZE] ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
    if(datainfo.getModeFlag ()==IM)
    {
    fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "***Error in reading the  key value of the NFILES ***" , infofile_in) ;
   }
    else
    {
          fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "***Error in reading the key value of the EVTFILES ***" , infofile_in) ;
    }
    fits_read_key (finfo_in , TSTRING , "DARKDIR" , darkdir , NULL , &status) ;
    printError (status , "***Error in reading the  key value of the NFILES ***" , infofile_in) ;

    char **sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    //reading frame names from information file into vector
    if(datainfo.getModeFlag ()==IM){
    fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
    }
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in reading the column value of the Input Signal List" , infofile_in) ;
    

    if(datainfo.getModeFlag ()==IM)
    {
    if (wtd_ac == 1)//Accumulation
    {
        status = setDirectoryStructure (moduleoutdir_ac , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue;
        }
        
    }
    }
//    else if(datainfo.getModeFlag ()==PC)
//    {
//        if(wtd_fi==1){
//             status = setDirectoryStructure (moduleoutdir_ac , "SignalFrames") ;
//        if (status)
//        {
//            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
//            continue;
//        }
//        status = setDirectoryStructure (moduleoutdir_ac , "ExposureFrames") ;
//        if (status)
//        {
//            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
//            continue;
//        }
//        }
//        
//    }
    
    if(datainfo.getModeFlag ()==IM)
    {
        float *sum_ac = new float[IMG_DIM_DI* IMG_DIM_DI] ;
    float *sum_ac_exp = new float[IMG_DIM_DI* IMG_DIM_DI] ;
     initArray (sum_ac,IMG_DIM_DI* IMG_DIM_DI,(float)0.0);
     initArray (sum_ac_exp,IMG_DIM_DI* IMG_DIM_DI,(float)0.0);
 
     double avg_time=0.0;
     int count_frame_ac=0;
          fitsfile *fptr ;
    
  
        int Accu_fra_no = Nacc;
        
        int No_discard=frames_toDiscard;
        double t1 , t2 , x1 , x2 ,factor;
        if(nframes<=Accu_fra_no)
        {
            cout<<"Total number of frames to be accumulated are greater  than total available frames,Accumulation will be done on individual frame"<<endl;
            Accu_fra_no=1;            
        }
        if(nframes<=No_discard)
        {
             cout<<"Total number of frames to be discarded are greater  than total available frames,No  frames will be discarded  "<<endl;
            No_discard=0;
        }
         int tempx=0,tempy=0;
        for (int k =0 ; k <IMG_DIM_DI ; k++)
            {
             for(int l=0;l<IMG_DIM_DI;l++)
             { 
                
                 track_invalidPix.push_back (0);
                
             }
            }
    
  //loop for total number of frames.  
         char infile[NAMESIZE] ;
    for (int i = 0 ; i <nframes ; i++)
    { 
        
        
        sprintf (errstr , "Error at iteration number %d" , i) ;    
        sprintf (infile , "%s/%s/%s" , moduleIndir , "SignalFrames" , sigframelist[i]) ;
        
        frame_Data= (float*)malloc (xsize*ysize*sizeof(float));    initArray (frame_Data , xsize*ysize , -9999.0f) ;
        //frame_ExpData= (float*)malloc (xsize*ysize*sizeof(float)); initArray (frame_ExpData , xsize*ysize , -9999.0f) ;
         
        status = readImage (infile , 1 , frame_Data,xsize,ysize) ;
       //opening DataIngest Signal frame
        fits_open_file (&fptr , infile , READONLY , &status) ; printError (status , "Error in opening the input file" , infile) ;
        fits_movabs_hdu (fptr , 1 , NULL , &status) ; printError (status , "Error in  moving to the 2nd HDU of the out information file" , infile) ;
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ; printError (status , "Error in  reading the FRAMENO keyword" , infile) ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ; printError (status , "Error in  reading the FRAMETIME keyword" , infile) ;
        fits_read_key (fptr , TDOUBLE , "INT_TIME" , &integrationtime , NULL , &status) ; printError (status , "Error in  reading the INT_TIME keyword" , infile) ;
        fits_close_file (fptr , &status) ;printError (status , "Error in  closing the file" , infile) ;
       
        status = copyAllheaderKeys (infile) ;     
        
     
        
      
    tempx=0,tempy=0;
      
 
     for (int k =tempx= 0 ; k < IMG_DIM_DI ; k++)
            {
             for(int l=tempy=0;l<IMG_DIM_DI;l++)
             {
                 if(frame_Data[k*IMG_DIM_DI+l]!=INVALID_PIX_VALUE  )
                 {
                     //LOG(INFO)<<"INDIIFIFI";exit(1);
                   sum_ac[k*IMG_DIM_DI+l] = sum_ac[k*IMG_DIM_DI+l] + frame_Data[k*IMG_DIM_DI+l] ;
                   //sum_ac_exp[k*IMG_DIM_DI+l] = sum_ac_exp[k*IMG_DIM_DI+l] + frame_ExpData[k*IMG_DIM_DI+l] ;
                 }
                 else {                              
                     track_invalidPix[k*IMG_DIM_DI+l]= track_invalidPix[k*IMG_DIM_DI+l]+1;
                        }
             }
            }
     
   
         avg_time = avg_time + frametime ;      //Accumulate the time  
       
         //check for accumulation.
          if ((i + 1) % Accu_fra_no == 0 )
          {
//           for (int i =0;i<xsize*ysize;i++)
//         LOG(INFO)<<"Frame Data"<<sum_ac[i]<<endl;
              count_frame_ac++ ;
              LOG(INFO)<<"performing Accumulation on frame "<<count_frame_ac;
              
              
             status= performAccOFFrame (frame_Data,sum_ac,IMG_DIM_DI,IMG_DIM_DI,Accu_fra_no);
              if (status)
            {
                LOG (ERROR) << "Error in performing Accumulation" ;
               break;
            }
//              status= performAccOFFrame (frame_ExpData,sum_ac_exp,IMG_DIM_DI,IMG_DIM_DI,Accu_fra_no);
//              if (status)
//            {
//                LOG (ERROR) << "Error in performing Accumulation" ;
//                 break;
//            }
              track_invalidPix.clear ();
              
              for (int k =0 ; k <IMG_DIM_DI ; k++)
            {
             for(int l=0;l<IMG_DIM_DI;l++)  track_invalidPix.push_back (0);
           
            }
              //vect_cnt_duplicates.clear ();
              frameno = count_frame_ac ;           

              frametime = avg_time / Accu_fra_no ;//calculating frametime for Accumulated frame
              avg_time = 0.0 ; 
           
             // cout<<frametime<<endl;exit(1);
               if (wtd_ac == 1)
              {
                status = writeOutputImageToDisk ("ac" , moduleoutdir_ac , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , IMG_DIM_DI, IMG_DIM_DI) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
//                status = writeOutputImageToDisk ("ac" , moduleoutdir_ac , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , IMG_DIM_DI , IMG_DIM_DI) ;
//                if (status)
//                {
//                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
//                    break;
//                }
            }
                     
          }
          
          
         //releasing memory
     free(frame_Data);
     free(frame_ExpData);
    
      }
        if(flatfieldFlag){
         delete[] flatfielddata;
     }   
        
        
    }
    else if(datainfo.getModeFlag ()==PC)
   {
        long nrows;
        char file_in[FLEN_FILENAME];
        sprintf (file_in , "%s/%s" , moduleIndir , eventfile) ; //taking event file full path
        fitsfile *fevt_in , *fout ;
        //reading the event file information
        fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
        printError (status , "Error in opening the input event file" , file_in) ;
        status = copyAllheaderKeys (file_in) ;     
        fits_movabs_hdu (fevt_in , 2 , NULL , &status) ;
        printError (status , "Error in moving to the 2nd HDU" , eventfile) ;
        fits_get_num_rows (fevt_in , &nrows , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
          unsigned short *psc = new unsigned short[nrows] ;
        unsigned short *frame_no = new unsigned short[nrows] ;
        unsigned short *Ix = new unsigned short[nrows] ;
        unsigned short *Iy = new unsigned short[nrows] ;

        double *time_frame = new double[nrows] ;
        float *fx = new float[nrows] ;
        float *fy = new float[nrows] ;
        unsigned  char *dmm = new unsigned char[nrows] ;
        unsigned char *Min = new  unsigned char[nrows] ;
        unsigned short *multflag = new unsigned short [nrows] ;
        unsigned short *badflag = new unsigned short[nrows] ;

        fits_read_col (fevt_in , TUSHORT , 1 , 1 , 1 , nrows , NULL , psc , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TUSHORT , 2 , 1 , 1 , nrows , NULL , frame_no , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TDOUBLE , 3 , 1 , 1 , nrows , NULL , time_frame , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TUSHORT , 4 , 1 , 1 , nrows , NULL , Ix , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TFLOAT , 5 , 1 , 1 , nrows , NULL , fx , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TUSHORT , 6 , 1 , 1 , nrows , NULL , Iy , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TFLOAT , 7 , 1 , 1 , nrows , NULL , fy , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TBYTE , 8 , 1 , 1 , nrows , NULL , dmm , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fevt_in , TBYTE , 9 , 1 , 1 , nrows , NULL , Min , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_close_file (fevt_in , &status) ;
        float *newXfract= new float [nrows];
        float *newYfract = new float[nrows];
        float *effective_NumPhotons = new float[nrows];
        float *one_dim_img,*one_dim_exp;
     //   unsigned char *badflag= new char [nrows];
     //   unsigned char *multflag = new char [nrows];
        xsize=600;
        ysize=600;
        
       for(int i=0;i<nrows;i++)
        {
           badflag[i]=0.0f;
           multflag[i]=0.0f;     
           effective_NumPhotons[i]=0.0;
        }
       
        for (int i=0;i<nrows;i++)
        {
            if(Ix[i]+fx[i]+44<xsize  && Iy[i]+fy[i]+44<ysize){
            newXfract[i]=Ix[i]+fx[i]+44;
            
            newYfract[i]=Iy[i]+fy[i]+44;
            badflag[i]=1;
            multflag[i]=1;
            effective_NumPhotons[i]=1.0;
            }
        }
//        FrameIntegration_Arr obj_1;
//          obj_1.pixels_sig_fi= new float [subDivision_size*subDivision_size];
//        obj_1.pixels_exp_fi=new float[subDivision_size*subDivision_size];
        
       // FrameIntegration_Arr.pixels_sig_fi= new float [subDivision_size*subDivision_size];
       // FrameIntegration_Arr.pixels_exp_fi=new float[subDivision_size*subDivision_size];
        vector<FrameIntegration_Arr> frameIntegration_track ;
        //   vector<obj_1> frameIntegration_track ;
        //frameIntegration_track
       
        sprintf (moduleoutdir_fi , "%s/%s_%s" , outputdir , moduleoutdir_FrameIntegration , VERSION) ;
      
        
        status = performFrameIntegration (nrows , frame_no , xsize , nFrameDiscard_fi , nFrameIntegrate_fi , newXfract , newYfract , multflag , badflag , time_frame , effective_NumPhotons , one_dim_img , one_dim_exp , frameIntegration_track) ;
        if (wtd_fi)
        {
            status = setDirectoryStructure (moduleoutdir_fi , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***"  ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_fi , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
                continue ;
            }
            for (int i = 0 ; i < frameIntegration_track.size () ; i ++)
            {
                status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "SignalFrames" , "sig" , frameIntegration_track[i].pixels_sig_fi.data () , nameprefix , frameIntegration_track[i].frame_time_fi , frameIntegration_track[i].frame_number_fi , subDivision_size , subDivision_size) ; //this is for the SignalFrame output
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" ;
                    continue ;
                }
                status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "ExposureFrames" , "exp" , frameIntegration_track[i].pixels_exp_fi.data () , nameprefix , frameIntegration_track[i].frame_time_fi , frameIntegration_track[i].frame_number_fi , subDivision_size , subDivision_size) ; //this is for the SignalFrame output
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" ;
                    continue ;
                }

            }
        }
        
        
    }
 
    
       
     //creating output information file to the reference frame calculation.
        
        
        
   }
       
   return (EXIT_SUCCESS) ;
}

int uvtImRa_commonArray::setDirectoryStructure (char *Dir , const char *subdir)
{
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s/" , Dir , subdir) ;
    //LOG(INFO)<<dir<<endl;
    string cmd ;
    system (cmd.c_str ()) ;
    if (DirExists (dir) && clobber == YES)
    {
        LOG(ERROR) << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) dir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (dir) && clobber == NO)
    {
        LOG(ERROR) << endl << dir << "  already exists " ;
        LOG(ERROR) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;

    return (EXIT_SUCCESS) ;

}

int uvtImRa_commonArray::writeOutputImageToDisk (char *id , char *outDir , char *dir , char *subscript , float *Array , char *namepre , double ftime , unsigned short fno , int sizex , int sizey)
{
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char outfile[NAMESIZE] ;
    long naxes[2] ;
    naxes[0] = naxes[1] = sizex ;
    int naxis = 2 ;
    int bitpix = FLOAT_IMG ;
    fitsfile *fout ;
    sprintf (outfile , "%s/%s/%s_t%.4f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , id , subscript) ;
    int numhdu = 0 ;
    fits_create_file (&fout , outfile , &status) ;
    printError(status , "Error in creating the output Signal File" , outfile) ;
    fits_create_img (fout , bitpix , naxis , naxes , &status) ;
    printError(status , "Error in Creating the image for Signal Fie" , outfile) ;
    fits_write_pix (fout , TFLOAT , fpixel , sizex*sizey , Array , &status) ;
    printError(status , "***Error in writing the pixels to output***" , outfile) ;
    fits_get_num_hdus (fout , &numhdu , &status) ;
    char temp[1000] ;
    writeCommonKeywords (fout , outDir) ;
    for (int i = 0 ; i < numhdu ; i++)
    {
        for (int p = 0 ; p < header_info.size () ; p++)
        {
            sprintf (temp , "%s" , header_info[p].c_str ()) ;
            fits_write_record (fout , (char*) &temp , &status) ;
        }
    }
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the  output Signal fits file" , outfile) ;
    return (EXIT_SUCCESS) ;
}

int uvtImRa_commonArray::copyAllheaderKeys (char* infile)
{

    int status = 0 , keyexist ;
    char record[FLEN_CARD] ;
    header_info.clear () ;
    fitsfile *fin ;
    fits_open_file (&fin , infile , READONLY , &status) ;
    if (status) return (EXIT_FAILURE) ;
    fits_get_hdrspace (fin , &keyexist , NULL , &status) ;
    if (status)
    {
        LOG(INFO) << endl << "***Could not find number of keywords in file " << infile << "***" ;
        fits_report_error (stderr , status) ;
        return (EXIT_FAILURE) ;
    }
    int keyclass ;
    //LOG(INFO)<<"\nNumber of keywords found:"<<keyexist;
    for (int i = 1 ; i <= keyexist ; i++)
    {
        fits_read_record (fin , i , record , &status) ;
        if (status)
        {
            LOG(INFO) << endl << "***Error in reading record number " << i << " in file " << infile << "***" ;
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }
        keyclass = fits_get_keyclass (record) ;
        if (keyclass == TYP_COMM_KEY)
            continue ;
        else if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY)
        {
            header_info.push_back (record) ;
        }
    }
   
    fits_close_file (fin , &status) ;
    // delete[] record;
    return (EXIT_SUCCESS) ;
}



int uvtImRa_commonArray::performAccOFFrame(float *frmsigdata,float *sumdata,int sizex,int sizey,int numoffrmAcc)
{
   
     if (numoffrmAcc==0)
    {
        LOG(INFO) << "Divide by Zero" << endl ;
        return (EXIT_FAILURE) ;
    }

     int tempx=0,tempy=0;
     int num_times_repeat = 0 ;

         for (int k = tempx=0 ; k <IMG_DIM_DI ; k++)
         {
            
             for(int l=tempy=0;l<IMG_DIM_DI;l++)
             {
                 
                   //  num_times_repeat = count (track_invalidPix.begin () , track_invalidPix.end () , k*sizex+l) ; //STL method for count
                    frmsigdata[k*sizex+l]=INVALID_PIX_VALUE;
                    if(track_invalidPix[k*sizex+l] == numoffrmAcc) 
                    {
                        frmsigdata[k*sizex+l]=INVALID_PIX_VALUE;
                    }else
                    {
                 //   frmsigdata[k*sizex+l] = (float) (sumdata[k*sizex+l] / (numoffrmAcc- num_times_repeat)) ;
                  frmsigdata[k*sizex+l] = (float) (sumdata[k*sizex+l] / (numoffrmAcc- track_invalidPix[k*sizex+l])) ;
                //  LOG(INFO)<<sumdata[k*sizex+l]<<" "<<frmsigdata[k*sizex+l]<<endl;
                    }
                    sumdata[k*sizex+l]=0.0;
             }
         }
     //exit(1);

    return(EXIT_SUCCESS);
}

int uvtImRa_commonArray:: performFrameIntegration (long nrows , unsigned short*frame_no , int xsize , int Ndiscard , int  Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float  *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr> &vect)
{
    FrameIntegration_Arr obj ;
    int jj = 0 ;
    int x_tem , y_tem ;
    double Start_Time , End_Time , Avg_Time = 0.0 ;
    
    one_dim_img = new float[subDivision_size * subDivision_size] ;
    one_dim_exp = new float[subDivision_size * subDivision_size] ;
    for (int i=0;i<subDivision_size * subDivision_size;i++)
    {
     one_dim_exp[i]=1.0f;
     one_dim_img[i]=0.0f;
    }
    long firstpix1[2] ;
    firstpix1[0] = firstpix1[1] = 1 ;
    int naxis = 2 ;
    long naxes[2] = {subDivision_size , subDivision_size} ;
    float Eff_NoPhtn = 0.0f ;
    float multi_factor = subDivision_size/ xsize ;
  //  LOG(INFO)<<"mult flG IS "<<multi_factor;exit(1);
    int file_num = 1 ;
    int frame_check ;
    LOG (INFO) << "MULTIPLICATION FACTOR  for converting image to "<<subDivision_size<<" is " << multi_factor  ;
    int frame_init = frame_check = Nacc ;
    double time_frame = 0.0 , time_frame_final = 0.0 ;
    int cnt_frame = 0 ;
     //   vector<string> vhistorystr ;
    //if(history==YES)
    // getHistory (vhistorystr) ;
    /****provision for if frames are not starting from 1****/
    int start_row = 0 ;
    int cmpr_term = 0 ;
    //Ndiscard=frame_no[0]+Ndiscard-1;
    if (frame_no[0] != 1)
    {
        cmpr_term = Ndiscard + (frame_no[0] - 1) ;
    }
    for (int j = 0 ; j < nrows ; j ++)
    {
        if (frame_no[j] > cmpr_term)
        {
            start_row = j ;
            break ;
        }
    }
    long int tot_num_frames = 0 ;
    //cout<<Nacc<<" "<<Ndiscard<<endl;exit(1) ;
    LOG (INFO) << "\033[1;34mPerforming Frame Integration...\033[0m." ;
    long p = 0 ;
   // LOG(INFO)<<start_row<<" "<<Nacc + Ndiscard<<endl;exit(1);
    Nacc= Nacc+frame_no[0]-1;
    for (int k = start_row ; k < nrows ; k ++)
    {
      
        /**in 'IF condition' accumulated array of  Nacc frames will be  generate   & In 'ELSE condition' array will be written to the frame.
         x and y location  for number of events are taken as a indexes for the array .For each Nacc events one array will be written as a frame. 
         **/
        //LOG(INFO)<<frame_no[k]<<" "<<Nacc + Ndiscard<<endl;exit(1);
        if (frame_no[k] <= (Nacc))
        {
            

                        if (Nacc> frame_no[nrows - 1])
                        {
                            frame_init = frame_no[nrows - 1 ] % frame_check ;
                        }
            while (jj == 0)
            {
                Start_Time = t[k] ;
                jj ++ ;
            }

            x_tem = (int) xFrac[k] ;
            y_tem = (int) yFrac[k] ;
            x_tem = x_tem * multi_factor ;
            y_tem = y_tem * multi_factor ;


            if (y_tem < subDivision_size && x_tem < subDivision_size)
            {

                float mult = (float) mult_phn[k] ;
                float badflg = (float) bad_Flag[k] ;

                Eff_NoPhtn = (float) ENP[k] ;
                one_dim_exp[y_tem * subDivision_size+ x_tem] = (one_dim_exp[y_tem * subDivision_size + x_tem] + 1.0f) * mult*badflg ;

                one_dim_img[y_tem * subDivision_size + x_tem] = (one_dim_img[y_tem * subDivision_size + x_tem] + Eff_NoPhtn) * mult*badflg ;
               

            }

            if (frame_no[k - 1] == frame_no[k])
            {
                time_frame = t[k] ;
            }
            else
            {
                time_frame_final = time_frame_final + time_frame ;
            }
                        if (k == nrows - 1)
                        {
                            goto label_else ;
                        }
        }
        else//frame will be written here.
        {
           // LOG(INFO)<<"Last frame"<<k<<;
label_else:
            jj = 0 ;
            End_Time = t[k - 1] ;

            Avg_Time = Start_Time + (End_Time - Start_Time) / 2 ;

            
            Nacc = Nacc + frame_init ;
            tot_num_frames ++ ;
            LOG (INFO) << "Frame created " << tot_num_frames ;
            cnt_frame = 0 ;
            time_frame = 0.0 ;
            time_frame_final = 0.0 ;
            if (k != nrows - 1)
            {
                k -- ;
            }
            //obj.pixels_sig_fi = new float[subDivision_size*subDivision_size];
            //obj.pixels_exp_fi= new float[subDivision_size*subDivision_size];
            for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
            {
                obj.pixels_sig_fi.push_back (one_dim_img[i]);
                obj.pixels_exp_fi.push_back (one_dim_exp[i])     ;   
//                obj.pixels_sig_fi[i] = 0.0f ;
//                obj.pixels_exp_fi[i] = 0.0f ;
//                obj.pixels_sig_fi[i] = one_dim_img[i] ;
//                obj.pixels_exp_fi[i] = one_dim_exp[i] ;
                one_dim_img[i] = 0.0f ;
                one_dim_exp[i] = 1.0f ;

            }

            obj.frame_time_fi = Avg_Time ;
            obj.frame_number_fi = tot_num_frames ;
            vect.push_back (obj) ;
            obj.pixels_sig_fi.clear ();
            obj.pixels_exp_fi.clear ();
         //   delete[] obj.pixels_sig_fi,obj.pixels_exp_fi;
        }

    }
    delete[] one_dim_img,one_dim_exp;
    return (EXIT_SUCCESS) ;
}



