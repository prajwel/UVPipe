/* 
 * File:   Directory.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */


#include "Directory.h"

void  Directory::setOrbitNo(string orbno)
{
    orbit_no=orbno;
}
/*
 Function to set up files for level-1 and also create level2 directory structure*/

//int Directory::setup (string level1dir , string level2dir , int channel)
//{
//    string level2dir_temp;
//    struct dirent **namelist ;
//    int n = 0 ;
//    bool uvitdirflag = false ;
//    string uvitdir ;
//  
//    n = scandir (level1dir.c_str () , &namelist , 0 , alphasort) ;
//    while (n--)
//    {
//     
//        if (namelist[n]->d_type == DT_DIR)
//        {
//            if (strcmp (namelist[n]->d_name , (const char *) UVIT_DIR) == 0)                 //look for UVIT directory in level1 directory
//            {
//                uvitdir.assign (namelist[n]->d_name) ;
//                uvitdirflag == true ;
//                break ;
//            }
//            
//        }
//        free (namelist[n]) ;
//    }
//    
//    free (namelist) ;
//
//    level1dir.append ((string) FILEPATHSEPARATOR) ;
////    string temp = level1dir + uvitdir ;                  //creating path upto uvit directory inside level-1 directory
//    string temp = level1dir + uvitdir+FILEPATHSEPARATOR+orbit_no ;                  //creating path upto uvit directory inside level-1 directory
//    //cout<<"PPPP"<<temp<<endl;exit(1);
//    uvitdir.clear () ;
//    uvitdir = temp ;
//
//    LOG(INFO) << endl << "UVIT directory :" << uvitdir << endl ;
//
//    if (uvitdirflag)
//    {
//        LOG(ERROR) << endl << "uvit directory not found in level1 directory" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//
//    vector<string> filelist ;
//
//    string channeldir ;
//    n = 0 ;
//    if (channel == FUV)
//    {
//        n = scandir (uvitdir.c_str () , &namelist , uvitF_filter , alphasort) ;
//        if (n <= 0)
//        {
//            return (EXIT_FAILURE) ;
//        }
//        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) FUV_DIR ;
//    }
//    else if (channel == NUV)
//    {
//        n = scandir (uvitdir.c_str () , &namelist , uvitN_filter , alphasort) ;
//        if (n <= 0)
//        {
//            return (EXIT_FAILURE) ;
//        }
//        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) NUV_DIR ;
//    }
//    else if (channel == VIS)
//    {
//        n = scandir (uvitdir.c_str () , &namelist , uvitV_filter , alphasort) ;
//        if (n <= 0)
//        {
//            return (EXIT_FAILURE) ;
//        }
//        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) VIS_DIR ;
//    }
//    else
//    {
//        LOG(ERROR) << endl << "Invalid channel value" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//
//    for (int i = 0 ; i < n ; i++) free (namelist[i]) ;
//    free (namelist) ;
//
//    int pos = 0 ;
//    if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtN") == 0){
//         darkDirectory=channeldir + (string) "/" + DARKDIRN;
//    }
//    else if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtF") == 0){
//        darkDirectory=channeldir + (string) "/" + DARKDIRF;
//    }
//    else if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtV") == 0){
//        darkDirectory=channeldir + (string) "/" + DARKDIRV;
//    }
//    else{
//        
//    }
//  //darkDirectory=channeldir + (string) "/" + DARKDIR;
//  if(!DirExists ((char*)darkDirectory.c_str ())){
//      LOG(ERROR)<<"***Level1 directory not containing Dark directory ,Expected directory = ***"<<" "<<darkDirectory<<endl;
//      exit(1);
//  }
//  
////    int status = getDirs (channeldir , dirlist) ;
////    for (int i = 0 ; i < dirlist.size () ; i++)
////    {
////        pos = dirlist[i].find ((string)DARKDIR) ;
////        if (pos > 0 && pos < dirlist[i].size ()){ darkDirectory=dirlist[i] ; break;}
////    }
//  
//  int  status = getFiles (uvitdir , filelist) ;
//  if(status){
//       LOG(ERROR)<<endl<<"Error in reading files from input directory  "<<channeldir;
//       return (EXIT_FAILURE);
//   }
//   for (int i = 0 ; i < filelist.size () ; i++)
//    {
//        pos = filelist[i].find ((string) TILDE) ;
//        if (pos > 0 && pos < filelist[i].size ()) continue ;
//        
////        pos = filelist[i].find ((string) SCIENCEDATA_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) sciencedatafile.push_back (filelist[i]) ;
//
//        pos = filelist[i].find ((string) MKF_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) mkffile = filelist[i] ;
//
//        pos = filelist[i].find ((string) ORBIT_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) orbitfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) ATTITUDE_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) attfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) GTI_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) gtifile.push_back (filelist[i]) ;
//
//        pos = filelist[i].find ((string) BTI_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) btifile.push_back (filelist[i]) ;
//
//        pos = filelist[i].find ((string) LBT_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) lbtfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) HK_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) hkfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) TCT_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) tctfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) GYRO_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) gyrofile = filelist[i] ;
//        //pos = filelist[i].find ((string)DARKDIR) ;
//      
//    }
//
//   
//  filelist.clear ();
//         
//   status = getFiles (channeldir , filelist) ;
//  if(status){
//       LOG(ERROR)<<endl<<"Error in reading files from input directory  "<<channeldir;
//       return (EXIT_FAILURE);
//   }
//      
//    //setting all filepaths in variable for the structure Directory
//    for (int i = 0 ; i < filelist.size () ; i++)
//    {
//        pos = filelist[i].find ((string) TILDE) ;
//        if (pos > 0 && pos < filelist[i].size ()) continue ;
//        
//        pos = filelist[i].find ((string) SCIENCEDATA_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) sciencedatafile.push_back (filelist[i]) ;
//
////        pos = filelist[i].find ((string) MKF_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) mkffile = filelist[i] ;
////
////        pos = filelist[i].find ((string) ORBIT_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) orbitfile = filelist[i] ;
////
////        pos = filelist[i].find ((string) ATTITUDE_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) attfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) GTI_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) gtifile.push_back (filelist[i]) ;
//
//        pos = filelist[i].find ((string) BTI_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) btifile.push_back (filelist[i]) ;
//
////        pos = filelist[i].find ((string) LBT_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) lbtfile = filelist[i] ;
//
////        pos = filelist[i].find ((string) HK_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) hkfile = filelist[i] ;
////
////        pos = filelist[i].find ((string) TCT_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) tctfile = filelist[i] ;
////
////        pos = filelist[i].find ((string) GYRO_EXT) ;
////        if (pos > 0 && pos < filelist[i].size ()) gyrofile = filelist[i] ;
////        pos = filelist[i].find ((string)DARKDIR) ;
//      
//    }
//
//    string dirpath ;
//    int pos1 = 0 , pos2 = 0 ;
//    
//   // string tempuvitdir = level1dir + (string) FILEPATHSEPARATOR+(string) UVIT_DIR+ (string) FILEPATHSEPARATOR ;
//    //tempuvitdir.find ();
//   // LOG(INFO)<<sciencedatafile.size ()<<endl;exit(1);
//   // strcpy((char*)level2dir_temp.c_str (),(char*)level2dir.c_str ());
//    level2dir_temp=level2dir;
//    for (int i = 0 ; i < sciencedatafile.size () ; i++)
//    {
//        string tempuvitdir = (string) FILEPATHSEPARATOR + (string) UVIT_DIR + (string) FILEPATHSEPARATOR  ;
//        pos1 = sciencedatafile[i].rfind ((string)tempuvitdir ,sciencedatafile[i].size ()-1) ;
//        pos2 = sciencedatafile[i].find_last_of ((string) FILEPATHSEPARATOR) ;
//        //LOG(INFO)<<pos1<<" "<<pos2<<endl;exit(1);
//        if (pos1 < sciencedatafile[i].size () && pos1 >= 0 && pos2 >= 0 && pos2 < sciencedatafile[i].size ())
//        {
//            dirpath = sciencedatafile[i].substr (pos1 + 1 , pos2 - pos1) ;
//            level2dir = level2dir_temp + (string) (FILEPATHSEPARATOR) + dirpath ;
//            level2path.push_back (level2dir) ;
//            string cmd = "mkdir -p " + level2dir ;
//           
//            system (cmd.c_str ()) ;
//
//        }
//    }
//
//        //  string tempuvitdir = (string) FILEPATHSEPARATOR + (string) UVIT_DIR ;
//      
//   return (EXIT_SUCCESS) ;
//}

int Directory::setup (string level1dir , string level2dir , int channel)
{
    string level2dir_temp;
    struct dirent **namelist ;
    int n = 0 ;
    bool uvitdirflag = false ;
    string uvitdir ;
  
    n = scandir (level1dir.c_str () , &namelist , 0 , alphasort) ;
    while (n--)
    {
        LOG(INFO)<<(int)namelist[n]->d_type<<" "<<namelist[n]->d_name;
        //if (namelist[n]->d_type == DT_DIR)
        //{
            if (strcmp (namelist[n]->d_name , UVIT_DIR) == 0)                 //look for UVIT directory in level1 directory
            {
                uvitdir.assign (namelist[n]->d_name) ;
                uvitdirflag = true ;
                break ;
            }
            
        //}
        free (namelist[n]) ;
    }
    
    free (namelist) ;

    level1dir.append ((string) FILEPATHSEPARATOR) ;
    //string temp = level1dir + uvitdir ;                  //creating path upto uvit directory inside level-1 directory
   string temp = level1dir + uvitdir+FILEPATHSEPARATOR+orbit_no ;                  //creating path upto uvit directory inside level-1 directory
    //cout<<"PPPP"<<temp<<endl;exit(1);
    uvitdir.clear () ;
    uvitdir = temp ;

    LOG(INFO) << endl << "UVIT directory :" << uvitdir << endl ;
    
    if (!uvitdirflag)
    {
        LOG(ERROR) << endl << "uvit directory not found in level1 directory" << endl ;
        return (EXIT_FAILURE) ;
    }

    vector<string> filelist ;

    string channeldir ;
    n = 0 ;
    if (channel == FUV)
    {
        n = scandir (uvitdir.c_str () , &namelist , uvitF_filter , alphasort) ;
        if (n <= 0)
        {
            return (EXIT_FAILURE) ;
        }
        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) FUV_DIR ;
    }
    else if (channel == NUV)
    {
        n = scandir (uvitdir.c_str () , &namelist , uvitN_filter , alphasort) ;
        if (n <= 0)
        {
            return (EXIT_FAILURE) ;
        }
        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) NUV_DIR ;
    }
    else if (channel == VIS)
    {
        n = scandir (uvitdir.c_str () , &namelist , uvitV_filter , alphasort) ;
        if (n <= 0)
        {
            return (EXIT_FAILURE) ;
        }
        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) VIS_DIR ;
    }
    else
    {
        LOG(ERROR) << endl << "Invalid channel value" << endl ;
        return (EXIT_FAILURE) ;
    }

    for (int i = 0 ; i < n ; i++) free (namelist[i]) ;
    free (namelist) ;

    int pos = 0 ;
    if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtN") == 0){
         darkDirectory=channeldir + (string) "/" + DARKDIRN;
    }
    else if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtF") == 0){
        darkDirectory=channeldir + (string) "/" + DARKDIRF;
    }
    else if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtV") == 0){
        darkDirectory=channeldir + (string) "/" + DARKDIRV;
    }
    else{
        
    }
  //darkDirectory=channeldir + (string) "/" + DARKDIR;
  if(!DirExists ((char*)darkDirectory.c_str ())){
      LOG(ERROR)<<"***Level1 directory not containing Dark directory ,Expected directory = ***"<<" "<<darkDirectory<<endl;
      exit(1);
  }
  
//    int status = getDirs (channeldir , dirlist) ;
//    for (int i = 0 ; i < dirlist.size () ; i++)
//    {
//        pos = dirlist[i].find ((string)DARKDIR) ;
//        if (pos > 0 && pos < dirlist[i].size ()){ darkDirectory=dirlist[i] ; break;}
//    }
    //LOG(INFO)<<"uvit_dir "<<uvitdir;exit(1);
  int  status = getFiles (uvitdir , filelist) ;
  if(status){
       LOG(ERROR)<<endl<<"Error in reading files from input directory  "<<channeldir;
       return (EXIT_FAILURE);
   }
   for (int i = 0 ; i < filelist.size () ; i++)
    {
       LOG(INFO)<<filelist[i];
   }
   for (int i = 0 ; i < filelist.size () ; i++)
    {
        pos = filelist[i].find ((string) TILDE) ;
        if (pos > 0 && pos < filelist[i].size ()) continue ;
        
//        pos = filelist[i].find ((string) SCIENCEDATA_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) sciencedatafile.push_back (filelist[i]) ;

        pos = filelist[i].find ((string) MKF_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) mkffile = filelist[i] ;

        pos = filelist[i].find ((string) ORBIT_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) orbitfile = filelist[i] ;

        pos = filelist[i].find ((string) ATTITUDE_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) attfile = filelist[i] ;

        pos = filelist[i].find ((string) GTI_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) gtifile.push_back (filelist[i]) ;

        pos = filelist[i].find ((string) BTI_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) btifile.push_back (filelist[i]) ;

        pos = filelist[i].find ((string) LBT_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) lbtfile = filelist[i] ;

        pos = filelist[i].find ((string) HK_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) hkfile = filelist[i] ;

        pos = filelist[i].find ((string) TCT_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) tctfile = filelist[i] ;

        pos = filelist[i].find ((string) GYRO_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) gyrofile = filelist[i] ;
        //pos = filelist[i].find ((string)DARKDIR) ;
      
    }

   
  filelist.clear ();
         
   status = getFiles (channeldir , filelist) ;
  if(status){
       LOG(ERROR)<<endl<<"Error in reading files from input directory  "<<channeldir;
       return (EXIT_FAILURE);
   }
      
    //setting all filepaths in variable for the structure Directory
    for (int i = 0 ; i < filelist.size () ; i++)
    {
        pos = filelist[i].find ((string) TILDE) ;
        if (pos > 0 && pos < filelist[i].size ()) continue ;
        
        pos = filelist[i].find ((string) SCIENCEDATA_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) sciencedatafile.push_back (filelist[i]) ;

//        pos = filelist[i].find ((string) MKF_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) mkffile = filelist[i] ;
//
//        pos = filelist[i].find ((string) ORBIT_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) orbitfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) ATTITUDE_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) attfile = filelist[i] ;

        pos = filelist[i].find ((string) GTI_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) gtifile.push_back (filelist[i]) ;

        pos = filelist[i].find ((string) BTI_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) btifile.push_back (filelist[i]) ;

//        pos = filelist[i].find ((string) LBT_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) lbtfile = filelist[i] ;

//        pos = filelist[i].find ((string) HK_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) hkfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) TCT_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) tctfile = filelist[i] ;
//
//        pos = filelist[i].find ((string) GYRO_EXT) ;
//        if (pos > 0 && pos < filelist[i].size ()) gyrofile = filelist[i] ;
//        pos = filelist[i].find ((string)DARKDIR) ;
      
    }

    string dirpath ;
    int pos1 = 0 , pos2 = 0 ;
    
   // string tempuvitdir = level1dir + (string) FILEPATHSEPARATOR+(string) UVIT_DIR+ (string) FILEPATHSEPARATOR ;
    //tempuvitdir.find ();
   // LOG(INFO)<<sciencedatafile.size ()<<endl;exit(1);
   // strcpy((char*)level2dir_temp.c_str (),(char*)level2dir.c_str ());
    level2dir_temp=level2dir;
    for (int i = 0 ; i < sciencedatafile.size () ; i++)
    {
        string tempuvitdir = (string) FILEPATHSEPARATOR + (string) UVIT_DIR + (string) FILEPATHSEPARATOR  ;
        pos1 = sciencedatafile[i].rfind ((string)tempuvitdir ,sciencedatafile[i].size ()-1) ;
        pos2 = sciencedatafile[i].find_last_of ((string) FILEPATHSEPARATOR) ;
        //LOG(INFO)<<pos1<<" "<<pos2<<endl;exit(1);
        if (pos1 < sciencedatafile[i].size () && pos1 >= 0 && pos2 >= 0 && pos2 < sciencedatafile[i].size ())
        {
            dirpath = sciencedatafile[i].substr (pos1 + 1 , pos2 - pos1) ;
            level2dir = level2dir_temp + (string) (FILEPATHSEPARATOR) + dirpath ;
            level2path.push_back (level2dir) ;
            string cmd = "mkdir -p " + level2dir ;
           
            system (cmd.c_str ()) ;

        }
    }

        //  string tempuvitdir = (string) FILEPATHSEPARATOR + (string) UVIT_DIR ;
      
   return (EXIT_SUCCESS) ;
}

/* function for finding uvtF directory*/
int uvitF_filter (const struct dirent *dptr)
{
    int retval = 0 ;
    if (strcmp (dptr->d_name , (const char *) FUV_DIR) == 0) retval = 1 ;
    return retval ;
}

/* function for finding uvtN directory*/
int uvitN_filter (const struct dirent *dptr)
{
    int retval = 0 ;
    if (strcmp (dptr->d_name , (const char *) NUV_DIR) == 0) retval = 1 ;
    return retval ;
}

/* function for finding uvtV directory*/
int uvitV_filter (const struct dirent *dptr)
{
    int retval = 0 ;
    if (strcmp (dptr->d_name , (const char *) VIS_DIR) == 0) retval = 1 ;
    return retval ;
}

/*Function to get all the files in vector filelist from a channel directory uvitF, uvitN or UVITV*/

int getFiles (string dirName , vector<string> &filelist)
{
    bool darkFlag=FALSE;
  
    struct dirent **namelist ;  
    int num = scandir (dirName.c_str () , &namelist , defaultfilter , alphasort) ;
    if (num <= 0)
    {
        //LOG(ERROR)<<endl<<"No files found in directory "<<dirName<<endl;
        return (EXIT_FAILURE) ;
    }

    bool fileflag = true ;
    struct stat filestat;
    while (num--)
    {
        
        string tempname = dirName + (string) "/" + (string) namelist[num]->d_name ;
        if (stat( tempname.c_str(), &filestat )) continue;
        if (S_ISDIR( filestat.st_mode ))      
        {
        
            fileflag = false ;
            getFiles (tempname , filelist) ;
         }
        else
        {
        filelist.push_back (tempname) ;
        }
//        if (namelist[num]->d_type == DT_REG)
//        {
//            filelist.push_back (tempname) ;
//        }
//        else if (namelist[num]->d_type == DT_DIR)
//        {
//            
//         
//            fileflag = false ;
//            getFiles (tempname , filelist) ;
//        }
    }


    if (fileflag == true)
    {
        return (EXIT_SUCCESS) ;
    }

    return (EXIT_SUCCESS) ;
}


int getDirs (string dirName , vector<string> &subDirs)
{
   
    struct dirent **namelist ;
    int num = scandir (dirName.c_str () , &namelist , defaultfilter , alphasort) ;
    if (num <= 0)
    {
        //LOG(ERROR)<<endl<<"No files found in directory "<<dirName<<endl;
        return (EXIT_FAILURE) ;
    }

    bool fileflag = true ;
    while (num--)
    {
        string tempname = dirName + (string) "/" + (string) namelist[num]->d_name ;
        if (namelist[num]->d_type == DT_REG)
        {
            
        }
        else if (namelist[num]->d_type == DT_DIR)
        {   fileflag = false ;
            subDirs.push_back (tempname) ;
            getDirs (tempname,subDirs);
        }
    }
      if (fileflag == true)
    {
        return (EXIT_SUCCESS) ;
    }

    return (EXIT_SUCCESS) ;
}

int Directory::CopyAuxDir(string level1dir , string level2dir , int channel)
{
    
   struct dirent **namelist ;
   int n = 0 ;
   bool uvitdirflag = false ;

    string uvitdir ;

    n = scandir (level1dir.c_str () , &namelist , 0 , alphasort) ;
    while (n--)
    {
       
        if (namelist[n]->d_type == DT_DIR)
        {
            if (strcmp (namelist[n]->d_name , (const char *) UVIT_DIR) == 0)                 //look for UVIT directory in level1 directory
            {
                uvitdir.assign (namelist[n]->d_name) ;
                uvitdirflag == true ;
                break ;
            }
        }
        free (namelist[n]) ;
    }
    
    free (namelist) ;

    level1dir.append ((string) FILEPATHSEPARATOR) ;
    string temp = level1dir + uvitdir ;                  //creating path upto uvit directory inside level-1 directory
    uvitdir.clear () ;
    uvitdir = temp ;

    LOG(INFO) << endl << "UVIT directory :" << uvitdir << endl ;

    if (uvitdirflag)
    {
        LOG(ERROR) << endl << "uvit directory not found in level1 directory" << endl ;
        return (EXIT_FAILURE) ;
    }
  string    cmd ="cp  -r "+uvitdir+"  "+level2dir ;
            system (cmd.c_str ()) ;

//    vector<string> dirlist ;
//
//    string channeldir ;
//    n = 0 ;
//    if (channel == FUV)
//    {
//        n = scandir (uvitdir.c_str () , &namelist , uvitF_filter , alphasort) ;
//        if (n <= 0)
//        {
//            return (EXIT_FAILURE) ;
//        }
//        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) FUV_DIR ;
//    }
//    else if (channel == NUV)
//    {
//        n = scandir (uvitdir.c_str () , &namelist , uvitN_filter , alphasort) ;
//        if (n <= 0)
//        {
//            return (EXIT_FAILURE) ;
//        }
//        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) NUV_DIR ;
//    }
//    else if (channel == VIS)
//    {
//        n = scandir (uvitdir.c_str () , &namelist , uvitV_filter , alphasort) ;
//        if (n <= 0)
//        {
//            return (EXIT_FAILURE) ;
//        }
//        channeldir = uvitdir + (string) FILEPATHSEPARATOR + (string) VIS_DIR ;
//    }
//    else
//    {
//        LOG(ERROR) << endl << "Invalid channel value" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//
//    for (int i = 0 ; i < n ; i++) free (namelist[i]) ;
//    free (namelist) ;
//
// int status = getDirs (uvitdir , dirlist) ;
// if(status)
// {
//       LOG(ERROR)<<endl<<"Error in reading files from input directory  "<<channeldir;
//       return (EXIT_FAILURE);
// }
//      
//    //setting all filepaths in variable for the structure Directory
//    vector<string> name_dir;
//    string dirpath,dirpath1 ;
//    int pos1 = 0 , pos2 = 0 ;
//    char tempname2[FLEN_FILENAME];
//    string  cmd;
//    for (int i = 0 ; i < dirlist.size () ; i++)
//    {
//        string tempuvitdir = (string) FILEPATHSEPARATOR + (string) UVIT_DIR +(string) FILEPATHSEPARATOR;
//        pos1 = dirlist[i].rfind ((string)tempuvitdir , dirlist[i].size ()-1) ;
//        pos2 = dirlist[i].find_last_of ((string) FILEPATHSEPARATOR) ;
//       // LOG(INFO)<<"888"<<dirlist[i]<<endl;
//        if (pos1 < dirlist[i].size () && pos1 >= 0 && pos2 >= 0 && pos2 < dirlist[i].size () )
//        {
//            dirpath = dirlist[i].substr (pos1 + 1 , pos2 - pos1) ;
//            dirpath1 = level2dir + (string) (FILEPATHSEPARATOR) + dirpath ;
//                
//            
//            if(strstr(basename(dirlist[i].c_str ()),".")!=NULL && strcmp(basename(dirlist[i].c_str()),".")!=0)
//            { 
//               sprintf(tempname2,"%s/%s",dirpath1.c_str (),basename(dirlist[i].c_str ()));                  
//               name_dir.push_back (tempname2);
//            }
//            else
//            {
//                 cmd ="cp  -r "+dirlist[i]+"  "+dirpath1 ;
//           
//            system (cmd.c_str ()) ;
//                
//            }
//         }       
//            
//    }
   // exit(1);
//    for(int i=0;i<name_dir.size ();i++)
//    {
//        
//            cmd ="mkdir  -p  "+(string)name_dir[i] ;
//            LOG(INFO)<<"Executing command .."<<cmd<<endl;
//            system (cmd.c_str ()) ;
//        
//    }
    
   
   return (EXIT_SUCCESS) ;
}
int readScienceDataFilesPath (string level1dir  , int channel,string orbit_no,string darkDirectory,vector<string> &sciencedatafile)
{
    string level2dir_temp;
    struct dirent **namelist ;
    int n = 0 ;
    bool uvitdirflag = false ;
    string uvitdir ;
  
    n = scandir (level1dir.c_str () , &namelist , 0 , alphasort) ;
    while (n--)
    {
        if (strcmp (namelist[n]->d_name , UVIT_DIR) == 0)                 //look for UVIT directory in level1 directory
            {
                uvitdir.assign (namelist[n]->d_name) ;
                uvitdirflag = true ;
                break ;
            }
            
        //}
        free (namelist[n]) ;
    }
    
    free (namelist) ;

    level1dir.append ((string) "/") ;
    //string temp = level1dir + uvitdir ;                  //creating path upto uvit directory inside level-1 directory
   string temp = level1dir + uvitdir+"/"+orbit_no ;                  //creating path upto uvit directory inside level-1 directory
    //cout<<"PPPP"<<temp<<endl;exit(1);
    uvitdir.clear () ;
    uvitdir = temp ;

   //cout << endl << "UVIT directory :" << uvitdir << endl ;
    
    if (!uvitdirflag)
    {
        cerr << endl << "uvit directory not found in level1 directory" << endl ;
        return (EXIT_FAILURE) ;
    }

    vector<string> filelist ;

    string channeldir ;
    n = 0 ;
    if (channel == 1)
    {
        n = scandir (uvitdir.c_str () , &namelist , uvitF_filter , alphasort) ;
        if (n <= 0)
        {
            return (EXIT_FAILURE) ;
        }
        channeldir = uvitdir + (string) "/" + (string) FUV_DIR ;
    }
    else if (channel == 2)
    {
        n = scandir (uvitdir.c_str () , &namelist , uvitN_filter , alphasort) ;
        if (n <= 0)
        {
            return (EXIT_FAILURE) ;
        }
        channeldir = uvitdir + (string) "/" + (string) NUV_DIR ;
    }
    else if (channel == 3)
    {
        n = scandir (uvitdir.c_str () , &namelist , uvitV_filter , alphasort) ;
        if (n <= 0)
        {
            return (EXIT_FAILURE) ;
        }
        channeldir = uvitdir + (string) "/" + (string) VIS_DIR ;
    }
    else
    {
        cerr << endl << "Invalid channel value" << endl ;
        return (EXIT_FAILURE) ;
    }

    for (int i = 0 ; i < n ; i++) free (namelist[i]) ;
    free (namelist) ;

    int pos = 0 ;
    if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtN") == 0){
         darkDirectory=channeldir + (string) "/" + DARKDIRN;
    }
    else if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtF") == 0){
        darkDirectory=channeldir + (string) "/" + DARKDIRF;
    }
    else if(strcmp (basename(channeldir.c_str ()) , (const char *)"uvtV") == 0){
        darkDirectory=channeldir + (string) "/" + DARKDIRV;
    }
    else{
        
    }
  //darkDirectory=channeldir + (string) "/" + DARKDIR;
  if(!DirExists ((char*)darkDirectory.c_str ())){
      cerr<<"***Level1 directory not containing Dark directory ,Expected directory = ***"<<" "<<darkDirectory<<endl;
      exit(1);
  }
  
  int  status = getFiles (uvitdir , filelist) ;
  if(status){
       cerr<<endl<<"Error in reading files from input directory  "<<channeldir;
       return (EXIT_FAILURE);
   }
      
  filelist.clear ();
         
   status = getFiles (channeldir , filelist) ;
  if(status){
       cerr<<endl<<"Error in reading files from input directory  "<<channeldir;
       return (EXIT_FAILURE);
   }
      
    //setting all filepaths in variable for the structure Directory
    for (int i = 0 ; i < filelist.size () ; i++)
    {
        pos = filelist[i].find ((string) TILDE) ;
        if (pos > 0 && pos < filelist[i].size ()) continue ;
        
        pos = filelist[i].find ((string) SCIENCEDATA_EXT) ;
        if (pos > 0 && pos < filelist[i].size ()) sciencedatafile.push_back (filelist[i]) ;


      
    }   
      
   return (EXIT_SUCCESS) ;
}
int get_firstLevelSubDirs (string dirName , vector<string> &subDirs)
{
   
    struct dirent **namelist ;
    int num = scandir (dirName.c_str () , &namelist , defaultfilter , alphasort) ;
    if (num <= 0)
    {
        //LOG(ERROR)<<endl<<"No files found in directory "<<dirName<<endl;
        return (EXIT_FAILURE) ;
    }

    bool fileflag = true ;
    struct stat filestat;
    while (num--)
    {
        string tempname = dirName + (string) "/" + (string) namelist[num]->d_name ;
        if (stat( tempname.c_str(), &filestat )) continue;
       if (S_ISDIR( filestat.st_mode ))
        {   fileflag = false ;
            subDirs.push_back (tempname) ;
            
        }
       else
        {
          LOG(INFO)  <<"FILE";
        }
        
    }
//      if (fileflag == true)
//    {
//        return (EXIT_SUCCESS) ;
//    }

    return (EXIT_SUCCESS) ;
}