
/* 
 * File:   uvtUtils.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include<uvtUtils.h>
#include <glog/logging.h>
#include<macro_def.h>
#include<spMatrix.h>
#include<algorithm>

#include "uvtUtils.h"
//--------------------------------------------------------------------


string ReturnYorN(int num){
    if (num==1){
        return "y";
    }
    else return "n";
    
}
bool FileExists (char *filename)
{
    //LOG(ERROR)<<"\n----------FileExists----------\n";
    struct stat filestat ;
    int filest ;
    filestat.st_mode = -1 ;
    bool exist = false ;
    filest = stat (filename , &filestat) ;
    if (S_ISREG (filestat.st_mode))
    {
        LOG(ERROR)<<endl<<filename<<endl;
        exist = true ;
    }
    else exist = false ;

    //LOG(ERROR)<<"\n----------FileExists End----------\n";
    return exist ;
}

//---------------------------------------------------------------------

bool DirExists (char *dirname)
{
    struct stat dirstat ;
    int t ;
    dirstat.st_mode = -1 ;
    t = stat (dirname , &dirstat) ;
    if (S_ISDIR (dirstat.st_mode))
    {
        return true ;
    }
    return false ;
}
//---------------------------------------------------------------------

int copyKeywords (fitsfile *fptr1 , fitsfile *fptr2 , int n , ...)
{
    //LOG(ERROR)<<"\n----------copyKeywords----------\n";
    va_list keywords ;
    va_start (keywords , n) ;

    int status = 0 ;
    char record[FLEN_CARD] , *key ;

    for (int i = 0 ; i < n ; i++)
    {
        key = va_arg (keywords , char *) ;

        fits_read_card (fptr1 , key , record , &status) ;
        if (status == 202)
        {
            LOG(ERROR) << endl << "***" << key << " not found in header***\n" ;
        }
        if (status)
        {
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }

        fits_write_record (fptr2 , record , &status) ;
        if (status)
        {
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }

    }
    va_end (keywords) ;
    return (EXIT_SUCCESS) ;
}
//------------------------------------------------------------------------------

int copyUserKeywords (fitsfile *fin , fitsfile *fout)
{
    int status = 0 , keyexist ;
    char record[FLEN_CARD] ;
    char str[3] ;
    char filename[FLEN_FILENAME] ;
    fits_file_name (fin , filename , &status) ;
    if (status) return (EXIT_FAILURE) ;
    fits_get_hdrspace (fin , &keyexist , NULL , &status) ;
    if (status)
    {
        LOG(ERROR)<< endl << "***Could not find number of keywords in file " << filename << "***" ;
        fits_report_error (stderr , status) ;
        return (EXIT_FAILURE) ;
    }
    int keyclass ;
    //LOG(ERROR)<<"\nNumber of keywords found:"<<keyexist;
    for (int i = 1 ; i <= keyexist ; i++)
    {
        fits_read_record (fin , i , record , &status) ;
        if (status)
        {
            LOG(ERROR) << endl << "***Error in reading record number " << i << " in file " << filename << "***" ;
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }
        keyclass = fits_get_keyclass (record) ;
//        if (keyclass == TYP_COMM_KEY)
//            continue ;
         if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY || keyclass ==TYP_COMM_KEY)
        {
            fits_write_record (fout , record , &status) ;
            if (status)
            {
                LOG(ERROR) << endl << "***Error in writing keyword " << record << "***" ;
                fits_report_error (stderr , status) ;
                return (EXIT_FAILURE) ;
            }
        }
    }
    return (EXIT_SUCCESS) ;
}

//-----------------------------------------------------------------------------

void joinStrings (char *outstring , int n , ...)
{
    va_list strings ;
    va_start (strings , n) ;
    char *str ;
    int len = 0 ;
    for (int i = 0 ; i < n ; i++)
    {
        str = va_arg (strings , char *) ;
        len = len + strlen (str) ;
        if (i == 0) strcpy (outstring , str) ;
        else strcat (outstring , str) ;
    }
    outstring[len] = '\0' ;
    va_end (strings) ;
}
//-------------------------------------------------------------------------

int deleteDir (char *dir)
{
    vector<string> directories ;
    //LOG(ERROR)<<"\nDir:"<<dir;
    struct dirent **namelist ;
    int n , n_copy ;
    char temp[2048] ;
    n = scandir (dir , &namelist , defaultfilter , alphasort) ;
    n_copy = n ;
    //LOG(ERROR)<<endl<<"n:"<<n;
    if (n < 0)
        return (EXIT_FAILURE) ;
    else if (n == 0)
    {
        LOG(INFO) << endl << "1.Removing " << dir ;
        if (remove (dir))
        {
            LOG(ERROR) << endl << "***1. Error in removing " << dir << "***" << endl ;
            return (EXIT_FAILURE) ;
        }
    }
    else
    {
        while (n--)
        {
            strcpy (temp , dir) ;
            strcat (temp , "/") ;
            strcat (temp , namelist[n]->d_name) ;
            //LOG(ERROR)<<endl<<temp;
            if (namelist[n]->d_type == DT_REG)
            {
                //LOG(ERROR)<<endl<<"Deleting "<<temp;
                if (unlink (temp))
                {
                    LOG(ERROR) << endl << "***Error in removing " << temp << "***" << endl ;
                    return (EXIT_FAILURE) ;
                }
            }
                /*else if((namelist[n]->d_type==DT_DIR) && 
                        (strcmp(namelist[n]->d_name,".")!=0) && 
                        (strcmp(namelist[n]->d_name,"..")!=0)){*/
            else if (namelist[n]->d_type == DT_DIR)
                deleteDir (temp) ;
            //LOG(ERROR)<<endl<<"Namelist N: "<<namelist[n]->d_name;

            free (namelist[n]) ;
        }
        free (namelist) ;
    }

    if (DirExists (dir))
    {
        LOG(INFO)<< endl << ".Removing " << dir ;
        if (remove (dir))
        {
            LOG(ERROR) << endl << "***Error in removing " << dir << "***" << endl ;
            return (EXIT_FAILURE) ;
        }
    }
    return (EXIT_SUCCESS) ;
}
//------------------------------------------------------------------------------

int defaultfilter (const struct dirent *dptr)
{
    int retval = 1 ;
    if (strcmp (dptr->d_name , ".") == 0) retval = 0 ;
    if (strcmp (dptr->d_name , "..") == 0) retval = 0 ;
    return retval ;
}
//------------------------------------------------------------------------------

void printError (int status , char *msg)
{

    if (status)
    {
        fits_report_error (stderr , status) ;
        LOG(ERROR) << endl << msg ;
        exit (status) ;
    }
    return ;
}

void printError (int status , char *errmsg , char *fitsname)
{

    if (status)
    {
        fits_report_error (stderr , status) ;
        LOG(ERROR)<< "***" << errmsg <<endl;
        LOG(ERROR)<<"*** File Name-" << fitsname <<" ***" ;
        LOG(ERROR)<< endl ;
        exit (status) ;
    }
    return ;
}

int parseString (string str , char delim , vector<string> &substr)
{
    int pos = 0 ;
    string temp ;
    int len = str.size () ;
    if (len <= 0)
    {
        LOG(ERROR) << "\n***Invalid input string" << str << "***\n" ;
        return (-1) ;
    }
    //LOG(ERROR)<<"\nString Passed is "<<str;
    //LOG(ERROR)<<"\nLength:"<<len;
    while (len > 0)
    {
        //LOG(ERROR)<<endl<<"Length:"<<len;
        pos = str.find (delim , 0) ;
        //LOG(ERROR)<<endl<<"pos:"<<pos;
        if (pos < 0)
        {
            substr.push_back (str) ;
            break ;
        }
        else
        {
            temp = str.substr (0 , pos) ;
            str = str.substr (pos + 1) ;
            substr.push_back (temp) ;
            len = len - pos ;
            //LOG(ERROR)<<"\nLength:"<<len;
        }
    }
    return (substr.size ()) ;
}
//------------------------------------------------------------------------------

int getbins (char *bins , double *start , double *end , int nebins)
{
    //LOG(ERROR)<<"\nInside getbins\n";
    int i = 0 ;
    int temp ;
    if (strcmp (bins , "-") == 0)
    {
        return 2 ;
    }
    else
    {
        string ebinranges = (string) bins ;
        //string sets[MAX_BINS];
        string temp , temp2 ;
        //LOG(ERROR)<<endl<<ebinranges<<endl;
        i = 0 ;
        int pos = 0 ;
        for (i = 0 ; i < nebins ; i++)
        {
            //LOG(ERROR)<<endl<<ebinranges;
            pos = ebinranges.find ('-' , 0) ;
            temp = ebinranges.substr (0 , pos) ;
            temp2 = ebinranges.substr (pos + 1) ;
            ebinranges.clear () ;
            ebinranges.assign (temp2) ;
            start[i] = atof (temp.c_str ()) ;
            temp.clear () ;
            temp2.clear () ;

            pos = ebinranges.find (',' , pos) ;
            temp = ebinranges.substr (0 , pos) ;
            temp2 = ebinranges.substr (pos + 1) ;
            ebinranges.clear () ;
            ebinranges.assign (temp2) ;
            end[i] = atof (temp.c_str ()) ;
            temp.clear () ;
            temp2.clear () ;
        }
    }
    return 0 ;
}
//------------------------------------------------------------------------------

int get_nbins (char *binstr)
{
    int nbins = 0 ;
    if (strcmp (binstr , "-") == 0)
    {
        nbins = 1 ;
    }
    else
    {
        char bins[2000] ;
        strcpy (bins , binstr) ;
        char *token = strtok (bins , ",") ;
        while (token != NULL)
        {
            token = strtok (NULL , ",") ;
            nbins++ ;
        }
    }
    return nbins ;
}

//--------------Quaternion related-----------

/*
 * 
 */

//Member Functions of struct Q

Q::Q ()
{
    q1 = 0 ;
    q2 = 0 ;
    q3 = 0 ;
    q4 = 0 ;
    norm = false ;
    AxisAngle = false ;
    qflag = false ;
    theta = 0 ;
    x = 0 ;
    y = 0 ;
    z = 0 ;
}

Q::Q (double a , double b , double c , double d)
{
    q1 = a ;
    q2 = b ;
    q3 = c ;
    q4 = d ;
    norm = false ;
    AxisAngle = false ;
}

Q::Q (const Q &q)
{
    //LOG(ERROR)<<endl<<"\nInside COPY CONSTRUCTOR\n";
    q1 = q.q1 ;
    q2 = q.q2 ;
    q3 = q.q3 ;
    q4 = q.q4 ;
    theta = q.theta ;
    x = q.x ;
    y = q.y ;
    z = q.z ;
    norm = false ;
    AxisAngle = false ;
}

void Q::readQ ()
{
    LOG(ERROR) << "\nEnter q1,q2,q3,q4 [q1 is scalar]:" ;
    cin >> q1 >> q2 >> q3 >> q4 ;
}

void Q::readAxisAngle ()
{
    LOG(ERROR) << "\nEnter theta, x, y and z:" ;
    cin >> theta >> x >> y >> z ;
}

int Q::getAxisAngle ()
{
    //LOG(ERROR)<<"\nInside getAxisAngle\n";
    if (AxisAngle == false)
    {
        double denominator = sqrt (1 - q1 * q1) ;
        theta = 2 * acos (q1) ;
        //theta=ceil(theta-0.5);  //rounding the theta angle
        if (denominator != 0)
        {
            x = q2 / denominator ;
            y = q3 / denominator ;
            z = q4 / denominator ;
        }
        else
        {
            x = 1 ;
            y = 0 ;
            z = 0 ;
        }
    }
    AxisAngle = true ;
    return 0 ;
}

int Q::getQuat ()
{
    if (qflag == false)
    {
        q1 = cos (theta / 2) ;
        q2 = x * sin (theta / 2) ;
        q3 = y * sin (theta / 2) ;
        q4 = z * sin (theta / 2) ;
        qflag = true ;
    }

}

void Q::display ()
{
    LOG(INFO) << "(" << q1 << ")+(" << q2 << ")i+(" <<
            q3 << ")j+(" << q4 << ")k" ;
    if (AxisAngle == true)
    {
        LOG(INFO) << "Angle=" << theta * (180 / M_PI) << " deg\t" ;
        LOG(INFO) << "Axis: x=" << x << "   y=" << y << "  z=" << z ;
    }
    //LOG(ERROR)<<endl;
}

void Q::normalize ()
{

    if (norm == false)
    {
        mod = sqrt (q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4) ;
        //LOG(ERROR)<<"Denominator:"<<denominator;
        if (mod != 0)
        {
            q1 = q1 / mod ;
            q2 = q2 / mod ;
            q3 = q3 / mod ;
            q4 = q4 / mod ;
            norm = true ;
        }
    }
}
//-------------------------------------------------------------------
//Member functions of axis

void Axis::normalize ()
{
    mod = sqrt (x * x + y * y + z * z) ;
    if (mod != 0)
    {
        x = x / mod ;
        y = y / mod ;
        z = z / mod ;
    }
    norm = true ;
}

void Axis::display ()
{
    //LOG(INFO) << "(" << x << ")i + (" << y << ")j + (" << z << ")k" ;
}

double Axis::getMod ()
{
    return (sqrt (x * x + y * y + z * z)) ;
}

void Axis::update(double a, double b, double c){
      x=a; y=b; z=c; norm=false;
}
//------------------------------------------------------------------------------

void matrix_product (float **A , float **B , float **C , int a , int b , int c)
{
    int i , j , k ;
    float sum ;
    for (i = 0 ; i < a ; i++)
    {
        for (j = 0 ; j < c ; j++)
        {
            sum = 0 ;
            for (k = 0 ; k < b ; k++)
            {
                sum = sum + (A[i][k] * B[k][j]) ;
                C[i][j] = sum ;
            }
        }
    }
}
//------------------------------------------------------------------------------

/**
 * Function to get product of two quaternions
 * @param Q1
 * @param Q2
 * @param Q3
 */
void quaternion_product (Q &Q1 , Q &Q2 , Q &Q3)
{
    Q3.q1 = Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3 - Q1.q4 * Q2.q4 ;
    Q3.q2 = Q1.q1 * Q2.q2 + Q1.q2 * Q2.q1 + Q1.q3 * Q2.q4 - Q1.q4 * Q2.q3 ;
    Q3.q3 = Q1.q1 * Q2.q3 - Q1.q2 * Q2.q4 + Q1.q3 * Q2.q1 + Q1.q4 * Q2.q2 ;
    Q3.q4 = Q1.q1 * Q2.q4 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2 + Q1.q4 * Q2.q1 ;
}
//------------------------------------------------------------------------------

Q &Inverse (Q &q)
{
    //LOG(ERROR)<<"\nInside inverse\n";
    Q inverseq ;
    inverseq.q1 = q.q1 ;
    inverseq.q2 = (-1) * q.q2 ;
    inverseq.q3 = (-1) * q.q3 ;
    inverseq.q4 = (-1) * q.q4 ;
    return inverseq ;
}

//---------------------------------------------------------------------------------

/**
 * Function to rotate one axis from one frame to another using quaternion
 * @param axis 
 * @param q      
 * @param axisrot 
 * @return 0
 */
int rotate (Axis &axis , Q &q , Axis &axisrot)
{
    //LOG(ERROR)<<"\nInside Rotation\n";
    Q qinv ;
    qinv = Inverse (q) ;
    Q temp ;
    Q v1 (0 , axis.x , axis.y , axis.z) , v2 ;
    //LOG(ERROR)<<"input vector:";  v1.display();
    quaternion_product (v1 , qinv , temp) ;
    quaternion_product (q , temp , v2) ;
    axisrot.x = v2.q2 ;
    axisrot.y = v2.q3 ;
    axisrot.z = v2.q4 ;
    //LOG(ERROR)<<"\noutput vector";  axisrot.display();

    return 0 ;
}
int compute_Rotmatrix(double &q1,double &q2,double &q3,double &q4,double &theeta,double &phi,double &psi )
{
  spMatrix RQ(3,3);
  RQ(0,0) = q1*q1 - q2*q2 - q3*q3 +q4*q4;
  RQ(0,1) = 2*(q1*q2+q3*q4);
  RQ(0,2) = 2*(q1*q3-q2*q4);
  RQ(1,0) = 2*(q1*q2-q3*q4);
  RQ(1,1) = -q1*q1 + q2*q2 - q3*q3 +q4*q4;
  RQ(1,2) = 2*(q2*q3+q1*q4);
  RQ(2,0) = 2*(q1*q3+q2*q4);
  RQ(2,1) = 2*(q2*q3-q1*q4);
  RQ(2,2) = -q1*q1 - q2*q2 + q3*q3 +q4*q4;
  
  // Roll.
   theeta = asin(- RQ(0,2));
// Pitch.
   phi=atan(RQ(0,1)/RQ(0,0));
// Yaw.
   psi=atan(RQ(1,2)/RQ(2,2));
   
//  LOG(INFO)<<theeta*180/M_PI<<" "<<phi<<" "<<psi;
  return (EXIT_SUCCESS);
  
}
//-----------------------------------------------------------------------------------------

int writeHistory (char *filename , vector<string> &vhistory)
{
   
    int status = 0 , numhdu ;
    fitsfile *fptr ;
    fits_open_file (&fptr , filename , READWRITE , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in opening file-writeHistory()\n" ;
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    fits_get_num_hdus (fptr , &numhdu , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in getting HDU number-writeHistory()\n" ;
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    for (int i = 1 ; i <= numhdu ; i++)
    {
        fits_movabs_hdu (fptr , i , NULL , &status) ;
        if (status)
        {
            LOG(ERROR) << "***Error moving in file " << filename << "-writeHistory()***" ;
            fits_report_error (stderr , status) ;
            return -1 ;
        }
        for (int j = 0 ; j < vhistory.size () ; j++)
        {
            fits_write_history (fptr , (char *) vhistory[j].c_str () , &status) ;
            if (status)
            {
                LOG(ERROR) << "\nError writing History-writeHistory()\n" ;
                fits_report_error (stderr , status) ;
                return -1 ;
            }
        }
    }
    fits_close_file (fptr , &status) ;
    if (status)
    {
        LOG(ERROR)<< "\nError closing file " << filename << " - writeHistory()\n" ;
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    return 0 ;
}
//------------------------------------------------------------------------------

int updateKeywords (char *filename , char *creator)
{
    int status = 0 , numhdu ;
    fitsfile *fptr ;
    fits_open_file (&fptr , filename , READWRITE , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in updateKeywords()\n" ;
        fits_report_error (stderr , status) ;
        return status ;
    }
    fits_get_num_hdus (fptr , &numhdu , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in updateKeywords()\n" ;
        fits_report_error (stderr , status) ;
        return status ;
    }
    for (int i = 1 ; i <= numhdu ; i++)
    {
        fits_movabs_hdu (fptr , i , NULL , &status) ;
        if (status)
        {
            LOG(ERROR) << "\nError in updateKeywords()\n" ;
            fits_report_error (stderr , status) ;
            return status ;
        }
        fits_write_date (fptr , &status) ;
        fits_update_key (fptr , TSTRING , "ORIGIN" , (char *) ORIGIN , NULL , &status) ;
        fits_update_key (fptr , TSTRING , "CREATOR" , creator , NULL , &status) ;
        fits_write_chksum (fptr , &status) ;
    }
    fits_close_file (fptr , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in updateKeywords()\n" ;
        fits_report_error (stderr , status) ;
        return status ;
    }
    return 0 ;
}
//------------------------------------------------------------------------------

int writeImg (char *file , long *img , int m , int n) //just to check;
{
    int status = 0 ;
    fitsfile *fptr ;
    fits_create_file (&fptr , file , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in creating file:" << file ;
        return -1 ;
    }
    int bitpix = LONG_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = m ;
    naxes[1] = n ;
    fits_create_img (fptr , bitpix , naxis , naxes , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in creating image:" << file ;
        return -1 ;
    }
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fits_write_pix (fptr , TLONG , fpixel , m*n , img , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in writing pixels:" << file ;
        return -1 ;
    }
    fits_close_file (fptr , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in closing file:" << file ;
        return -1 ;
    }
    return 0 ;
}

int writeImg (char *file , float *img , int m , int n) //just to check;
{
    int status = 0 ;
    fitsfile *fptr ;
    fits_create_file (&fptr , file , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in creating file:" << file ;
        return -1 ;
    }
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = m ;
    naxes[1] = n ;
    fits_create_img (fptr , bitpix , naxis , naxes , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in creating image:" << file ;
        return -1 ;
    }
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fits_write_pix (fptr , TFLOAT , fpixel , m*n , img , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in writing pixels:" << file ;
        return -1 ;
    }
    fits_close_file (fptr , &status) ;
    if (status)
    {
        LOG(ERROR) << "\nError in closing file:" << file ;
        return -1 ;
    }
    return 0 ;
}

int writeArray (char *file , float **array, int m , int n)
{
    fitsfile *fptr ;
    int status = 0 ;
    float data[m * n] ;
    int index = 0 ;
    for (int i = 0 ; i < m ; i++)
        for (int j = 0 ; j < n ; j++)
            data[index++] = array[i][j] ;
    fits_create_file (&fptr , file , &status) ;
    if (status)
    {
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = n ;
    naxes[1] = m ;
    fits_create_img (fptr , bitpix , naxis , naxes , &status) ;
    if (status)
    {
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    long fpixel[2] ;
    fpixel[0] = 1 ;
    fpixel[1] = 1 ;
    fits_write_pix (fptr , TFLOAT , fpixel , m*n , data , &status) ;
    if (status)
    {
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    fits_close_file (fptr , &status) ;
    if (status)
    {
        fits_report_error (stderr , status) ;
        return -1 ;
    }
    return (EXIT_SUCCESS) ;
}

void checkPFILESenv ()
{
    if (getenv (ENV_PFILES) == NULL)
    {
        LOG(ERROR) << "\n***Environment variable PFILES not set***\n" ;
        exit (EXIT_FAILURE) ;
    }
}

void checkParFile (char *modulename)
{
    char *pardir , parfilename[NAMESIZE] ;
    pardir = getenv (ENV_PFILES) ;
    if (pardir == NULL)
    {
        LOG(ERROR) << "\n***Environment variable PFILES not set***\n" ;
        exit (EXIT_FAILURE) ;
    }
    
    char *token;
    bool flag=false;	
	
   do{
	token=strtok(pardir,":");
	if(token==NULL) break;
    	strcpy (parfilename , token) ;
  	strcat (parfilename , "/") ;
    	strcat (parfilename , modulename) ;
    	strcat (parfilename , ".par") ;
    	//LOG(ERROR)<<"\n"<<parfilename<<endl;
    	if(FileExists(parfilename))
	{
		flag=true;
		break;
	}
	pardir=NULL;
       }while(token!=NULL);
    if (flag==false)
    {
        LOG(ERROR) << endl << "***Parameter file not found for the module '" << modulename << "'***" << endl ;
        exit (EXIT_FAILURE) ;
    }
}

int updateFilenameExtension (string *file , string oldExt , string newExt)
{
    //LOG(ERROR)<<endl<<*file<<"   "<<oldExt<<"   "<<newExt;
    int pos = (*file).find (oldExt) ;
    int pos_expected = (*file).size () - oldExt.size () ;
    if (pos != pos_expected)
    {
        LOG(ERROR) << endl << "Error in updating extension of file  " << *file << " from " << oldExt << " to " << newExt ;
        exit (EXIT_FAILURE) ;
    }
    string temp = (*file).substr (0 , pos) ;
    *file = temp + newExt ;
    //LOG(ERROR)<<endl<<*file;
    return (EXIT_SUCCESS) ;
}

void writeCommonKeywords (fitsfile *fptr , char *creator)
{
    int status = 0 ,numhdu=0;
    char filename[FLEN_FILENAME] ;
    fits_file_name (fptr , filename , &status) ;
   
//    if (status) return (EXIT_FAILURE) ;
    fits_get_num_hdus (fptr , &numhdu , &status) ;
  
    for (int i = 1 ; i <= numhdu ; i++)
    {
     fits_movabs_hdu (fptr ,i , NULL , &status) ;
    printError (status , "***Error in moving to perticuler HDU***") ;
    //fits_write_date (fptr , &status) ;
   // printError (status , "***Error writing date***") ;
    fits_update_key (fptr , TSTRING , "ORIGIN" , (char *) ORIGIN , NULL , &status) ;
    printError (status , "***Error writing Origin***") ;
    fits_update_key (fptr , TSTRING , "CREATOR" , creator , NULL , &status) ;
    printError (status , "***Error writing Creator***") ;
    fits_write_chksum (fptr , &status) ;
    printError (status , "***Error writing checksum***") ;
    }
//    fits_close_file(fptr,&status);
//    printError (status , "555") ;
    
}

void linear_fit (double * x_val , double * y_val , int n , double *cf)
{
    int ind ;
    double sum_x , sum_y , sum_xy , sum_x2 ;

    sum_x = sum_y = sum_xy = 0.0 ;

    for (ind = 0 ; ind < n ; ind++)
    {
        sum_x += x_val[ind] ;
        sum_x2 += x_val[ind]
                * x_val[ind] ;
        sum_y += y_val[ind] ;
        sum_xy += x_val[ind]
                * y_val[ind] ;
    }

    cf[1] = ((double) n * sum_xy - sum_x * sum_y) //cf[1] is gain
            / ((double) n * sum_x2 - sum_x * sum_x) ; //cf[0] is bias

    cf[0] = (sum_y - cf[1] * sum_x) / (double) n ;
}

string searchFile (char* dirname , char* pattern)
{
    //LOG(INFO)<<endl<<"Dirname :"<<dirname<<endl;
    DIR *dp ;
    char **sub ;
//LOG(INFO)<<(string)pattern;
    char *fitname =new char[100];
  fitname = pattern ;
//LOG(INFO)<<(string)fitname;
    string fitfile = " " ;
    sub = new char* [100] ;
    for (int i = 0 ; i < 100 ; i++)
        sub[i] = new char[30] ;
    int p = 0 ;
    int count = 0 ;
    struct dirent *dirp ;
    if ((dp = opendir ((const char *)dirname) )== NULL)
   {  LOG(ERROR) << endl << "***Error in opening " << dirname << " directory***" << endl ;
   return NULL;
    }
    while (dirp = readdir (dp))
    {
        p++ ;
        if (dirp->d_name == " . " || dirp->d_name == " .. "){
           
	
 continue ;
}
        sub[p - 1] = dirp->d_name ;
//	LOG(INFO)
        // LOG(ERROR)<<"\nthe sub"<<sub[p]<<endl;
        string temp = (string) sub[p - 1] ;
        //LOG(ERROR)<<"\nthe temp"<<temp<<endl;
        int size_temp = temp.size () ;
        //  LOG(ERROR)<<"\nthe fits name"<<fitname;
       
        int npos = temp.find (fitname , 0) ;
      
        if (npos > 0 && npos < size_temp)
        {

            fitfile = temp ;
            LOG(INFO)<<fitfile;
            //  exit(1);
            break ;
           // count++ ;
        }
        //LOG(ERROR)<<"\n"<<"the sub "<< fitfile<<"  position"<<npos<<endl;
        //  p++;

    }
    closedir (dp) ;
//const char *strout =fitfile.c_str();
	//return strout;;
//return (const char*)fitfile;
LOG(INFO)<<fitfile.data();
    return fitfile.data() ;

}

void getKeywordVal (char *file , char *keyname , int hdunum , char *val)
{
    fitsfile *fptr ;
    int status = 0 ;
    fits_open_file (&fptr , file , READONLY , &status) ;
    printError (status , "") ;
    fits_movabs_hdu (fptr , hdunum , NULL , &status) ;
    printError (status , "") ;
    fits_read_key (fptr , TSTRING , keyname , val , NULL , &status) ;
    printError (status , "") ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file") ;
}

int Applypadding (float *inputArray , int input_xsize,int input_ysize,
        float *outputArray,int padding_xsize ,int padding_ysize)
{
    int count = 0 ;
//    for(int i=0;i<padding_xsize*padding_ysize;i++){
//        outputArray[i]=0.0;
//    }
    int col_padd = padding_xsize - input_xsize ;
    int eachside_col = col_padd / 2 ;
  int row_padd = padding_ysize - input_ysize ;
  int eachside_row=row_padd/2;
    int p=0 ;
    for (long r = 0 ; r <input_xsize ; r++)
    {
        for (long j = p = (eachside_row + r) * padding_xsize + eachside_col ; j < p + input_xsize ; j++)
        {
            outputArray[j] = inputArray[count] ;
            count++ ;

        }
    }
    return (EXIT_SUCCESS) ;
}
int  copyUsrkeywrdsTovect(fitsfile *fin , vector<string> &strvect){
   
     strvect.clear();
     int status = 0 , keyexist ;
    char record[FLEN_CARD] ;
    char filename[FLEN_FILENAME] ;
    fits_file_name (fin , filename , &status) ;
    if (status) return (EXIT_FAILURE) ;
    fits_get_hdrspace (fin , &keyexist , NULL , &status) ;
    if (status)
    {
        LOG(ERROR)<< endl << "***Could not find number of keywords in file " << filename << "***" ;
        fits_report_error (stderr , status) ;
        return (EXIT_FAILURE) ;
    }
    int keyclass ;
    //LOG(ERROR)<<"\nNumber of keywords found:"<<keyexist;
    for (int i = 1 ; i <= keyexist ; i++)
    {
        fits_read_record (fin , i , record , &status) ;
        if (status)
        {
            LOG(ERROR) << endl << "***Error in reading record number " << i << " in file " << filename << "***" ;
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }
     //   keyclass = fits_get_keyclass (record) ;

       //  if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY || keyclass ==TYP_COMM_KEY)
       // {
       
            strvect.push_back (record);
       // }
    }
//    fits_close_file(fin,&status);
//    printError (status , "Error in closing the input File" , filename) ;
    return (EXIT_SUCCESS) ;
     
     }
     
int writeUsrkeywordsFrmvect(char *filename,vector<string> &strvect){
            
    fitsfile *fout;
            int status=0;
            int numhdu=0;
             fits_open_file (&fout , filename , READWRITE , &status) ;
              printError (status , "Error in opening the input File" , filename) ;
           // LOG(INFO)<<"Successfully opened"<<endl;
               fits_get_num_hdus (fout , &numhdu , &status) ;
                               
               for(int hno=1;hno<=1;hno++){
                    fits_movabs_hdu (fout , hno , NULL , &status) ;

        for (int i = 0; i < strvect.size(); i++) {
            if (strstr(strvect[i].c_str(), "NAXIS") == NULL && strstr(strvect[i].c_str () ,"BITPIX")==NULL && strstr(strvect[i].c_str () ,"NAXIS1")==NULL 
                     && strstr(strvect[i].c_str () ,"NAXES2")==NULL   && strstr(strvect[i].c_str(),"SIMPLE")==NULL  && strstr(strvect[i].c_str(),"EXTEND")==NULL   && strstr(strvect[i].c_str(),"DATASUM")==NULL   && strstr(strvect[i].c_str(),"CREATOR")==NULL && strstr(strvect[i].c_str(),"XSIZE")==NULL && strstr(strvect[i].c_str(),"YSIZE")==NULL && strstr(strvect[i].c_str(),"EXP_TIME")==NULL && strstr(strvect[i].c_str(),"NAMEPRFX")==NULL && strstr(strvect[i].c_str(),"CHECKSUM")==NULL) {
                fits_write_record(fout, strvect[i].c_str(), &status);


                //LOG(INFO)<<strvect[i].c_str ()<<endl;
                if (status) {
                    LOG(ERROR) << endl << "***Error in writing keyword " << strvect[i].c_str() << "***";
                    fits_report_error(stderr, status);
                    fits_close_file(fout, &status);
                    printError(status, "Error in closing the input File", filename);
                    //     LOG(INFO)<<"Successfully closed with error"<<endl;
                    return (EXIT_FAILURE);
                }
            }
        }
                           
               }
               
            fits_close_file(fout,&status);
            printError (status , "Error in closing the input File" , filename) ;
           // LOG(INFO)<<"Successfully closed"<<endl;
//            fits_open_file (&fout , filename , READWRITE , &status) ;
//              printError (status , "Error in opening the input File" , filename) ;
//                fits_close_file(fout,&status);
//            printError (status , "Error in closing the input File" , filename) ;
            
            return(EXIT_SUCCESS);
        }

int createOutputDirectory(int clobber, char *outpath)
{
    if (DirExists (outpath) && !(clobber))  
    {
        LOG (ERROR) << "***Output Directory " << outpath << "    already exists....***" ;
        LOG (ERROR) << "***Use clobber = y for overwriting***" << endl ;
        return (EXIT_FAILURE) ;
    }
    else
    {
      if(DirExists(outpath)){
                string command1 = "rm -rf  " + (string) outpath;
                LOG(INFO)<<"Executing......... "<<command1;
                system (command1.c_str ()) ;    //Executing the shell command for removing the output directory 
                LOG(INFO)<<outpath<<" directory deleted";
        }
        string command2 = "mkdir -p " + (string) outpath ;
        system (command2.c_str ()) ; //Executing the shell command for creating the output directory
        LOG (INFO) << outpath << " directory created" ;
    }
    return (EXIT_SUCCESS);
}

int writeColumnsToFITS(char *filename,int hduno,int n, ...)
{
   va_list arglist;
   va_start(arglist,n*3);
    
    int status=0;   
    fitsfile *fptr;
    
    fits_open_file(&fptr, filename, READWRITE, &status);
    printError(status, " Error in opening file " , filename);
    
    fits_movabs_hdu(fptr,hduno,NULL,&status);
    printError(status, "Error in moving to HDU ",filename);
        
    int datatype,coln,totalEle;
    char *str;
    int *intval;
    float *fval;
    double *dval;
    unsigned short *usval;
    short *sval ;
    unsigned char *ucval;
    unsigned int *uintval;
    long *longval ;
    unsigned long *ulongval;
        
    for(int i=0;i<n;i++)
    {
        datatype = va_arg(arglist,int);
       
         coln=va_arg(arglist,int);
         
         switch(datatype){
             
             case TSTRING:
                str = va_arg(arglist,char *);
                 totalEle=va_arg(arglist,int);
                fits_write_col (fptr , TSTRING , coln , 1 , 1 , totalEle , str , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                
            case TLOGICAL:
                LOG(ERROR)<<endl<<"TLOGICAL not supported by this function";
                return (EXIT_FAILURE);
                break;
                
            case TBYTE:
                 ucval = va_arg(arglist, unsigned char *);
                  totalEle=va_arg(arglist,int);  
                 fits_write_col (fptr , TBYTE , coln , 1 , 1 , totalEle , ucval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                 
            case TSHORT:
                 sval = va_arg(arglist, short*);
                  totalEle=va_arg(arglist,int);
                 fits_write_col (fptr , TSHORT , coln , 1 , 1 , totalEle , sval , &status) ;
                 printError (status , "Error in writing the column of  time ",filename) ;
                break;
                 
            case TUSHORT:
                 usval = va_arg(arglist, unsigned short *);
                  totalEle=va_arg(arglist,int);
                  fits_write_col (fptr , TUSHORT , coln , 1 , 1 , totalEle , usval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                  
            case TINT:
                 intval = va_arg(arglist, int *);
                  totalEle=va_arg(arglist,int);
                 fits_write_col (fptr , TINT , coln , 1 , 1 , totalEle , intval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                 break;
                  
            case TUINT:
                 uintval = va_arg(arglist, unsigned int *);
                  totalEle=va_arg(arglist,int);
                 fits_write_col (fptr , TUINT , coln , 1 , 1 , totalEle , uintval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                
            case TLONG:
                 longval = va_arg(arglist, long*);
                  totalEle=va_arg(arglist,int);
                 fits_write_col (fptr , TLONG , coln , 1 , 1 , totalEle , longval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                 
            case TULONG:
                  ulongval = va_arg(arglist, unsigned long *);
                  totalEle=va_arg(arglist,int);
                  fits_write_col (fptr , TULONG , coln , 1 , 1 , totalEle ,ulongval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                  
            case TFLOAT:
                fval = va_arg(arglist, float *);
                 totalEle=va_arg(arglist,int);
                fits_write_col (fptr , TFLOAT , coln , 1 , 1 , totalEle , fval , &status) ;
                printError (status , "Error in writing the column of  time ",filename) ;
                 break;
                 
            case TDOUBLE:
                  dval = va_arg(arglist, double *);
                   totalEle=va_arg(arglist,int);
                  fits_write_col (fptr , TDOUBLE , coln , 1 , 1 , totalEle , dval , &status) ;
                  printError (status , "Error in writing the column of  time ",filename) ;
                break;
                 
            case TCOMPLEX:
                LOG(ERROR) << endl << "TCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            case TDBLCOMPLEX:
                LOG(ERROR) << endl << "TDBLCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            default: 
                LOG(ERROR) << endl << "Unsupported datatype : "<<datatype << endl;
                return (EXIT_FAILURE);
         }
    }
    fits_close_file(fptr,&status);
    printError(status, "Error in closing file  ",filename);
    
    va_end(arglist);
    return(EXIT_SUCCESS);
}


int readColumnsFromFITS(char *filename,int hduno,int n, ...)
{
   va_list arglist;
   va_start(arglist,n*4);
    
    int status=0;
    
    fitsfile *fptr;
    
    fits_open_file(&fptr, filename, READONLY, &status);
    printError(status, " Error in opening file " , filename);
    
    fits_movabs_hdu(fptr,hduno,NULL,&status);
    printError(status, "Error in moving to HDU ",filename);
        
    int datatype,coln,totalEle;
   
    char *str;
    int *intval;
    float *fval ;
    double *dval;
    unsigned short *usval;
    short *sval;
    unsigned char *ucval;
    unsigned int *uintval;
    long *longval;
    unsigned long *ulongval;
        
    for(int i=0;i<n;i++)
    {
        datatype = va_arg(arglist,int);
        coln=va_arg(arglist,int);
                
         switch(datatype){
             
            case TSTRING:
                str = va_arg(arglist,char *);
                totalEle=va_arg(arglist,int);
                fits_read_col (fptr , TSTRING , coln , 1 , 1 , totalEle , NULL , str , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                
            case TLOGICAL:
                LOG(ERROR)<<endl<<"TLOGICAL not supported by this function";
                return (EXIT_FAILURE);
                break;
                
            case TBYTE:
                 ucval = va_arg(arglist, unsigned char *);
                   totalEle=va_arg(arglist,int);
                 fits_read_col (fptr , TBYTE , coln , 1 , 1 , totalEle , NULL , ucval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                 
            case TSHORT:
                 sval = va_arg(arglist, short*);
                   totalEle=va_arg(arglist,int);
                 fits_read_col (fptr , TSHORT , coln , 1 , 1 , totalEle , NULL , sval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                 
            case TUSHORT:
                  usval = va_arg(arglist, unsigned short *);
                  totalEle=va_arg(arglist,int);
                 
                  fits_read_col (fptr , TUSHORT , coln , 1 , 1 , totalEle , NULL ,usval , NULL , &status) ;
                 printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                  
            case TINT:
                 intval = va_arg(arglist, int *);
                   totalEle=va_arg(arglist,int);
                   fits_read_col (fptr , TINT , coln , 1 , 1 , totalEle , NULL , intval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                  
            case TUINT:
                 uintval = va_arg(arglist, unsigned int *);
                   totalEle=va_arg(arglist,int);
                  fits_read_col (fptr , TUINT , coln , 1 , 1 , totalEle , NULL , uintval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                
            case TLONG:
                 longval = va_arg(arglist, long*);
                   totalEle=va_arg(arglist,int);
                  fits_read_col (fptr , TLONG , coln , 1 , 1 , totalEle , NULL , longval , NULL , &status) ;
                   printError (status , "Reading a column Fails in caldb",filename) ;
                   break;
                 
            case TULONG:
                  ulongval = va_arg(arglist, unsigned long *);
                    totalEle=va_arg(arglist,int);
                   fits_read_col (fptr , TULONG , coln , 1 , 1 , totalEle , NULL , ulongval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                  
            case TFLOAT:
                 fval = va_arg(arglist, float *);
                   totalEle=va_arg(arglist,int);
                fits_read_col (fptr , TFLOAT , coln , 1 , 1 , totalEle , NULL , fval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                 
            case TDOUBLE:
                  dval = va_arg(arglist, double *);
                    totalEle=va_arg(arglist,int);
                 fits_read_col (fptr , TDOUBLE , coln , 1 , 1 , totalEle , NULL , dval , NULL , &status) ;
                printError (status , "Reading a column Fails in caldb",filename) ;
                 break;
                 
            case TCOMPLEX:
                LOG(ERROR) << endl << "TCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            case TDBLCOMPLEX:
                LOG(ERROR) << endl << "TDBLCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            default: 
                LOG(ERROR) << endl << "Unsupported datatype : "<<datatype << endl;
                return (EXIT_FAILURE);
         }
        
    }
    fits_close_file(fptr,&status);
    printError(status, "Error in closing file  ",filename);
    
    va_end(arglist);
    return(EXIT_SUCCESS);
}



int readParams(int argc, char **argv, int n, ...)
 {
     va_list params_list;
     va_start(params_list,n*3);
        
    int status = 0 ;
    if (PIL_OK != (status = PILInit (argc , argv)))            //Initialize PIL
    {
        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
     
     int datatype;
     char *param_name;
     int *intname, *boolname;
     float *real4name;
     double *realname;
     char *fname,*stringname;
     int tempval;
     float fval;
    double dval;
     
    for(int i=0;i<n;i++)
    {               //loop for number of parameters
        datatype=va_arg(params_list,int);
        param_name=va_arg(params_list,char *);
            
        switch(datatype){
            case BOOL:  
                boolname = va_arg(params_list, int*);
                if (PIL_OK != (status = PILGetBool (param_name , &tempval)))
                {
                    LOG(ERROR) << endl << "***Error reading parameter : "<<param_name<<"***" ;
                    return status ;
                }
                *boolname = tempval;
                break;
                
            case INT: 
                intname = va_arg(params_list, int*);
                if (PIL_OK != (status = PILGetInt (param_name ,&tempval)))
                {
                    LOG(ERROR) << endl <<"***Error reading parameter : "<<param_name<<"***" ;
                    return status ;
                }
                *intname = tempval;
                break;
                
            case REAL: 
                realname = va_arg(params_list, double*);
                if (PIL_OK != (status = PILGetReal (param_name , &dval)))
                {
                    LOG(ERROR) << endl << "***Error reading parameter : "<<param_name<<"***" ;
                    return status ;
                }
                *realname = dval;
                break;
                
            case REAL4: 
                real4name = va_arg(params_list,float*);
                if (PIL_OK != (status = PILGetReal4(param_name , &fval)))
                {
                    LOG(ERROR) << endl << "***Error reading parameter : "<<param_name<<"***" ;
                    return status ;
                }
                *real4name = fval;
                break;
                
            case STRING: 
                stringname = va_arg(params_list,char*);
                //LOG(INFO)<<"HEAR"<<endl;
                if (PIL_OK != (status = PILGetString(param_name , stringname)))
                {
                    LOG(ERROR) << endl << "***Error reading parameter : "<<param_name<<"***" ;
                    return status ;
                }
                break;
                
            case FNAME:  
                fname = va_arg(params_list,char *);
                if (PIL_OK != (status = PILGetFname (param_name , fname)))
                {
                    LOG(ERROR) << endl << "***Error reading parameter : "<<param_name<<"***" ;
                    return status ;
                }
                     
        }
        
    }
    
     PILClose (status) ; 
     va_end(params_list);
     
     return (EXIT_SUCCESS);
  }

// int updateKeywords(char *filename,char *origin){
//
//    fitsfile *fptr ;
//    int status = 0 ;
//    int numhdu = 0 ;
//    fits_open_file (&fptr , filename , READWRITE , &status) ;
//    printError (status , "Error in opening the input File" , filename) ;
//    fits_get_num_hdus (fptr , &numhdu , &status) ;
//
//    for (int hno = 1 ; hno <= numhdu ; hno++)
//    {
//        fits_movabs_hdu (fptr , hno , NULL , &status) ;
//        fits_write_date (fptr , &status) ;
//        printError (status , "***Error writing date***") ;
//        fits_update_key (fptr , TSTRING , "ORIGIN" , (char *) ORIGIN , NULL , &status) ;
//        printError (status , "***Error writing Origin***") ;
//        fits_update_key (fptr , TSTRING , "CREATOR" , creator , NULL , &status) ;
//        printError (status , "***Error writing Creator***") ;
//        fits_write_chksum (fptr , &status) ;
//        printError (status , "***Error writing checksum***") ;
//
//
//    }
//      fits_close_file(fptr,&status);
//      printError (status , "***Error in closing the file ***") ;
//   
//    
//    return(EXIT_SUCCESS);
// }

//------------------------------------------------------------------------------------------------------------------------------------------------
/**
 * Function to read multiple keywords fitsfile
 * @param fitsfilename
 * @param hdunumber
 * @param n
 * @param ... - Argument list in the sets of datatype (of keyword), keyword name and variable to store 
 *                      the keyword value
 * @return 
 */

int readKeywords(char *fitsfilename, int hdunumber,int n,...)
{
    va_list arglist;
    va_start(arglist,n*3);
    
    int status=0;
    
    fitsfile *fptr;
    
    fits_open_file(&fptr, fitsfilename, READONLY, &status);
    printError(status, " Error in opening file " , fitsfilename);
   // LOG(INFO)<<"Successfully opened in PC mode"<<endl;
    fits_movabs_hdu(fptr,hdunumber,NULL,&status);
    printError(status, "Error in moving to HDU ",fitsfilename);
        
    int datatype;
    char *keyname;
   char *str;
     //string str;
    int *intval;
    float *fval;
    double *dval;
    unsigned short *usval ;
    short *sval;
    unsigned char *ucval;
    unsigned int *uintval;
    long *longval;
    unsigned long *ulongval;
        
    for(int i=0;i<n;i++)
    {
        datatype = va_arg(arglist,int);
         keyname = va_arg(arglist,char *);
       
         string errmsg = "***Error in reading "+(string)keyname+"  keyword***";
         
         switch(datatype){
            case TSTRING:
                str = va_arg(arglist,char* );
            
                fits_read_key(fptr, datatype, keyname, str, NULL, &status);
                break;
       
            case TLOGICAL:
                LOG(ERROR)<<endl<<"TLOGICAL not supported by this function";
                return (EXIT_FAILURE);
                break;
                
            case TBYTE:
                 ucval = va_arg(arglist, unsigned char *);
                 fits_read_key(fptr, datatype, keyname,ucval, NULL, &status);
                 break;
                 
            case TSHORT:
                 sval = va_arg(arglist, short*);
                 fits_read_key(fptr, datatype, keyname,sval,NULL, &status);
                 break;
                 
            case TUSHORT:
                 usval = va_arg(arglist, unsigned short *);
                 fits_read_key(fptr, datatype, keyname,usval, NULL, &status);
                  break;
                  
            case TINT:
                 intval = va_arg(arglist, int *);
                 fits_read_key(fptr, datatype, keyname,intval, NULL, &status);
                  break;
                  
            case TUINT:
                 uintval = va_arg(arglist, unsigned int *);
                 fits_read_key(fptr, datatype, keyname,uintval, NULL, &status);
                break;
                
            case TLONG:
                 longval = va_arg(arglist, long*);
                 fits_read_key(fptr, datatype, keyname,longval, NULL, &status);
                 break;
                 
            case TULONG:
                  ulongval = va_arg(arglist, unsigned long *);
                  fits_read_key(fptr, datatype, keyname,ulongval, NULL, &status);
                  break;
                  
            case TFLOAT:
                 fval = va_arg(arglist, float *);
                 fits_read_key(fptr, datatype, keyname,fval, NULL, &status);
                 break;
                 
            case TDOUBLE:
                  dval = va_arg(arglist, double *);
                 fits_read_key(fptr, datatype, keyname,dval, NULL, &status);
                 break;
                 
            case TCOMPLEX:
                LOG(ERROR) << endl << "TCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            case TDBLCOMPLEX:
                LOG(ERROR) << endl << "TDBLCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            default: 
                LOG(ERROR) << endl << "Unsupported datatype : "<<datatype << endl;
                return (EXIT_FAILURE);
         }
         printError(status,(char*)errmsg.c_str (),fitsfilename);
    }
    fits_close_file(fptr,&status);
    printError(status, "Error in closing file  ",fitsfilename);
  //  LOG(INFO)<<"Successfully closed in PC mode"<<endl;
    va_end(arglist);
    
    return (EXIT_SUCCESS);
}
//----------------------------------------------------------------------------------------------------------------------------------------------


/**
 * Function to update keywords in fitsfile
 * @param fitsfilename
 * @param hdunumber
 * @param n
 * @param ... - Argument list in the sets of datatype (of keyword), keyword name and variable to store 
 *                      the keyword value
 * @return 
 */
int updateKeywords(char *fitsfilename, int hdunumber,int n,...){
    va_list arglist;
    va_start(arglist,n*3);
    
    int status=0;
    fitsfile *fptr;
    fits_open_file(&fptr, fitsfilename, READWRITE, &status);
    printError(status, " Error in opening file " , fitsfilename);
    fits_movabs_hdu(fptr,hdunumber,NULL,&status);
    printError(status, "Error in moving to other hdu ",fitsfilename);
        
    int datatype;
    char *keyname;
    char *str;
    int *intval;
    float *fval;
    double *dval;
    unsigned short *usval ;
    short *sval;
    unsigned char *ucval;
    unsigned int *uintval;
    long *longval;
    unsigned long *ulongval;
        
    for(int i=0;i<n;i++)
    {
         datatype = va_arg(arglist,int);
         keyname = va_arg(arglist,char *);
               
         string errmsg = "***Error in updating "+(string)keyname+"  keyword***";
         
         switch(datatype){
            case TSTRING:
                str = va_arg(arglist,char *);
                fits_update_key(fptr, datatype, keyname, str, NULL, &status);
                break;
                
            case TLOGICAL:
                LOG(ERROR)<<endl<<"TLOGICAL not supported by this function";
                return (EXIT_FAILURE);
                break;
                
            case TBYTE:
                 ucval = va_arg(arglist, unsigned char *);
                 fits_update_key(fptr, datatype, keyname,ucval, NULL, &status);
                 break;
                 
            case TSHORT:
                 sval = va_arg(arglist, short*);
                 fits_update_key(fptr, datatype, keyname,sval,NULL, &status);
                 break;
                 
            case TUSHORT:
                 usval = va_arg(arglist, unsigned short *);
                 fits_update_key(fptr, datatype, keyname,usval, NULL, &status);
                  break;
                  
            case TINT:
                 intval = va_arg(arglist, int *);
                 fits_update_key(fptr, datatype, keyname,intval, NULL, &status);
                  break;
                  
            case TUINT:
                 uintval = va_arg(arglist, unsigned int *);
                 fits_update_key(fptr, datatype, keyname,uintval, NULL, &status);
                break;
                
            case TLONG:
                 longval = va_arg(arglist, long*);
                 fits_update_key(fptr, datatype, keyname,longval, NULL, &status);
                 break;
                 
            case TULONG:
                  ulongval = va_arg(arglist, unsigned long *);
                  fits_update_key(fptr, datatype, keyname,ulongval, NULL, &status);
                  break;
                  
            case TFLOAT:
                 fval = va_arg(arglist, float *);
                 fits_update_key(fptr, datatype, keyname,fval, NULL, &status);
                 break;
                 
            case TDOUBLE:
                  dval = va_arg(arglist, double *);
                 fits_update_key(fptr, datatype, keyname,dval, NULL, &status);
                 break;
                 
            case TCOMPLEX:
                LOG(ERROR) << endl << "TCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            case TDBLCOMPLEX:
                LOG(ERROR) << endl << "TDBLCOMPLEX not supported by this function" << endl;
                return (EXIT_FAILURE);
                break;
                
            default: 
                LOG(ERROR) << endl << "Unsupported datatype : "<<datatype << endl;
                return (EXIT_FAILURE);
         }
         printError(status,(char *)errmsg.c_str(),fitsfilename);
    }
    fits_close_file(fptr,&status);
    printError(status, "Error in closing file  ",fitsfilename);
    va_end(arglist);
    return (EXIT_SUCCESS);
}


int  getFileNamesfrmDir(char* dirname, char* pattern,vector<string> &filenames)
{
    DIR *dp ;
    char **sub ;
    char *fitname = pattern ;
    string fitfile = "" ;
    sub = new char* [100] ;
    for (int i = 0 ; i < 100 ; i++)
        sub[i] = new char[30] ;
    int p = 0 ;
    int count = 0 ;
    struct dirent *dirp ;
   
    if ((dp = opendir ((const char *)dirname) )== NULL)
   //  if ((dp = opendir (dirname))= NULL)
    {  LOG(ERROR) << endl << "***Error in opening " << dirname << " directory***" << endl ;
        exit (EXIT_FAILURE) ;
    }
    
    while (dirp = readdir (dp))
    {
        p++ ;
        if (dirp->d_name == " . " || dirp->d_name == " .. ")
            continue ;
        sub[p - 1] = dirp->d_name ;
        // LOG(ERROR)<<"\nthe sub"<<sub[p]<<endl;
        string temp = (string) sub[p - 1] ;
        //LOG(ERROR)<<"\nthe temp"<<temp<<endl;
        int size_temp = temp.size () ;
        //  LOG(ERROR)<<"\nthe fits name"<<fitname;
        int npos = temp.find (fitname , 0) ;
        if (npos > 0 && npos < size_temp)
        {
           
            fitfile = temp ;
            filenames.push_back (fitfile);
            count++ ;
            
            
        }  

    }

   // closedir (dp) ;
  //exit(1);
    return (EXIT_SUCCESS);
}

string  getSerialNo(int &p){
    string  strnum;
    p=p+1;
    sprintf((char*)strnum.c_str(),"P%d",p);
    cout<<strnum;
   
    
    return strnum;
}


//----------------------------------------------------------------------------------------------------------------------------------------------
int extractTars(string tar_name,string &fileoutname,string &orbitnum)
{
    string tarname=tar_name;
    string strtemp (tarname) ;
     int pos = strtemp.find_last_of ("/") ;
    string tempFilename="";
    string cmd;
    //cmd.repla
     string tar_date_data, tar_OBS_id_data;
    try{
    tar_date_data=tarname.substr (pos+1+11,8);
     tar_OBS_id_data=tarname.substr (pos+1+19,21);
    orbitnum=tarname.substr(pos+1+41,5);
    }
    catch(...){
        LOG(ERROR)<<"***Level-1 Tar file is not in a proper format.***";
        return(EXIT_FAILURE);
    }
    // LOG(INFO)<<"The orbit num is "<<orbitnum;exit(1);
     if (pos > 0 && pos < strtemp.length ())
    {
      tempFilename =strtemp.substr (0 , pos) ;  
       cmd="tar -xvf  "+tarname+" --directory="+tempFilename;
       fileoutname=tempFilename+"/"+tar_date_data+"_"+tar_OBS_id_data+"_level1";
   
    }
     else
     {
         cmd ="tar -xvf  "+tarname;
         fileoutname=tar_date_data+"_"+tar_OBS_id_data+"_level1";
     }
    
    LOG(INFO)<<"Executing "<<cmd;
    system(cmd.c_str ());
    return(EXIT_SUCCESS);
    
}

int extractZip(string tar_name,string &fileoutname,string &orbitnum)
{
    string tarname=tar_name;
    string strtemp (tarname) ;
     int pos = strtemp.find_last_of ("/") ;
    string tempFilename="";
    string cmd;
    //cmd.repla
     string tar_date_data, tar_OBS_id_data;
    try{
    tar_date_data=tarname.substr (pos+1+11,8);
     tar_OBS_id_data=tarname.substr (pos+1+19,21);
    orbitnum=tarname.substr(pos+1+41,5);
    }
    catch(...){
        LOG(ERROR)<<"***Level-1 Tar file is not in a proper format.***";
        return(EXIT_FAILURE);
    }
    // LOG(INFO)<<"The orbit num is "<<orbitnum;exit(1);
     if (pos > 0 && pos < strtemp.length ())
    {
      tempFilename =strtemp.substr (0 , pos) ;  
       cmd="unzip -o "+tarname+" -d "+tempFilename;
       fileoutname=tempFilename+"/"+tar_date_data+"_"+tar_OBS_id_data+"_level1";
   
    }
     else
     {
         cmd ="unzip -o "+tarname;
         fileoutname=tar_date_data+"_"+tar_OBS_id_data+"_level1";
     }
    
    LOG(INFO)<<"Executing "<<cmd;
    system(cmd.c_str ());
    return(EXIT_SUCCESS);
    
}

int extractTars(string tar_name,string &fileoutname)
{
    string tarname=tar_name;
   
     string strtemp (tarname) ;
    int pos = strtemp.find_last_of ("/") ;
    string tempFilename="";
    string cmd;
     string tar_date_data=tarname.substr (pos+1+11,8);
    string tar_OBS_id_data=tarname.substr (pos+1+19,21);
    //orbitnum=tarname.substr(pos+1+41,5);
   // LOG(INFO)<<"The orbit num is "<<orbitnum;exit(1);
     if (pos > 0 && pos < strtemp.length ())
    {
      tempFilename =strtemp.substr (0 , pos) ;  
       cmd="tar -xvf  "+tarname+" --directory="+tempFilename;
       fileoutname=tempFilename+"/"+tar_date_data+"_"+tar_OBS_id_data+"_level1";
   
    }
     else
     {
         cmd ="tar -xvf  "+tarname;
         fileoutname=tar_date_data+"_"+tar_OBS_id_data+"_level1";
     }
    
    LOG(INFO)<<"Executing "<<cmd;

    system(cmd.c_str ());
        
    
    return(EXIT_SUCCESS);
}

int getTemp (char *lbtfile,char *detector,double *time_lbt,float *insideTemp,float *outsideTemp,long nrows_lbt)
{
    int status = 0 ;
    //reading observation id from data
   
    
    fitsfile *flbt ;
    fits_open_file (&flbt , lbtfile , READONLY , &status) ;
    printError (status , "Error opening LBT file",lbtfile) ;
   
    char obsid[FLEN_KEYWORD];
    
    //reading observation id from header of LBT file
    fits_read_key(flbt, TSTRING, "OBS_ID",obsid, NULL, &status);
    printError(status, " Error reading OBS_ID from header " ,lbtfile);
    
    fits_movabs_hdu (flbt , 2 , NULL , &status) ;
    printError (status , "Error moving to HDU 2 of lbtfile",lbtfile) ;
    //check when  number of Rows in lbt file > 1(TBD)
//    fits_get_num_rows (flbt , &nrows_lbt , &status) ;
//    printError (status , "Error reading the number of rows in lbt file",lbtfile) ;
     int colinside = 0 , coloutside = 0 ;             //variables to store column numbers to be used from LBT file
    char *md = detector;      //read channel for data
    
    LOG(INFO)<<"Channel :" <<md<<endl;
          
    if (md = (char *) "NUV")
    {
        colinside = INSIDE_TEMP_NUV ;
        coloutside = OUTSIDE_TEMP_NUV ;
    }
    else if (md = (char *) "FUV")
    {
        colinside =INSIDE_TEMP_FUV ;
        coloutside = OUTSIDE_TEMP_FUV ;
    }
    else if (md = (char *) "VIS")
    {
        colinside = INSIDE_TEMP_VIS;
        coloutside = OUTSIDE_TEMP_VIS;
    }
  // time_lbt=new double[nrows_lbt];
//   insideTemp = new float [nrows_lbt];
   //outsideTemp= new float[nrows_lbt];
//    int rowno;              
    fits_read_col (flbt , TDOUBLE , 1 , 1 , 1 ,nrows_lbt , NULL , time_lbt , NULL , &status) ;
    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
    fits_read_col (flbt , TFLOAT , colinside , 1 , 1 ,nrows_lbt , NULL , insideTemp , NULL , &status) ;
    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
    fits_read_col (flbt , TFLOAT , coloutside , 1 , 1 ,nrows_lbt , NULL , outsideTemp , NULL , &status) ;
    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
                   //row number in LBT file from which the temperatures will be read
     fits_close_file (flbt , &status) ;
   printError (status , "Error in closing file",lbtfile) ;
   
  
    return (EXIT_SUCCESS) ;
}
int readNumRowsFromFits(char *filename,int hduno,long &numrows)
{
    fitsfile *fptr;
    int status=0;
    fits_open_file (&fptr , filename , READONLY , &status) ;
    printError (status , "Error opening LBT file",filename) ;
   
     
    fits_movabs_hdu (fptr , hduno , NULL , &status) ;
    printError (status , "Error moving to HDU 2 of lbtfile",filename) ;
    //check when  number of Rows in lbt file > 1(TBD)
    fits_get_num_rows (fptr , &numrows , &status) ;
    printError (status , "Error reading the number of in lbt file",filename) ;

    fits_close_file (fptr , &status) ;
   printError (status , "Error in closing file",filename) ;
    return (EXIT_SUCCESS);
}
int readQEMCPFile (char *qeFile,char *detector, long nrows,float *temp,float *f0,float *f1,float *f2,float *f3,float *f4,float*f5,float*f6,float *f7)
{
    LOG(INFO) << "Reading QE MCP  Temperature vs filter file from calDB........" ;
    fitsfile *fqemcp ;
    int status = 0 ;
    fits_open_file (&fqemcp , qeFile , READONLY , &status) ;
    printError (status , "Error in opening the qeFile ",qeFile) ;
    int hdutype ;
    fits_movabs_hdu (fqemcp , 2 , &hdutype , &status) ;
    printError (status , "Error in moving to 2nd HDU in qeFile  ",qeFile) ;
    if (hdutype != BINARY_TBL)
    {
        LOG(ERROR) << endl << "***Expected binary table at hdu 2 of temperature vs filter file*** " << endl ;
        return (EXIT_FAILURE) ;
    }
//    long nrows ;
//    fits_get_num_rows (fqemcp , &nrows , &status) ;
//    printError (status , "Error in readQEMCP()",qeFile) ;
//    nCalDBTempValues = nrows ;
//    temp = new float[nrows] ;
//    f0 = new float[nrows] ;
//    f1 = new float[nrows] ;
//    f2 = new float[nrows] ;
//    f3 = new float[nrows] ;
//    f5 = new float[nrows] ;
//    f6 = new float[nrows] ;
//    f7 = new float[nrows] ;
//    f4 = new float[nrows] ;
    
    fits_read_col (fqemcp , TFLOAT , 1 , 1 , 1 , nrows , NULL , (void*) temp , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 2 , 1 , 1 , nrows , NULL , (void*) f0 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;   
    fits_read_col (fqemcp , TFLOAT , 3 , 1 , 1 , nrows , NULL , (void*) f1 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 4 , 1 , 1 , nrows , NULL , (void*) f2 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 5 , 1 , 1 , nrows , NULL , (void*) f3 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 6 , 1 , 1 , nrows , NULL , (void*) f4 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 7 , 1 , 1 , nrows , NULL , (void*) f5 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    
    //tobe removed
   fits_read_col (fqemcp , TFLOAT , 8 , 1 , 1 , nrows , NULL , (void*) f6 , NULL , &status) ;
          printError (status , "Reading a column Fails in caldb",qeFile) ;
    if(strcmp(detector,"FUV")==0 || strcmp(detector,"NUV")==0)         
    {     
         fits_read_col (fqemcp , TFLOAT , 8 , 1 , 1 , nrows , NULL , (void*) f6 , NULL , &status) ;
          printError (status , "Reading a column Fails in caldb",qeFile) ;
        fits_read_col (fqemcp , TFLOAT , 9 , 1 , 1 , nrows , NULL , (void*) f7 , NULL , &status) ;
            printError (status , "Reading a column Fails in caldb",qeFile) ;
    }
    LOG(INFO) <<"Reading QEMCP Temperature vs filter file from caldb Finished.." ;
    fits_close_file (fqemcp , &status) ;
    printError (status , "Error in closing qeFile",qeFile) ;
    return (EXIT_SUCCESS) ;
} 
int performUnitConversionIM(float *frmdata, float *expdata,double  intgrntime,int sizex,int sizey){
    
      if (intgrntime == 0)
        {
            LOG(ERROR) << "***Divide by Zero***" << endl ;
            return (EXIT_FAILURE) ;

        }
     
      for(int i=0;i<sizex*sizey;i++){
          if(frmdata[i]!=INVALID_PIX_VALUE){
          frmdata[i]=frmdata[i]/intgrntime;
          }
          else if(frmdata[i]==0.0f){
              frmdata[i]=0.0f;
          }
      }
    return(EXIT_SUCCESS);
}
int performCorrectionIM(float *frmsigdata, float *frmexpdata, float *badpixarry,int sizex,int sizey,double  intgrntime)
{
    float  temp_value=0.0f;
   
    for (int pixno = 0 ; pixno <sizex*sizey  ; pixno++)
     {
         temp_value=frmsigdata[pixno];
         if(frmsigdata[pixno]!=INVALID_PIX_VALUE && temp_value!=0.0f  ){
            frmsigdata[pixno] = frmsigdata[pixno] * badpixarry[pixno] ;
            if(frmsigdata[pixno]==0 )
            { 
               
                frmsigdata[pixno]=INVALID_PIX_VALUE;
            }
         }
         else if(temp_value==0.0f  && badpixarry[pixno]==1)
         {
            frmsigdata[pixno]=0.0f;
              //LOG(INFO)<< (pixno+1)%512<<" "<<(pixno+1)/512<<temp_value;
         }
         else if (temp_value==0.0f && badpixarry[pixno]==0)
         {
             frmsigdata[pixno]=INVALID_PIX_VALUE;
         }
     }
   
        for (int pix = 0 ; pix < sizex * sizey; pix++)
        {
            frmexpdata[pix] = INVALID_PIX_VALUE ;
            frmexpdata[pix] = badpixarry[pix] * intgrntime ;
        }
    return(EXIT_SUCCESS);
}
int performCosmicRayCorrIM(float *frmsigdata, float *frmexpdata,int sizex,int sizey,float threshold_cr)
{
    vector<int> xpix,ypix;
    int cnt_cosmicAffected = 0 ;
     for (int pixno = 0 ; pixno <sizex*sizey ; pixno++)
        {
         if(frmsigdata[pixno]!=INVALID_PIX_VALUE){
           
             if (frmsigdata[pixno] > threshold_cr)
            {
                 xpix.push_back (pixno % sizey+1) ;
                 ypix.push_back (pixno / sizex+1) ;
                //frmsigdata[pixno] = frmexpdata[pixno] = 0.0 ;
            }
         }
            
        }    
      for (int i = 0 ; i < xpix.size () ; i ++)
        {
            cnt_cosmicAffected = 0 ;

            for (int j = xpix[i] - 3/2 ; j <= xpix[i] +3/2 ; j ++)
            {
                for (int k = ypix[i] -3/2 ; k <= ypix[i] +3/2 ; k ++)
                {

                    if (j < sizex && j > 0 && k < sizey && k > 0)
                    {
                        if (frmsigdata[k * sizey + j] > threshold_cr)
                        {
                            cnt_cosmicAffected ++ ;
                        }

                    }

                }

            }
             if (cnt_cosmicAffected == 1)
            {
                frmexpdata[ypix[i] * sizex + xpix[i]] = frmsigdata[ypix[i] * sizey + xpix[i]] = -9999 ; //cr effected
                //cr_effected.push_back (frameno) ;
               // x_crFailed.push_back (X_pixel[i]) ;
                //y_crFailed.push_back (Y_pixel[i]) ;

            }
       }  
    
    return(EXIT_SUCCESS);
}
int performSubDivisionIM(float *frmsigdata,int sizex,int sizey,float *subdivideddata,int size_subdivx,int size_subdivy)
{
    
      if (sizex == 0 || sizey == 0)
    {
        LOG(INFO) << "Divide by Zero" << endl ;
        return (EXIT_FAILURE) ;
    }
    int xfactor = size_subdivx / sizex ;
    int yfactor = size_subdivy / sizey ;
    for (int i = 0 ; i < size_subdivx ; i++)
    {
        for (int j = 0 ; j < size_subdivx ; j++)
        {
//            subdivideddata[i * size_subdivx + j] = frmsigdata[(int)(round(i / yfactor) * sizey + round(j / xfactor))] ;
             subdivideddata[i * size_subdivx + j] = frmsigdata[((int)(i / yfactor) * sizex + (int)(j / xfactor))] ;
        }
    }
    return(EXIT_SUCCESS);
}
//int  writeOutputTblToDisk (char *id , char *outDir , char *dir , char *subscript  , char *namepre , double ftime , unsigned short fno ,char **type1,char**tform1,int tfields1,char *extname,vector<float> &X ,vector<float> &Y,vector<float> &val)
//{
//    char out_file[FLEN_FILENAME];
//    sprintf (out_file , "%s/%s/%s_t%f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , id,subscript ) ;
//    if(strcmp (id,"rf")==0){
//        name_track.push_back (basename(out_file));
//    }
//                fitsfile *fout1 ;
//                int status=0;
//                fits_create_file (&fout1 , out_file , &status) ;
//                printError (status , "Error creating the output file ",out_file) ;
//                fits_create_tbl (fout1 , BINARY_TBL , 0 , tfields1 , type1,tform1 ,NULL,extname , &status) ;
//                printError (status , "Error in creating the table",out_file) ;
//               fits_update_key (fout1 , TDOUBLE , "FRMTIME" , &ftime , " Average Frame time" , &status) ;
//               printError (status , "Error in updating the key value of the FRMTIME",out_file) ;
//                fits_write_col (fout1 , TFLOAT , 1 , 1 , 1 , X.size () , X.data () , &status) ;
//                printError (status , "Error in writing  column",out_file) ;
//                fits_write_col (fout1 , TFLOAT , 2 , 1 , 1 , Y.size () , Y.data () , &status) ;
//                printError (status , "Error in writing the column",out_file) ;
//                fits_write_col (fout1 , TFLOAT , 3 , 1 , 1 , val.size () , val.data () , &status) ;
//                printError (status , "Error in writing the column",out_file) ;
//                fits_close_file(fout1,&status);                
//                return (EXIT_SUCCESS) ;
//
//}
int readImage(char * caldb_file,int hduno,float *frm_data ,int xsize,int ysize)
{
    
    fitsfile *fptr;
    int status=0;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fits_open_file (&fptr , caldb_file , READONLY , &status) ;
    printError (status , "Error in readFlatField()" , caldb_file) ;
   fits_movabs_hdu (fptr , hduno , NULL , &status) ;
    printError (status , "Error in  moving to the 2nd HDU of the out information file" , caldb_file) ;
    for (int q = 0 ; q < xsize * ysize ; q++) frm_data[q] = 0.0 ;
    fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , frm_data , NULL , &status) ;
    printError (status , "Error in reading the pixels from caldb file" , caldb_file) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file" , caldb_file) ;
    return(EXIT_SUCCESS);
}

int darkFrameSubtraction (float *Array , float *frame_data,int xsize,int ysize)
{
  
  
 // vector<float> ::iterator it=min_element (min_value.begin (),min_value.end ());
 //;exit(1);
// int loc=distance (min_value.begin (),it);
 
  
  for (int p = 0 ; p < xsize * ysize ; p++)
    {
        if(frame_data[p]!=INVALID_PIX_VALUE ){
//            LOG(INFO)<<frame_data[p]<<" "<<Array[p];
           
        frame_data[p] = (frame_data[p] - Array[p]);
        
        }
    }

      vector<float> min_value;
       for(int i=0;i<xsize*ysize;i++)
       {
           if(frame_data[i]<0.0f && frame_data[i]!=INVALID_PIX_VALUE) 
           {
            //   LOG(INFO)<<frame_data[i]<<" "<<i;
               min_value.push_back (frame_data[i]);
           }
       }
    
 float min_value_final=0.0f;
  for (int i=0;i<min_value.size ();i++)
  {
    if(min_value[i]<min_value_final)
    {
        min_value_final=min_value[i];
    }        
  }
// min_value_final=0.0;
//LOG(INFO)<<min_value_final;
 
 for (int p = 0 ; p < xsize * ysize ; p++)
    {
//       if(frame_data[p]!=INVALID_PIX_VALUE )
//        frame_data[p] = frame_data[p] + abs (min_value_final)*MIN_MULTNO;
      
    }
    return (EXIT_SUCCESS) ;


}
int performFlatFieldCorrIM(float *frmsigdata,float *flatfieldarry,int sizex,int sizey)
{
for (int pixno = 0 ; pixno <sizex*sizey  ; pixno++)
     { 
         if(frmsigdata[pixno]!=INVALID_PIX_VALUE){
            frmsigdata[pixno] = frmsigdata[pixno] * flatfieldarry[pixno] ;
         }
     }
  
    return(EXIT_SUCCESS);
}
 int  performQEMCPcorrection(float *sigdata,int xsize,int ysize,double fact)
 {
     for(int i=0;i<xsize*ysize;i++){
         if(sigdata[i]!=INVALID_PIX_VALUE){
         sigdata[i]=sigdata[i]*fact;
         }
     }
     
     
     return(EXIT_SUCCESS);
 }
 int readpathOfRASfile(string sciencedataFile,char *channel,string &output_RasFile_Dir,string &outputrasFile)
 {
     char* temp_ch = new char[100];
     if(strcmp (channel,"NUV")==0)
     {
         strcpy(temp_ch,"/uvtN");
     }
     else if(strcmp (channel,"FUV")==0){
         strcpy (temp_ch,"/uvtF");
     }
     else if(strcmp (channel,"VIS")==0){
         strcpy (temp_ch,"/uvtV");
     }
     else {
         cerr<<"***Invalid channel value*** ";
     }
      
     int pos=sciencedataFile.find (temp_ch);
     LOG(INFO)<<pos;
     string lengh=sciencedataFile.substr (pos-6,20);
  
     string final_str_rasfile=output_RasFile_Dir+"/uvit/"+lengh+"/uvtComputeDrift_"+VERSION+"/";   
    // LOG(INFO)<<"ooo"<<final_str_rasfile;
     string new_name=searchFile((char*)final_str_rasfile.c_str (),".fits");
     
      outputrasFile=final_str_rasfile+"/"+new_name;
      cout<<"The Ras file used is "<<outputrasFile;
      //pos = sciencedataFile.find (temp_ch) ;
     //string tar_date_data=sciencedataFile.substr (pos+10,3);
     
    //rasFilepath=temp_ch+(string)tar_date_data; 
     
    //rasFilepath=output_RasFile_Dir+rasFilepath;
//    pos = sciencedataFile.find (rasFilepath) ;
//    int pos2=sciencedata   tar_date_data=sciencedataFile.substr (pos2+1,pos2+23);File.find ("/");
//    tar_date_data=sciencedataFile.substr (pos2+1,pos2+23);
//    cout<<"tar_data "<<tar_date_data<<" "<<rasFilepath;exit(1);
//   tar_date_data=tar_date_data+(string)"uvtComputeDrift_1.0";
//   tar_date_data=output_RasFile_Dir+"/"+tar_date_data;
//    cout<<tar_date_data;exit(1);
//    const char * new_name=searchFile((char*)tar_date_data.c_str (),".fits");
//    rasFilepath=tar_date_data+"/"+new_name;//exit(1)   ;
//    pos=rasFilepath.find ("/");
//    string outStr=rasFilepath.substr (pos,rasFilepath.size ());
//    rasFilepath=output_Dir+outStr;
//    cout<<"the output Dir "<<outStr;
    return(EXIT_SUCCESS);    
 }
 int  checkMasterClock(char *filename  ,string &clockMaster,int &cnt_VIS,int &cnt_FUV,int &cnt_NUV,long &cnt_Tot){
     fitsfile *fptr ;
     int status=0;
     long numrows=0;
     long felement=1;
     long nelements=1;
     int tot_bytes=48;
     unsigned char *Moniter_Data;
     fits_open_file (&fptr , filename , READONLY , &status) ;
    printError (status , "Error opening data ingest output file" ,  filename ) ;
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU in input information File" ,   filename) ;

    fits_get_num_rows (fptr , &numrows , &status) ;
    printError (status , "Error in getting number of rows" ,   filename) ;
    nelements=numrows;
   // cout<<numrows<<endl;
    cnt_VIS=0;cnt_FUV=0;cnt_Tot=0,cnt_NUV=0;
    cnt_Tot=nelements;
    Moniter_Data= new unsigned char[nelements*tot_bytes];
      fits_read_col (fptr , TBYTE , 15 , 1, felement , nelements*tot_bytes , NULL , Moniter_Data , NULL , &status) ;
      printError (status , "Error in reading the column of the centroids in inputFile" , filename) ;
     
      for (int i=0;i<nelements;i++){
          if(((Moniter_Data[i*tot_bytes+15]>>4)&3)==0) cnt_FUV++;
          else if(((Moniter_Data[i*tot_bytes+15]>>4)&3)==1) cnt_NUV++;
          else if(((Moniter_Data[i*tot_bytes+15]>>4)&3)==2) cnt_VIS++;
      }
      for (int i=1;i<nelements;i++){
          if(Moniter_Data[i*tot_bytes+15]!=Moniter_Data[(i-1)*tot_bytes+15])
          {
              LOG(ERROR)<<"***Not Continuous clock master !!!! ***";
              break;
              
          }
          
      }
      
      int  clock_Master=(Moniter_Data[15]>>4) & 3;
      if(clock_Master==0) clockMaster="FUV";
      else if(clock_Master==1) clockMaster="NUV";
      else if(clock_Master==2) clockMaster="VIS";
      
    
    fits_close_file(fptr,&status);
     
     return(EXIT_SUCCESS);
 }
 //int ReadAtd_status(char * filename ,vector<double> &list_time_att,int &val_att_flag,int &att_bit_0,int &att_bit_1,int &att_bit_2)
 int ReadAtd_status(char * filename ,double &frmTime,int &att_bit_0,int &att_bit_1,int &att_bit_2)
 {
     
     //list_time_att.clear ();
    // int valid_flag;
  //   if(val_att_flag==1) valid_flag=1;
    // else if(val_att_flag==2) valid_flag=0;
     long numrows=0;
     long felement=1;
   
     int tot_bytes=8;
     int status=0;
     fitsfile *fptr;
    fits_open_file (&fptr , filename , READONLY , &status) ;
    printError (status , "Error opening data ingest output file" ,  filename ) ;
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU in input information File" ,   filename) ;
    fits_get_num_rows (fptr , &numrows , &status) ;
    printError (status , "Error in getting number of rows" ,   filename) ;
  
      unsigned  char* atd_status= new unsigned char[numrows*tot_bytes];
      double  *time_attFile = new double [numrows];
      fits_read_col (fptr , TDOUBLE, 1 , 1, felement , numrows , NULL ,time_attFile , NULL , &status) ;
      printError (status , "Error in reading the column of the centroids in inputFile" , filename) ;
      fits_read_col (fptr , TBIT, 6 , 1, felement , numrows*tot_bytes , NULL ,atd_status , NULL , &status) ;
      printError (status , "Error in reading the column of the centroids in inputFile" , filename) ;
     
      //LOG(INFO)<<" ATD time "<<setprecision (20)<<time_attFile[0]<<" "<<setprecision (20)<<time_attFile[numrows-1];
     
      //added
      for (int i =0;i<numrows-1;i++)
      {
          if(frmTime >time_attFile[i]  &&  frmTime <= time_attFile[i+1])
          {
              att_bit_0=atd_status[i*8];
              att_bit_1=atd_status[i*8+1];
              att_bit_2=atd_status[i*8+2];
              //LOG(INFO)<<att_bit_0<<" "<<att_bit_1<<" "<<att_bit_2;exit(1);
              //to be return
              break;
           }
          
      }
      
      //till this
//      for (int i=0;i<numrows;i++)
//      {
//          if(atd_status[i*8]==1 )
//          {
//              if(atd_status[i*8+1]==valid_flag || atd_status[i*8+2]==valid_flag)
//              {
//                  list_time_att.push_back (time_attFile[i]);                  
//              }             
//          }          
//      }
//      
      fits_close_file(fptr,&status); 
      printError (status , "Error in closing the file" , filename) ;
      
      return(EXIT_SUCCESS);
 }
 
 
 string convertIntToStr(int val)
 {
     char temp[PIL_LINESIZE];
     sprintf(temp,"%d",val);
     
      return  (string)temp;
 }
 string convertFloatToStr(float val)
 {
     char temp[PIL_LINESIZE];
     sprintf(temp,"%f",val);
     
      return  (string)temp;
 }
 int editframe(float *frmdata,int height,int width,float thrval,bool &flag_All,int colnum){
     flag_All=FALSE;

      vector<int> xlocations_counter,ylocation_counter;
      int cnt=0;
      int cnt_x=0,cnt_y=0;
  
float *tempArray= new float[height*width];
for(int i=0;i<height*width;i++){
tempArray[i]=frmdata[i];
}
//logic for stripe removal
     for (int i=0;i<height;i++)
     {
         
        
         for(int j=0;j<width;j++)
         {
             cnt_x=0;
             cnt_y=0;
             
             if(frmdata[i*width+j]>thrval){
                   tempArray[i*width+j]=INVALID_PIX_VALUE;
//                
                
                 
                 
             }
             
             
         }
                 
     }


//Addition check for saturated stars and handling
bool flag_Saturated=FALSE;
for(int i=0;i<width;i++)
{
   for (int j=0;j<height;j++)
    {
       flag_Saturated=FALSE;
       if((j*height+i)> 0 && (j*height+i)<height*width)
       {
           if (frmdata[j * height + i] > SATURATED_PIXEL_THR && tempArray[j*height+i] ==INVALID_PIX_VALUE) {

               //check for neighboring pixels
               if(frmdata[(j-1)*height+i]>SATURATED_NEIGHBOUR_PIXEL_THR  && frmdata[(j+1)*height+i]>SATURATED_NEIGHBOUR_PIXEL_THR
                       && frmdata[(j)*height+(i-1)]>SATURATED_NEIGHBOUR_PIXEL_THR && frmdata[(j)*height+(i+1)]>SATURATED_NEIGHBOUR_PIXEL_THR)
               {
               flag_Saturated=TRUE;    
               
               }
               

                }
       }
       if(flag_Saturated==TRUE){
           //flag all 17X17 pixel around this pix
                for (int index_col = i - SATURATED_PIX_ENVELOPESIZE; index_col <= i + SATURATED_PIX_ENVELOPESIZE; index_col++) {
                    for (int index_row = j - SATURATED_PIX_ENVELOPESIZE; index_row <= j + SATURATED_PIX_ENVELOPESIZE; index_row++) {
                        if((index_row*height+index_col)> 0 && (index_row*height+index_col)<height*width)
                        {
                            
                     tempArray[index_row*height+index_col]=INVALID_PIX_VALUE;
                        
                        }
                }
                
           
           
           }
           
           
       }
       
        
    }

}

//finding star like thing

//till this
for(int i=0;i<height*width;i++){
frmdata[i]=tempArray[i];
}

delete[] tempArray;
     
//      flag_All=FALSE;
//
//      vector<int> xlocations_counter,ylocation_counter;
//      int cnt=0;
//      int cnt_x=0,cnt_y=0;
//  
//float *tempArray= new float[height*width];
//for(int i=0;i<height*width;i++){
//tempArray[i]=frmdata[i];
//}
//     for (int i=0;i<height;i++)
//     {       
//        
//         for(int j=0;j<width;j++)
//         {
//             cnt_x=0;
//             cnt_y=0;
//             
//             if(frmdata[i*width+j]>thrval)
//             {
//                     flag_All=TRUE;
//                     tempArray[i*width+j]=INVALID_PIX_VALUE;                
//                 
//             }
//             
//             
//         }
//                 
//     }
//
////taking mean after all higher intensity pixels are removed
//float meanImage=getmean(tempArray,width*height);
//float sdImage=getSD(tempArray,width*height);
//int cnt_Window=0;
//int totalPixFoundenvelope=0;
//
//float thrVal=meanImage+2.5*sdImage;
//
//LOG(INFO)<<meanImage<<" "<<sdImage;
//float finalWindow=8;
//float *FinalArray = new float[height*width];
//for (int i =0;i<height*width;i++)
//{
//FinalArray[i]=tempArray[i];    
//}
//for (int i=0;i<width;i++)
//{
//    for (int j=0;j<height;j++)
//    {
//        cnt_Window=0;
//        if(tempArray[i*width+j]==INVALID_PIX_VALUE){
////           do{
////                cnt_Window++;
////                //method call
////                totalPixFoundenvelope=envelopeCheckForBrightPix(tempArray,height,width,j,i,cnt_Window,thrVal);
////                //LOG(INFO)<<totalPixFoundenvelope<<" "<<i<<" "<<j;
////                if(cnt_Window==6) break; //maximum window size 17X17
////            }while(totalPixFoundenvelope>0);
////            
////             
////            //now That whole window to be marked as a INVALID.
////            finalWindow=cnt_Window-1;
//              for (int i_index=j-finalWindow;i_index<=j+finalWindow;i_index++)
//     {
//         for (int j_index=i-finalWindow;j_index<=i+finalWindow;j_index++)
//         {
//             if((j_index*height+i_index)> 0 && (j_index*width+i_index) <height*width)
//             {
////                 FinalArray[j_index*height+i_index]=INVALID_PIX_VALUE;
//                  tempArray[j_index*height+i_index]=-99999;
//             }
//             
//             
//         }
//     
//     }
//      
//        }
//    }
//    
//}
//
//
//for(int i=0;i<height*width;i++)
//{
//    if(tempArray[i]==-99999){
//        frmdata[i]=INVALID_PIX_VALUE;
//    }
//    else{
//    frmdata[i]=tempArray[i];
//    }
//    //frmdata[i]=FinalArray[i];
//}
//
//delete[] tempArray;
//     
     
    
     
     return(EXIT_SUCCESS);
 }

 int envelopeCheckForBrightPix(float *Array,int XSIZE,int YSIZE,int Xloc,int Yloc,int window,float thrval){
     int cntpix=0;
     for (int i=Xloc-window;i<=Xloc+window;i++)
     {
         for (int j=Yloc-window;j<=Yloc+window;j++)
         {
             if((j*XSIZE+i)> 0 && (j*XSIZE+i) <XSIZE*YSIZE)
             {
                 if((i==(Xloc-window) || i == (Xloc+window) ) || (j==(Yloc-window) || j== (Yloc+window) )){
                 if(Array[j*XSIZE+i]>thrval && i!=Xloc && j !=Yloc)
                 {
                     cntpix++;
                     break;
                 }
                 }
             }
             
             
         }
     
     }
     
     
     
     
     return cntpix;
 }
int compressTars(string tar_name,string directoryname,string pathtar)
{
//    cout<<sizeof(tar_name);exit(1);
   string fileoutname;
string trgfilename;
    string tarname=(string)basename(tar_name.c_str());
    string strtemp (tarname) ;
     int pos = strtemp.find_last_of ("/") ;
    string tempFilename="";
    string cmd;
    //cmd.repla
     string tar_date_data, tar_OBS_id_data;
     string tar_nameprfx,ver_tar;
    try{
        tar_nameprfx=tarname.substr (pos+1,46);
        int index1=tar_nameprfx.find ("LEVL1",0);
        if(index1==std::string::npos) {
          LOG(INFO)<<"Error";
            return(EXIT_FAILURE);
        }
     tar_nameprfx.replace (index1,5,"LEVL2");
        //tar_nameprfx.replace ("LEVL1","LEVL2");
        ver_tar=tarname.substr (52,3);
        //cout<<"Output tar name "<<tar_nameprfx<<" yy: "<<ver_tar;
        fileoutname=pathtar+"/"+tar_nameprfx+".tar_V"+ver_tar;
        trgfilename=pathtar+"/"+tar_nameprfx+".trg_V"+ver_tar;
        //cout<<fileoutname;exit(1);
    //tar_date_data=tarname.substr (pos+1+11,8);
     //tar_OBS_id_data=tarname.substr (pos+1+19,21);
    //orbitnum=tarname.substr(pos+1+41,5);
    }
    catch(...){
        LOG(ERROR)<<"***Level-1 Tar file is not in a proper format.***";
        return(EXIT_FAILURE);
    }
    // LOG(INFO)<<"The orbit num is "<<orbitnum;exit(1);
     try{
   cmd="tar -cvf  "+fileoutname+" "+directoryname;
       
    LOG(INFO)<<"Executing.. "<<cmd;

    system(cmd.c_str ());
     }
     catch(...){
         LOG(INFO)<<"Error in compressing the directory";
         return(EXIT_FAILURE);
     }
    
     FILE *ptr=NULL;
    ptr=fopen(fileoutname.c_str (),"rb");
    fseek(ptr,0,SEEK_END);
    int size=ftell(ptr);
    fclose(ptr);
   //cout<<size;exit(1);
    ofstream out1(trgfilename.c_str ());
    out1<<size;
    out1.close ();
    return(EXIT_SUCCESS);
}

int compressZip(string tar_name,string directoryname,string pathtar)
{
//    cout<<sizeof(tar_name);exit(1);
   string fileoutname;
string trgfilename;
    string tarname=(string)basename(tar_name.c_str());
    string strtemp (tarname) ;
     int pos = strtemp.find_last_of ("/") ;
    string tempFilename="";
    string cmd;
    //cmd.repla
     string tar_date_data, tar_OBS_id_data;
     string tar_nameprfx,ver_tar;
    try{
        tar_nameprfx=tarname.substr (pos+1,46);
        int index1=tar_nameprfx.find ("LEVL1",0);
        if(index1==std::string::npos) {
          LOG(INFO)<<"Error";
            return(EXIT_FAILURE);
        }
     tar_nameprfx.replace (index1,5,"LEVL2");
        //tar_nameprfx.replace ("LEVL1","LEVL2");
        //ver_tar=tarname.substr (52,3);
        //cout<<"Output tar name "<<tar_nameprfx<<" yy: "<<ver_tar;
        fileoutname=pathtar+"/"+tar_nameprfx+".zip";//+ver_tar;
        trgfilename=pathtar+"/"+tar_nameprfx+".trg";//+ver_tar;
        //cout<<fileoutname;exit(1);
    //tar_date_data=tarname.substr (pos+1+11,8);
     //tar_OBS_id_data=tarname.substr (pos+1+19,21);
    //orbitnum=tarname.substr(pos+1+41,5);
    }
    catch(...){
        LOG(ERROR)<<"***Level-1 Tar file is not in a proper format.***";
        return(EXIT_FAILURE);
    }
    // LOG(INFO)<<"The orbit num is "<<orbitnum;exit(1);
     try{
   cmd="zip -r  "+fileoutname+" "+directoryname;
       
    LOG(INFO)<<"Executing.. "<<cmd;

    system(cmd.c_str ());
     }
     catch(...){
         LOG(INFO)<<"Error in compressing the directory";
         return(EXIT_FAILURE);
     }
    
     FILE *ptr=NULL;
    ptr=fopen(fileoutname.c_str (),"rb");
    fseek(ptr,0,SEEK_END);
    int size=ftell(ptr);
    fclose(ptr);
   //cout<<size;exit(1);
    ofstream out1(trgfilename.c_str ());
    out1<<size;
    out1.close ();
    return(EXIT_SUCCESS);
}


int   getNearByValues(double *arr, double *arr_output,long size,float thrval)
{
    
 double number;
 double ** diffrence_map= new double*[size*size];
for (int i=0;i<size;i++)
{
diffrence_map[i]=new double[size];    
}
 for(int i=0;i<size;i++){
     for (int j=0;j<size;j++)
     {
         diffrence_map[i][j]=-9999;
     }
 }
 for (int i=0;i<size;i++)
 {
     arr_output[i]=-9999;
 }
 double diff_arr;
 for (int i=0;i<size;i++)
 {
     for(int j=i+1;j<size;j++)
     {
             diff_arr=arr[j]-arr[i];
             diffrence_map[i][j]=abs(diff_arr);
             //cout<<diffrence_map[i][j]<<" "<<i<<" "<<j<<endl;
         
         
     }    
 }
 int cnt_Forvalidmatch=0;
 vector<int> cnt_track;
 vector<int> validStar_indices;
 for(int i=0;i<size;i++)
 {
     cnt_Forvalidmatch=0;
     for(int j=0;j<size;j++)
     {
     if(diffrence_map[i][j]!=-9999)
     {
         if(diffrence_map[i][j]<thrval)
         {
             cnt_Forvalidmatch++;
             
            // arr_output[i]=1;
            // arr_output[j]=1;
         }
     }
     }
     cnt_track.push_back(cnt_Forvalidmatch);
//     LOG(INFO)<<cnt_Forvalidmatch;
//     if(cnt_Forvalidmatch>1){
//         validStar_indices.push_back(i);
//     }
//     else validStar_indices.push_back(-9999);
 }
 int max_val=0;
 int indices_ForSearch=0;
 for(int i=0;i<cnt_track.size();i++)
 {
     if(cnt_track[i]>max_val){
         indices_ForSearch=i;
         max_val=cnt_track[i];
     }     
 }
 
 float FinalValToCompare=arr[indices_ForSearch];
 for(int i=0;i<size;i++){
     if(abs(FinalValToCompare - arr[i])<thrval){
          arr_output[i]=1;
     }
 }
 
 

 
// for (int i=0;i<validStar_indices.size();i++)
// {
//     arr_output[i]=0;
//     if(validStar_indices[i]!=-9999){
//         arr_output[i]=1;
//     }
//     else arr_output[i]=0;   
// }
    
return 0;
}
//long checkAvailableSpace(const char *path){
//    struct statvfs sf1;
//    if(!DirExists((char*)path)){
//        cout<<"***Path not exist***,Path -> "<<path<<endl;
//        exit(1);
//    }
//   if (statvfs(path, &sf1) != 0) {
//    // error happens, just quits here
//       cerr<<"***Something went wrong for calculating the available space***"<<endl;
//       
//       exit(1);
//  }
//    return sf1.f_bsize*sf1.f_bavail<<endl;
//    
//    
//}
const char* searchFile_RecDir (char* dirname , char* pattern)
{
    //LOG(INFO)<<endl<<"Dirname :"<<dirname<<endl;
    DIR *dp ;
    char **sub ;
    char *fitname = pattern ;
    string fitfile = "" ;
    sub = new char* [100] ;
    for (int i = 0 ; i < 100 ; i++)
        sub[i] = new char[30] ;
    int p = 0 ;
    int count = 0 ;
    struct dirent *dirp ;
    if ((dp = opendir ((const char *)dirname) )== NULL)
   {  LOG(ERROR) << endl << "***Error in opening " << dirname << " directory***" << endl ;
        exit (EXIT_FAILURE) ;
    }
    while (dirp = readdir (dp))
    {
        p++ ;
        if (dirp->d_name == " . " || dirp->d_name == " .. ")
            continue ;
        sub[p - 1] = dirp->d_name ;
        // LOG(ERROR)<<"\nthe sub"<<sub[p]<<endl;
        string temp = (string) sub[p - 1] ;
        //LOG(ERROR)<<"\nthe temp"<<temp<<endl;
        int size_temp = temp.size () ;
        //  LOG(ERROR)<<"\nthe fits name"<<fitname;
        int npos = temp.find (fitname , 0) ;
        if (npos > 0 && npos < size_temp)
        {

            fitfile = temp ;
            //  exit(1);
            break ;
           // count++ ;
        }
        //LOG(ERROR)<<"\n"<<"the sub "<< fitfile<<"  position"<<npos<<endl;
        //  p++;

    }
    closedir (dp) ;
    return fitfile.c_str () ;
}
int Convert_J2000_To_ObervationRADECVal(float dayRef,float RA_orig,float DEC_orig,float &RAdiff,float &DECdiff)
{
//     float m = (3.07496+0.00186*fract_Year);
//    float n = (1.33621-0.00057*fract_Year);
//    
//    RAdiff=(m+n*sin(RA_orig*M_PI/180)*tan(DEC_orig*M_PI/180))/240;//convert Seconds to Degree.
//    DECdiff=(n*cos(RA_orig*M_PI/180))/240;
//    
    
    
    
    RAdiff=(3.075+1.336*sin(RA_orig*M_PI/180)*tan(DEC_orig*M_PI/180));//convert Seconds to Degree.
    RAdiff=(0.0000114077)*RAdiff;
    DECdiff=(0.0000152407)*cos(RA_orig*M_PI/180);;
    float RAmid =RA_orig+(RAdiff*dayRef)/2;
    float DECmid=DEC_orig+(DECdiff*dayRef)/2;
    RAdiff=3.075+1.336*sin(RAmid*M_PI/180)*tan(DECmid*M_PI/180);
    RAdiff=(0.0000114077)*RAdiff;
    DECdiff=(0.0000152407)*cos(RAmid*M_PI/180);
    
    return EXIT_SUCCESS;
}
int ConvertDateToMJD(float day,float  month ,float year){
    float yearp,monthp,B,A,C,D,jd;
    if (month == 1 || month == 2){
yearp = year - 1;
monthp = month + 12;
    }
else{
yearp = year;
monthp = month;
}

if ((year < 1582) or
(year == 1582 and month < 10) or
(year == 1582 and month == 10 and day < 15)){

B = 0;
}
else{

A = round(yearp / 100.0);
B = 2 - A + round(A / 4.0);
}
    
if (yearp < 0){
C = round((365.25 * yearp) - 0.75);
}
    
else{
C = round(365.25 * yearp);
}
    
D = round(30.6001 * (monthp + 1));
        
jd = B + C + D + day + 1720994.5;

return (jd-2400000.5);
   
}

string  get_process_name_by_pid(const int pid)
{
    char* name = (char*)calloc(256,sizeof(char));
    if(name){
        sprintf(name, "/proc/%d/cmdline",pid);
        FILE* f = fopen(name,"r");
        if(f){
            size_t size = fread(name, sizeof(char), 256, f);
            if(size>0){
                if('\n'==name[size-1])
                    name[size-1]='\0';
            }
            fclose(f);
        }
    }
    string temp=name;
    return temp;
}

int ApplyBinning (float* inputarray , int in_xsize , int in_ysize , float* outputarray , int out_xsize , int out_ysize,float *counterArray)
{
    
    if(in_xsize < out_xsize){
        LOG (ERROR) << "***Input Array size is LESS-THAN output Array size***" ;
        return (EXIT_FAILURE) ;
    }
    if (in_xsize % out_xsize != 0)
    {
        LOG (ERROR) << "Can't sub sampling with  this input array ,input Array size is not matching" ;
        return (EXIT_FAILURE) ;
    }
    float devide_term_x = in_xsize / out_xsize ;
    float devide_term_y = in_ysize / out_ysize ;
 //   LOG(INFO)<<devide_term_x<<" "<<devide_term_y<<" "<<in_xsize<<" "<<in_ysize<<" "<<out_xsize<<" "<<out_ysize;exit(1);
   
    int cnt_win = 0 ;
    int cnt_nonZeroPix=0;
    double sum_win = 0.0f ;
    int index_finalArray = 0 ;
    LOG(INFO)<<"Entered binning";
    
    for (int temp_y = 0 ; temp_y < in_ysize ; temp_y = temp_y + devide_term_y)
    {
        for (int temp_x = 0 ; temp_x < in_xsize ; temp_x = temp_x + devide_term_x)
        {
            if (index_finalArray > out_xsize * out_ysize)
            {
                LOG (ERROR) << "Array is out of bound, EXEED to " << out_xsize << " !!!" ;
                return (EXIT_FAILURE) ;
            }
            cnt_win = 0 ;
            sum_win = 0.0f ;
            cnt_nonZeroPix=0;
            for (int i = temp_y ; i < temp_y + devide_term_y ; i ++)
            {
                for (int j = temp_x ; j < temp_x + devide_term_x ; j ++)
                {
                    if (inputarray[i * in_xsize + j] != INVALID_PIX_VALUE)
                    {
//                        if(counterArray[j * in_xsize + i] !=0){
//                            cnt_nonZeroPix=cnt_nonZeroPix+counterArray[j * in_xsize + i];
//                        }
                        sum_win = sum_win + inputarray[i * in_xsize + j] ;
                        cnt_win ++ ;
                    }

                }
            }
            if (cnt_win != 0 )
            {
               // if(cnt_nonZeroPix==(devide_term_x*devide_term_y)-1)
                   
                //else
                  //  outputarray[index_finalArray ++] = (float) (sum_win);// / cnt_win) ;
                   // if(cnt_nonZeroPix==0){
                     //   outputarray[index_finalArray ++]=0.0f;
                   // }
                  //  else{
//                         outputarray[index_finalArray ++] = (float)((((devide_term_x*devide_term_y)/(cnt_nonZeroPix))*(sum_win)));
outputarray[index_finalArray ++] = (float)((sum_win));  
                //  }
            }
            else
            {
                outputarray[index_finalArray ++] = INVALID_PIX_VALUE ;
            }


        }

    }
    LOG(INFO)<<"EXITING "<<index_finalArray;
    return (EXIT_SUCCESS) ;

}
 void RemoveDuplicateStringsFromVECT(vector<string> &strvect){
     vector<string> temp_copy;
    //temp_copy= windows;
    bool flag_Dupli=0;
    for (int i=0;i<strvect.size();i++){
        flag_Dupli=0;
        for (int j=i+1;j<strvect.size();j++){
            if(strcmp(strvect[i].c_str(),strvect[j].c_str())==0){
                flag_Dupli=1;
            }
        }
        if(flag_Dupli==0){
            temp_copy.push_back(strvect[i]);
        }
    }
    strvect=temp_copy;
     
 }
 
 int ApplySubSampling_Addition (float* inputarray , int in_xsize , int in_ysize , float* outputarray , int out_xsize , int out_ysize)
{

    if (in_xsize % out_xsize != 0)
    {
        LOG (ERROR) << "Can't sub sampling with  this input array ,input Array size is not matching" ;
        return (EXIT_FAILURE) ;
    }
    float devide_term_x = in_xsize / out_xsize ;
    float devide_term_y = in_ysize / out_ysize ;
    int cnt_win = 0 ;
    float sum_win = 0.0f ;
    int index_finalArray = 0 ;
    for (int temp_x = 0 ; temp_x < in_xsize ; temp_x = temp_x + devide_term_x)
    {
        for (int temp_y = 0 ; temp_y < in_ysize ; temp_y = temp_y + devide_term_y)
        {
            if (index_finalArray > out_xsize * out_ysize)
            {
                LOG (ERROR) << "Array is out of bound, EXEED to " << out_xsize << " !!!" ;
                return (EXIT_FAILURE) ;
            }
            cnt_win = 0 ;
            sum_win = 0.0f ;
            for (int i = temp_y ; i < temp_y + devide_term_y ; i ++)
            {
                for (int j = temp_x ; j < temp_x + devide_term_x ; j ++)
                {
                    if (inputarray[j * in_xsize + i] != INVALID_PIX_VALUE)
                    {
                        sum_win = sum_win + inputarray[j * in_xsize + i] ;
                        cnt_win ++ ;
                    }

                }
            }
            if (cnt_win != 0)
            {
                outputarray[index_finalArray ++] = (float) (sum_win) ;
            }
            else
            {
                outputarray[index_finalArray ++] = INVALID_PIX_VALUE ;
            }


        }

    }

    return (EXIT_SUCCESS) ;

}
