
#include<Attitude.h>
#include<fitsio.h>
#include<uvtUtils.h>
#include<glog/logging.h>


void Attitude::display(){
    LOG(INFO)<<time<<"  "<<q1<<"  "<<q2<<"   "<<q3<<"   "<<q4;
}

int readAttitude(char *attitudefile,char *timecol,char *qcol,double tstart,double tstop,Attitude &att){
    
    LOG(INFO)<<"Reading attitude data from "<<attitudefile;
    
    int status=0;
    
    fitsfile *fatt;
    fits_open_file(&fatt,attitudefile,READONLY,&status);
    printError(status,"Error opening attitude file",attitudefile);
    
    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
    printError(status,"Error in moving to 2nd HDU ", attitudefile);
    
    long nrows;
    fits_get_num_rows(fatt,&nrows,&status);
    printError(status,"Error reading number of rows from attitude file  ",attitudefile);
        
    long firstrow=1;
    long firstelem=1;
    long nelements=nrows;                  //for time
    double *time = new double[nelements];
    int tcoln,qcoln;               //column number for time and quaternion
    
    fits_get_colnum(fatt,CASEINSEN,timecol,&tcoln,&status);
     printError(status,"Error reading time column number from attitude file  ",attitudefile);   
    
     fits_read_col(fatt,TDOUBLE,tcoln, firstrow,firstelem, nelements, NULL,time,NULL,&status); 
     printError(status,"Error reading time column from attitude file  ",attitudefile);     
    
    long startindex=0, lastindex=nrows-1;
    LOG(INFO)<<time[0]<<" "<<time[nrows-1]<<" "<<tstart<<" "<<tstop;
    if(tstart>time[nrows-1] || tstop<time[0])
    {
        LOG(ERROR)<<"\033[1;31m***Time information in data file and attitude file do not match***\033[0m";
        return (EXIT_FAILURE);
    }
    double time_mid =( tstart+tstop)/2;
   long index_attFile=0;
    for (int i=1;i<nelements;i++)
    {
        if(time_mid>time[i-1] && time_mid<time[i]){
            index_attFile=i-1;
            break;
        }      
        
    }
    //for start time
//    if(tstart<time[0])                         //if starttime for data is earlier than attitude data starttime
//        startindex=0;                 
//    else if(tstart>time[0]){                     //start time of data is after starttime for attitude data
//        while(time[startindex]<tstart){
//            startindex++;
//            if(startindex>nrows-1)  break;
//        }
//        startindex--;                                  //this will be the start index for attitude data to be read
//    }
//    
//    //for stop time
//    if(tstop>time[nrows-1])                     //if data end time is after attitude data stop time
//        lastindex=nrows-1;
//    else if(tstop<time[lastindex]){
//        while(time[lastindex]>tstop){
//            lastindex--;
//            if(lastindex<0)  break;
//        }
//        lastindex++;                                 //this will be last index for attitude data to be read
//    }
    
   // LOG(INFO)<<"Start Index :"<<startindex<<" Last Index : "<<lastindex;
    
   
   //added 
   //index_attFile=1;
     //  long tempo=nelements;
   //to be removed
    nelements = 4;                    //number of elements to be read for quaternions
   double *qs = new double[nelements];
   //  double *qs = new double[nelements*tempo];
    fits_get_colnum(fatt,CASEINSEN,qcol,&qcoln,&status);
    printError(status,"Error reading quaternion column number from attitude file  ",attitudefile);
    
    fits_read_col(fatt,TDOUBLE,qcoln,index_attFile,1, nelements, NULL,(void*)qs,NULL,&status); 
    printError(status,"Error reading quaternions from attitude file  ",attitudefile);
    
    fits_close_file(fatt,&status);
    
    att.q1=qs[0];att.q2=qs[1];att.q3=qs[2];att.q4=qs[3];
    //long i,j;
    
    
   
    LOG(INFO)<<"-------Attitude data------";
    
//    for(i=startindex,j=0; i<=lastindex; i++,j=j+4){               //
//        Attitude att_t(time[i],qs[j+3],qs[j],qs[j+1],qs[j+2]);
//        att.push_back(att_t);
//        att_t.display();
//     }
    
    LOG(INFO)<<"----------------------------------";
        
    delete[] time, qs;
    return (EXIT_SUCCESS);
}


//int readAttitude(char *attitudefile,char *timecol,char *qcol,double tstart,double tstop,Attitude &att,vector<Attitude> &vect){
//    
//    LOG(INFO)<<"Reading attitude data from "<<attitudefile;
//    
//    int status=0;
//    
//    fitsfile *fatt;
//    fits_open_file(&fatt,attitudefile,READONLY,&status);
//    printError(status,"Error opening attitude file",attitudefile);
//    
//    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
//    printError(status,"Error in moving to 2nd HDU ", attitudefile);
//    
//    long nrows;
//    fits_get_num_rows(fatt,&nrows,&status);
//    printError(status,"Error reading number of rows from attitude file  ",attitudefile);
//        
//    long firstrow=1;
//    long firstelem=1;
//    long nelements=nrows;                  //for time
//    double *time = new double[nelements];
//    int tcoln,qcoln;               //column number for time and quaternion
//    
//    fits_get_colnum(fatt,CASEINSEN,timecol,&tcoln,&status);
//     printError(status,"Error reading time column number from attitude file  ",attitudefile);   
//    
//     fits_read_col(fatt,TDOUBLE,tcoln, firstrow,firstelem, nelements, NULL,time,NULL,&status); 
//     printError(status,"Error reading time column from attitude file  ",attitudefile);     
//    
//    long startindex=0, lastindex=nrows-1;
//    
////    if(tstart>time[nrows-1] || tstop<time[0]){
////        LOG(ERROR)<<"\033[1;31m***Time information in data file and attitude file do not match***\033[0m";
////        return (EXIT_FAILURE);
////    }
//    double time_mid =( tstart+tstop)/2;
//   long index_attFile=0;
//    for (int i=1;i<nelements;i++)
//    {
//        if(time_mid>time[i-1] && time_mid<time[i]){
//            index_attFile=i-1;
//            break;
//        }      
//        
//    }
//    //for start time
////    if(tstart<time[0])                         //if starttime for data is earlier than attitude data starttime
////        startindex=0;                 
////    else if(tstart>time[0]){                     //start time of data is after starttime for attitude data
////        while(time[startindex]<tstart){
////            startindex++;
////            if(startindex>nrows-1)  break;
////        }
////        startindex--;                                  //this will be the start index for attitude data to be read
////    }
////    
////    //for stop time
////    if(tstop>time[nrows-1])                     //if data end time is after attitude data stop time
////        lastindex=nrows-1;
////    else if(tstop<time[lastindex]){
////        while(time[lastindex]>tstop){
////            lastindex--;
////            if(lastindex<0)  break;
////        }
////        lastindex++;                                 //this will be last index for attitude data to be read
////    }
//    
//   // LOG(INFO)<<"Start Index :"<<startindex<<" Last Index : "<<lastindex;
//    
//   
//   //added 
//   index_attFile=1;
//       long tempo=nelements;
//   //to be removed
//    nelements = 4;                    //number of elements to be read for quaternions
//   // double *qs = new double[nelements];
//     double *qs = new double[nelements*tempo];
//    fits_get_colnum(fatt,CASEINSEN,qcol,&qcoln,&status);
//    printError(status,"Error reading quaternion column number from attitude file  ",attitudefile);
//    
//    fits_read_col(fatt,TDOUBLE,qcoln,index_attFile,1, nelements*tempo, NULL,(void*)qs,NULL,&status); 
//    printError(status,"Error reading quaternions from attitude file  ",attitudefile);
//    
//    fits_close_file(fatt,&status);
//    
//    att.q1=qs[0];att.q2=qs[1];att.q3=qs[2];att.q4=qs[3];
//    //long i,j;
//    
//    
//    //added block 
//    
//    for(int i=0;i<tempo*nelements;i=i+4){
//        //LOG(INFO)<<"ADDED"<<qs[i]<<" "<<qs[i+1]<<" "<<qs[i+2]<<" "<<qs[i+3];
//                
//        att.q1=qs[i];att.q2=qs[i+1];att.q3=qs[i+2];att.q4=qs[i+3];
//        vect.push_back (att);
//        //temp_vect_To_delete.push_back (att);
//    }
//    
//    
//   // exit(1);
//    //to be removed
//    LOG(INFO)<<"-------Attitude data------";
//    
////    for(i=startindex,j=0; i<=lastindex; i++,j=j+4){               //
////        Attitude att_t(time[i],qs[j+3],qs[j],qs[j+1],qs[j+2]);
////        att.push_back(att_t);
////        att_t.display();
////     }
//    
//    LOG(INFO)<<"----------------------------------";
//        
//    delete[] time, qs;
//    return (EXIT_SUCCESS);
//}
