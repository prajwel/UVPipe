/* 
 * File:   Attitude.h
 * Author: uvit
 *
 * Created on September 30, 2014, 11:10 AM
 */

#ifndef ATTITUDE_H
#define	ATTITUDE_H

#include<vector>

using namespace std;

typedef struct Attitude{
    double time;
    double q1,q2,q3,q4;                   //q1 is scaler
    
    Attitude(){}
    Attitude(double time, double q1,double q2,double q3, double q4){
        this->time=time;
        this->q1=q1;
        this->q2=q2;
        this->q3=q3;
        this->q4=q4;
    }
    void display();
};

 //int readAttitude(char *attitudefile,char *timecol,char *qcol,double tstart,double tstop,Attitude &att,vector<Attitude> &vect);
int readAttitude(char *attitudefile,char *timecol,char *qcol,double tstart,double tstop,Attitude &att);
    


#endif	/* ATTITUDE_H */

