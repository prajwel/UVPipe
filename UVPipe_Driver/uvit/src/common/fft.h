/* 
 * File:   fft.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on July 30, 2014, 2:23 PM
 */
#ifndef FFT_H
#define	FFT_H

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function

 double* DataFitting(double *t,double *X,int order,int nRecords);
  void  doFft(double data[],  long nn, int isign);


#endif	/* FFT1_H */

