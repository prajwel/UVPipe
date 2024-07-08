/* 
 * File:   fft.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on July 30, 2014, 2:23 PM
 */
#include "fft.h"
#include "spMatrix.h"
#include "spGeneral.h"
#include<stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include<math.h>

using namespace std;
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void doFft (double data[] , long nn , int isign)
{
    long n , mmax , m , j , istep , i ;
    double wtemp , wr , wpr , wpi , wi , theta ;
    float tempr , tempi ;

    n = nn << 1 ;
    j = 1 ;
    for (i = 1 ; i < n ; i += 2)
    {
        if (j > i)
        {
            SWAP (data[j] , data[i]) ;
            SWAP (data[j + 1] , data[i + 1]) ;
        }
        m = nn ;
        while (m >= 2 && j > m)
        {
            j -= m ;
            m >>= 1 ;
        }
        j += m ;
    }
    mmax = 2 ;
    while (n > mmax)
    {
        istep = mmax << 1 ;
        theta = isign * (6.28318530717959 / mmax) ;
        wtemp = sin (0.5 * theta) ;
        wpr = -2.0 * wtemp*wtemp ;
        wpi = sin (theta) ;
        wr = 1.0 ;
        wi = 0.0 ;
        for (m = 1 ; m < mmax ; m += 2)
        {
            for (i = m ; i <= n ; i += istep)
            {
                j = i + mmax ;
                tempr = wr * data[j] - wi * data[j + 1] ;
                tempi = wr * data[j + 1] + wi * data[j] ;
                data[j] = data[i] - tempr ;
                data[j + 1] = data[i + 1] - tempi ;
                data[i] += tempr ;
                data[i + 1] += tempi ;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr ;
            wi = wi * wpr + wtemp * wpi + wi ;
        }
        mmax = istep ;
    }
}


double* DataFitting (double *t , double *X , int order , int nRecords)
{

    double *Y ;
    Y = new double[nRecords] ;
    for (int i = 0 ; i < nRecords ; i++)
        Y[i] = 0.0 ;

    spMatrix Aa (nRecords , order + 1) ;
    spMatrix Bb (nRecords , 1) ;
    spMatrix coeff_mat (order + 1 , 1) ;

    for (int i = 0 ; i < nRecords ; i++)
    {
        for (int j = 0 ; j < order + 1 ; j++)
            Aa (i , j) = pow (t[i]-t[0] , j) ;
            Bb (i , 0) = X[i] ;
        
    }
    //cout<<Aa;
   // cout<<Bb;
    coeff_mat.ApplyLeastSquare (Aa , Bb) ;

    for (int i = 0 ; i < nRecords ; i++)
    {
        for (int j = 0 ; j < order + 1 ; j++)
        {
            Y[i] = Y[i] + coeff_mat (j , 0) * pow (t[i]-t[0] , j) ;

        }
    }

    return Y ;
}
