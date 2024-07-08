#include "spGeneral.h"
// File  Name 	 	:	spMatrix.h
//
// CORE SCIENTIFIC 
// CONTRIBUTOR          	: P. K. Srivastava      
//				: T. P. Srinivasan
//
// ANALYST/CONCEPTUAL     
// DESIGNER             :       Amit Gupta
//
// TECHNICAL DESIGNER   :       Amit Gupta      

// IMPLEMENTER          :       Amit Gupta       

// VALIDATORS             
// Code Walkthrough     :        
// Test case Results    :       T. P. Srinivasan
//
// Date			:	14-07-1999
//
// Description		:	This program file contains declarations for
//				class "spMatrix" & "spMatrix33" 
//
// Version		:	1.1

#ifndef 	_SP_MATRIX_H_
#define	_SP_MATRIX_H_ 

/*--------------------------- INCLUDE FILES -------------------------------*/
#include<iostream>
#include<iomanip> 
#include<stdlib.h>
#include<fstream>

using namespace std;

/*-------------------------- CONSTANTS    ---------------------------------*/ 
#define OUT_OF_RANGE_MATRIX_ACCESS   		(ERROR_NO_BIAS + 10)
#define MATRIX_CONSTRUCTION_ERROR    		(ERROR_NO_BIAS + 11)
#define INVALID_MATRIX_ORDER_ASSIGN		(ERROR_NO_BIAS + 12)
#define INVALID_MATRIX_ORDER_ADD		(ERROR_NO_BIAS + 13)
#define INVALID_MATRIX_ORDER_SUBTRACT		(ERROR_NO_BIAS + 14)
#define INVALID_MATRIX_ORDER_MULTIPLY		(ERROR_NO_BIAS + 15)
#define INVALID_MATRIX_ORDER_TRACE		(ERROR_NO_BIAS + 16)
#define INVALID_MATRIX_ORDER_SOLVE_EQ		(ERROR_NO_BIAS + 17)
#define INVALID_MATRIX_ORDER_LEASE_SQ		(ERROR_NO_BIAS + 18)

#define SINGULAR_MATRIX  (ERROR_NO_BIAS + 19)
#define ILL_CONDITIONED_EQUATIONS_MATRIX     (ERROR_NO_BIAS + 20)

//#pragma once

/*--------------------- CLASS "spMatrix" DECLARATIONS -----------------------*/
/******************************************************************************
*	This program file contains implementation code
*	for all the methods of class "spMatrix" as
*	declared in "spMatrix.h" file

*	Program Name 	:	spMatrix.h </br>

*	CORE SCIENTIFIC CONTRIBUTOR          	:   P. K. Srivastava</br>
*	CORE SCIENTIFIC CONTRIBUTOR          	:	T. P. Srinivasan</br>

*	ANALYST/CONCEPTUAL DESIGNER            :       Amit Gupta</br>

*	TECHNICAL DESIGNER   :       Amit Gupta</br>

*	IMPLEMENTER          :       Amit Gupta</br>

*	VALIDATORS
*	Code Walkthrough     :       </br>
*	Test case Results    :       T. P. Srinivasan</br>

*	Created Date			:	14-07-1999</br>


*	@version		:	1.0</br>

*/
class  spMatrix 
{

/*---------------------    TYPEDEFS   --------xx--------------*/
public:
/******************************************************************************
*	spMatrixOrder is int.
*/
typedef  int  		spMatrixOrder;
/******************************************************************************
*	spMatrixIndex is int.
*/
typedef  int 		spMatrixIndex;
/******************************************************************************
*	spScalar is double.
*/
typedef  double		spScalar;
/******************************************************************************
*	spMatrixElement is double.
*/
typedef  double		spMatrixElement;
/******************************************************************************


/*---------------------  DATA MEMBERS ------------------------*/
protected:
/******************************************************************************
*	 Stores pointer to a matrix.
*/
spMatrixElement  **p;          		     // matrix pointer
/******************************************************************************
*	Stores number of rows.
*/
spMatrixOrder    no_of_rows;
/******************************************************************************
*	Stores number of columns.
*/
spMatrixOrder no_of_columns;  // dimensions of matrix 
/******************************************************************************
*	Stores status.
*/
int              status;

/*---------------------  METHODS  ---------------------------*/


public:

// CONSTRUCTORS
/**
* Default Constructor.
* This constructor initializes the data members to its default values.
*/
spMatrix();
spMatrix(spMatrixOrder  rows, spMatrixOrder columns = 1);
spMatrix(const spMatrix &); 			// COPY CONSTRUCTOR

// DESTRUCTOR
virtual ~spMatrix();

// OPERATORS

friend spMatrix    operator +(const spMatrix & A,const spMatrix & B);
friend spMatrix    operator -(const spMatrix & A,const spMatrix & B);
friend spMatrix    operator *(const spMatrix & A,const spMatrix  & B);
friend spMatrix    operator *(const spMatrix & A,const spScalar & val);
friend spMatrix    operator *(const spScalar & vaL,const spMatrix  & A );
friend istream &   operator >>(istream & in_stream,spMatrix  & A);
friend ostream &   operator <<(ostream & out_stream,const spMatrix  & A);

spMatrix & 	   operator =(const spMatrix &);
spMatrixElement &  operator() (spMatrixIndex row,spMatrixIndex column=0) const ;

// OPERATIONS
void     	 InterchangeRows(spMatrixIndex row1,spMatrixIndex row2);
spMatrix         Transpose() const;
spMatrix         Inverse();                // Gauss - Jordan Method
void		 SolveEquations(spMatrix   A,spMatrix   B); 
void 		 SolveEqsUsingGaussElim(spMatrix);
void		 ApplyLeastSquare(const spMatrix   A,const spMatrix   B);
spMatrix	 ExtractRow(spMatrixIndex row_number);
spMatrix	 ExtractColumn(spMatrixIndex column_number);
double		 Trace() const;
void             Resize(spMatrixOrder  d1, spMatrixOrder d2);
// DATA MEMBER ACCESSORS

/******************************************************************************
*	Returns the Value of no_of_rows.
*/
spMatrixOrder	 NO_OF_ROWS(void) const
{
return no_of_rows;
}

/******************************************************************************
*	Returns the Value of no_of_columns.
*/
spMatrixOrder 	 NO_OF_COLUMNS(void) const
{
return no_of_columns;
}

/******************************************************************************
*	Returns the Value of status.
*/
int              STATUS() const
{
return status;
}

// INQUIRY
spMatrixIndex   MaximumElementRow(spMatrixIndex  column, 
				   spMatrixIndex  offset_row = 0);

};

// INLINE METHOD

/*------------------- CLASS "spMatrix33" DECLARATIONS -----------------------*/

/******************************************************************************
*	This class creates a three by three matrix.
*/
class spMatrix33 : public spMatrix
{
 public:
// CONSTRUCTOR

spMatrix33();

spMatrix33(const spMatrix & A);    //copy constructor

// DESTRUCTOR

/******************************************************************************
*	Destructor 
*/
virtual  ~spMatrix33(){}

// OPERATOR

spMatrix33 & operator= (const spMatrix & A);

};

#endif // _SP_MATRIX_H_
