#include "spGeneral.h"

/******************************************************************
*	This program file contains implementation code
*	for all the methods of class "spMatrix" as
*	declared in "spMatrix.h" file

*	@Program Name 	:	spMatrix.cpp

*	@CORE SCIENTIFIC CONTRIBUTOR          	:   P. K. Srivastava
*	@CORE SCIENTIFIC CONTRIBUTOR          	:	T. P. Srinivasan

*	@ANALYST/CONCEPTUAL DESIGNER            :       Amit Gupta

*	@TECHNICAL DESIGNER   :       Amit Gupta

*	@IMPLEMENTER          :       Amit Gupta

*	VALIDATORS
*	@Code Walkthrough     :       
*	@Test case Results    :       T. P. Srinivasan

*	@Created Date			:	14-07-1999


*	@version		:	1.0
******************************************************************/
/**----------------------------- INCLUDE FILES ------------------------------*/

#include <math.h>
#include "spMatrix.h"

/*------------------------------ CONSTANTS ---------------------------------*/



#define NORMAL_MATRIX                           0
#define DEFAULT_ROW_ORDER          		3
#define DEFAULT_COLUMN_ORDER       		3
#define PIVOT_EPSILON		   		1.0e-16
#define MATRIX_INITIALIZATION_VALUE		0

/*----------------------------    TYPEDEFS   ---------------xx--------------*/
typedef  spMatrix :: spMatrixOrder	spMatrixOrder;
typedef  spMatrix :: spMatrixIndex	spMatrixIndex;
typedef  spMatrix :: spScalar		spScalar;
typedef  spMatrix :: spMatrixElement	spMatrixElement;

//MatrixInverseLcmap(double *InvPtr,int n,double*tmp1,double*tmp2)


/**---------------------------------------------------------------------------*/
// CONSTRUCTORS
/**
* Default Constructor.
* This constructor initializes the data members to its default values.
*/
spMatrix :: spMatrix()
{
 no_of_rows    = DEFAULT_ROW_ORDER;
 no_of_columns = DEFAULT_COLUMN_ORDER;
 status        = NORMAL_MATRIX;

 if( (p = new spMatrixElement *[no_of_rows]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;
    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;
	cerr << endl << "Unable to allocate memory for Matrix construction" << endl;
    cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;

    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

 spMatrixIndex i;
 for(i = 0;i < no_of_rows;i++)
   if(( p[i] = new spMatrixElement [no_of_columns]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;
    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;
    cerr << endl << "Unable to allocate memory for Matrix construction" << endl;
    cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;

    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

 for (i=0;i<no_of_rows;i++)
    for (spMatrixIndex j=0;j<no_of_columns;j++)
	p[i][j] = MATRIX_INITIALIZATION_VALUE;
}

/**
* 	Constructor with dimensions.
* 	This constructor initializes the data members to values given as *	parameters
*	and allocates memory required for storing values.
*	@param spMatrixOrder d1 (no of rows)
*	@param spMatrixOrder d2 (no of columns)
*/
spMatrix  :: spMatrix(spMatrixOrder d1,spMatrixOrder d2)
{
 no_of_rows    = d1;
 no_of_columns = d2;
 status        = NORMAL_MATRIX;

 if( (p = new spMatrixElement *[no_of_rows]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;
    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;
    cerr << endl << "Unable to allocate memory for Matrix construction" << endl;
    cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;

    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

 spMatrixIndex i;
 for(i = 0;i < no_of_rows;i++)
   if(( p[i] = new spMatrixElement [no_of_columns]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;
    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;

	cerr << endl << "Unable to allocate memory for Matrix construction" << endl;
    cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;

    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

 for (i=0;i<no_of_rows;i++)
    for (spMatrixIndex j=0;j<no_of_columns;j++)
	p[i][j] = MATRIX_INITIALIZATION_VALUE;
}
/**
*	Copy Constructor.
* 	Assigns data member values of one matrix to another.
*	@param spMatrix  &A
*/
spMatrix :: spMatrix(const spMatrix  &A)		// COPY CONSTRUCTOR
{
 status     = A.status;
 no_of_rows = A.no_of_rows; no_of_columns= A.no_of_columns;

 if( (p = new spMatrixElement *[no_of_rows]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;
    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;

	cerr << endl << "Unable to allocate memory for Matrix construction" << endl;
    cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;

    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

 spMatrixIndex i;
 for(i = 0;i < no_of_rows;i++)
   if(( p[i] = new spMatrixElement [no_of_columns]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;
    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;

    cerr << endl << "Unable to allocate memory for Matrix construction" << endl;
    cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;

    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

 for ( i=0;i<no_of_rows;i++)
    for (spMatrixIndex j=0;j<no_of_columns;j++)
	p[i][j] = A.p[i][j];
}

// DESTRUCTOR
/**
*	Destructor.
*	Deallocate the memory.
*/
spMatrix  :: ~spMatrix()
{
 for(spMatrixIndex i = 0;i < no_of_rows;i++)
    delete [] p[i];

 delete [] p;
}

/**---------------------------------------------------------------------------*/
// OVERLOADED OPERATORS FOR MATRIX OPERATIONS
/**
*	= operator.
* 	Assigns data member values of one matrix to another.
*	If Left hand sid matrix is not of the same order that of
*	the right hand side matrix then this method exits with
*	error code INVALID_MATRIX_ORDER_ASSIGN.
*
*	@param spMatrix  &A
*	@return an instance of spMatrix.
*/

spMatrix  & spMatrix  :: operator =(const spMatrix  &A)
{
 status     = A.status;
 if((no_of_rows == A.no_of_rows) && (no_of_columns == A.no_of_columns))
   for ( spMatrixIndex i=0;i<no_of_rows;i++)
     for (spMatrixIndex j=0;j<no_of_columns;j++)
	p[i][j] = A.p[i][j];
 else
  {
   clog << endl
        << "Error: Assigning a Matrix to another Matrix of different order."
	<< endl;

   cerr << "Error: Assigning a Matrix to another Matrix of different order."
		<< endl;

    status = INVALID_MATRIX_ORDER_ASSIGN;
   exit(INVALID_MATRIX_ORDER_ASSIGN);
  }

return *this;
}

spMatrix  operator +(const spMatrix  &A,const spMatrix &B)
/**
*	This method defines addition operator('+') for adding matrices of appropriate orders.
*	If the matrices are not of appropriate orders then this method
*	exits with errorcode INVALID_MATRIX_ORDER_ADD.
*
*	@param 				spMatrix  A
*	@param 				spMatrix  B </br>
*
*	<b>Pre-condition</b>     :  Matrix 'A' & Matrix 'B' must be of same order</br>
*
*	<b>Post-condition</b>    :  Resultant matrix must be of same order as Matrix 'A'
*	@return instance of spMatrix containig values of A+B
*/
{
spMatrix  C(A.no_of_rows,A.no_of_columns);

if((A.no_of_rows == B.no_of_rows)  &&  (A.no_of_columns == B.no_of_columns))
{
for (spMatrixIndex i=0;i<C.no_of_rows;i++)
  for (spMatrixIndex j=0;j<C.no_of_columns;j++)
    C.p[i][j]=A.p[i][j]+B.p[i][j];

}

else
{
clog << endl << "Error! An attempt to add order-incompatible matrices." << endl;

cerr << endl << "Error! An attempt to add order-incompatible matrices." << endl;
C.status = INVALID_MATRIX_ORDER_ADD;
exit(INVALID_MATRIX_ORDER_ADD);
}

return(C);
}

spMatrix   operator -(const spMatrix &A,const spMatrix &B)

/**
*	This method defines subtraction operator('-')
*	for subtracting matrices of appropriate orders.
*	If the matrices are not of appropriate orders then this method
*	exits with errorcode INVALID_MATRIX_ORDER_SUBTRACT.
*
*	@param 				spMatrix  A
*	@param 				spMatrix  B </br>
*
*	<b>Pre-condition</b>     :  Matrix 'A' & Matrix 'B' must be of same order</br>
*
*	<b>Post-condition</b>    :  Resultant matrix must be of same order as Matrix 'A'
*	@return instance of spMatrix containig values of A-B
*/
{
spMatrix  C(A.no_of_rows,A.no_of_columns);

if((A.no_of_rows == B.no_of_rows)  &&  (A.no_of_columns == B.no_of_columns))
{
for (spMatrixIndex i=0;i<C.no_of_rows;i++)
  for (spMatrixIndex j=0;j<C.no_of_columns;j++)
    C.p[i][j]=A.p[i][j] - B.p[i][j];

}

else
{
clog << endl << "Error! An attempt to subtract order-incompatible matrices."
     << endl;

cerr << "Error! An attempt to subtract order-incompatible matrices."
     << endl;

C.status = INVALID_MATRIX_ORDER_SUBTRACT;
exit(INVALID_MATRIX_ORDER_SUBTRACT);
}

return(C);
}

spMatrix   operator *(const spMatrix &A,const spMatrix &B)
/**
*	This method defines multiplication operator('*')
*	for multiplying matrices of appropriate orders.
*	If the matrices are not of appropriate orders then this method
*	exits with errorcode INVALID_MATRIX_ORDER_MULTIPLY.
*
*	@param 				spMatrix  A
*	@param 				spMatrix  B </br>
*
*	<b>Pre-condition</b>     :  No. of columns of Matrix 'A' must be equal to no. of
*				        rows of matrix 'B'.</br>
*	<b>Post-condition</b>    :  Let A has order m*n & B has order n*l. Resultant matrix
*				        must be of order m*l.
*	@return instance of spMatrix containig values of A*B
*/
{
spMatrix   C(A.no_of_rows,B.no_of_columns);

if(A.no_of_columns == B.no_of_rows)
{
 for (spMatrixIndex i=0;i<A.no_of_rows;i++)
   for (spMatrixIndex j=0;j<B.no_of_columns;j++)
     for (spMatrixIndex k=0;k<A.no_of_columns;k++)
     C.p[i][j] += A.p[i][k]*B.p[k][j];

}

else
{
clog <<endl<< "Error: An attempt to multiply matrices of invalid orders" <<endl;

cerr << "Error: An attempt to multiply matrices of invalid orders" <<endl;
C.status = INVALID_MATRIX_ORDER_MULTIPLY;
exit(INVALID_MATRIX_ORDER_MULTIPLY);
}

return(C);
}

spMatrix   operator *(const  spMatrix & A,const spScalar & val)
/**
*	This method defines multiplication operator('*')
*				for post-multiplication of a scalar with a matrix
*
*	@param 		spMatrix A
*	@param 		a scalar value val </br>
*
*	<b>Pre-condition</b>     :  None</br>
*
*	<b>Post-condition</b>    :  Resultant matrix must be of same order as Matrix 'A'
*
*	@return instance of spMatrix containig values of A*val
*/
{
spMatrix  R(A.no_of_rows,A.no_of_columns);

 for (spMatrixIndex i=0;i<A.no_of_rows;i++)
   for (spMatrixIndex j=0;j<A.no_of_columns;j++)
	R.p[i][j] = A.p[i][j] * val;

return R;
}

spMatrix   operator *(const spScalar &val ,const spMatrix   &A)
/**
*	This method defines multiplication operator('*')
*				for pre-multiplication of a scalar with a matrix
*
*	@param		a scalar value val
*	@param 		spMatrix A</br>
*
*	<b>Pre-condition</b>     :  None</br>
*
*	<b>Post-condition</b>    :  Resultant matrix must be of same order as Matrix 'A'
*
*	@return instance of spMatrix containig values of val*A
*/
{
spMatrix   R(A.no_of_rows,A.no_of_columns);

 for (spMatrixIndex i=0;i<A.no_of_rows;i++)
   for (spMatrixIndex j=0;j<A.no_of_columns;j++)
	R.p[i][j] = A.p[i][j] * val;

return R;
}

/**
* This method defines operator '>>' to support
* Matrix object
*/
istream & operator >>(istream & in_stream,spMatrix  & A)
{

for (spMatrixIndex i=0;i < A.no_of_rows;i++)
     for (spMatrixIndex j=0;j < A.no_of_columns;j++)
       in_stream >> A.p[i][j];

return in_stream;
}

/**
* This method defines operator '<<' to support
*			Matrix object
*/
ostream & operator <<(ostream & out_stream,const spMatrix  & A)
{

int width = out_stream.width();

for (spMatrixIndex i=0;i < A.no_of_rows;i++)
     {
     for (spMatrixIndex j=0;j < A.no_of_columns;j++)
       out_stream << setw(width) << A.p[i][j] << " ";
       out_stream << endl;
     }


 return out_stream;
}

/*--------------------------------------------------------------------------*/
// OTHER MATRIX OPERATIONS
/**
* This method interchanges specified rows of the
*                      matrix which calls this method
*
* @param	      spMatrixIndex row1
* @param	      spMatrixIndex row2
*	Two rows which are to be interchanged
*/
void spMatrix :: InterchangeRows(spMatrixIndex row1,spMatrixIndex row2)
  {
   spMatrixElement *temporary_pointer;

   if((row1 >= 0) && (row1 < no_of_rows) && (row2 >= 0) && (row2 < no_of_rows))
     {
	temporary_pointer = p[row1];
	p[row1]           = p[row2];
	p[row2]           = temporary_pointer;
     }
   else
     {
	clog << "Error: An out-of-range row has been specified in a call";
        clog << endl <<  "to function InterchangeRows()" << endl;
	clog << "row1 = " <<  row1 << ", row2 = " << row2 << endl;

	cerr << "Error: An out-of-range row has been specified in a call"
 		<<  "to function InterchangeRows()"
		<< "row1 = " <<  row1 << ", row2 = " << row2 << endl;
        status = OUT_OF_RANGE_MATRIX_ACCESS;
        exit(OUT_OF_RANGE_MATRIX_ACCESS);
     }
  }

	 
/**
* Computes the transpose of the matrix
*			which calls this method </br>
*
* <b>Post-condition</b>     :	Resultant matrix must be of reverse order.
*								This means if calling matrix is *	of order m*n then
*								resultant matrix must be of *	order n*m.
*/
spMatrix   spMatrix :: Transpose() const
{
spMatrix  X(no_of_columns,no_of_rows);

for (spMatrixIndex i=0;i<X.no_of_rows;i++)
   for (spMatrixIndex j=0;j<X.no_of_columns;j++)
         X.p[i][j] = p[j][i];


return(X);
}

/**
* Computes the trace  of the matrix
*			which calls this method.</br>
*
* <b>Pre-condition</b>      :	Calling matrix must be a square matrix </br>
*
* <b>Post-condition</b>      :	None
*/
double  spMatrix :: Trace() const
{
 spMatrix   A = *this;

 if(A.no_of_rows != A.no_of_columns)
   {
    clog << endl << "Error: A non-square matrix has been attempted for trace" << endl;
    cerr << "Error: A non-square matrix has been attempted for trace" << endl;

    A.status = INVALID_MATRIX_ORDER_TRACE;
    exit(INVALID_MATRIX_ORDER_TRACE);
   }

 double trace=0.0;
 for ( int i=0; i<A.no_of_rows; i++ )
   {
    trace += A(i,i);
   }

 return trace;
}

/**
* Computes the inverse  of the matrix
*			which calls this method. "Gauss-Jordan" method has been
*			used in the following code. </br>
*
* <b>Pre-condition</b>      :	Calling matrix must be a square matrix </br>
*
* <b>Post-condition</b>      :	Resultant matrix must be of same  order as the 
*								calling matrix.
*/
spMatrix  spMatrix :: Inverse()
{

 spMatrix   A = *this;
 spMatrixOrder N = A.no_of_rows;
 spMatrix INVERSE(N,N);

 if(A.no_of_rows != A.no_of_columns)
   {
    clog << endl << "Error: A non-square matrix has been attempted for inversion" << endl;
    cerr << "Error: A non-square matrix has been attempted for inversion" << endl;
    INVERSE.status = SINGULAR_MATRIX;
    exit(SINGULAR_MATRIX);
   }



 spMatrixIndex n;
 for(n = 0;n < N;n++)   // Preparing INVERSE as an identity matrix
    INVERSE(n,n) = 1.0;

 for(n = 0;n < N;n++)
    {
     spMatrixIndex r = A.MaximumElementRow(n,n);
     A.InterchangeRows(n,r);
     INVERSE.InterchangeRows(n,r);

     if(fabs(A(n,n))  <= PIVOT_EPSILON)
	{
	 clog << endl << "Error: Matrix not invertible" << endl;
	 cerr << "Error: Matrix not invertible" << endl;
         INVERSE.status = SINGULAR_MATRIX;
	 exit(SINGULAR_MATRIX);
	}


     spMatrixElement     pivot = A(n,n);

     spMatrixIndex j;
     for(j = 0;j < N;j++)
	{

	 A(n,j)    	 = A(n,j)  / pivot;
	 INVERSE(n,j)    = INVERSE(n,j)  / pivot;
	}
     for(spMatrixIndex i = 0;i < N;i++)
	if(i != n)
          {
           spMatrixElement temp = A(i,n);
	   for(j = 0;j < N;j++)
	      {
		A(i,j) 	     = A(i,j) - A(n,j) * temp;
		INVERSE(i,j) = INVERSE(i,j) - INVERSE(n,j) * temp;
	      }
          }
     }
 return(INVERSE);

}


/**
* Solves the simultaneous equations	represented in the following form: </br>
*			A*X = B
*    				where A is coefficient matrix; X is a vector
*			of unknown variables & B is resultant vector.
*
* Arguments	      :	Matrix 'A' & 'B'</br>
*
* <b>Pre-condition</b>      : No. of columns in matrix 'A' must be equal to no. 
*					of rows in matrix 'B'. Matrix 'A' must be a square 
*					matrix.</br>
*
* <b>Post-condition</b>    : No. of rows in resultant vector 'X' must be equal to
*			no. of rows in Matrix 'A'
*/
void spMatrix :: SolveEquations(spMatrix A, spMatrix  B)
{
 spMatrix & X = *this;

 if(A.NO_OF_ROWS() !=  A.NO_OF_COLUMNS())
   {
    clog << endl << "Error: Coefficients matrix 'A' is not a square matrix ";
    clog << "inside solve_equations" << endl;
    
	cerr << "Error: Coefficients matrix 'A' is not a square matrix "
    	<< "inside solve_equations" << endl;
    X.status = SINGULAR_MATRIX;
    exit(SINGULAR_MATRIX);
   }


 if((A.NO_OF_COLUMNS() == X.NO_OF_ROWS()) && (A.NO_OF_ROWS() == B.NO_OF_ROWS())
    && (X.NO_OF_COLUMNS() == 1) && (B.NO_OF_COLUMNS() == 1))
	{

 	spMatrixOrder N = A.NO_OF_ROWS();

         spMatrixIndex n;
	 for(n = 0;n < B.NO_OF_ROWS();n++)
    	     X(n) = B(n);

	 for(n = 0;n < N;n++)
    	    {
     		spMatrixIndex r = A.MaximumElementRow(n,n);
     		A.InterchangeRows(n,r);
     		X.InterchangeRows(n,r);

     		if(fabs(A(n,n)) <= PIVOT_EPSILON)
        	  {
               
        	  	clog << endl
                	 << "Error: Ill-conditioned equations, solution ";
		  		clog << "not possible in solve_equations" << endl;
        	  	
				cerr << "Error: Ill-conditioned equations, solution "
		  			 << "not possible in solve_equations" << endl;

                X.status = ILL_CONDITIONED_EQUATIONS_MATRIX;
	      	  	exit(ILL_CONDITIONED_EQUATIONS_MATRIX);
        	  }


     		spMatrixElement     pivot = A(n,n);
                spMatrixIndex j;
     		for(j = 0;j < N;j++)
	           A(n,j)    = A(n,j)  / pivot;

	        X(n)    = X(n)  / pivot;

	     	for(spMatrixIndex i = 0;i < N;i++)
        	   if(i != n)
	             {
           	     spMatrixElement temp = A(i,n);
           	     for(j = 0;j < N;j++)
               		 A(i,j) = A(i,j) - A(n,j) * temp;
                     X(i) = X(i) - X(n) * temp;
	             }
     	    }
	 }
 else
	 {
	 clog << "Error: Please check the orders of Matrices A, X, & B that ";
	 clog << "were passed to spMatrix :: SolveEquations()";
	 
	 cerr << "Error: Please check the orders of Matrices A, X, & B that "
	  	  << "were passed to spMatrix :: SolveEquations()"<<endl;
         X.status = INVALID_MATRIX_ORDER_SOLVE_EQ;
	 exit(INVALID_MATRIX_ORDER_SOLVE_EQ);
	 }
}


/**
* Solves the simultaneous equations
*			using least-squares method. Equations are represented
*			in the following form:</br>
*
*			<b>A*X = B</b>
*    				where A is coefficient matrix; X is a vector
*			of unknown variables & B is resultant vector.
*
* @param      	spMatrix A
* @param		spMatrix B</br>
*
* <b>Pre-condition</b>      : No. of rows in matrix 'A' must be equal to no. of
* 			rows in matrix 'B'.</br>
*
* <b>Post-condition</b>     : No. of rows in resultant vector 'X' must be equal 
*								to no. of rows in Matrix 'A'.
*/
void spMatrix :: ApplyLeastSquare(const spMatrix A, const spMatrix  B)
{
 spMatrix & X = *this;

	spMatrix temp = X;
 if((A.NO_OF_COLUMNS() == X.NO_OF_ROWS()) && (A.NO_OF_ROWS() == B.NO_OF_ROWS())
    && (X.NO_OF_COLUMNS() == 1) && (B.NO_OF_COLUMNS() == 1))
	{
 	spMatrixOrder n = A.NO_OF_COLUMNS();
 	spMatrix A1(n,n);
 	spMatrix B1(n,1);


 	A1 = A.Transpose() * A;
 	B1 = A.Transpose() * B;
 	X.SolveEquations(A1,B1);
	if(X.status == ILL_CONDITIONED_EQUATIONS_MATRIX)
	{
		X = temp;
	}
	}
 else
	{
	clog << "Error: Please check the orders of Matrices A, X, & B that ";
	clog << "were passed to spMatrix :: ApplyLeastSquare()";
	
	cerr << "Error: Please check the orders of Matrices A, X, & B that ";
	cerr << "were passed to spMatrix :: ApplyLeastSquare()"<<endl;
        X.status = INVALID_MATRIX_ORDER_LEASE_SQ;
	exit(INVALID_MATRIX_ORDER_LEASE_SQ);
	}
}


//MatrixInverseLcmap(InvPtr,n,tmp1,tmp2)


//	int	n; /* input */
//	double	*InvPtr; /* Input & Output */
//	double	tmp1[3],tmp2[3]; /* Input */
/*{

	int	i,j,k,l,m,p,q,r,s,t,u,v,w,x,y,z,a;
	double	d,biga,Hold,rabs1,rabs2;


	d = 1.0;
	p = -n;

	for(k=1; k<n+1; k++){
	  p = p + n;
	  tmp1[k-1] = (double)k;
	  tmp2[k-1] = (double)k;
	  q = p + k;
	  
	  biga = *(InvPtr+q-1);
	  for(j=k; j<n+1; j++){
	    r = n * (j-1);
	    for(i=k; i<n+1; i++){
		
	 	s = r + i;
		if(fabs(biga) < 0)
		     rabs1 = -(biga);
		else
		     rabs1 = biga;
		if(fabs(*(InvPtr+s-1)) < 0)
		     rabs2 = -(*(InvPtr+s-1));
		else
		     rabs2 = (*(InvPtr+s-1));

		if(rabs1 < rabs2){

		   biga = *(InvPtr + s -1);
		   tmp1[k-1] = (double)i;
		   tmp2[k-1] = (double)j;
		}
	    }
	  }
	  j = (int)tmp1[k-1];
	  if((j-k) > 0){
	    t = k - n;
	    for(i=1; i<n+1; i++){
	      t +=n;
	      Hold = - (*(InvPtr + t -1));
	      u = t - k + j;
	      *(InvPtr+t-1) = *(InvPtr+u-1);
	      *(InvPtr+u-1) = Hold;
	    }
	  }
	  i = (int)tmp2[k-1];
	  if((i-k) > 0){
	    v = n * (i - 1);
	    for(j=1; j<n+1; j++){
		w = p + j;
		u = v + j;
		Hold = -(*(InvPtr+w-1));
		*(InvPtr+w-1) = *(InvPtr+u-1);
		*(InvPtr+u-1) = Hold;
	    }
	  }
	  if(biga != 0.0){
	    for(i=1; i<n+1; i++){
		if((i-k) != 0){
		   x = p + i;

		   *(InvPtr+x-1) = (*(InvPtr+x-1))/(-(biga));
		}
	    }
	    for(i=1; i<n+1; i++){
	     	x = p + i;
		Hold = *(InvPtr+x-1);
		s = i - n;
		for(j=1; j<n+1; j++){
		   s = s + n;
		   if((i-k) != 0){
			if((j-k) != 0){
			   y = s - i + k;
			   *(InvPtr+s-1) = Hold * (*(InvPtr+y-1)) + *(InvPtr+s-1);
		        }
		   }
	        }
	     }
	     y = k - n;
	     for(j=1; j<n+1; j++){
		y = y + n;
		if((j-k) != 0){
		   *(InvPtr+y-1) = (*(InvPtr+y-1)) / biga ;
	        }
	     }
	     d = d * biga;
	     *(InvPtr+q-1) = 1.0 / biga;
	  }
	}

	k = n;
	k -=1;
	
	 while( k > 0){
	     i = (int)tmp1[k-1];
	     if((i-k) > 0){
		z = n * (k-1);
		a = n * (i-1);
  		for(j=1; j<n+1; j++){
		   w = z + j;
		   u = a + j;

		   Hold = *(InvPtr+w-1);
		   *(InvPtr+w-1) = -(*(InvPtr+u-1));
		   *(InvPtr+u-1) = Hold;

		}
	     }
	     j = (int)tmp2[k-1];
	     if((j-k) > 0){
		t = k - n;
		for(i=1; i<n+1; i++){
		   t = t + n;
		   Hold = *(InvPtr+t-1);
		   u = t - k + j;
		   *(InvPtr+t-1) = - (*(InvPtr+u-1));
		   *(InvPtr+u-1) = Hold;
		}
	     }
	     k = k - 1; 
	  }

} 
 


*/
/**
* Extracts the specified row, stores 
*			it in a Matrix object, & returns that object
*
* @param	      Matrix index for the row to be extracted </br>
*
* <b>Pre-condition</b>      : Row index should be a valid index, i.e.
* 			0 <= row_number < no_of_rows. </br>
*
* <b>Post-condition</b>     : No. of columns in resultant vector 'X' must be 
*					equal to no. of columns in Matrix A & no of rows 
*					in X should be 1.
*/
spMatrix spMatrix :: ExtractRow(spMatrixIndex	row_number)
{
 spMatrix X(1,no_of_columns);

 if((0 <= row_number) && (row_number < no_of_rows))
   {
	for(spMatrixOrder j = 0; j < no_of_columns; j++)
	   X(0,j) = p[row_number][j];
   }
 else
   {
	clog << "Error: An out of range access has been attempted"
	     << " in spMatrix :: ExtractRow(" << row_number << ")" << endl;
	
	cerr << "Error: An out of range access has been attempted"
	     << " in spMatrix :: ExtractRow(" << row_number << ")" << endl;
        X.status = OUT_OF_RANGE_MATRIX_ACCESS;
	exit(OUT_OF_RANGE_MATRIX_ACCESS);
   }

 return X;
}

/**
* Extracts the specified column,
*			stores it in a Matrix object, & returns that object
*
* @param      Matrix index for the column to be extracted </br>
*
* <b>Pre-condition</b>      : Column index should be a valid index, i.e.
* 			0 <= column_number < no_of_columns . </br>
*
* <b>Post-condition</b>     : No. of rows in resultant vector 'X' must be equal 
* 					to no. of rows in Matrix A & no of columns in X 
*					should be 1
*/
spMatrix spMatrix :: ExtractColumn(spMatrixIndex	column_number)
{
 spMatrix X(no_of_rows,1);

 if(( column_number >= 0 ) && (column_number < no_of_columns))
   {
	for(spMatrixOrder j = 0; j < no_of_rows; j++)
	   X(j,0) = p[j][column_number];
   }
 else
   {
	clog << "Error: An out of range access has been attempted"
	     << " in spMatrix :: ExtractColumn(" << column_number << ")" << endl;
	
	cerr << "Error: An out of range access has been attempted"
	     << " in spMatrix :: ExtractColumn(" << column_number << ")" << endl;
        X.status = OUT_OF_RANGE_MATRIX_ACCESS;
	exit(OUT_OF_RANGE_MATRIX_ACCESS);
   } 
 
 return X;
}

/*-------------------------------------------------------------------------*/

// INQUIRY OPERATION

/**
*	Find the row which has the maximum value in the specified column
*	@param spMatrixIndex column
*	@param spMatrixIndex offset_row (default = 0)
*	
*	@return Index of the row which has maximum value.
*/
spMatrixIndex spMatrix :: MaximumElementRow(spMatrixIndex column,
					    spMatrixIndex offset_row)
 {
  spMatrixIndex r = offset_row;
  if(r < (no_of_rows - 1))
    { 
     spMatrixElement maximum;

     maximum = fabs(p[offset_row][column]);
	
     for(spMatrixIndex i = offset_row + 1;i < no_of_rows;i++)
        if(maximum < fabs(p[i][column]))
	  {
	   maximum =  fabs(p[i][column]);
           r = i;
	  }
     }

  return(r);
 }

/*-------------- METHODS : spMatrix33 ----------------------------------------*/

// CONSTRUCTORS

/**
* 	Default Constructor.
*	Sets the dimension as 3X3.
*/
spMatrix33 :: spMatrix33():spMatrix(3,3){}

/**
*	Copy Constructor.
* 	@param spMatrix A
*/
spMatrix33 :: spMatrix33(const spMatrix & A): spMatrix(A)
{}

// OPERATOR

/**
*	This method defines = operator the class spMatrix33.
*	@param spMatrix A
*	@return Instance of spMatrix33 with same data values as A.
*/
spMatrix33 & spMatrix33 :: operator= (const spMatrix & A)
{
 spMatrix :: operator =(A); 

 return (*this);
}


/**
*	Resizes the matrix dimensions with maintaining data values.
*	@param spMatrixOrder d1
*	@param spMatrixOrder d2
*	New dimensions of matrix.
*/
void spMatrix :: Resize(spMatrixOrder  d1, spMatrixOrder d2)
{
// keeping old reference and values.
 spMatrixElement **p_old = p;
 int no_of_rows_old = no_of_rows;
 int no_of_columns_old = no_of_columns;
 status = NORMAL_MATRIX;

// Creating resized matrix.
 no_of_rows = d1;
 no_of_columns = d2;

 if( (p = new spMatrixElement *[no_of_rows]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;
    
	cerr << "Unable to allocate memory for Matrix construction" << endl;
	cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;
    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }
 spMatrixIndex i;
 for(i = 0;i < no_of_rows;i++)
   if(( p[i] = new spMatrixElement [no_of_columns]) == NULL)
    {
    clog << endl << "Unable to allocate memory for Matrix construction" << endl;    clog << "no_of_rows = " << no_of_rows << endl;
    clog << "no_of_columns = " << no_of_columns << endl;
	cerr << "Unable to allocate memory for Matrix construction" << endl;
	cerr << "no_of_rows = " << no_of_rows << endl;
    cerr << "no_of_columns = " << no_of_columns << endl;
    status = MATRIX_CONSTRUCTION_ERROR;
    exit(MATRIX_CONSTRUCTION_ERROR);
    }

// Copying the content from previous matrix;
 for ( i=0;i<no_of_rows && i<no_of_rows_old;i++)
    for (spMatrixIndex j=0;j<no_of_columns && j<no_of_columns_old;j++)
        p[i][j] = p_old[i][j];

// Deallocating previous memory;
 for(i = 0;i < no_of_rows_old;i++)
    delete [] p_old[i];

 delete [] p_old;
}
/**
*
*	Method to access an element of the matrix at the row,column position.
*	If the position of row or column is out of range or less then 0 then 
*	this method exits with error code : OUT_OF_RANGE_MATRIX_ACCESS.
*	
*	@param  spMatrixIndex row
*	@param  spMatrixIndex column
*	@return Returns the value at (row,column) position in Matrix
*	
*/
spMatrix :: spMatrixElement & spMatrix:: operator()(spMatrixIndex row, 
					spMatrixIndex column) const 
{

#ifndef OPTIMIZATION

 if(!( (0 <= row) && (row < no_of_rows) && (0 <= column) && 
       (column < no_of_columns)))
   {
    clog << "Error: Accessed matrix element is out of range" << endl;
    clog << "Maximum allowed row index    = " << no_of_rows - 1 << endl;
    clog << "Maximum allowed column index = " << no_of_columns - 1 << endl;
    clog << "Accessed row index           = " << row << endl;
    clog << "Accessed column index        = " << column << endl;
    
	cerr << "Error: Accessed matrix element is out of range" << endl;
    cerr << "Maximum allowed row index    = " << no_of_rows - 1 << endl;
    cerr << "Maximum allowed column index = " << no_of_columns - 1 << endl;
    cerr << "Accessed row index           = " << row << endl;
    cerr << "Accessed column index        = " << column << endl;
    exit(OUT_OF_RANGE_MATRIX_ACCESS);
   }

#endif // OPTIMIZATION
 return(p[row][column]);

}


