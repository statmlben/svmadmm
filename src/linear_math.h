#ifndef LB_LINEAR_MATH_H
#define LB_LINEAR_MATH_H

#include <stdlib.h>	/* for default allocator/deallocator malloc/free */
#ifdef __c
extern "C" {
#endif


#ifdef LB_SINGLE_PRECISION
	typedef float lbScalar;
#else
	#define LB_DOUBLE_PRECISION
	typedef double lbScalar;
#endif

#ifndef LB_SIZE
	#define LB_SIZE
	typedef unsigned int lbSize;
#endif

#ifndef LB_BOOL
	#define LB_BOOL
	typedef int lbBool;
#endif

#ifndef LB_TRUE_FALSE
	#define LB_TRUE_FALSE
	#define LB_TRUE 1
	#define LB_FALSE 0
#endif

#ifndef LB_USE_ROW_MAJOR_MATRIX
	#define LB_USE_COLUMN_MAJOR_MATRIX
#endif

#ifndef LB_LM_DISABLE_ARG_CHECK
	#define LB_LM_ENABLE_ARG_CHECK
#endif

typedef struct _lbVector
{
	lbSize size;
	lbScalar *data;
} lbVector;

typedef struct _lbMatrix
{
	lbSize nRow;
	lbSize nCol;
	lbScalar *data;
} lbMatrix;

typedef void * (*lbAllocator) (size_t size);
typedef void (*lbDeallocator) (void *ptr);
typedef void (*lbErrorReportFunc) (const char *errMsg);

typedef lbScalar (*lbUnaryOperation) (lbScalar value, void *params);
typedef lbScalar (*lbBinaryOperation) (lbScalar left, lbScalar right, void *params);

/* the default error reporter calls
 * fprintf(stderr, "LinearMath: %s\n", errMsg) */
void lbSetLinearMathErrorReporter(lbErrorReportFunc errorReporter);

/**** vector functions ****/

lbVector lbAllocVector(int size, lbAllocator allocator);
void lbFreeVector(lbVector *vec, lbDeallocator deallocator);
lbVector *lbAllocVectors(int size, int count, lbAllocator allocator);
void lbFreeVectors(lbVector *vecs, lbDeallocator deallocator);

lbBool lbVectorEmpty(const lbVector vec);

lbSize lbVectorSize(const lbVector vec);
lbScalar lbGetVectorElem(const lbVector vec, int index);
void lbSetVectorElem(lbVector *vec, int index, lbScalar value);
void lbVectorClear(lbVector *vec, lbScalar value);

/* vector operations */

void lbVectorCopy(const lbVector in, lbVector *out);
void lbVectorAddition(const lbVector left, const lbVector right, lbVector *out);
void lbVectorSubtraction(const lbVector left, const lbVector right, lbVector *out);
lbScalar lbVectorDotProduct(const lbVector left, const lbVector right);
void lbVectorScalarProduct(const lbVector vec, lbScalar scalar, lbVector *out);
lbScalar lbVectorSumElem(const lbVector vec);
void lbVectorProduct(const lbVector left, const lbVector right, lbVector *out);
lbScalar lbVectorNorm(const lbVector vec);
/**
 * Use each element of [in] and [params] as the input to operation [uOper], and put the
 * result into vector [out].
 * the prototype of lbUnaryOperation is: lbScalar (*lbUnaryOperation) (lbScalar value, void *params)
 * note: [in] and [out] must have the same size
 */
void lbVectorForeach(const lbVector in, lbVector *out, lbUnaryOperation uOper, void *params);

/**
 * Use each element of [vec1] and [vec2] and [params] as the input to operation [biOper], and
 * put the result of the operation into [out].
 * the prototype of lbBinaryOperation is:
 * 		lbScalar (*lbBinaryOperation) (lbScalar left, lbScalar right, void *params);
 * note: [vec1], [vec2] and [out] must have the same size
 */
void lbVectorBiForeach(const lbVector vec1, const lbVector vec2, lbVector *out, lbBinaryOperation biOper, void *params);


/**** matrix functions ****/

lbBool lbMatrixRowMajor();
lbMatrix lbAllocMatrix(int nRow, int nCol, lbAllocator allocator);
void lbFreeMatrix(lbMatrix *mat, lbDeallocator deallocator);

lbBool lbMatrixEmpty(const lbMatrix mat);

lbSize lbMatrixRowCount(const lbMatrix mat);
lbSize lbMatrixColumnCount(const lbMatrix mat);

lbScalar lbGetMatrixElem(const lbMatrix mat, int rowIndex, int columnIndex);
void lbSetMatrixElem(lbMatrix *mat, int rowIndex, int columnIndex, lbScalar value);

lbBool lbMatrixSquare(const lbMatrix mat);
lbBool lbMatrixSameSize(const lbMatrix mat1, const lbMatrix mat2);

/* matrix operations */

void lbMatrixCopy(const lbMatrix in, lbMatrix *out);
void lbMatrixTranspose(const lbMatrix in, lbMatrix *out);
void lbMatrixAddition(const lbMatrix matL, const lbMatrix matR, lbMatrix *out);
void lbMatrixSubtraction(const lbMatrix matL, const lbMatrix matR, lbMatrix *out);
void lbVectorMatrixMultiply(const lbVector vec, const lbMatrix mat, lbVector *out);
void lbMatrixVectorMultiply(const lbMatrix mat, const lbVector vec, lbVector *out);
void lbVectorMatrixProduct(const lbVector vec, const lbMatrix mat, lbMatrix *out);
void lbMatrixMultiply(const lbMatrix left, const lbMatrix right, lbMatrix *out);
void lbMatrixColSum(const lbMatrix mat, lbVector *out);
void lbMatrixBiForeach(const lbMatrix mat1, const lbMatrix mat2, lbMatrix *out, lbBinaryOperation biOper, void *params);

void lbMatrixSetDiag(lbMatrix *mat, const lbVector values);
void lbMatrixClearDiag(lbMatrix *mat, lbScalar value);

#ifdef __c
}
#endif

#endif  /* LB_LINEAR_MATH_H */
