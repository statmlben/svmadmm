#include "linear_math.h"
#include <stdio.h>	/* for fprintf */
#include <string.h>	/* for memset and memcpy */
#include <math.h>


#ifdef NDEBUG
#define HANDLE_ERROR
#else
#define HANDLE_ERROR /* if it's debugging, abort the program */
#endif

#ifdef LB_USE_ROW_MAJOR_MATRIX
#define MATRIX_INDEX(mat, i, j) i * mat.nCol + j
#define PMATRIX_INDEX(pMat, i, j) i * pMat->nCol + j
#else 	/* use column major matrix */
#define MATRIX_INDEX(mat, i, j) j * mat.nRow + i
#define PMATRIX_INDEX(pMat, i, j) j * pMat->nRow + i
#endif	/* LB_USE_ROW_MAJOR_MATRIX */

static void _defaultErrorReporter(const char *errMsg)
{
	// fprintf(stderr, "LinearMath: %s\n", errMsg);
}

static lbErrorReportFunc lm_reportError = _defaultErrorReporter;

void lbSetLinearMathErrorReporter(lbErrorReportFunc errorReporter)
{
	lm_reportError = errorReporter;
}

/**** vector functions ****/

lbVector lbAllocVector(int size, lbAllocator alloc)
{
	lbVector vec;
	if (size > 0)
	{
		vec.size = size;
		vec.data = (lbScalar *) alloc( size * sizeof(lbScalar) );
		if (vec.data == 0) {	/* allocation failure */
			lm_reportError("lbAllocVector(size, allocator): failed to allocate memory");
			vec.size = 0;
			HANDLE_ERROR
		}
	} else if (size == 0) {
		vec.size = 0;
		vec.data = 0;
	} else {
		lm_reportError("lbAllocVector(size, allocator): invalid size");
		HANDLE_ERROR
	}
	return vec;
}

void lbFreeVector(lbVector *vec, lbDeallocator dealloc) {
	if ( !lbVectorEmpty(*vec) ) {
		vec->size = 0;
		dealloc(vec->data);
		vec->data = 0;
	}
}

lbVector *lbAllocVectors(int size, int count, lbAllocator alloc)
{
	lbVector *vecs;
	lbSize i;
	vecs = 0;
	if (size < 0)
	{
		lm_reportError("lbAllocVectors(size, count, allocator): invalid size");
		HANDLE_ERROR
		return 0;
	}
	if (count <= 0)
	{
		lm_reportError("lbAllocVectors(size, count, allocator): invalid count");
		HANDLE_ERROR
		return 0;
	}
	if (size == 0)
	{
		vecs = (lbVector *) alloc(count * sizeof(lbVector));
		if (vecs == 0)
		{
			lm_reportError("lbAllocVectors(size, count, allocator): failed to allocate memory");
			HANDLE_ERROR
			return 0;
		}
		for (i = 0; i < count; ++i) {
			vecs[i].size = 0;
			vecs[i].data = 0;
		}
	}
	else
	{
		/* alloc for both vectors and their data */
		vecs = (lbVector *) alloc(count * ( sizeof(lbVector) + size * sizeof(lbScalar)) );
		if (vecs == 0)
		{
			lm_reportError("lbAllocVectors(size, count, allocator): failed to allocate memory");
			HANDLE_ERROR
			return 0;
		}
		vecs[0].data = (lbScalar *)(vecs + count);
		for (i = 0; i < count; ++i) {
			vecs[i].size = size;
			vecs[i].data = vecs[0].data + i * size;
		}
	}
	return vecs;
}

void lbFreeVectors(lbVector *vecs, lbDeallocator dealloc)
{
	if ( !lbVectorEmpty(vecs[0]) )
	{
		dealloc(vecs);
	}
}

lbBool lbVectorEmpty(const lbVector vec)
{
	if (vec.size > 0 && vec.data != 0)
	{
		return LB_FALSE;
	}
	else if (vec.size == 0 && vec.data == 0)
	{
		return LB_TRUE;
	}
	else
	{
		lm_reportError("lbVectorEmpty(vec): invalid vector");
		HANDLE_ERROR
		return LB_TRUE;
	}
}

lbSize lbVectorSize(const lbVector vec) {
#ifndef LB_LM_DISABLE_ARG_CHECK
	if (vec.size > 0 && vec.data != 0)
	{
		return vec.size;
	}
	else if (vec.size == 0 && vec.data == 0)
	{
		return 0;
	}
	else
	{
		lm_reportError("lbGetVectorSize(vec): invalid vector");
		HANDLE_ERROR
		return 0;
	}
#else
	return vec.size;
#endif
}

lbScalar lbGetVectorElem(const lbVector vec, int index)
{
#ifndef LB_LM_DISABLE_ARG_CHECK
	if (index >= 0 && (lbSize)index < lbVectorSize(vec)) {
		return vec.data[index];
	} else {
		lm_reportError("lbGetVectorElem(vec, index): index out of bound");
		HANDLE_ERROR
		return 0;
	}
#else
	return vec.data[index];
#endif
}

void lbSetVectorElem(lbVector *vec, int index, lbScalar value)
{
#ifndef LB_LM_DISABLE_ARG_CHECK
	if (index >= 0 && (lbSize)index < lbVectorSize(*vec)) {
		vec->data[index] = value;
	}
	else {
		lm_reportError("lbSetVectorElem(vec, index, value): index out of bound");
		HANDLE_ERROR
	}
#else
	vec->data[index] = value;
#endif
}

void lbVectorClear(lbVector *vec, lbScalar value)
{
	lbSize i, sz;
	sz = lbVectorSize(*vec);
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		vec->data[i] = value;
	}
}

void lbVectorCopy(const lbVector in, lbVector *out)
{
	lbSize sz;
	sz = lbVectorSize(in);
	if (lbVectorSize(*out) != sz) {
		lm_reportError("lbVectorCopy(in, out): in and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	memcpy(out->data, in.data, sz * sizeof(lbScalar));
}

void lbVectorAddition(const lbVector left, const lbVector right, lbVector *out)
{
	lbSize i, sz;
	sz = left.size;
	if (right.size != sz || out->size != sz)
	{
		lm_reportError("lbVectorAddition(left, right, out): left, right and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		out->data[i] = left.data[i] + right.data[i];
	}
}

void lbVectorSubtraction(const lbVector left, const lbVector right, lbVector *out)
{
	lbSize i, sz;
	sz = left.size;
	if (right.size != sz || out->size != sz)
	{
		lm_reportError("lbVectorSubtraction(left, right, out): left, right and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		out->data[i] = left.data[i] - right.data[i];
	}
}

lbScalar lbVectorDotProduct(const lbVector left, const lbVector right)
{
	lbScalar sum = 0;
	lbSize i, sz = left.size;
	if (right.size != sz) {
		lm_reportError("lbVectorDotProduct(left, right): left and right do not have the same size");
		HANDLE_ERROR
		return 0;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		sum += left.data[i] * right.data[i];
	}
	return sum;
}

lbScalar lbVectorNorm(const lbVector vec)
{
	lbScalar sum = 0;
	lbSize i, sz = vec.size;
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		sum += vec.data[i] * vec.data[i];
	}
	return sqrt(sum);
}


void lbVectorProduct(const lbVector left, const lbVector right, lbVector *out)
{
	lbSize i, sz = left.size;
	if (right.size != out->size) {
		lm_reportError("lbVectorProduct(vec, scalar, out): vec and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	if (right.size != sz) {
		lm_reportError("lbVectorProduct(left, right): left and right do not have the same size");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < right.size; ++i) {
		out->data[i] = right.data[i] * left.data[i];
	}
}


void lbVectorScalarProduct(const lbVector vec, lbScalar scalar, lbVector *out)
{
	lbSize i;
	if (vec.size != out->size) {
		lm_reportError("lbVectorScalarProduct(vec, scalar, out): vec and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < vec.size; ++i) {
		out->data[i] = vec.data[i] * scalar;
	}
}

lbScalar lbVectorSumElem(const lbVector vec)
{
	lbSize i;
	lbScalar sum;
	sum = 0;
	// #pragma omp parallel for private(i)
	for (i = 0; i < vec.size; ++i) {
		sum += vec.data[i];
	}
	return sum;
}

void lbVectorForeach(const lbVector in, lbVector *out, lbUnaryOperation uOper, void *params)
{
	lbSize i;
	if (in.size != out->size) {
		lm_reportError("lbVectorForeach(in, out, uOper, params): in and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < in.size; ++i) {
		out->data[i] = uOper(in.data[i], params);
	}
}

void lbVectorBiForeach(const lbVector vec1, const lbVector vec2, lbVector *out, lbBinaryOperation biOper, void *params)
{
	lbSize i, sz;
	sz = vec1.size;
	if (vec2.size != sz || out->size != sz) {
		lm_reportError("lbVectorBiForeach(vec1, vec2, out, biOper, params): vec1, vec2 and out"
			" do not have the same size");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		out->data[i] = biOper(vec1.data[i], vec2.data[i], params);
	}
}

/**** matrix functions ****/

lbBool lbMatrixRowMajor() {
#ifdef LB_USE_ROW_MAJOR_MATRIX
	return LB_TRUE;
#else
	return LB_FALSE;
#endif
}

lbMatrix lbAllocMatrix(int nRow, int nCol, lbAllocator alloc)
{
	lbMatrix mat;
	if (nRow > 0 && nCol > 0)
	{
		mat.nRow = nRow;
		mat.nCol = nCol;
		mat.data = (lbScalar* ) alloc( nRow * nCol * sizeof(lbScalar) );
		if (mat.data == 0) {	/* allocation failure */
			mat.nRow = mat.nCol = 0;
		}
	}
	else if (nRow == 0 && nCol == 0) {
		mat.nRow = mat.nCol = 0;
		mat.data = 0;
	}
	else
	{
		mat.nRow = mat.nCol = 0;
		mat.data = 0;
		lm_reportError("lbAllocMatrix(nRow, nCol, allocator): nRow and nCol invalid");
		HANDLE_ERROR
	}
	return mat;
}

void lbFreeMatrix(lbMatrix *mat, lbDeallocator dealloc) {
	if ( !lbMatrixEmpty(*mat) ) {
		mat->nRow = mat->nCol = 0;
		dealloc(mat->data);
		mat->data = 0;
	}
}

lbBool lbMatrixEmpty(const lbMatrix mat)
{
	if (mat.nRow > 0 && mat.nCol > 0 && mat.data != 0)
	{
		return LB_FALSE;
	}
	else if (mat.nRow == 0 && mat.nCol == 0 && mat.data == 0)
	{
		return LB_TRUE;
	}
	else
	{
		lm_reportError("lbMatrixEmpty(mat): invalid matrix");
		HANDLE_ERROR
		return LB_TRUE;
	}
}

lbSize lbMatrixRowCount(const lbMatrix mat)
{
	if (mat.nRow > 0 && mat.nCol > 0 && mat.data != 0)
	{
		return mat.nRow;
	}
	else if (mat.nRow == 0 && mat.nCol == 0 && mat.data == 0)
	{
		return 0;
	}
	else
	{
		lm_reportError("lbMatrixRowCount(mat): invalid matrix");
		HANDLE_ERROR
		return 0;
	}
}

lbSize lbMatrixColumnCount(const lbMatrix mat)
{
	if (mat.nRow > 0 && mat.nCol > 0 && mat.data != 0)
	{
		return mat.nCol;
	}
	else if (mat.nRow == 0 && mat.nCol == 0 && mat.data == 0)
	{
		return 0;
	}
	else
	{
		lm_reportError("lbMatrixColumnCount(mat): invalid matrix");
		HANDLE_ERROR
		return 0;
	}
}

lbScalar lbGetMatrixElem(const lbMatrix mat, int row, int col)
{
	if ((lbSize)row < mat.nRow && (lbSize)col < mat.nCol) {
		return mat.data[MATRIX_INDEX(mat, row, col)];
	} else {
		lm_reportError("lbGetMatrixElem(mat, row, col): row or col out of bound");
		HANDLE_ERROR
		return 0;
	}
}

void lbSetMatrixElem(lbMatrix *mat, int row, int col, lbScalar value)
{
	if ((lbSize)row < mat->nRow && (lbSize)col < mat->nCol) {
		mat->data[PMATRIX_INDEX(mat, row, col)] = value;
	} else {
		lm_reportError("lbSetMatrixElem(mat, row, col, value): row or col out of bound");
		HANDLE_ERROR
	}
}

lbBool lbMatrixSquare(const lbMatrix mat)
{
	if (mat.nRow == mat.nCol && mat.nRow != 0) {
		return LB_TRUE;
	} else {
		return LB_FALSE;
	}
}

lbBool lbMatrixSameSize(const lbMatrix mat1, const lbMatrix mat2)
{
	if (mat1.nRow == mat2.nRow && mat1.nCol == mat2.nCol)
	{
		return LB_TRUE;
	}
	else
	{
		return LB_FALSE;
	}
}

void lbMatrixCopy(const lbMatrix in, lbMatrix *out)
{
	lbSize sz;
	if (!lbMatrixSameSize(in, *out)) {
		lm_reportError("lbMatrixCopy(in, out): in and out do not have the same size");
		HANDLE_ERROR
	}
	sz = in.nRow * in.nCol;
	memcpy(out->data, in.data, sz * sizeof(lbScalar));
}

void lbMatrixTranspose(const lbMatrix in, lbMatrix *out)
{
	lbSize i, j;
	
	if (in.nRow != out->nCol || in.nCol != out->nRow) {
		lm_reportError("lbMatrixTranspose(in, out): the row or col of in and out matrix are not match");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i,j)
    for (i = 0; i < out->nRow; ++i)
	{
	    for (j = 0; j < out->nCol; ++j)
		{
			out->data[PMATRIX_INDEX(out, i, j)] = in.data[MATRIX_INDEX(in, j, i)];
		}
	} 
}


void lbMatrixAddition(const lbMatrix matL, const lbMatrix matR, lbMatrix *out)
{
	lbSize i, sz;
	if (!lbMatrixSameSize(matL, matR) || !lbMatrixSameSize(matL, *out)) {
		lm_reportError("lbMatrixAddition(matL, matR, out): matL, matR and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	sz = matL.nRow * matL.nCol;
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		out->data[i] = matL.data[i] + matR.data[i];
	}
}

void lbMatrixSubtraction(const lbMatrix matL, const lbMatrix matR, lbMatrix *out)
{
	lbSize i, sz;
	if (!lbMatrixSameSize(matL, matR) || !lbMatrixSameSize(matL, *out)) {
		lm_reportError("lbMatrixSubtraction(matL, matR, out): matL, matR and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	sz = matL.nRow * matL.nCol;
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		out->data[i] = matL.data[i] - matR.data[i];
	}
}

void lbVectorMatrixMultiply(const lbVector vec, const lbMatrix mat, lbVector *out)
{
	lbSize i, j;
	lbScalar sum;
	if (vec.size != mat.nRow) {
		lm_reportError("lbVectorMatrixMultiply(vec, mat, out): the size of vec does not "
			"match the number of rows of the matrix");
		HANDLE_ERROR
		return;
	}
	if (out->size != mat.nCol) {
		lm_reportError("lbVectorMatrixMultiply(vec, mat, out): size of out does not match "
			"the number of columns of the matrix");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i,j)
	for (i = 0; i < out->size; ++i) {
		sum = 0;
		for (j = 0; j < vec.size; ++j) {
			sum += vec.data[j] * mat.data[MATRIX_INDEX(mat, j, i)];
		}
		out->data[i] = sum;
	}
}

void lbVectorMatrixProduct(const lbVector vec, const lbMatrix mat, lbMatrix *out)
{
	lbSize i, j;
	if (vec.size != mat.nRow) {
		lm_reportError("lbVectorMatrixMultiply(vec, mat, out): the size of vec does not "
			"match the number of rows of the matrix");
		HANDLE_ERROR
		return;
	}
	if (!lbMatrixSameSize(mat, *out)) {
		lm_reportError("lbVectorMatrixProduct(vec, mat, *out): mat and out matrix do not have same size");
		HANDLE_ERROR
	}
	// #pragma omp parallel for private(i,j)

	for (i = 0; i < mat.nRow; ++i) 
	{
		for (j = 0; j < mat.nCol; ++j) 
		{
			out->data[PMATRIX_INDEX(out, i, j)] = vec.data[i] * mat.data[MATRIX_INDEX(mat, i, j)];
		}
	}
}


void lbMatrixVectorMultiply(const lbMatrix mat, const lbVector vec, lbVector *out)
{
	lbSize i, j;
	lbScalar sum;
	if (mat.nCol != vec.size) {
		lm_reportError("lbMatrixVectorMultiply(vec, mat, out): the size of vec does not "
			"match the number of columns of the matrix");
		HANDLE_ERROR
		return;
	}
	if (out->size != mat.nRow) {
		lm_reportError("lbVectorMatrixMultiply(vec, mat, out): size of out does not match "
			"the number of rows of the matrix");
		HANDLE_ERROR
		return;
	}
	// #pragma omp parallel for private(i,j)
	for (i = 0; i < out->size; ++i) {
		sum = 0;
		for (j = 0; j < mat.nCol; ++j) {
			sum += mat.data[MATRIX_INDEX(mat, i, j)] * vec.data[j];
		}
		out->data[i] = sum;
	}
}

void lbMatrixMultiply(const lbMatrix left, const lbMatrix right, lbMatrix *out)
{
	lbSize i, j, k;
	lbScalar sum;
	if (left.nCol != right.nRow) {
		lm_reportError("lbMatrixMultiply(left, right, out): the number of columns of left does not match"
			"the number of rows of right");
		HANDLE_ERROR
		return;
	}
	if (out->nRow != left.nRow) {
		lm_reportError("lbMatrixMultiply(left, right, out): the number of rows of out does not match "
			"the number of rows of left");
		HANDLE_ERROR;
		return;
	}
	if (out->nCol != right.nCol) {
		lm_reportError("lbMatrixMultiply(left, right, out): the number of columns of out does not match "
			"the number of columns of right");
		HANDLE_ERROR;
		return;
	}
	// #pragma omp parallel for private(i,j,k)
	for (i = 0; i < left.nRow; ++i)
	{
		for (j = 0; j < right.nCol; ++j)
		{
			sum = 0;
			for (k = 0; k < left.nCol; ++k)
			{
				sum += left.data[MATRIX_INDEX(left, i, k)] * right.data[MATRIX_INDEX(right, k, j)];
			}
			out->data[PMATRIX_INDEX(out, i, j)] = sum;
		}
	}
}

void lbMatrixColSum(const lbMatrix mat, lbVector *out)
{
	lbSize i, j;
	lbScalar sum;
	if (out->size != mat.nCol)
	{
		lm_reportError("lbMatrixColSum(mat, out): size of out does not match "
			"the number of cols of the matrix");
		HANDLE_ERROR;
		return;
	}
	// #pragma omp parallel for private(i,j)
	for (i = 0; i < mat.nCol; ++i)
	{
		sum = 0;
		for (j = 0; j < mat.nRow; ++j)
		{
			sum += mat.data[MATRIX_INDEX(mat, j, i)];
		}
	out->data[i] = sum;
	} 
}


void lbMatrixBiForeach(const lbMatrix mat1, const lbMatrix mat2, lbMatrix *out, lbBinaryOperation biOper, void *params)
{
	lbSize i, sz;
	if (!lbMatrixSameSize(mat1, mat2) || !lbMatrixSameSize(mat1, *out)) {
		lm_reportError("lbMatrixBiForeach(mat1, mat2, out): mat1, mat2 and out do not have the same size");
		HANDLE_ERROR
		return;
	}
	sz = mat1.nRow * mat1.nCol;
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		out->data[i] = biOper(mat1.data[i], mat2.data[i], params);
	}
}

void lbMatrixSetDiag(lbMatrix *mat, const lbVector values)
{
	lbSize i, sz;
	if ( !lbMatrixSquare(*mat) ) {
		lm_reportError("lbMatrixSetDiag(mat, value): mat should be a square matrix");
		HANDLE_ERROR
		return;
	}
	sz = lbVectorSize(values);
	if (sz != mat->nRow) {
		lm_reportError("lbMatrixSetDiag(mat, values): size of values does match the "
			"demension of mat");
		HANDLE_ERROR
		return;
	}
	memset(mat->data, 0, mat->nRow * mat->nCol * sizeof(lbScalar));
	// #pragma omp parallel for private(i)
	for (i = 0; i < sz; ++i) {
		mat->data[PMATRIX_INDEX(mat, i, i)] = values.data[i];
	}
}

void lbMatrixClearDiag(lbMatrix *mat, lbScalar value)
{
	lbSize i;
	if ( !lbMatrixSquare(*mat) ) {
		lm_reportError("lbMatrixSetDiag(mat, value): mat should be a square matrix");
		HANDLE_ERROR
		return;
	}
	memset(mat->data, 0, mat->nRow * mat->nCol * sizeof(lbScalar));
	// #pragma omp parallel for private(i)
	for (i = 0; i < mat->nRow; ++i) {
		mat->data[PMATRIX_INDEX(mat, i, i)] = value;
	}
}
