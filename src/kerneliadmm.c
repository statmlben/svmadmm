#include "linear_math.h"

#ifdef LB_USE_ROW_MAJOR_MATRIX
#define MATRIX_INDEX(mat, i, j) i * mat.nCol + j
#define PMATRIX_INDEX(pMat, i, j) i * pMat->nCol + j
#else 	/* use column major matrix */
#define MATRIX_INDEX(mat, i, j) j * mat.nRow + i
#define PMATRIX_INDEX(pMat, i, j) j * pMat->nRow + i
#endif	/* LB_USE_ROW_MAJOR_MATRIX */

static lbScalar _abs(lbScalar scalar, void *ptr)
{
	return scalar < 0 ? -scalar : scalar;
}

static lbScalar _bigger(lbScalar scalar1, void *ptr)
{
	lbScalar *pScalar2 = (lbScalar *)ptr;
	return scalar1 >= *pScalar2 ? scalar1 : *pScalar2;
}

lbBool kerneliadmm(lbSize *size, double *kdata, lbScalar *zeta, lbScalar *lambda,
	lbScalar *rho, lbScalar *eps1, double *outData)
{
#define ALPHA	0
#define ETA		1
#define MU		2

	lbScalar sz;
	lbScalar eps, scalar, sum1, sum2, ratio;
	lbVector *vecs, vec, *out;
	lbVector *iadmm, *old, *diff, *alphaAbs, *vecC, *vecOne, *vecRho, *tmpVec;
	lbMatrix kmat, tmpMat;
	sz = *size;
	kmat = lbAllocMatrix(sz, sz, malloc);
	kmat.nRow = kmat.nCol = sz;
	kmat.data = kdata;
	vec.size = sz;
	vec.data = outData;
	out = &vec;

	if (*lambda == 0) {
		/* will cause divide by 0 */
		return LB_FALSE;
	}
	vecs = lbAllocVectors(sz, 12, malloc);
	if (vecs == 0) {
		/* out of memory */
		return LB_FALSE;
	}
	iadmm = &vecs[0];
	old = &vecs[3];
	diff = &vecs[6];
	alphaAbs = &vecs[7];
	vecC = &vecs[8];
	vecOne = &vecs[9];
	vecRho = &vecs[10];
	tmpVec = &vecs[11];

	tmpMat = lbAllocMatrix(sz, sz, malloc);
	if ( lbMatrixEmpty(tmpMat) ) {
		/* out of memory */
		lbFreeVectors(vecs, free);
		return LB_FALSE;
	}

	lbVectorClear(&iadmm[ALPHA], 0);
	lbVectorClear(&iadmm[ETA], 0);
	lbVectorClear(&iadmm[MU], 0);

	lbVectorClear(&old[ALPHA], 1);
	lbVectorClear(&old[ETA], 1);
	lbVectorClear(&old[MU], 1);

	eps = *eps1;
	lbVectorClear(vecC, 1.0 / *lambda);
	lbVectorClear(vecOne, 1);
	lbVectorClear(vecRho, *rho);
	lbVectorClear(tmpVec, 1);
	scalar = 0;

	/* temp = diag(zeta, n) - k */
	lbMatrixClearDiag(&tmpMat, *zeta);
	lbMatrixSubtraction(tmpMat, kmat, &tmpMat);

	/* abs(iadmm$alpha - old$alpha) */
	lbVectorSubtraction(iadmm[ALPHA], old[ALPHA], diff);
	lbVectorForeach(*diff, diff, _abs, NULL);
	/* abs(iadmm$alpha) */
	lbVectorForeach(iadmm[ALPHA], alphaAbs, _abs, NULL);

	sum1 = lbVectorSumElem(*diff);
	sum2 = lbVectorSumElem(*alphaAbs);
	ratio = sum1 / sum2;

	while ( ratio > eps )
	{
		/* old = iadmm */
		lbVectorCopy(iadmm[ALPHA], &old[ALPHA]);
		lbVectorCopy(iadmm[ETA], &old[ETA]);
		lbVectorCopy(iadmm[MU], &old[MU]);
		lbMatrixVectorMultiply(kmat, iadmm[ALPHA], tmpVec);

		/* iadmm$eta = pmax(c - iadmm$alpha - iadmm$mu, 0) */
		lbVectorSubtraction(*vecC, iadmm[ALPHA], &iadmm[ETA]);
		lbVectorSubtraction(iadmm[ETA], iadmm[MU], &iadmm[ETA]);
		lbVectorForeach(iadmm[ETA], &iadmm[ETA], _bigger, &scalar);

		/* iadmm$mu = old$mu + ( - c + iadmm$alpha + iadmm$eta) */
		lbVectorSubtraction(old[MU], *vecC, &iadmm[MU]);
		lbVectorAddition(iadmm[MU], iadmm[ALPHA], &iadmm[MU]);
		lbVectorAddition(iadmm[MU], iadmm[ETA], &iadmm[MU]);

		/* iadmm$alpha = as.vector( 1 / (zeta + rho) * (zeta * iadmm$alpha - k %*% old$alpha + 1 - rho * (iadmm$eta - c + iadmm$mu))) */
		lbMatrixVectorMultiply(tmpMat, old[ALPHA], &iadmm[ALPHA]);
		lbVectorAddition(iadmm[ALPHA], *vecOne, &iadmm[ALPHA]);
		lbVectorSubtraction(iadmm[ETA], *vecC, tmpVec);
		lbVectorAddition(*tmpVec, iadmm[MU], tmpVec);
		lbVectorScalarProduct(*tmpVec, *rho, tmpVec);
		lbVectorSubtraction(iadmm[ALPHA], *tmpVec, &iadmm[ALPHA]);
		lbVectorScalarProduct(iadmm[ALPHA], 1.0f / (*zeta + *rho), &iadmm[ALPHA]);

		/* abs(iadmm$alpha - old$alpha) */
		lbVectorSubtraction(iadmm[ALPHA], old[ALPHA], diff);
		lbVectorForeach(*diff, diff, _abs, NULL);
		/* abs(iadmm$alpha) */
		lbVectorForeach(iadmm[ALPHA], alphaAbs, _abs, NULL);

		sum1 = lbVectorSumElem(*diff);
		sum2 = lbVectorSumElem(*alphaAbs);
		ratio = sum1 / sum2;
	}
	lbVectorCopy(iadmm[ALPHA], out);
	lbFreeVectors(vecs, free);
	lbFreeMatrix(&tmpMat, free);
	return LB_TRUE;

#undef ALPHA
#undef ETA
#undef MU
}
