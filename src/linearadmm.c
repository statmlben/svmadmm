#include "linear_math.h"
#include <math.h>

static lbScalar _bigger(lbScalar scalar1, void *ptr)
{
	lbScalar *pScalar2 = (lbScalar *)ptr;
	return scalar1 >= *pScalar2 ? scalar1 : *pScalar2;
}

lbBool linearadmm(lbSize *Nsize, lbSize *Psize, double *xdata, double *ydata, double *tyxdata,
 double *invdata, lbScalar *lambda, lbScalar *rho, lbSize *q, lbScalar *eps1, double *outData)
{
#define xi	0
#define gamma	1
#define mu	2
#define update 0
#define wait 1

	lbScalar n, p;
	lbScalar epsR, epsS, scalar, s, r, com1, com2, com3, eps, temp;
	lbVector *Nvecs, *Pvecs, vec1, vec2, *out, *y;
	lbVector *admm, *old, *vecOne, *vecConst, *fitted, *err, *d, *tmpVec, *pvec, *beta;
	lbMatrix x, tyx, tmpMat, inv;

	n = *Nsize;
	p = *Psize;

	x.nRow = n;
	x.nCol = p;
    x.data = xdata;

	tyx.nRow = p;
	tyx.nCol = n;
    tyx.data = tyxdata;

	inv.nRow = inv.nCol = p;
	inv.data = invdata;

	vec1.size = p;
	vec1.data = outData;
	out = &vec1;

	vec2.size = n;
	vec2.data = ydata;
	y = &vec2;

	if (*lambda == 0) {
		/* will cause divide by 0 */
		return LB_FALSE;
	}
	Nvecs = lbAllocVectors(n, 12, malloc);
	if (Nvecs == 0) {
		/* out of memory */
		return LB_FALSE;
	}
	Pvecs = lbAllocVectors(p, 3, malloc);
	if (Pvecs == 0) {
		/* out of memory */
		return LB_FALSE;
	}
	beta = &Pvecs[0];
	pvec = &Pvecs[2];
	admm = &Nvecs[0];
	old = &Nvecs[3];
	fitted = &Nvecs[6];
	err = &Nvecs[7];
	d = &Nvecs[8];
	vecOne = &Nvecs[9];
	vecConst = &Nvecs[10];
	tmpVec = &Nvecs[11];


	tmpMat = lbAllocMatrix(p, n, malloc);

	if ( lbMatrixEmpty(tmpMat) ) {
		/* out of memory */
		lbFreeVectors(Nvecs, free);
		lbFreeVectors(Pvecs, free);

		return LB_FALSE;
	}

	lbVectorClear(&beta[update], 0);
	lbVectorClear(&beta[wait], 1);

	lbVectorClear(&admm[xi], 0);
	lbVectorClear(&admm[gamma], 0);
	lbVectorClear(&admm[mu], 0);

	lbVectorClear(&old[xi], 1);
	lbVectorClear(&old[gamma], 1);
	lbVectorClear(&old[mu], 1);

	lbVectorClear(fitted, 0);
	lbVectorClear(err, 0);
	lbVectorClear(d, 0);

    eps = *eps1;
	epsR = eps;
    epsS = eps;
	r = (lbScalar)10;
	s = (lbScalar)10;
	lbVectorClear(vecConst, - 1.0 / *rho);
	lbVectorClear(vecOne, 1);
	scalar = 0;
	com1 = 0;
	com2 = 0;
	com3 = 0;
	

	while ( r > epsR || s > epsS )
	{
		/* old = admm */
		lbVectorCopy(admm[xi], &old[xi]);
		lbVectorCopy(admm[gamma], &old[gamma]);
		lbVectorCopy(admm[mu], &old[mu]);
		lbVectorCopy(beta[update], &beta[wait]);

		if (*q == 2)
		{
			/* admm $ beta = as.vector(solve(t(y.tr * x.tr) %*% (y.tr * x.tr) + diag(2 * lambda / rho, p + 1)) %*% t(y.tr * x.tr) %*% (1 - admm$xi + admm$gamma - admm$mu)) */
			lbMatrixMultiply(inv, tyx, &tmpMat);
			lbVectorSubtraction(*vecOne, admm[xi], tmpVec);
			lbVectorAddition(*tmpVec, admm[gamma], tmpVec);
			lbVectorSubtraction(*tmpVec, admm[mu], tmpVec);
			lbMatrixVectorMultiply(tmpMat, *tmpVec, &beta[update]);
		}
		if (*q == 1)
		{
		    lbMatrixMultiply(inv, tyx, &tmpMat);
			lbVectorSubtraction(*vecOne, admm[xi], tmpVec);
			lbVectorAddition(*tmpVec, admm[gamma], tmpVec);
			lbVectorSubtraction(*tmpVec, admm[mu], tmpVec);
			lbMatrixVectorMultiply(tmpMat, *tmpVec, &beta[update]);
		}


		/* fitted = y.tr * (x.tr %*% admm$beta) */
		lbMatrixVectorMultiply(x, beta[update], tmpVec);
		lbVectorProduct(*tmpVec, *y, fitted);

		/* d = fitted - 1 + admm$mu */
		lbVectorSubtraction(*fitted, *vecOne, d);
		lbVectorAddition(*d, admm[mu], d);

		/* admm$xi = as.vector(pmax( -1 / rho - d, 0)) */
		lbVectorSubtraction(*vecConst, *d, &admm[xi]);
		lbVectorForeach(admm[xi], &admm[xi], _bigger, &scalar);

		/* admm$gamma = as.vector(pmax(d, 0)) */
		lbVectorForeach(*d, &admm[gamma], _bigger, &scalar);

		/* err = as.vector(fitted - 1 + admm$xi - admm$gamma) */
		lbVectorSubtraction(*fitted, *vecOne, err);
		lbVectorAddition(*err, admm[xi], err);
		lbVectorSubtraction(*err, admm[gamma], err);

		/* admm$mu = admm$mu + err */
		lbVectorAddition(old[mu], *err, &admm[mu]);

		/* r = norm( as.matrix(fitted - 1 + admm$xi - admm$gamma), type = "f") */
		lbVectorSubtraction(*fitted, *vecOne, tmpVec);
		lbVectorAddition(*tmpVec, admm[xi], tmpVec);
		lbVectorSubtraction(*tmpVec, admm[gamma], tmpVec);
		r = lbVectorNorm(*tmpVec);


		/* s = rho * norm( as.matrix( t(y.tr * x.tr) %*% (admm$xi - admm$gamma - (old$xi - old$gamma)) ), type = "f") */
        lbVectorSubtraction(admm[xi], admm[gamma], tmpVec);
		lbVectorSubtraction(*tmpVec, old[xi], tmpVec);
		lbVectorAddition(*tmpVec, old[gamma], tmpVec);
		lbMatrixVectorMultiply(tyx, *tmpVec, pvec);
		s = *rho * lbVectorNorm(*pvec);
		/* com1 = norm( as.matrix(fitted), type = "f") */
		com1 = lbVectorNorm(*fitted);

		/* com2 = norm( as.matrix(admm$xi - admm$gamma), type = "f") */
        lbVectorSubtraction(admm[xi], admm[gamma], tmpVec);
		com2 = lbVectorNorm(*tmpVec);

		/* com3 = norm( as.matrix(rep(1, n)), type = "f") */
        com3 = lbVectorNorm(*vecOne);

		/* e.r = sqrt(n) * eps + eps * max(c(com1, com2, com3)) */
        temp = com1;
		temp = temp > com2 ? temp : com2;
		temp = temp > com3 ? temp : com3;

		epsR = sqrt(n) * eps + eps * temp;
		/* e.s = sqrt(p) * eps + eps * rho * norm( as.matrix( t(y.tr * x.tr) %*% admm$mu), type = "f") */
		lbMatrixVectorMultiply(tyx, admm[mu], pvec);
		temp = lbVectorNorm(*pvec);
		epsS = sqrt(p) * eps + eps * (*rho) * temp;

	}
	lbVectorCopy(beta[update], out);

	lbFreeMatrix(&tmpMat, free);
	lbFreeVectors(Nvecs, free);
	lbFreeVectors(Pvecs, free);
	return LB_TRUE;

#undef xi
#undef gamma
#undef mu
#undef update
#undef wait
}
