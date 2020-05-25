#include "linear_math.h"
#include <math.h>

static lbScalar _bigger(lbScalar scalar1, void *ptr)
{
	lbScalar *pScalar2 = (lbScalar *)ptr;
	return scalar1 >= *pScalar2 ? scalar1 : *pScalar2;
}

lbBool lineariadmm(lbSize *Nsize, lbSize *Psize, double *xdata, double *ydata, double *tyxdata, 
	lbScalar *zeta, lbScalar *lambda, lbScalar *rho, lbSize *q, lbScalar *eps1, double *outData)
{
#define xi	0
#define gamma	1
#define mu	2
#define update 0
#define wait 1

	lbScalar n, p;
	lbScalar epsR, epsS, scalar, s, r, com1, com2, com3, eps, temp;
	lbVector *Nvecs, *Pvecs, vec1, vec2, *out, *y;
	lbVector *iadmm, *old, *vecOne, *vecConst, *fitted, *err, *d, *tmpVec, *pvec, *beta, *nu;
	lbMatrix x, tyx, tmpMat;

	n = *Nsize;
	p = *Psize;

	x.nRow = n;
	x.nCol = p;
    x.data = xdata;

	tyx.nRow = p;
	tyx.nCol = n;
    tyx.data = tyxdata;

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
	Pvecs = lbAllocVectors(p, 4, malloc);
	if (Pvecs == 0) {
		/* out of memory */
		return LB_FALSE;
	}
	beta = &Pvecs[0];
	pvec = &Pvecs[2];
	nu = &Pvecs[3];
	iadmm = &Nvecs[0];
	old = &Nvecs[3];
	fitted = &Nvecs[6];
	err = &Nvecs[7];
	d = &Nvecs[8];
	vecOne = &Nvecs[9];
	vecConst = &Nvecs[10];
	tmpVec = &Nvecs[11];


	tmpMat = lbAllocMatrix(n, p, malloc);

	if ( lbMatrixEmpty(tmpMat) ) {
		/* out of memory */
		lbFreeVectors(Nvecs, free);
		lbFreeVectors(Pvecs, free);

		return LB_FALSE;
	}

	lbVectorClear(&beta[update], 0);
	lbVectorClear(&beta[wait], 1);

	lbVectorClear(&iadmm[xi], 0);
	lbVectorClear(&iadmm[gamma], 0);
	lbVectorClear(&iadmm[mu], 0);

	lbVectorClear(&old[xi], 1);
	lbVectorClear(&old[gamma], 1);
	lbVectorClear(&old[mu], 1);

	lbVectorClear(fitted, 0);
	lbVectorClear(err, 0);
	lbVectorClear(d, 0);
	lbVectorClear(nu, 0);
	lbVectorClear(pvec, 0);

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
		/* old = iadmm */
		lbVectorCopy(iadmm[xi], &old[xi]);
		lbVectorCopy(iadmm[gamma], &old[gamma]);
		lbVectorCopy(iadmm[mu], &old[mu]);
		lbVectorCopy(beta[update], &beta[wait]);

		/* fitted = y.tr * (x.tr %*% iadmm$beta) */
		lbMatrixVectorMultiply(x, beta[update], tmpVec);
		lbVectorProduct(*tmpVec, *y, fitted);

		/* d = fitted - 1 + iadmm$mu */
		lbVectorSubtraction(*fitted, *vecOne, d);
		lbVectorAddition(*d, iadmm[mu], d);

	    /* iadmm$xi = as.vector(pmax( -1 / rho - d, 0)) */
		lbVectorSubtraction(*vecConst, *d, &iadmm[xi]);
		lbVectorForeach(iadmm[xi], &iadmm[xi], _bigger, &scalar);

		/* iadmm$gamma = as.vector(pmax(d, 0)) */
		lbVectorForeach(*d, &iadmm[gamma], _bigger, &scalar);

		/* err = as.vector(fitted - 1 + iadmm$xi - iadmm$gamma) */
		lbVectorSubtraction(*fitted, *vecOne, err);
		lbVectorAddition(*err, iadmm[xi], err);
		lbVectorSubtraction(*err, iadmm[gamma], err);

		/* iadmm$mu = iadmm$mu + err */
		lbVectorAddition(old[mu], *err, &iadmm[mu]);

		/* nu = apply((err + iadmm$mu) * y.tr * x.tr, 2, sum)*/
		lbMatrixTranspose(tyx, &tmpMat);
		lbVectorAddition(*err, iadmm[mu], tmpVec);
		lbVectorMatrixProduct(*tmpVec, tmpMat, &tmpMat);
		lbMatrixColSum(tmpMat, nu);

        if (*q == 2){
		    /* iadmm $ beta = (zeta * iadmm $ beta - rho * nu) / (2 * lambda + zeta) */
		    lbVectorScalarProduct(beta[update], *zeta, &beta[update]);
		    lbVectorScalarProduct(*nu, *rho, pvec);
		    lbVectorSubtraction(beta[update], *pvec, &beta[update]);
		    temp = 1 / (2 * *lambda + *zeta);
		    lbVectorScalarProduct(beta[update], temp, &beta[update]);
		}

		if (*q == 1){
			/* beta = soft(zeta * beta - rho * nu, lambda) / zeta */
			lbVectorScalarProduct(beta[update], *zeta, &beta[update]);
			lbVectorScalarProduct(*nu, *rho, pvec);
		    lbVectorSubtraction(beta[update], *pvec, &beta[update]);
			lbVectorForeach(beta[update], &beta[update], _bigger, lambda);
			lbVectorScalarProduct(beta[update], 1 / *zeta, &beta[update]);
		}

		/* r = norm( as.matrix(fitted - 1 + iadmm$xi - iadmm$gamma), type = "f") */
		lbVectorSubtraction(*fitted, *vecOne, tmpVec);
		lbVectorAddition(*tmpVec, iadmm[xi], tmpVec);
		lbVectorSubtraction(*tmpVec, iadmm[gamma], tmpVec);
		r = lbVectorNorm(*tmpVec);

		/* s = rho * norm( as.matrix( t(y.tr * x.tr) %*% (iadmm$xi - iadmm$gamma - (old$xi - old$gamma)) ), type = "f") */
        lbVectorSubtraction(iadmm[xi], iadmm[gamma], tmpVec);
		lbVectorSubtraction(*tmpVec, old[xi], tmpVec);
		lbVectorAddition(*tmpVec, old[gamma], tmpVec);
		lbMatrixVectorMultiply(tyx, *tmpVec, pvec);
		s = *rho * lbVectorNorm(*pvec);
		/* com1 = norm( as.matrix(fitted), type = "f") */
		com1 = lbVectorNorm(*fitted);

		/* com2 = norm( as.matrix(iadmm$xi - iadmm$gamma), type = "f") */
        lbVectorSubtraction(iadmm[xi], iadmm[gamma], tmpVec);
		com2 = lbVectorNorm(*tmpVec);

		/* com3 = norm( as.matrix(rep(1, n)), type = "f") */
        com3 = lbVectorNorm(*vecOne);

		/* e.r = sqrt(n) * eps + eps * max(c(com1, com2, com3)) */
        temp = com1;
		(temp < com2) && (temp = com2);
		(temp < com3) && (temp = com3);
    	epsR = sqrt(n) * eps + eps * temp;

		/* e.s = sqrt(p) * eps + eps * rho * norm( as.matrix( t(y.tr * x.tr) %*% iadmm$mu), type = "f") */
		lbMatrixVectorMultiply(tyx, iadmm[mu], pvec);
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


