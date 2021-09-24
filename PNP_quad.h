#ifndef PNP_QUAD_H
#define PNP_QUAD_H

#include<stdarg.h>
#include"PNP.h"
#include"PNP_func.h"
#include"PNP_coefficient.h"

/* include element quad and face quad */

//begin element quad

//1 bas

//----without gradient
/* FLOAT PNP_Quad_1_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int n, int order); */
#define	PNP_Quad_1_D_Bas(e, func, dof1, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, NULL, NULL, NULL, NULL, u, n, order)
/* FLOAT PNP_Quad_2_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int n, int order); */
#define	PNP_Quad_2_D_Bas(e, func, dof1, dof2, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, dof2, NULL, NULL, NULL, u, n, order)
/* FLOAT PNP_Quad_3_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int n, int order); */
#define	PNP_Quad_3_D_Bas(e, func, dof1, dof2, dof3, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, dof2, dof3, NULL, NULL, u, n, order)
/* FLOAT PNP_Quad_4_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int n, int order); */
#define	PNP_Quad_4_D_Bas(e, func, dof1, dof2, dof3, dof4, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, dof2, dof3, dof4, NULL, u, n, order)
FLOAT PNP_Quad_5_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order);

//----with gradient
/* FLOAT PNP_Quad_1_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_G_D_G_B(e, func, dof1, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, NULL, NULL, NULL, NULL, u, v, n, order)
/* FLOAT PNP_Quad_2_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_G_D_G_B(e, func, dof1, dof2, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, dof2, NULL, NULL, NULL, u, v, n, order)
/* FLOAT PNP_Quad_3_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_G_D_G_B(e, func, dof1, dof2, dof3, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, dof2, dof3, NULL, NULL, u, v, n, order)
/* FLOAT PNP_Quad_4_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_4_D_G_D_G_B(e, func, dof1, dof2, dof3, dof4, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, dof2, dof3, dof4, NULL, u, v, n, order)
FLOAT PNP_Quad_5_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, DOF *v, int n, int order);

//2 bas

//----without gradient bas
/* FLOAT PNP_Quad_1_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_Bas_Bas(e, func, dof1, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, NULL, NULL, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_2_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_Bas_Bas(e, func, dof1, dof2, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, dof2, NULL, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_3_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_Bas_Bas(e, func, dof1, dof2, dof3, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, dof2, dof3, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_4_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_4_D_Bas_Bas(e, func, dof1, dof2, dof3, dof4, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, dof2, dof3, dof4, NULL, u, m, v, n, order)
FLOAT PNP_Quad_5_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order);

//----with 2 gradient bas
/* FLOAT PNP_Quad_1_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_G_B_G_B(e, func, dof1, u, m, v, n, order, P) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, NULL, NULL, NULL, NULL, u, m, v, n, order, P)
/* FLOAT PNP_Quad_2_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_G_B_G_B(e, func, dof1, dof2, u, m, v, n, order, P) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, dof2, NULL, NULL, NULL, u, m, v, n, order, P)
/* FLOAT PNP_Quad_3_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_G_B_G_B(e, func, dof1, dof2, dof3, u, m, v, n, order, P) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, dof2, dof3, NULL, NULL, u, m, v, n, order, P)
/* FLOAT PNP_Quad_4_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_4_D_G_B_G_B(e, func, dof1, dof2, dof3, dof4, u, m, v, n, order, P) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, dof2, dof3, dof4, NULL, u, m, v, n, order, P)
FLOAT PNP_Quad_5_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order, int P);

//----with 1 gradient bas
/* FLOAT PNP_Quad_1_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int m, DOF *gdof, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_Bas_G_D_G_B(e, func, dof1, u, m, gdof, v, n, order) \
	PNP_Quad_4_D_Bas_G_D_G_B(e, func, dof1, NULL, NULL, NULL, u, m, gdof, v, n, order)
/* FLOAT PNP_Quad_2_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *gdof, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_Bas_G_D_G_B(e, func, dof1, dof2, u, m, gdof, v, n, order) \
	PNP_Quad_4_D_Bas_G_D_G_B(e, func, dof1, dof2, NULL, NULL, u, m, gdof, v, n, order)
/* FLOAT PNP_Quad_3_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int m, DOF *gdof, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_Bas_G_D_G_B(e, func, dof1, dof2, dof3, u, m, gdof, v, n, order) \
	PNP_Quad_4_D_Bas_G_D_G_B(e, func, dof1, dof2, dof3, NULL, u, m, gdof, v, n, order)
FLOAT PNP_Quad_4_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *gdof, DOF *v, int n, int order);

//----for SUPG
FLOAT PNP_Quad_SUPG_P1(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order);
//FLOAT PNP_Quad_SUPG_Pn(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order);
FLOAT PNP_Quad_SUPG_Pn2(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order);
FLOAT PNP_Quad_SUPG_P1_2_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *gdof1, DOF *u, int m, DOF * dof2, DOF *gdof2, DOF *v, int n, int order);
FLOAT PNP_Quad_SUPG_P1_2_D_2_G_D_2_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *gdof1, DOF *gdof2, DOF *u, int m , DOF *v, int n, int order);
FLOAT PNP_Quad_SUPG_P1_2_D_2_G_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *gdof1, DOF *gdof2, DOF *gdof3, DOF *u, int m, int order);
FLOAT PNP_Quad_SUPG_P1_2_D_2_G_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *gdof1, DOF *gdof2, DOF *u, int m, DOF *gdof3, DOF *v, int n, int order);
//dim of every dof is 1
//at least 1 dof
FLOAT PNP_Quad_5_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1;
	int i;
	FLOAT tmp;
	FLOAT d;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		order = i1 + i2 + i3 + i4 + i5 + j1;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetBasisValues(e, u, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(x1++) * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(x1++) * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(x1++) * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(x1++) * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * *(x1++) * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_5_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1 ,j2;
	int i, j;
	FLOAT tmp;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *x2, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetDofValues(e, u, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * d0 * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * d0 * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * d0 * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * d0 * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * d0 * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_5_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1, j2;
	int i;
	FLOAT tmp;
	FLOAT d;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *x2, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetBasisValues(e, u, m, quad);
	x2 = phgQuadGetBasisValues(e, v, n, quad);
	
	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(x1++) * *(x2++) * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(x1++) * *(x2++) * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(x1++) * *(x2++) * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(x1++) * *(x2++) * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * *(x1++) * *(x2++) * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_5_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order, int P) {
	int i1, i2 ,i3, i4, i5, j1 ,j2;
	int i, j;
	FLOAT tmp;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *x2, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetBasisGradient(e, u, m, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * d0 * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
         if(P)
         phgPrintf("\nphi=%le ", tmp);
			func(&tmp);
         if(P)
         phgPrintf("exp(phi)=%le\n", tmp);
			d += tmp * *(v2++) * d0 * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * d0 * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * d0 * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * d0 * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_4_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *gdof, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, j1 ,j2, j3;
	int i, j;
	FLOAT tmp;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *v3, *v4, *x1, *x2, *x3, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
				}
				else {
					i4 = 0;
				}
			}
			else {
				i3 = i4 = 0;
			}
		}
		else {
			i2 = i3 = i4 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(gdof, e);
		j3 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + j1 + j2 + j3;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
			}
		}
	}
	x1 = phgQuadGetBasisValues(e, u, m, quad);
	x2 = phgQuadGetDofValues(e, gdof, quad);
	x3 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * d0 * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * *(v2++) * d0 * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * *(v2++) * *(v3++) * d0 * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * *(v2++) * *(v3++) * *(v4++) * d0 * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

//----SUPG quad
//SUPG for P1
/*
FLOAT PNP_Quad_SUPG_P1(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2, j1, j2;
	int i, j;
	FLOAT d, d0, d1;
	const FLOAT *v1, *v2, *x1, *x2, *w; 
	QUAD *quad;
	if(order<0){
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e) - 1;
		j2 = DofTypeOrder(v, e) - 1;
		order = i1 + 2 * i2 + j1 + j2;	
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	v2 = phgQuadGetDofValues(e, dof2, quad);
	x1 = phgQuadGetBasisGradient(e, u, m, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	d = 0.;
	for(i = 0; i < quad->npoints; i++, v2 += 3, x1 += 3, x2 += 3){
		d0 = x1[0] * v2[0] + x1[1] * v2[1] + x1[2] * v2[2];
		d1 = v2[0] * x2[0] + v2[1] * x2[1] + v2[2] * x2[2];
		d += *(v1++) * d0 * d1 * *(w++);
	}

	return d * phgGeomGetVolume(u->g, e);
}
*/

FLOAT PNP_Quad_SUPG_P1(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2, j1, j2;
	int i, j;
	FLOAT d, d0, d1;
	const FLOAT *v1, *v2, *x1, *x2, *w; 
	QUAD *quad;
	if(order<0){
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e) - 1;
		j2 = DofTypeOrder(v, e) - 1;
		order = 2 * i1 + 2 * i2 + j1 + j2;	
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	v2 = phgQuadGetDofValues(e, dof2, quad);
	x1 = phgQuadGetBasisGradient(e, u, m, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	d = 0.;
	for(i = 0; i < quad->npoints; i++, v2 += 3, x1 += 3, x2 += 3){
		d0 = x1[0] * v2[0] + x1[1] * v2[1] + x1[2] * v2[2];
		d1 = v2[0] * x2[0] + v2[1] * x2[1] + v2[2] * x2[2];
		d += (*v1) * d0 * (*v1) * d1 * *(w++);
		v1 += 1;
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_SUPG_P1_2_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *gdof1, DOF *u, int m ,DOF *dof2, DOF *gdof2, DOF *v, int n, int order) {
int i1, i2, j1, j2, k1, k2;
int i, j;
FLOAT d, d0, d1;
FLOAT tmp;
const FLOAT *v1, *v2, *gv1, *gv2, *x1, *x2, *w;
QUAD *quad;
if(order < 0)
{
   i1 = DofTypeOrder(dof1, e);
   i2 = DofTypeOrder(dof2, e);
   j1 = DofTypeOrder(gdof1, e);
   j2 = DofTypeOrder(gdof2, e);
   k1 = DofTypeOrder(u, e) - 1;
   k2 = DofTypeOrder(v, e) - 1;
   order = i1 + i2 + j1 + j2 + k1 + k2;
}
quad = phgQuadGetQuad3D(order);

w = quad->weights;

v1 = phgQuadGetDofValues(e, dof1, quad);
v2 = phgQuadGetDofValues(e, dof2, quad);
gv1 = phgQuadGetDofValues(e, gdof1, quad);
gv2 = phgQuadGetDofValues(e, gdof2, quad);
x1 = phgQuadGetBasisGradient(e, u, m, quad);
x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

d = 0.0;

for(i = 0; i < quad->npoints; i++, gv1 += 3, gv2 += 3, x1 += 3, x2 += 3)
{
   tmp = *(v1++);
   func(&tmp);
   d0 = gv1[0] * x1[0] + gv1[1] * x1[1] + gv1[2] * x1[2];
   d1 = gv2[0] * x2[0] + gv2[1] * x2[1] + gv2[2] * x2[2];
   d += tmp * *(v2++) * d0 * d1 * *(w++);
}

return d * phgGeomGetVolume(u->g, e);
}
FLOAT PNP_Quad_SUPG_P1_2_D_2_G_D_2_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *gdof1, DOF *gdof2, DOF *u, int m , DOF *v, int n, int order) {
int i1, i2, j1, j2, k1, k2;
int i, j;
FLOAT d, d0, d1;
FLOAT tmp;
const FLOAT *v1, *v2, *gv1, *gv2, *x1, *x2, *w;
QUAD *quad;
if(order < 0)
{
   i1 = DofTypeOrder(dof1, e);
   i2 = DofTypeOrder(dof2, e);
   j1 = DofTypeOrder(gdof1, e);
   j2 = DofTypeOrder(gdof2, e);
   k1 = DofTypeOrder(u, e) - 1;
   k2 = DofTypeOrder(v, e) - 1;
   order = i1 + i2 + j1 + j2 + k1 + k2;
}
quad = phgQuadGetQuad3D(order);

w = quad->weights;

v1 = phgQuadGetDofValues(e, dof1, quad);
v2 = phgQuadGetDofValues(e, dof2, quad);
gv1 = phgQuadGetDofValues(e, gdof1, quad);
gv2 = phgQuadGetDofValues(e, gdof2, quad);
x1 = phgQuadGetBasisGradient(e, u, m, quad);
x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

d = 0.0;

for(i = 0; i < quad->npoints; i++, gv1 += 3, gv2 += 3, x1 += 3, x2 += 3)
{
   tmp = *(v1++);
   func(&tmp);
   d0 = gv1[0] * gv2[0] + gv1[1] * gv2[1] + gv1[2] * gv2[2];
   d1 = x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2];
   d += tmp * d0 * *(v2++) * d1 * *(w++);
}

return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_SUPG_P1_2_D_2_G_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *gdof1, DOF *gdof2, DOF *gdof3, DOF *u, int m, int order) {
int i1, i2, j1, j2, j3, k1;
int i, j;
FLOAT d, d0, d1;
FLOAT tmp;
const FLOAT *v1, *v2, *gv1, *gv2, *gv3, *x1, *w;
QUAD *quad;
if(order < 0)
{
   i1 = DofTypeOrder(dof1, e);
   i2 = DofTypeOrder(dof2, e);
   j1 = DofTypeOrder(gdof1, e);
   j2 = DofTypeOrder(gdof2, e);
   j3 = DofTypeOrder(gdof3, e);
   k1 = DofTypeOrder(u, e) - 1;
   order = i1 + i2 + j1 + j2 + j3 + k1;
}
quad = phgQuadGetQuad3D(order);

w = quad->weights;

v1 = phgQuadGetDofValues(e, dof1, quad);
v2 = phgQuadGetDofValues(e, dof2, quad);
gv1 = phgQuadGetDofValues(e, gdof1, quad);
gv2 = phgQuadGetDofValues(e, gdof2, quad);
gv3 = phgQuadGetDofValues(e, gdof3, quad);
x1 = phgQuadGetBasisGradient(e, u, m, quad);

	if(func == NULL)
		func = func_id;

d = 0.0;

for(i = 0; i < quad->npoints; i++, gv1 += 3, gv2 += 3, gv3 += 3, x1 += 3)
{
   tmp = *(v1++);
   func(&tmp);
   d0 = gv1[0] * gv2[0] + gv1[1] * gv2[1] + gv1[2] * gv2[2];
   d1 = gv3[0] * x1[0] + gv3[1] * x1[1] + gv3[2] * x1[2];
   d += tmp * d0 * *(v2++) * d1 * *(w++);
}

return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_SUPG_P1_2_D_2_G_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *gdof1, DOF *gdof2, DOF *u, int m, DOF *gdof3, DOF *v, int n, int order) {
int i1, i2, j1, j2, j3, k1, k2;
int i, j;
FLOAT d, d0, d1;
FLOAT tmp;
const FLOAT *v1, *v2, *gv1, *gv2, *gv3, *x1, *x2,*w;
QUAD *quad;
if(order < 0)
{
   i1 = DofTypeOrder(dof1, e);
   i2 = DofTypeOrder(dof2, e);
   j1 = DofTypeOrder(gdof1, e);
   j2 = DofTypeOrder(gdof2, e);
   j3 = DofTypeOrder(gdof3, e);
   k1 = DofTypeOrder(u, e);
   k1 = DofTypeOrder(v, e) - 1;
   order = i1 + i2 + j1 + j2 + j3 + k1;
}
quad = phgQuadGetQuad3D(order);

w = quad->weights;

v1 = phgQuadGetDofValues(e, dof1, quad);
v2 = phgQuadGetDofValues(e, dof2, quad);
gv1 = phgQuadGetDofValues(e, gdof1, quad);
gv2 = phgQuadGetDofValues(e, gdof2, quad);
gv3 = phgQuadGetDofValues(e, gdof3, quad);
x1 = phgQuadGetBasisValues(e, u, m, quad);
x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

d = 0.0;

for(i = 0; i < quad->npoints; i++, gv1 += 3, gv2 += 3, gv3 += 3, x2 += 3)
{
   tmp = *(v1++);
   func(&tmp);
   d0 = gv1[0] * gv2[0] + gv1[1] * gv2[1] + gv1[2] * gv2[2];
   d1 = gv3[0] * x2[0] + gv3[1] * x2[1] + gv3[2] * x2[2];
   d += tmp * *(v2++) * d0 * *(x1++) * d1 * *(w++);
}
return d * phgGeomGetVolume(u->g, e);
}
//SUPG for Pn(n>1)
/*
FLOAT PNP_Quad_SUPG_Pn(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2, j1, j2;
	int i;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *x1, *x2, *w; 
	QUAD *quad;
	if(order<0){
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e) - 2;
		j2 = DofTypeOrder(v, e) - 1;
		order = 2 * i1 + i2 + j1 + j2;
		if(order < 0)
			order = 0;	
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	v2 = phgQuadGetDofValues(e, dof2, quad);
	x1 = phgQuadGetBasisHessian(e, u, m, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);
 
	d = 0.;
	for(i = 0; i < quad->npoints; i++, v1 +=1, v2 += 3, x1 += 6, x2 += 3) {
		d0 = v2[0] * x2[0] + v2[1] * x2[1] +v2[2] * x2[2];
		d += *(v1 + 1) * (x1[0] + x1[2] + x1[5]) * *(v1 + 1) * d0 * *(w++);
	}

	return d * phgGeomGetVolume(u->g, e);
}
*/
FLOAT PNP_Quad_SUPG_Pn2(ELEMENT *e, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2, j1, j2;
	int i, j;
	FLOAT d, d0, lambda[Dim + 1], v3;
	const FLOAT *v1, *v2, *x1, *x2, *w, *p; 
	QUAD *quad;
	if(order<0){
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e) - 1;
		order = 2 * i1 + (2 * i2 -1) + j1 + j2;	
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;
	p = quad->points;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	v2 = phgQuadGetDofValues(e, dof2, quad);
	x1 = phgQuadGetBasisValues(e, u, m, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	d = 0.;
	for(i = 0; i < quad->npoints; i++, v1 += 1, v2 += 3, x2 += 3) {
		for(j = 0; j < Dim + 1; j++) {
			lambda[j] = *(p++);
		}
		phgDofEvalDivergence(dof2, e, lambda, NULL, &v3);
		d0 = v2[0] * x2[0] + v2[1] * x2[1] +v2[2] * x2[2];
		d += *(v1 + 1) * *(x1++) * v3 * *(v1 + 1) * d0 * *(w++);
	}

	return d * phgGeomGetVolume(u->g, e);
}
	

//end element quad

/*---------------------------------------------------------------------------*/

//begin face quad

//2 bas

/* FLOAT PNP_Quad_Face_1_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_1_D_Bas_Bas(e, face, func, dof1, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, NULL, NULL, NULL, NULL, u, n, v, m, order) 
/* FLOAT PNP_Quad_Face_2_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_2_D_Bas_Bas(e, face, func, dof1, dof2, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, dof2, NULL, NULL, NULL, u, n, v, m, order) 
/* FLOAT PNP_Quad_Face_3_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_3_D_Bas_Bas(e, face, func, dof1, dof2, dof3, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, dof2, dof3, NULL, NULL, u, n, v, m, order) 
/* FLOAT PNP_Quad_Face_4_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_4_D_Bas_Bas(e, face, func, dof1, dof2, dof3, dof4, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, dof2, dof3, dof4, NULL, u, n, v, m, order) 
FLOAT PNP_Quad_Face_5_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, DOF *v, int m, int order);

//1 bas

/* FLOAT PNP_Quad_Face_1_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *u, int n, int order); */
#define PNP_Quad_Face_1_D_Bas(e, face, func, dof1, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, NULL, NULL, NULL, NULL, u, n, order) 
/* FLOAT PNP_Quad_Face_2_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int n, int order); */
#define PNP_Quad_Face_2_D_Bas(e, face, func, dof1, dof2, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, dof2, NULL, NULL, NULL, u, n, order) 
/* FLOAT PNP_Quad_Face_3_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int n, int order); */
#define PNP_Quad_Face_3_D_Bas(e, face, func, dof1, dof2, dof3, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, dof2, dof3, NULL, NULL, u, n, order) 
/* FLOAT PNP_Quad_Face_4_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int n, int order); */
#define PNP_Quad_Face_4_D_Bas(e, face, func, dof1, dof2, dof3, dof4, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, dof2, dof3, dof4, NULL, u, n, order) 
FLOAT PNP_Quad_Face_5_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order);

//for Face Stabilization
FLOAT PNP_Quad_Face_Stab_Jump(ELEMENT *e, int face, DOF *u, int m, DOF *v, int n, int order);
FLOAT PNP_Quad_Face_Stab_P1(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order);
//FLOAT PNP_Quad_Face_Stab_Pn1(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order);
FLOAT PNP_Quad_Face_Stab_Pn2(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order);
FLOAT PNP_Quad_Face_5_D_G_D_G_B(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, DOF *v, int n, int order);
#define PNP_Quad_Face_1_D_G_D_G_B(e, face, dof1, u, v, n, order) \
        PNP_Quad_Face_5_D_G_D_G_B(e, face, dof1, NULL, NULL, NULL, NULL, u, v, n, order)
#define PNP_Quad_Face_2_D_G_D_G_B(e, face, dof1, dof2, u, v, n, order) \
        PNP_Quad_Face_5_D_G_D_G_B(e, face, dof1, dof2, NULL, NULL, NULL, u, v, n, order)
#define PNP_Quad_Face_3_D_G_D_G_B(e, face, dof1, dof2, dof3, u, v, n, order) \
        PNP_Quad_Face_5_D_G_D_G_B(e, face, dof1, dof2, dof3, NULL, NULL, u, v, n, order)
#define PNP_Quad_Face_4_D_G_D_G_B(e, face, dof1, dof2, dof3, dof4, u, v, n, order) \
        PNP_Quad_Face_5_D_G_D_G_B(e, face, dof1, dof2, dof3, dof4, NULL, u, v, n, order)


//dim of every dof is 1
//at least 1 dof
FLOAT PNP_Quad_Face_5_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, DOF *v, int m, int order) {
	DOF *dof[5] = {dof1, dof2, dof3, dof4, dof5};
	int i[5], count, j, j1, j2, v0, v1, v2;
	FLOAT x[5];
	FLOAT tmp, d, lambda[Dim + 1];
	FLOAT *buffer;
	const FLOAT *bas0, *bas1, *p, *w;
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);

	type_u = u->type;
	type_v = v->type;

	if (order < 0) {
		order = 0;
		for(count = 0; count < 5; count++) {
			if(dof[count] != NULL) {
				i[count] = DofTypeOrder(dof[count], e);
				order += i[count];
			}
			else {
				break;
			}
		}
		order += DofTypeOrder(u, e);
		order += DofTypeOrder(v, e);
	}
	quad = phgQuadGetQuad2D(order);

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	if (n != m || u != v)
		buffer = phgAlloc(sizeof(*buffer));
	else
		buffer = NULL;
	p = quad->points;
	w = quad->weights;

	if(func == NULL)
		func = func_id;


	d = 0.;
	for (count = 0; count < quad->npoints; count++) {
	        lambda[v0] = *(p++);
	        lambda[v1] = *(p++);
	        lambda[v2] = *(p++);
		/* Note: type->BasFuncs returns an internal static buffer, its contents
		 * are changed by the next call if type_u==type_v */
		for(j = 0; j < 5; j++) {
			if(dof[j] != NULL) {
				phgDofEval(dof[j], e, lambda, &x[j]);
			}
			else x[j] = 1;
		}
		bas0 = type_u->BasFuncs(u, e, n, n + 1, lambda);
		if (n == m && u == v) {
			bas1 = bas0;
		}
		else {
			memcpy(buffer, bas0, sizeof(*buffer));
			bas0 = buffer;
			bas1 = type_v->BasFuncs(v, e, m, m + 1, lambda);
		}
		tmp = x[0];
		func(&tmp);
		d += tmp * x[1] * x[2] * x[3] * x[4] * *(bas0++) * *(bas1++) * *(w++);
	}

	if (n != m || u != v)
		phgFree(buffer);

	return d * phgGeomGetFaceArea(u->g, e, face);
}

FLOAT PNP_Quad_Face_5_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order) {
	DOF *dof[5] = {dof1, dof2, dof3, dof4, dof5};
	int count, i[5], j, v0, v1, v2;
	FLOAT x[5];
	FLOAT tmp, d, lambda[Dim + 1];
	const FLOAT *bas, *p, *w;
	QUAD *quad;
	DOF_TYPE *type_u;

	assert(face >= 0 && face < NFace);

	type_u = u->type;

	if (order < 0) {
		order = 0;
		for(count = 0; count < 5; count++) {
			if(dof[count] != NULL) {
				i[count] = DofTypeOrder(dof[count], e);
				order += i[count];
			}
			else {
				break;
			}
		}
		order += DofTypeOrder(u, e);
	}
	quad = phgQuadGetQuad2D(order);

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	p = quad->points;
	w = quad->weights;

	if(func == NULL)
		func = func_id;


	d = 0.;
	for (count = 0; count < quad->npoints; count++) {
	        lambda[v0] = *(p++);
	        lambda[v1] = *(p++);
	        lambda[v2] = *(p++);
		/* Note: type->BasFuncs returns an internal static buffer, its contents
		 * are changed by the next call if type_u==type_v */
		for(j = 0; j < 5; j++) {
			if(dof[j] != NULL) {
				phgDofEval(dof[j], e, lambda, &x[j]);
			}
			else x[j] = 1;
		}
		bas = type_u->BasFuncs(u, e, n, n + 1, lambda);

		tmp = x[0];
		func(&tmp);
		d += tmp * x[1] * x[2] * x[3] * x[4] * *bas * *(w++);
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}


//for face stabilization
FLOAT PNP_Quad_Face_Stab_Jump(ELEMENT *e, int face, DOF *u, int m, DOF *v, int n, int order) {
	GRID *g = u->g;
	SIMPLEX *e_neigh;
	int j1, j2;
	int i, j, k, face_neigh, v0, v1, v2, v3 = -1, v0_n, v1_n, v2_n, v3_n = -1;
	FLOAT d, d0, d1;
	FLOAT gph1[3], gph1_neigh[3], gph2[3], gph2_neigh[3], lambda[Dim + 1], lambda_neigh[Dim + 1];
	FLOAT *gbas1, *gbas1_neigh, *gbas2, *gbas2_neigh;
	const FLOAT *p, *w, *n0, (*J)[Dim + 1], (*J_neigh)[Dim + 1];
	NEIGHBOUR_DATA *neigh_u, *neigh_v;
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);

	type_u = u->type;	
	type_v = v->type;

//	neigh_u = phgDofInitNeighbourData(u, NULL);
//	neigh_v = phgDofInitNeighbourData(v, NULL);

	e_neigh = e->neighbours[face];
//	face_neigh = phgOppositeFace(g, e, face, e_neigh);
	for(i = 0; i < NFace; i++) {
		if(e_neigh->neighbours[i] == e) {
			face_neigh = i; 
			break;
		}
	}
	J = (void *)(phgGeomGetJacobian(g, e));
	J_neigh = (void *)(phgGeomGetJacobian(g, e_neigh));

	if(order < 0) {
		j1 = DofTypeOrder(u, e) - 1;
		j2 = DofTypeOrder(v, e) - 1;
		order = j1 + j2;
	}
	quad = phgQuadGetQuad2D(order);

	p = quad->points;
	w = quad->weights;

	n0 = phgGeomGetFaceOutNormal(g, e, face);

	GetFaceVertices(e, face, v0, v1, v2, v3);
	assert(v3 == -1);
	lambda[face] = 0.;

	GetFaceVertices(e_neigh, face_neigh, v0_n, v1_n, v2_n, v3_n);
	assert(v3_n == -1);
	lambda_neigh[face_neigh] = 0.;

	d = 0.;
	for(i = 0; i < quad->npoints; i++) {
		lambda[v0] = *p; lambda_neigh[v0_n] = *p;
		lambda[v1] = *(p+1); lambda_neigh[v1_n] = *(p+1);
		lambda[v2] = *(p+2); lambda_neigh[v2_n] = *(p+2);
		p = p + 3;
		gbas1 = type_u->BasGrads(u, e, m, m + 1, lambda);
		for(j = 0; j < 3; j++) {
			gph1[j] = gbas1[0] * J[0][j] + gbas1[1] * J[1][j] + gbas1[2] * J[2][j] + gbas1[3] * J[3][j];
		}
		gbas1_neigh = type_u->BasGrads(u, e_neigh, m, m + 1, lambda_neigh);
		for(j = 0; j < 3; j++) {
			gph1_neigh[j] = gbas1_neigh[0] * J_neigh[0][j] + gbas1_neigh[1] * J_neigh[1][j] + gbas1_neigh[2] * J_neigh[2][j] + gbas1_neigh[3] * J_neigh[3][j];
		}
		gbas2 = type_v->BasGrads(v, e, n, n + 1, lambda);
		for(j = 0; j < 3; j++) {
			gph2[j] = gbas2[0] * J[0][j] + gbas2[1] * J[1][j] + gbas2[2] * J[2][j] + gbas2[3] * J[3][j];
		}	
		gbas2_neigh = type_v->BasGrads(v, e_neigh, n, n + 1, lambda_neigh);
		d0 = 0.; d1 = 0.; 
		for(j = 0; j < 3; j++) {
			gph2_neigh[j] = gbas2_neigh[0] * J_neigh[0][j] + gbas2_neigh[1] * J_neigh[1][j] + gbas2_neigh[2] * J_neigh[2][j] + gbas2_neigh[3] * J_neigh[3][j];
			d0 += (gph1[j] - gph1_neigh[j]) * n0[j];
			d1 += (gph2[j] - gph2_neigh[j]) * n0[j];
		}
		d += d0 * d1 * *(w++);
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}

FLOAT PNP_Quad_Face_Stab_P1(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order) {
	GRID *g = u->g;
	int i1, i2, j1, j2;
	int i, j, k, v0, v1, v2;
	FLOAT d, d0, d1;
	FLOAT gph1[3], gph2[3], lambda[Dim + 1];
	FLOAT dv1, dv2[3], *gbas1, *gbas2;
	const FLOAT *p, *w, (*J)[Dim + 1];
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);

	type_u = u->type;	
	type_v = v->type;

	J = (void *)(phgGeomGetJacobian(g, e));

	if(order < 0) {
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e) - 1;
		j2 = DofTypeOrder(v, e) - 1;
		order = 2 * i1 + 2 * i2 + j1 + j2;
	}
	quad = phgQuadGetQuad2D(order);

	p = quad->points;
	w = quad->weights;

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	d = 0.;
	for(i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		phgDofEval(dof1, e, lambda, &dv1);
		phgDofEval(dof2, e, lambda, dv2);
		gbas1 = type_u->BasGrads(u, e, m, m + 1, lambda);
		for(j = 0; j < 3; j++) {
			gph1[j] = gbas1[0] * J[0][j] + gbas1[1] * J[1][j] + gbas1[2] * J[2][j] + gbas1[3] * J[3][j];
		}
		if (n == m && u == v) {
			gbas2 = gbas1;
		}
		else {
			gbas2 = type_v->BasGrads(v, e, n, n + 1, lambda);
		}
		d0 = 0.; d1 = 0.; 
		for(j = 0; j < Dim; j++) {
			gph2[j] = gbas2[0] * J[0][j] + gbas2[1] * J[1][j] + gbas2[2] * J[2][j] + gbas2[3] * J[3][j];
			d0 += dv2[j] * gph1[j];
			d1 += dv2[j] * gph2[j];
		}
		d += dv1 * d0 * dv1 * d1 * *(w++);
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}

/*
FLOAT PNP_Quad_Face_Stab_P1(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order) {
	GRID *g = u->g;
	int i1, i2, j1, j2;
	int i, j, k, v0, v1, v2;
	FLOAT d, d0, d1;
	FLOAT gph1[3], gph2[3], lambda[Dim + 1];
	FLOAT dv1, dv2[3], *gbas1, *gbas2;
	const FLOAT *p, *w, (*J)[Dim + 1];
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);

	type_u = u->type;	
	type_v = v->type;

	J = (void *)(phgGeomGetJacobian(g, e));

	if(order < 0) {
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e) - 1;
		j2 = DofTypeOrder(v, e) - 1;
		order = i1 + 2 * i2 + j1 + j2;
	}
	quad = phgQuadGetQuad2D(order);

	p = quad->points;
	w = quad->weights;

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	d = 0.;
	for(i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		phgDofEval(dof1, e, lambda, &dv1);
		phgDofEval(dof2, e, lambda, dv2);
		gbas1 = type_u->BasGrads(u, e, m, m + 1, lambda);
		for(j = 0; j < 3; j++) {
			gph1[j] = gbas1[0] * J[0][j] + gbas1[1] * J[1][j] + gbas1[2] * J[2][j] + gbas1[3] * J[3][j];
		}
		if (n == m && u == v) {
			gbas2 = gbas1;
		}
		else {
			gbas2 = type_v->BasGrads(v, e, n, n + 1, lambda);
		}
		d0 = 0.; d1 = 0.; 
		for(j = 0; j < Dim; j++) {
			gph2[j] = gbas2[0] * J[0][j] + gbas2[1] * J[1][j] + gbas2[2] * J[2][j] + gbas2[3] * J[3][j];
			d0 += dv2[j] * gph1[j];
			d1 += dv2[j] * gph2[j];
		}
		d += dv1 * d0 * d1 * *(w++);
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}
*/
/*
FLOAT PNP_Quad_Face_Stab_Pn1(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order)
{
	GRID *g = u->g;
	int i1, i2, j1, j2;
	int i, j, k, v0, v1, v2;
	int N = u->type->nbas;
	FLOAT d, d0;
	FLOAT dv1, dv2[3], hbas[3 * N], gph[3], *gbas, lambda[Dim + 1];
	const FLOAT *w, *p, (*J)[Dim + 1];
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);

	type_u = u->type;
	type_v = v->type;

	if(order < 0) {
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e) - 2;
		j2 = DofTypeOrder(v, e) - 1;
		order = 2 * i1 + i2 + j1 + j2;
		if(order < 0)
			order = 0;
	}

	J = phgGeomGetJacobian(g, e);

	quad = phgQuadGetQuad2D(order);

	w = quad->weights;
	p = quad->points;

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	d = 0.;	
	for(i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		phgDofEval(dof1, e, lambda, &dv1);
		phgDofEval(dof2, e, lambda, dv2); 
		phgDofEvalBasisHessian(u, e, lambda, "D", hbas);
		gbas = type_v->BasGrads(v, e, n, n + 1, lambda);
		d0 = 0.;
		for(j = 0; j < Dim; j++) {
			gph[j] = 0.;
			for(k = 0; k < Dim + 1; k++)
				gph[j] += J[k][j] * gbas[k];
			d0 += dv2[j] * gph[j];
		}
		d += dv1 * (hbas[3*m] + hbas[3*m+1] + hbas[3*m+2]) * dv1 * d0 * *(w++);
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}
*/

FLOAT PNP_Quad_Face_Stab_Pn2(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order)
{
	GRID *g = u->g;
	int i1, i2, j1, j2;
	int i, j, k, v0, v1, v2;
	FLOAT d, d0;
	FLOAT dv1, dv2[3], dv3, bv1, gph[3], *bas, *gbas, lambda[Dim +1];
	const FLOAT *p, *w, (*J)[Dim + 1];
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;
	
	if(order < 0) {
		i1 = DofTypeOrder(dof1, e);
		i2 = DofTypeOrder(dof2, e);
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e) - 1;
		order = 2 * i1 + (2 * i2 - 1) + j1 + j2;
		if(order < 0)
			order = 0;
	}

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace); 

	type_u = u->type;
	type_v = v->type;

	J = phgGeomGetJacobian(g, e);

	quad = phgQuadGetQuad2D(order);	
	p = quad->points;
	w = quad->weights;

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	d = 0.;		
	for(i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		phgDofEval(dof1, e, lambda, &dv1);
		phgDofEval(dof2, e, lambda, dv2);
		phgDofEvalDivergence(dof2, e, lambda, NULL, &dv3);
		bas = type_u->BasFuncs(u, e, m, m + 1, lambda);
		bv1 = *bas;
		gbas = type_v->BasGrads(v, e, n, n + 1, lambda);
		d0 = 0.;
		for(j = 0; j < Dim; j++) {
			gph[j] = 0.;
			for(k = 0; k < Dim + 1; k++)
				gph[j] += J[k][j] * gbas[k];
			d0 += dv2[j] * gph[j];
		}
		d += dv1 * bv1 * dv3 * dv1 * d0 * *(w++);
	}
	
	return d * phgGeomGetFaceArea(u->g, e, face);
}
		

FLOAT PNP_Quad_Face_5_D_G_D_G_B(ELEMENT *e, int face, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, DOF *v, int n, int order) 
{
	GRID *g = u->g;
	DOF *dof[5] = {dof1, dof2, dof3, dof4, dof5};
	int i, j, k, count, v0, v1, v2;
	FLOAT d, d0;
	FLOAT dv[5], gv1[3], gv2[3], *gbas, lambda[Dim + 1];
	const FLOAT *p, *w, (*J)[Dim + 1];
	QUAD *quad;
	DOF_TYPE *type_v;

	assert(!SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);
	
	type_v = v->type;

	if(order < 0) {
		for(count = 0; count < 5; count++) {
			if(dof[count] != NULL)
				order += DofTypeOrder(dof[count], e);
			else
				break;
		}
		order += DofTypeOrder(u, e);
		order += DofTypeOrder(v, e) - 1;
	}

	J = phgGeomGetJacobian(g, e);

	quad = phgQuadGetQuad2D(order);
	p = quad->points;
	w = quad->weights;

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	d = 0.;
	for(i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		for(count = 0; count < 5; count++) {
			if(dof[count] != NULL)
				phgDofEval(dof[count], e, lambda, &dv[count]);
			else
				dv[count] = 1;
		}
		phgDofEval(u, e, lambda, gv1);
		gbas = type_v->BasGrads(v, e, n, n + 1, lambda);
		d0 = 0.;
		for(j = 0; j < Dim; j++) {
			gv2[j] = 0.;
			for(k = 0; k < Dim + 1; k++)
				gv2[j] += J[k][j] * gbas[k];
			d0 += gv1[j] * gv2[j];
		}
		d += dv[0] * dv[1] * dv[2] * dv[3] * dv[4] * d0 * *(w++);
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}
	

//end face quad

////////////////////////////////////////////////////////////////////////

FLOAT PNP_quad_channel_section(FLOAT *value, int d, ELEMENT *e, COORD *points, INT func_dof, PNP_FUNC_1 func, INT grad_dofs, DOF *dof, ...) {
	//points are coordinate of A, B, C
	//(the three points of the triangle ABC)
	//points[3][ ] for A B C
	//points[ ][3] for x y z
	//grad_dofs is used to calculate gradient of dof and select z axis

	GRID *g = dof->g;

	//get dofs
	DOF **dofs;
	int i, j, k, ndof;
	dofs = phgAlloc(8 * sizeof(*dofs));

	va_list ap;
	va_start(ap, dof);
	for(ndof = 0; ndof < 256; ndof++) {
		if(dof == NULL) break;
		dofs[ndof] = dof;
		dof = va_arg(ap, DOF *);
	}
	dofs[ndof] = NULL;
	va_end(ap);

	//change coordinate
	FLOAT point_lambda[3][4];
	const FLOAT *tmp_point;
	for(i = 0; i < 3; i++) {
		tmp_point = phgGeomXYZ2Lambda(g, e, points[i][0], points[i][1], points[i][2]);
		for(j = 0; j < Dim + 1; j++)
			point_lambda[i][j] = tmp_point[j];
	}

	FLOAT area;
	FLOAT *X = (FLOAT *)malloc(3 * sizeof(FLOAT));
	FLOAT *Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
	for(i = 0; i < 3; i++) {
		X[i] = points[i][(d+1)%3];
		Y[i] = points[i][(d+2)%3];
	}
	triangle_area(X, Y, &area);

	//face quad
	int order;
	QUAD *quad;
	FLOAT res, sum, lambda[Dim + 1];
	FLOAT dof_value[1], dof_grad[Dim];
	const FLOAT *p, *w;
	
	order = 0;
	for(i = 0; i < ndof; i++) {
		if(DofTypeOrder(dofs[i], e) > 0)
			order += DofTypeOrder(dofs[i], e);
	}
	order -= grad_dofs;
	//if(order < 1) order = 1;
	quad = phgQuadGetQuad2D(order);
	p = quad->points;
	w = quad->weights;
	
	sum = 0.0;
	for(i = 0; i < quad->npoints; i++) {
		for(j = 0; j < Dim + 1; j++) {
		    lambda[j] = 0.0;
		    for(k = 0; k < 3; k++) {
			lambda[j] += point_lambda[k][j] * p[k + i * 3];
			if(lambda[j] > 1.1) phgError(1, "lambda > 1");
		    }
		}
		res = 1.0;
		for(j = 0; j < ndof; j++) {
		    if(grad_dofs && !j) {
			phgDofEvalGradient(dofs[j], e, lambda, NULL, dof_grad);
			res *= dof_grad[d];
		    }
		    else {
			phgDofEval(dofs[j], e, lambda, dof_value);
			if(j + 1 == func_dof && func != NULL) {
				func(dof_value);
				res *= dof_value[0];
			}
			else {
				res *= dof_value[0];
			}
		    }
		}
		sum += res * *(w++);
	}
	
	phgFree(dofs);

	FLOAT tmp_value;
	if(value != NULL) {
		*value = sum * area;
		tmp_value = *value;
	}
	else {
		tmp_value = sum * area;
	}
	return tmp_value;
}

/****************************************************************************************************************************/
FLOAT PNP_Quad_Face_Current(ELEMENT *e, int d, int face, int which_ion, DOF *dof1, DOF *dof2, DOF *grad_dof)
{
	int order1, order2;
	int i, v0, v1, v2;
	FLOAT x1, x2, x3, x4, x_tmp[3];
	FLOAT d1, d2, lambda[Dim + 1];
	const FLOAT *p1, *p2, *w1, *w2;
	QUAD *quad1, *quad2;

	order1 = DofTypeOrder(dof1, e) + DofTypeOrder(dof2, e) - 1;
	order2 = DofTypeOrder(dof1, e) + DofTypeOrder(grad_dof, e) + DofTypeOrder(dof2, e);

	quad1 = phgQuadGetQuad2D(order1);
	quad2 = phgQuadGetQuad2D(order2);

	p1 = quad1->points;
	p2 = quad2->points;
	w1 = quad1->weights;
	w2 = quad2->weights;

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;
	
	d1 = 0.;
	for(i = 0; i < quad1->npoints; i++) {
		lambda[v0] = *(p1++);
		lambda[v1] = *(p1++);
		lambda[v2] = *(p1++);
		phgDofEval(dof1, e, lambda, &x1);
		phgDofEvalGradient(dof2, e, lambda, NULL, x_tmp);
		x4 = x_tmp[d];
		d1 += x1 * x4 * *(w1++);
	}

	d2 = 0.;
	for(i = 0; i < quad2->npoints; i++) {
		lambda[v0] = *(p2++);
		lambda[v1] = *(p2++);
		lambda[v2] = *(p2++);
		phgDofEval(dof1, e, lambda, &x1);
		phgDofEval(dof2, e, lambda, &x2); 
		phgDofEval(grad_dof, e, lambda, x_tmp);
		x3 = x_tmp[d];
		d2 += x1 * x2 * x3 * *(w2++);
	}
	d2 *= ion[which_ion].Z;
		
	return (d1 + d2) * phgGeomGetFaceArea(dof2->g, e, face);	
}	
	

#endif
