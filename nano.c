#include "phg.h"
#include "PNP_coefficient.h"
#include "PNP_analytic.c"
#include "PNP_more.c"

/* NOTE:
   1) The folloing codes implement the FVSG method,
      Zlamal FEM, Inverse Average Mixed FEM, and Exoponential 
      Basis FEM.
   2) The variables that we use are carrier concentrations.
   3) Only Gummel iteration is implemented.
   4) SRH Recombination model is implemented.
   5) The devices we simulate include PN diode, NPN bipolar junction
      transistor, and nMOSFET. The NPN transistor case is not carefully
      test, and may have problems.
*/
static BOOLEAN
BiasSweepPN(void)
{
    /* cathode grounded, anode ramps up/down */
    static BOOLEAN initilized = FALSE;

    if (!initilized) {
	initilized = TRUE;
    }

    if (Fabs(anode) + 1.0e-6 > Fabs(bias_anode))
	return TRUE;

    if (bias_anode >= bias_cathode)
	anode += h;
    else
	anode -= h;

    phgPrintf("\n  ================ BIAS: %0.4lf ================\n",
	      (double)(anode - bias_cathode));

    return FALSE;
}
static BOOLEAN
BiasSweepPN_new(void)
{
    /* cathode grounded, anode ramps up/down */
    static BOOLEAN stage1 = TRUE;
    static BOOLEAN stage2 = FALSE;
    static BOOLEAN initilized = FALSE;
    FLOAT eps = 1.0e-12;

    if(stage1){
    if (!initilized) {
	initilized = TRUE;
    }

    if (Fabs(anode) + eps > Fabs(bias_anode))
    {
         stage1 = FALSE;
         stage2 = TRUE;
         initilized = FALSE;
    }
    else{
         if (bias_anode >= bias_cathode)
	         anode += h;
         else
	         anode -= h;
      }
    }

#if 1
    if(stage2)
    {
         if (!initilized) {
	         initilized = TRUE;
         }

    if (Fabs(surface_charge) + eps > Fabs(bias_surface_charge))
	      return TRUE;
    else
       surface_charge *= h_surface_charge;
      
    }
    phgPrintf("\n  ================ anode: %0.4lf, surface charge: %0.4le, analytic_density: %0.4le ================\n",
	      (double)(anode - bias_cathode), (double)surface_charge, (double)analytic_density);
#else
    if(stage2)
    {
         if(!initilized){
            initilized = TRUE;
         }

         if(analytic_density + eps > bias_analytic_density)
            return TRUE;
         else
            analytic_density *= h_analytic_density;
    }
    phgPrintf("\n  ================ anode: %0.4lf, analytic_density: %0.4le ================\n",
	      (double)(anode - bias_cathode), (double)analytic_density);
#endif

    return FALSE;
}
static BOOLEAN
SurfaceSweepNANO(void)
{
    /* cathode grounded, anode ramps up/down */
    static BOOLEAN initilized = FALSE;

    if (!initilized) {
	initilized = TRUE;
    }

    if (Fabs(surface_charge) + 1.0e-10 > Fabs(bias_surface_charge))
	return TRUE;
   
    surface_charge += h_surface_charge;

    phgPrintf("\n  ================ BIAS: %0.8e ================\n", surface_charge);

    return FALSE;
}


static BOOLEAN
BiasSweep(void)
{
    if (device == PN)
	return BiasSweepPN_new();

    return TRUE;
}

static BOOLEAN SurfaceSweep(void)
{
   if(device == NANO || device == PN)
      return SurfaceSweepNANO();

   return TRUE;
}
static void
SetupVertsRegion(GRID *g, BOOLEAN flag)
{
    int i;

    if (flag == FALSE) {
	phgFree(verts_region);
	verts_region = NULL;
	return;
    }

    assert(verts_region == NULL);
    verts_region = phgCalloc(g->nvert, sizeof(*verts_region));

    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++) {
		*(verts_region + e->verts[i]) |= MAT_SEMICONDUCTOR;
	}
    }
}

#define B(x) (Fabs(x)>1.0e-4? (x)/(Exp(x)-1): \
	(((-1./720 *(x)*(x) + 1./12)*(x) -1./2)*(x) +1.))

FLOAT Ber(FLOAT x)
{         
   FLOAT y;
   if(Fabs(x) > 1.0e4)
   {
      if(x > 0)
      {
      y = ((-1./720 /x /x +1./12)/x -1./2)/x;
      }
      else
      y = -x;
   }
   else if(Fabs(x) < 1.0e-4)
   {
      y = 1.+((-1./720 *x *x +1./12)*x -1./2)*x;
   }
   else
   {
      y = x / (Exp(x)-1.0);
   }
   return y;
}



static void
exponential(GRID *g, SIMPLEX *e, DOF *u, EQUATION eqn, QUAD *quad, int n,
	    const FLOAT *a, FLOAT bas[], FLOAT flux[][Dim])
/* cf: Three-dimentional exponentially fitted conforming tetrahedral 
   finite elements for the semiconductor continuity equations */
{
    int i, j;
    FLOAT x[Dim];
    FLOAT lm[NVert][Dim];
    FLOAT sm[NVert], bs[NVert], nbs[NVert];
    FLOAT bas_val[NVert];
    FLOAT A[Dim][Dim], b[Dim][NVert];
    BOOLEAN flag;
    const FLOAT *lambda = quad->points + n * (Dim + 1);
    FLOAT sum;

    phgGeomLambda2XYZ(g, e, quad->points + n * (Dim + 1), &(x[0]), &(x[1]),&(x[2]));

    for (i = 0; i < NVert; i++) {
	sm[i] = 0;
	for (j = 0; j < Dim; j++) {
	    lm[i][j] = x[j] - g->verts[e->verts[i]][j];
	    if (eqn == ELECTRON)	/* electrons */
		sm[i] += a[j] * lm[i][j];
	    else if (eqn == HOLE)	/* holes */
		sm[i] += -a[j] * lm[i][j];
	}
	bs[i] = B(sm[i]);
	nbs[i] = B(-sm[i]);
//	bs[i] = Ber(sm[i]);
//	nbs[i] = Ber(-sm[i]);
    }

    sum = 0.;
    for (i = 0; i < NVert; i++)
	sum += lambda[i] * bs[i];

    for (i = 0; i < NVert; i++)
	bas_val[i] = lambda[i] * nbs[i] / sum;

    if (bas != NULL) {
	for (i = 0; i < NVert; i++)
	    bas[i] = bas_val[i];
    }

    if (flux == NULL)
	return;

    bzero(&(b[0][0]), Dim * NVert * sizeof(b[0][0]));
    for (i = 0; i < Dim; i++) {
	for (j = 0; j < Dim; j++)
	    A[i][j] = lm[i][j];
	for (j = 0; j < NVert; j++)
	    b[i][j] = bs[i] * bas_val[j];
	b[i][i] += -nbs[i];
    }

    flag = phgSolverDenseSolver(Dim, NVert, &(A[0][0]), &(b[0][0]));
    if (!flag)
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

    for (i = 0; i < NVert; i++)
	for (j = 0; j < Dim; j++)
	    flux[i][j] = b[j][i];
}

static int
cmp_float(const void *p0, const void *p1)
{
    FLOAT f = *((FLOAT *)p0) - *((FLOAT *)p1);
    return (f < 0.) ? -1 : ((f > 0.) ? 1 : 0);
}



static FLOAT *
JBJ(DOF *V, SIMPLEX *e, int n, EQUATION eqn)
{
    GRID *g = V->g;
    static FLOAT JT[Dim][Dim];
    FLOAT IJT[Dim][Dim];
    FLOAT psi0, psi1, b;
    int order[NVert];
    int i, j;

    for (i = 0; i < NVert; i++)
	order[i] = (n + i) % NVert;

    psi0 = *(DofVertexData(V, e->verts[order[0]]));
    for (i = 1; i < NVert; i++) {
	psi1 = *(DofVertexData(V, e->verts[order[i]]));
	b = eqn == ELECTRON ? B(psi0 - psi1) : B(psi1 - psi0);
//	b = eqn == HOLE ? B(psi0 - psi1) : B(psi1 - psi0);
//	b = eqn == HOLE ? Ber(psi0 - psi1) : Ber(psi1 - psi0);
//	b = eqn == ELECTRON ? Ber(psi0 - psi1) : Ber(psi1 - psi0);
	for (j = 0; j < Dim; j++) {
	    JT[i - 1][j] =
		g->verts[e->verts[order[i]]][j] -
		g->verts[e->verts[order[0]]][j];
	    IJT[i - 1][j] = JT[i - 1][j];
	    JT[i - 1][j] *= b;	/* B x JT */
	}
    }

    /* JT^-1 x B x JT */
    phgSolverDenseSolver(Dim, Dim, &(IJT[0][0]), &(JT[0][0]));

    return JT[0];
}

static void
exp_quad_dof_eval(SIMPLEX *e, DOF *u, EQUATION eqn, DOF *gradv,
		  QUAD *quad, FLOAT *val)
{
    int i, n;
    FLOAT bas_val[NVert];
    const FLOAT *a;

    a = DofElementData(gradv, e->index);
    for (n = 0; n < quad->npoints; n++) {
	exponential(u->g, e, u, eqn, quad, n, a, bas_val, NULL);
	val[n] = 0.;
	for (i = 0; i < NVert; i++)
	    val[n] += *(DofVertexData(u, e->verts[i])) * bas_val[i];
    }
}

static void
InitialGuess(DOF *V, DOF *P, DOF *N)
/* use values of potential and carrier concentrations at thermal
   equilibrium states as initial values */
{
    GRID *g = V->g;
    SIMPLEX *e;
    FLOAT  phi, v, p, n;
    BTYPE interface, btype;
    int i;
    INT ind;
    VEC *vecV, *vecP, *vecN;
    VEC *diagV, *diagP, *diagN;
    MAP *map;
//
    map = phgMapCreate(V, NULL);
    vecV = phgMapCreateVec(map, 1);
    vecP = phgMapCreateVec(map, 1);
    vecN = phgMapCreateVec(map, 1);
    diagV = phgMapCreateVec(map, 1);
    diagP = phgMapCreateVec(map, 1);
    diagN = phgMapCreateVec(map, 1);
    phgVecDisassemble(vecV);
    phgVecDisassemble(vecP);
    phgVecDisassemble(vecN);
    phgVecDisassemble(diagV);
    phgVecDisassemble(diagP);
    phgVecDisassemble(diagN);

    if (device == PN) {
	assert(cathode == anode);
	phi = anode / (KB * T);
	interface = JUNCTION_PN;
    }
    if(device == NANO)
    {
       //cathode = 0.0;
	    assert(cathode == anode);
       phi = anode / (KB * T);
       interface = JUNCTION_PN;
    }

    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++) {
	    ind = phgMapE2L(map, 0, e, i);
		btype = phgDofGetElementBoundaryType(V, e, i);
		if (btype & interface) {
		    v = phi * (g->verts[e->verts[i]][2] + analytic_period / 2) / analytic_period;
		    p = n = analytic_density * NA / Pow(DM2CM, 3) / SCALE;
		}
		else if (btype & CATHODE){
          v = cathode / (KB * T);
		    p = n = analytic_density * NA / Pow(DM2CM, 3) / SCALE;
		}
      else if (btype & ANODE)
      {
          v = phi;
		    p = n = analytic_density * NA / Pow(DM2CM, 3) / SCALE;
      }
      else
      {
//         v = Kappa_2 * analytic_density * NA /Pow(DM2CM, 3) * Pow(analytic_period / M_PI, 2) / Pow(CM2MICRON, 2) / 3 / EPSILONE;
       v = 0.0;
        p = 5.0 * nie / SCALE;
        n = 5.0 * nie / SCALE;
          //  p = 0.0;
          //  n = 0.0;
      }

		phgVecAddEntry(vecV, 0, ind, v);
		phgVecAddEntry(diagV, 0, ind, SCALE);
		phgVecAddEntry(vecP, 0, ind, p);
		phgVecAddEntry(diagP, 0, ind, 1.);
		phgVecAddEntry(vecN, 0, ind, n);
		phgVecAddEntry(diagN, 0, ind, 1.);
	}
    }

    phgVecAssemble(vecV);
    phgVecAssemble(vecP);
    phgVecAssemble(vecN);
    phgVecAssemble(diagV);
    phgVecAssemble(diagP);
    phgVecAssemble(diagN);

    for (i = 0; i < map->nlocal; i++) {
	*(vecV->data + i) /= *(diagV->data + i);
	if (*(diagP->data + i) > 0.5)
	    *(vecP->data + i) /= *(diagP->data + i);
	if (*(diagN->data + i) > 0.5)
	    *(vecN->data + i) /= *(diagN->data + i);
    }

    phgVecDestroy(&diagV);
    phgVecDestroy(&diagP);
    phgVecDestroy(&diagN);

    phgMapVecToDofArrays(map, vecV, FALSE, &V, NULL);
    phgMapVecToDofArrays(map, vecP, FALSE, &P, NULL);
    phgMapVecToDofArrays(map, vecN, FALSE, &N, NULL);

    phgMapDestroy(&map);
    phgVecDestroy(&vecV);
    phgVecDestroy(&vecP);
    phgVecDestroy(&vecN);
}

static void
TakeBoundary(DOF *V, DOF *P, DOF *N, int who)
    /* who == 0 ==> take bdry for V,
     * who == 1 ==> take bdry for P and N,
     * who == 2 ==> take bdry for V, P, and N */
{
    GRID *g = V->g;
    SIMPLEX *e;
    BTYPE btype;
    FLOAT phi, bias;
    int i;

    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++) {
	    if (phgDofGetElementBoundaryType(V, e, i) & ohmic_contact)
		break;
       if(analytic_test)
       {
	    if (phgDofGetElementBoundaryType(V, e, i) & JUNCTION_PN)
		break;
       }
	}
	if (i >= NVert)		/* no bdry d.o.f */
	    continue;


	for (i = 0; i < NVert; i++) {
	    btype = phgDofGetElementBoundaryType(V, e, i);
	    if (btype & ohmic_contact)
       {
	    if (device == PN)
   		bias = (btype & CATHODE) ? cathode : anode;
       if(device == NANO)
          bias = anode;
       }
       else{
         if(analytic_test && (btype & JUNCTION_PN))
         {
	       if (device == PN)
		    bias = anode * (g->verts[e->verts[i]][2] + analytic_period / 2) / analytic_period;
	       if (device == NANO)
		    bias = anode * (g->verts[e->verts[i]][2] + analytic_period / 2) / analytic_period;
         }
         else{
		         continue;
         }
       }

	    phi = bias / (KB * T);
       
	    if (who == 1 || who == 2) {
		*(P->data + phgDofMapE2D(P, e, i)) = analytic_density * NA / Pow(DM2CM, 3) / SCALE;
		*(N->data + phgDofMapE2D(N, e, i)) = analytic_density * NA / Pow(DM2CM, 3) / SCALE;
	    }
	    if (who == 0 || who == 2)
		*(V->data + phgDofMapE2D(V, e, i)) = phi;
	}
    }
}

FLOAT Average2(FLOAT x)
{
   if(Fabs(x)>1.0e4)
   {
      if(x < 0.0)
         return 0.0;
      else 
         return 1.0;
   }
      return Exp(x) / (1 + Exp(x));
}

FLOAT Average3_0(FLOAT x)
{
   FLOAT alpha = 0.5;
#if 1
   if(Fabs(x)<1.0e-4)
      alpha = 0.5;
   else
   {
      if(x > 0)
         alpha = 1. / x;
      else
         alpha = 1 + 1. /x;
   }
#endif
      return Exp(x) / (1 - alpha + alpha * Exp(x));
}
FLOAT Average3_1(FLOAT x)
{
   FLOAT alpha = 0.5;
#if 1
   if(Fabs(x)<1.0e-4)
      alpha = 0.5;
   else
   {
      if(x > 0)
         alpha = 1. / x;
      else
         alpha = 1 + 1. / x;
   }
#endif
      return 1.0 / (1 - alpha + alpha * Exp(x));
}

static void
average(FLOAT local_v[4], EQUATION eqn, FLOAT average_coeff[2], int psi0, int psi1)
{
    FLOAT s0, s1;
    FLOAT p[4];
    FLOAT p1, p2, p3, p4;
    FLOAT d41, d34, d31, d23, d24, d21;
    FLOAT v0, v1;
    FLOAT eps = 1.0e-3;
    FLOAT r;
    int i;

    for (i = 0; i < 4; i++) {
	if (eqn == ELECTRON)
	    p[i] = -local_v[i];
	else
	    p[i] = +local_v[i];
    }

    s0 = p[psi0];
    s1 = p[psi1];

    qsort(p, 4, sizeof(p[0]), cmp_float);
    p1 = p[0];
    p2 = p[1];
    p3 = p[2];
    p4 = p[3];

    d41 = p4 - p1;
    d34 = p3 - p4;
    d31 = p3 - p1;
    d23 = p2 - p3;
    d24 = p2 - p4;
    d21 = p2 - p1;

#define IB(x) (Fabs(x)>eps ? (Exp(x)-1.)/(x): \
	(((1./24 *(x) + 1./6)*(x) + 1./2)*(x) +1.))
#define W(x,y) (1./2 + 1./6*(x+y)+1./24*(x*x+x*y+y*y)+ \
	1./120*(x+y)*(x*x+y*y))
#define T(x,y,z) (1./6 + 1./24*(x+y+z) + 1./120*(x*x+y*y+z*z+x*y+x*z+y*z) + \
	1./720*(x*x*x+y*y*y+z*z*z+(y+z)*x*x+(y*y+z*y+z*z)*x+z*y*y+z*z*y))

    if (Fabs(d34) > eps)
	v0 = 1. / d34 * (Exp(d34) * IB(d23) - IB(d24));
    else
	v0 = IB(d23) * IB(d34) - W(d23, d24);

    if (Fabs(d31) > eps)
	v1 = 1. / d31 * (Exp(d31) * IB(d23) - IB(d21));
    else
	v1 = IB(d23) * IB(d31) - W(d23, d21);

    if (Fabs(d41) > eps) {
	average_coeff[0] = 6. / d41 * (Exp(p4 - s0) * v0 - Exp(p1 - s0) * v1);
	average_coeff[1] = 6. / d41 * (Exp(p4 - s1) * v0 - Exp(p1 - s1) * v1);
    }
    else {
	r = IB(d41) * v0;
	r -= IB(d23) * W(d34, d31);
	r += T(d24, d21, d23);
	average_coeff[0] = 6. * r * Exp(p1 - s0);
	average_coeff[1] = 6. * r * Exp(p1 - s1);
    }

#undef IB
#undef W
#undef T
}

static void
current_Average(DOF *V, DOF *P, DOF *N, CONTACT cont[], int ncont,
	       FLOAT currents[][3])
    /* currents[][0]: total currents,
       currents[][1]: electron currents, 
       currents[][2]: hole currents */
{
    GRID *g = V->g;
    SIMPLEX *e;
    BYTE flag[ncont];
    int i, j, k, c, order = 2;
    FLOAT psi0, psi1, vp[NVert], vn[NVert], local_v[NVert];
    FLOAT coeff_p, coeff_n;
    FLOAT rp, rn;
    FLOAT Bp_ij, Bp_ii, Bn_ij, Bn_ii;
    FLOAT Ap[NVert][NVert], An[NVert][NVert];
    FLOAT Ap_bas, An_bas;
    FLOAT averagep_coeff[2], averagen_coeff[2];


    bzero(currents[0], ncont * 3 * sizeof(currents[0][0]));
    
    ForAllElements(g, e) {
	bzero(flag, ncont * sizeof(flag[0]));
	for (i = 0; i < NVert; i++) {
	    for (j = 0; j < ncont; j++)
		if (g->types_vert[e->verts[i]] & cont[j])
		    flag[j] = 1;
	}
	for (i = 0; i < NVert; i++) {
	    local_v[i] = *(DofVertexData(V, e->verts[i]));
	    vp[i] = *(DofVertexData(P, e->verts[i]));
	    vn[i] = *(DofVertexData(N, e->verts[i]));
	}
	coeff_n = Dn * SCALE;
	coeff_p = Dp * SCALE;
	for (c = 0; c < ncont; c++) {
	    if (flag[c] == 0)
		continue;
	    for (i = 0; i < NVert; i++) {
		if (!(g->types_vert[e->verts[i]] & cont[c]))
		    continue;
      psi0 = local_v[i];
      Ap[i][i] = An[i][i] = 0.0;
   #if 1
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ij = B(psi0 - psi1);
            Bp_ii = B(psi1 - psi0);
            Bn_ij = B(psi1 - psi0);
            Bn_ii = B(psi0 - psi1);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p * Bp_ij * Ap_bas;
            Ap[i][i] -= coeff_p * Bp_ii * Ap_bas;
            An[i][j] = coeff_n * Bn_ij * An_bas;
            An[i][i] -= coeff_n * Bn_ii * An_bas;
         }
      }
   #endif
   #if 0
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ii = Average2(psi0 - psi1);
            Bp_ij = 1.0 - Bp_ii;
            Bn_ii = Average2(psi1 - psi0);
            Bn_ij = 1.0 - Bn_ii;
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p * 2 * Bp_ij * Ap_bas;
            Ap[i][i] -= coeff_p * 2 * Bp_ii * Ap_bas;
            An[i][j] = coeff_n * 2 * Bn_ij * An_bas;
            An[i][i] -= coeff_n * 2 * Bn_ii * An_bas;
         }
      }
   #endif
   #if 0
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ii = Average3_0(psi0 - psi1);
            Bp_ij = Average3_1(psi0 - psi1);
            Bn_ii = Average3_0(psi1 - psi0);
            Bn_ij = Average3_1(psi1 - psi0);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p * Bp_ij * Ap_bas;
            Ap[i][i] -= coeff_p * Bp_ii * Ap_bas;
            An[i][j] = coeff_n * Bn_ij * An_bas;
            An[i][i] -= coeff_n * Bn_ii * An_bas;
         }
      }
   #endif
   #if 0
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
            psi1 = local_v[j];
            Bp_ij = Exp(psi1 - psi0);
            Bp_ii = Exp(psi0 - psi1);
            Bn_ij = Exp(psi0 - psi1);
            Bn_ii = Exp(psi1 - psi0);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p / 2 * (1 + Bp_ij) * Ap_bas;
            Ap[i][i] -= coeff_p / 2 * (1 + Bp_ii) * Ap_bas;
            An[i][j] = coeff_n / 2 * (1 + Bn_ij) * An_bas;
            An[i][i] -= coeff_n / 2 * (1 + Bn_ii) * An_bas;
         }
      }
   #endif
   #if 0
      for(j = 0; j < NVert; j++)
      {
         if(j != i)
         {
	         bzero(averagep_coeff, 2 * sizeof(averagep_coeff[0]));
	         bzero(averagen_coeff, 2 * sizeof(averagen_coeff[0]));
            average(local_v, HOLE, averagep_coeff, i, j);
            average(local_v, ELECTRON, averagen_coeff, i, j);
            Ap_bas = phgQuadGradBasDotGradBas(e, P, j, P, i, 3);
            An_bas = phgQuadGradBasDotGradBas(e, N, j, N, i, 3);
            Ap[i][j] = coeff_p / averagep_coeff[1] * Ap_bas;
            Ap[i][i] -= coeff_p / averagep_coeff[0] * Ap_bas;
            An[i][j] = coeff_n / averagen_coeff[1] * An_bas;
            An[i][i] -= coeff_n / averagen_coeff[0] * An_bas;
         }
      }

   #endif
      rp = rn = 0.0;
      for(j = 0; j < NVert; j++)
      {
         rp += Ap[i][j] * vp[j];
         rn += An[i][j] * vn[j];
      }
		    currents[c][1] += rn;
		    currents[c][2] -= rp;
	    }/*loop on NVert*/
	}/*loop on ncont*/
    }/*loop on element*/

    for (c = 0; c < ncont; c++) {
	/* handle units */
	for (i = 1; i < 3; i++)
	    currents[c][i] /= CM2MICRON;
	currents[c][0] = currents[c][1] + currents[c][2];
    }

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT currents0[ncont][3];
	memcpy(currents0[0], currents[0],
	       ncont * 3 * sizeof(currents0[0][0]));
	MPI_Allreduce(currents0[0], currents[0], ncont * 3, PHG_MPI_FLOAT,
		      MPI_SUM, g->comm);
    }
#endif
}

static void
current_zlamal(DOF *V, DOF *P, DOF *N, CONTACT cont[], int ncont,
	       FLOAT currents[][3])
    /* currents[][0]: total currents,
       currents[][1]: electron currents, 
       currents[][2]: hole currents */
{
    GRID *g = V->g;
    SIMPLEX *e;
    BYTE flag[ncont];
    int i, j, k, n, c, order = 2;
    QUAD *quad = phgQuadGetQuad3D(order);
    int v0, v1;
    const FLOAT *w, *dphi[NVert];
    FLOAT psi[NVert], vp[NVert], vn[NVert];
    FLOAT expdp[NVert], expdn[NVert];
    FLOAT coeff_p, coeff_n;
    FLOAT dp[Dim], dn[Dim], fp[Dim], fn[Dim];
    FLOAT *jbj_p, *jbj_n;
    FLOAT rp, rn;


    bzero(currents[0], ncont * 3 * sizeof(currents[0][0]));
    jbj_p = (FLOAT *)malloc(Dim * Dim * sizeof(*jbj_p));
    jbj_n = (FLOAT *)malloc(Dim * Dim * sizeof(*jbj_n));
    
    ForAllElements(g, e) {
	FLOAT vol = phgGeomGetVolume(g, e);


	bzero(flag, ncont * sizeof(flag[0]));
	for (i = 0; i < NVert; i++) {
	    for (j = 0; j < ncont; j++)
		if (g->types_vert[e->verts[i]] & cont[j])
		    flag[j] = 1;
	}


	for (i = 0; i < NVert; i++) {
	    dphi[i] = phgQuadGetBasisGradient(e, N, i, quad);
	    psi[i] = *(DofVertexData(V, e->verts[i]));
	    vp[i] = *(DofVertexData(P, e->verts[i]));
	    vn[i] = *(DofVertexData(N, e->verts[i]));
	}


	coeff_n = Dn * SCALE;
	coeff_p = Dp * SCALE;


	for (c = 0; c < ncont; c++) {
	    if (flag[c] == 0)
		continue;

	    for (i = 0; i < NVert; i++) {
		if (!(g->types_vert[e->verts[i]] & cont[c]))
		    continue;


		memcpy(jbj_p, JBJ(V, e, i, HOLE), Dim * Dim * sizeof(*jbj_p));
		memcpy(jbj_n, JBJ(V, e, i, ELECTRON),Dim * Dim * sizeof(*jbj_n));
		FLOAT (*Jp)[Dim] = (void *)jbj_p;
		FLOAT (*Jn)[Dim] = (void *)jbj_n;


		for (j = 0; j < NVert; j++) {
		    expdn[j] = Exp(psi[i] - psi[j]);
		    expdp[j] = Exp(psi[j] - psi[i]);
		}


		w = quad->weights;
		for (n = 0; n < quad->npoints; n++) {
		    dp[0] = dp[1] = dp[2] = 0.;
		    dn[0] = dn[1] = dn[2] = 0.;
		    for (j = 0; j < NVert; j++) {
			for (k = 0; k < Dim; k++) {
			    dp[k] += vp[j] * expdp[j] * dphi[j][n * Dim + k];
			    dn[k] += vn[j] * expdn[j] * dphi[j][n * Dim + k];
			}
		    }

		    rp = rn = 0.;
		    for (k = 0; k < Dim; k++) {
			fp[k] =
			    Jp[k][0] * dp[0] + Jp[k][1] * dp[1] +
			    Jp[k][2] * dp[2];
			rp += fp[k] * dphi[i][n * Dim + k];
			fn[k] =
			    Jn[k][0] * dn[0] + Jn[k][1] * dn[1] +
			    Jn[k][2] * dn[2];
			rn += fn[k] * dphi[i][n * Dim + k];
		    }
		    rp *= *w * vol;
		    rn *= *w * vol;

		    currents[c][1] += coeff_n * rn;
		    currents[c][2] -= coeff_p * rp;

		    w++;
		}
	    }
	}
    }

    for (c = 0; c < ncont; c++) {
	/* handle units */
	for (i = 1; i < 3; i++)
	    currents[c][i] /= CM2MICRON;
	currents[c][0] = currents[c][1] + currents[c][2];
    }

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT currents0[ncont][3];
	memcpy(currents0[0], currents[0],
	       ncont * 3 * sizeof(currents0[0][0]));
	MPI_Allreduce(currents0[0], currents[0], ncont * 3, PHG_MPI_FLOAT,
		      MPI_SUM, g->comm);
    }
#endif

    free(jbj_p);
    free(jbj_n);
}
//static void
//func_NOT(SIMPLEX *e,FLOAT *values)
//{
//    swich(e->region_mark){



//	}




//}
static void
build_linear_system_poisson_fem(SOLVER *solver, DOF *V, DOF *P, DOF *N)
{
    int i, j, n, f, order;
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *gradv = phgDofGradient(V, NULL, NULL, NULL);
//    DOF *gradv = phgDofCopy(analytic_grad_u, NULL, NULL, NULL);
    FLOAT A[NVert][NVert], rhs[NVert];
    FLOAT vol, ck, cr, qc, coeff, expp, expn;
    FLOAT *valp, *valn, *valu;
    BTYPE btype;
    QUAD *quad;
    const FLOAT *w, *phi[NVert];
    INT I[NVert];
    DOF *NOT;
    order = 2;			/* order 1 too less, order 3 too more, 2011.5.7 */
    quad = phgQuadGetQuad3D(order);

    valp = phgAlloc(quad->npoints * 3 * sizeof(*valp));
    valn = valp + quad->npoints;
    valu = valp + 2 * quad->npoints;
    
    NOT = phgDofNew(g,DOF_DEFAULT,1,"NOT",DofInterpolation);
     phgDofSetDataByValue(NOT,1.0e12);

    ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(rhs, NVert * sizeof(rhs[0]));


	    for (i = 0; i < NVert; i++)
		phi[i] = phgQuadGetBasisValues(e, V, i, quad);

	    exp_quad_dof_eval(e, P, HOLE, gradv, quad, valp);
	    exp_quad_dof_eval(e, N, ELECTRON, gradv, quad, valn);

	    memcpy(valu, phgQuadGetDofValues(e, V, quad),
		   quad->npoints * sizeof(*valu));
	    vol = phgGeomGetVolume(g, e);

	    qc = 100.0 * Q / (EPSILONE0 * KB * T) /* cm */ ;
	    qc /= CM2MICRON * CM2MICRON * CM2MICRON;	/* handle units */

	    w = quad->weights;
	    for (n = 0; n < quad->npoints; n++) {
		ck = qc * (valp[n] + valn[n]) * SCALE;
		cr = qc * (valp[n] - valn[n]) * SCALE;

		for (i = 0; i < NVert; i++) {
		    rhs[i] += cr * phi[i][n] * (*w);
		    for (j = 0; j < NVert; j++)
			A[i][j] += ck * phi[i][n] * phi[j][n] * (*w);
		}
		w++;
	    }

	    for (i = 0; i < NVert; i++) {
		rhs[i] *= vol;
		for (j = 0; j < NVert; j++)
		    A[i][j] *= vol;
	    }


	coeff =  EPSILONE;
	coeff /= CM2MICRON;	/* handle units */
	for (i = 0; i < NVert; i++) {
	    rhs[i] -=
		coeff * phgQuadDofDotGradBas(e, gradv, V, i, QUAD_DEFAULT);
	    for (j = 0; j < NVert; j++)
		A[i][j] += coeff * phgQuadGradBasDotGradBas
		    (e, V, i, V, j, QUAD_DEFAULT);
	}

	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    btype = phgDofGetElementBoundaryType(V, e, i);
	    if (btype & ohmic_contact) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
      if(analytic_test)
    {
	    if (btype & JUNCTION_PN) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
   }
	}
//for(i = 0; i < NVert; i++)
  // phgPrintf("\n%.8e\n", rhs[i]);
	/* add entries */
	for (i = 0; i < NVert; i++)
	I[i] = phgSolverMapE2L(solver, 0, e, i);
	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
    }

    phgFree(valp);
    phgDofFree(&gradv);
    phgDofFree(&NOT);
}

static void
build_linear_system_poisson_fem_new(SOLVER *solver, DOF *V, DOF *Log_P, DOF *Log_N)
{
    int i, j, n, f, order;
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *gradv = phgDofGradient(V, NULL, NULL, NULL);
    FLOAT A[NVert][NVert], rhs[NVert];
    FLOAT vol, ck, cr, qc, coeff, expp, expn;
    FLOAT *valp, *valn, *valu;
    BTYPE btype;
    QUAD *quad;
    const FLOAT *w, *phi[NVert];
    INT I[NVert];
   // DOF *NOT;
       DOF *E;
       BTYPE interface;
    if(surface)
    {
       E = phgDofNew(g, DOF_CONSTANT, 1, "unit", DofNoAction);
       phgDofSetDataByValue(E, 1.0);
       if(device == PN)
          interface = JUNCTION_PN;
       if(device == NANO)
          interface = ANODE;
    }
    order = 2;			/* order 1 too less, order 3 too more, 2011.5.7 */
    quad = phgQuadGetQuad3D(order);

    valp = phgAlloc(quad->npoints * 3 * sizeof(*valp));
    valn = valp + quad->npoints;
    valu = valp + 2 * quad->npoints;
    
   // NOT = phgDofNew(g,DOF_DEFAULT,1,"NOT",DofInterpolation);
   //  phgDofSetDataByValue(NOT,1.0e12);

    ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(rhs, NVert * sizeof(rhs[0]));


	    for (i = 0; i < NVert; i++)
		phi[i] = phgQuadGetBasisValues(e, V, i, quad);

//	    exp_quad_dof_eval(e, P, HOLE, gradv, quad, valp);
//	    exp_quad_dof_eval(e, N, ELECTRON, gradv, quad, valn);

	    memcpy(valp, phgQuadGetDofValues(e, Log_P, quad),
		   quad->npoints * sizeof(*valp));
	    memcpy(valn, phgQuadGetDofValues(e, Log_N, quad),
		   quad->npoints * sizeof(*valn));
	    memcpy(valu, phgQuadGetDofValues(e, V, quad),
		   quad->npoints * sizeof(*valu));
	    vol = phgGeomGetVolume(g, e);

	    qc = 100.0 * Q / (EPSILONE0 * KB * T) /* cm */ ;
	    qc /= CM2MICRON * CM2MICRON * CM2MICRON;	/* handle units */

	    w = quad->weights;
	    for (n = 0; n < quad->npoints; n++) {
		ck = qc * nie * (Exp(valp[n]) + Exp(valn[n]));
		cr = qc * nie * (Exp(valp[n]) - Exp(valn[n]));

		for (i = 0; i < NVert; i++) {
		    rhs[i] += cr * phi[i][n] * (*w);
		    for (j = 0; j < NVert; j++)
			A[i][j] += ck * phi[i][n] * phi[j][n] * (*w);
		}
		w++;
	    }

	    for (i = 0; i < NVert; i++) {
		rhs[i] *= vol;
		for (j = 0; j < NVert; j++)
		    A[i][j] *= vol;
	    }


	coeff =  EPSILONE;
	coeff /= CM2MICRON;	/* handle units */
	for (i = 0; i < NVert; i++) {
	    rhs[i] -=
		coeff * phgQuadDofDotGradBas(e, gradv, V, i, QUAD_DEFAULT);
	    for (j = 0; j < NVert; j++)
		A[i][j] += coeff * phgQuadGradBasDotGradBas
		    (e, V, i, V, j, QUAD_DEFAULT);
	}

	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    btype = phgDofGetElementBoundaryType(V, e, i);
	    if (btype & ohmic_contact) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
      if(analytic_test)
    {
	    if (btype & JUNCTION_PN) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
   }
   if(surface)
   {
       for(j = 0; j < NFace; j++)
       {
          if(e->bound_type[j] & interface)
          {
            rhs[i] += 100 / (KB * T * EPSILONE0) *surface_charge /Pow(CM2MICRON, 2) * phgQuadFaceDofDotBas(e, j, E, DOF_PROJ_NONE, V, i, QUAD_DEFAULT); 
          }
      }
   }
	}
//for(i = 0; i < NVert; i++)
  // phgPrintf("\n%.8e\n", rhs[i]);
	/* add entries */
	for (i = 0; i < NVert; i++)
	I[i] = phgSolverMapE2L(solver, 0, e, i);
	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
    }

    phgFree(valp);
    phgDofFree(&gradv);
    if(surface)
    phgDofFree(&E);
    //phgDofFree(&NOT);
}

static void
build_linear_system_continuity_Average(SOLVER *solver, DOF *V, DOF *P, DOF *N, EQUATION eqn)
{
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *dof = eqn == ELECTRON ? N : P;
    FLOAT A[NVert][NVert], rhs[NVert], A_bas;
    FLOAT u[NVert];
    INT I[NVert];
    int i, j, k, n;
    int order = 2;
    QUAD *quad = phgQuadGetQuad3D(order);
    FLOAT *fr;
    const FLOAT* bas_val[NVert];
    FLOAT vol, r, coeff, psi0, psi1;
    const FLOAT *w, *dphi[NVert];
    FLOAT B_ij, B_ii;
    FLOAT local_v[NVert];
    FLOAT average_coeff[2];
    FLOAT alph_upw_p = 0.0, alph_upw_n = 0.0;/* control upwind*/

   ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(local_v, NVert * sizeof(local_v[0]));
	bzero(rhs, NVert * sizeof(rhs[0]));

	vol = phgGeomGetVolume(g, e);
      for (i = 0; i < NVert; i++)
      {
         u[i] = *(DofVertexData(dof, e->verts[i]));
	      bas_val[i] = phgQuadGetBasisValues(e, dof, i, quad);
         local_v[i] = *(DofVertexData(V, e->verts[i]));     
      }

   if(analytic_test)
   {
      if(eqn == ELECTRON)
         fr = phgQuadGetDofValues(e, analytic_f_n, quad);
      else
         fr = phgQuadGetDofValues(e, analytic_f_p, quad);
   }


	coeff = eqn == ELECTRON ? Dn : Dp;
	coeff /= CM2MICRON;	/* handle units */

#if 1
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
            //   psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ij = eqn == HOLE ? B(psi0 - psi1) : B(psi1 - psi0);
               B_ii = eqn == HOLE ? B(psi1 - psi0) : B(psi0 - psi1);
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff * B_ij * A_bas;
               A[i][i] -= coeff * B_ii * A_bas;
            }
         }
      }
#endif
#if 0
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
               //psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ii = eqn == HOLE ? Average2(psi0 - psi1) : Average2(psi1 - psi0);
               B_ij = 1.0 - B_ii;
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff * 2 * B_ij * A_bas;
               A[i][i] -= coeff * 2 * B_ii * A_bas;
            }
         }
      }
#endif
#if 0
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
               //psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ii = eqn == HOLE ? Average3_0(psi0 - psi1) : Average3_0(psi1 - psi0);
               B_ij = eqn == HOLE ? Average3_1(psi0 - psi1) : Average3_1(psi1 - psi0);
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff *  B_ij * A_bas;
               A[i][i] -= coeff *  B_ii * A_bas;
            }
         }
      }

#endif
#if 0
      for(i = 0; i < NVert; i++)
      {
         //psi0 = *(DofVertexData(V, e->verts[i]));
         psi0 = local_v[i];
         A[i][i] = 0.0;
         for(j = 0; j < NVert; j++)
         {
            if(j != i)
            {
              // psi1 = *(DofVertexData(V, e->verts[j]));
               psi1 = local_v[j];
               B_ij = eqn == HOLE ? Exp(psi1 - psi0) : Exp(psi0 - psi1);
               B_ii = eqn == HOLE ? Exp(psi0 - psi1) : Exp(psi1 - psi0);
               A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
               A[i][j] =  coeff / 2 * (1 + B_ij) * A_bas;
               A[i][i] -= coeff / 2 * (1 + B_ii) * A_bas;
            }
         }
      }

#endif
#if 0
      for(i = 0; i < NVert; i++)
      {
         A[i][i] = 0.0;
            for(j = 0; j < NVert; j++)
            {
               if(j != i)
               {
	               bzero(average_coeff, 2 * sizeof(average_coeff[0]));
                  average(local_v, eqn, average_coeff, i, j);
                  A_bas = phgQuadGradBasDotGradBas(e, dof, j, dof, i, 3);
                  A[i][j] = coeff / average_coeff[1] * A_bas;
                  A[i][i] -= coeff / average_coeff[0] * A_bas;
               }
            }
      }
#endif

      for(i = 0; i < NVert; i++)
      {
         r = 0.0;
         for(j = 0; j < NVert; j++)
            r += A[i][j] * u[j]; 
         rhs[i] = -1.0 * r;
         if(analytic_test)
         {
	         w = quad->weights;
            for(n = 0; n < quad->npoints; n++)
            {
               rhs[i] += coeff *  fr[n] * bas_val[i][n] * (*w) * vol; 
		         w++;
            }
         }
      }
	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    BTYPE btype = phgDofGetElementBoundaryType(dof, e, i);
	    if (btype & ohmic_contact) {
		   for (j = 0; j < NVert; j++)
		       A[i][j] =  A[j][i] = 0.0;
		   rhs[i] = 0.0;
		   A[i][i] = 1.0;
	    }
       if(analytic_test) 
       {
          if(btype & JUNCTION_PN)
         {
		      for (j = 0; j < NVert; j++)
		         A[i][j] = A[j][i] = 0.0;
		   rhs[i] = 0.0;
		   A[i][i] = 1.0;
         }
       }
	}
	/* add entries */
	for (i = 0; i < NVert; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
   }
}

static void
build_linear_system_continuity_zlamal(SOLVER *solver, DOF *V, DOF *P, DOF *N,
				      EQUATION eqn)
{
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *dof = eqn == ELECTRON ? N : P;
    DOF *gradv = phgDofGradient(V, NULL, NULL, NULL);
    FLOAT A[NVert][NVert], rhs[NVert], expd[NVert];
    FLOAT u[NVert], psi[NVert], du[Dim], f[Dim], du2[Dim];
    INT I[NVert];
    int i, j, k, n;
    int order = 1;
    QUAD *quad = phgQuadGetQuad3D(order);
    FLOAT *fr;
    const FLOAT* bas_val[NVert];
    FLOAT ND, NA, vol, r, coeff, psi0, psi1;
    const FLOAT *w, *dphi[NVert];

    ForAllElements(g, e) {
	bzero(A[0], NVert * NVert * sizeof(A[0][0]));
	bzero(rhs, NVert * sizeof(rhs[0]));


	vol = phgGeomGetVolume(g, e);
	for (i = 0; i < NVert; i++) {
	    dphi[i] = phgQuadGetBasisGradient(e, dof, i, quad);
	    u[i] = *(DofVertexData(dof, e->verts[i]));
	    psi[i] = *(DofVertexData(V, e->verts[i]));
	    bas_val[i] = phgQuadGetBasisValues(e, dof, i, quad);
	}

   if(analytic_test)
   {
      if(eqn == ELECTRON)
      {
         fr = phgQuadGetDofValues(e, analytic_f_n, quad);
      }
      else
      {
         fr = phgQuadGetDofValues(e, analytic_f_p, quad);
      }
   }

	coeff = eqn == ELECTRON ? Dn : Dp;
	coeff /= CM2MICRON;	/* handle units */

	for (i = 0; i < NVert; i++) {
	    FLOAT (*J)[Dim] = (void *)JBJ(V, e, i, eqn);

	    for (j = 0; j < NVert; j++) {
		expd[j] = eqn == ELECTRON ?
		    Exp(psi[i] - psi[j]) : Exp(psi[j] - psi[i]);
	    }

	    w = quad->weights;
	    for (n = 0; n < quad->npoints; n++) {
		du[0] = du[1] = du[2] = 0.;
//		du2[0] = du2[1] = du2[2] = 0.;
		for (j = 0; j < NVert; j++) {
		    for (k = 0; k < Dim; k++) {
			du[k] += u[j] * expd[j] * dphi[j][n * Dim + k];
//			du[k] += u[j] * expd[k+1] * dphi[j][n * Dim + k];
	//		du[k] += u[j] *  dphi[j][n * Dim + k];
  /*       if(eqn == ELECTRON)
         {
            du2[k] -= u[i] * psi[j] * dphi[j][n * Dim + k];
         }
         else if(eqn == HOLE)
         {
            du2[k] += u[i] * psi[j] * dphi[j][n * Dim + k];
         }*/
		    }
		}

		r = 0.;
		for (k = 0; k < Dim; k++) {
		    f[k] =
			J[k][0] * du[0] + J[k][1] * du[1] + J[k][2] * du[2];
		    r += f[k] * dphi[i][n * Dim + k];
       //   r += du2[k] * dphi[i][n * Dim + k];
		}
		r *= *w * vol;
		rhs[i] -= coeff * r;	/* residual */

		for (j = 0; j < NVert; j++) {
		    for (k = 0; k < Dim; k++) {
			f[k] =
			    J[k][0] * dphi[j][n * Dim + 0] +
			    J[k][1] * dphi[j][n * Dim + 1] +
			    J[k][2] * dphi[j][n * Dim + 2];
/*		f[k] =
			    J[k][0] * expd[k+1] * dphi[j][n * Dim + 0] +
			    J[k][1] * expd[k+1] * dphi[j][n * Dim + 1] +
			    J[k][2] * expd[k+1] * dphi[j][n * Dim + 2];
         if(eqn == ELECTRON)
         f[k] -= u[i] * dphi[j][n * Dim + k];
         else if(eqn == HOLE)
         f[k] += u[i] * dphi[j][n * Dim + k];*/
		    }

		    r = f[0] * dphi[i][n * Dim + 0] +
			f[1] * dphi[i][n * Dim + 1] +
			f[2] * dphi[i][n * Dim + 2];
		    r *= *w * vol;
		   A[i][j] += coeff * expd[j] * r;
		 //   A[i][j] += coeff  * r;
		}

      if(analytic_test)
      {
         rhs[i] += coeff *  fr[n] * bas_val[i][n] * (*w) * vol; 
      }

		w++;
	    }
	}

	/* handle bdry */
	for (i = 0; i < NVert; i++) {
	    BTYPE btype = phgDofGetElementBoundaryType(V, e, i);
	    if (btype & ohmic_contact) {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
	    }
       if(analytic_test) 
       {
       if(btype & JUNCTION_PN)
       {
		for (j = 0; j < NVert; j++)
		    A[i][j] = A[j][i] = 0.0;
		rhs[i] = 0.0;
		A[i][i] = 1.0;
       }
       }

	}

      add_entries:
	/* add entries */
	for (i = 0; i < NVert; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);

	phgSolverAddMatrixEntries(solver, NVert, I, NVert, I, A[0]);
	phgSolverAddRHSEntries(solver, NVert, I, rhs);
    }

    phgDofFree(&gradv);
}

RelativeError(DOF *V, DOF *P, DOF *N, BOOLEAN flag)
{
    GRID *g = V->g;
    SIMPLEX *e;
    static DOF *V_old, *P_old, *N_old;
    static FLOAT err[3];
    INT i;

    if (!flag) {
	V_old = phgDofCopy(V, NULL, NULL, NULL);
	P_old = phgDofCopy(P, NULL, NULL, NULL);
	N_old = phgDofCopy(N, NULL, NULL, NULL);
	return NULL;
    }
    else {
	FLOAT v_old, v_new, p_old, p_new, n_old, n_new, r, rv, rp, rn;
	FLOAT C;

#if 0
	C = 1.e14 / SCALE;
#else
	C = 1.;
#endif

#define max(a,b) ((a)>(b)? (a):(b))

	rv = rp = rn = 0.;
	for (i = 0; i < g->nvert; i++) {
	    if (g->types_vert[i] == UNREFERENCED)
		continue;

	    v_old = *(V_old->data + i);
	    v_new = *(V->data + i);
	    r = Fabs(v_new - v_old) / max(Fabs(v_old), 1.);
	    if (rv < r)
		rv = r;

	    if (!(verts_region[i] & MAT_SEMICONDUCTOR))
		continue;

	    p_old = *(P_old->data + i);
	    p_new = *(P->data + i);
	    n_old = *(N_old->data + i);
	    n_new = *(N->data + i);

	    if (p_old < 1. / SCALE)
		p_old = 1. / SCALE;
	    if (p_new < 1. / SCALE)
		p_new = 1. / SCALE;
	    if (n_old < 1. / SCALE)
		n_old = 1. / SCALE;
	    if (n_new < 1. / SCALE)
		n_new = 1. / SCALE;

	    /* concentration to fermi level */
       p_old = v_old + Log(p_old * SCALE / nie);
	    n_old = v_old - Log(n_old * SCALE / nie);
	    p_new = v_new + Log(p_new * SCALE / nie);
	    n_new = v_new - Log(n_new * SCALE / nie);

	    r = Fabs(p_new - p_old) / max(Fabs(p_old), C);
	    if (rp < r)
		rp = r;
	    r = Fabs(n_new - n_old) / max(Fabs(n_old), C);
	    if (rn < r)
		rn = r;
	}
#undef max

	err[0] = rv;
	err[1] = rp;
	err[2] = rn;

#if USE_MPI
	if (g->nprocs > 1) {
	    FLOAT err0[3];
	    memcpy(err0, err, 3 * sizeof(err0[0]));
	    MPI_Allreduce(err0, err, 3, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	}
#endif

	phgDofFree(&V_old);
	phgDofFree(&P_old);
	phgDofFree(&N_old);

	return err;
    }
}

RelativeError_new(DOF *V, DOF *P, DOF *N, BOOLEAN flag)
{
    GRID *g = V->g;
    SIMPLEX *e;
    static DOF *V_old, *P_old, *N_old;
    static FLOAT err[3];
    INT i;

    if (!flag) {
	V_old = phgDofCopy(V, NULL, NULL, NULL);
	P_old = phgDofCopy(P, NULL, NULL, NULL);
	N_old = phgDofCopy(N, NULL, NULL, NULL);
	return NULL;
    }
    else {
	FLOAT v_old, v_new, p_old, p_new, n_old, n_new, r, rv, rp, rn;
	FLOAT C;

#if 0
	C = 1.e14 / SCALE;
#else
	C = 1.;
#endif

#define max(a,b) ((a)>(b)? (a):(b))

	rv = rp = rn = 0.;
/*	for (i = 0; i < g->nvert; i++) {
	    if (g->types_vert[i] == UNREFERENCED)
		continue;

	    v_old = *(V_old->data + i);
	    v_new = *(V->data + i);
	    r = Fabs(v_new - v_old) / max(Fabs(v_old), 1.);
	    if (rv < r)
		rv = r;

	    if (!(verts_region[i] & MAT_SEMICONDUCTOR))
		continue;

	    p_old = *(P_old->data + i);
	    p_new = *(P->data + i);
	    n_old = *(N_old->data + i);
	    n_new = *(N->data + i);

	    if (p_old < 1. / SCALE)
		p_old = 1. / SCALE;
	    if (p_new < 1. / SCALE)
		p_new = 1. / SCALE;
	    if (n_old < 1. / SCALE)
		n_old = 1. / SCALE;
	    if (n_new < 1. / SCALE)
		n_new = 1. / SCALE;
*/
	    /* concentration to fermi level */
  /*     p_old = v_old + Log(p_old * SCALE / nie);
	    n_old = v_old - Log(n_old * SCALE / nie);
	    p_new = v_new + Log(p_new * SCALE / nie);
	    n_new = v_new - Log(n_new * SCALE / nie);

	    r = Fabs(p_new - p_old) / max(Fabs(p_old), C);
	    if (rp < r)
		rp = r;
	    r = Fabs(n_new - n_old) / max(Fabs(n_old), C);
	    if (rn < r)
		rn = r;
      }*/
   phgDofAXPBY(-1.0, V, 1.0, &V_old);
   phgDofAXPBY(-1.0, P, 1.0, &P_old);
   phgDofAXPBY(-1.0, N, 1.0, &N_old);
   r = max(phgDofNormL2(V), C);
   rv = phgDofNormL2(V_old)/r;
   r = max(phgDofNormL2(P), C);
   rp = phgDofNormL2(P_old)/r;
   r = max(phgDofNormL2(N), C);
   rn = phgDofNormL2(N_old)/r;
#undef max

	err[0] = rv;
	err[1] = rp;
	err[2] = rn;

#if USE_MPI
	if (g->nprocs > 1) {
	    FLOAT err0[3];
	    memcpy(err0, err, 3 * sizeof(err0[0]));
	    MPI_Allreduce(err0, err, 3, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	}
#endif

	phgDofFree(&V_old);
	phgDofFree(&P_old);
	phgDofFree(&N_old);

	return err;
    }
}
static double
elapsed_time(GRID *g, BOOLEAN flag, SOLVER *solver)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (solver != NULL) {
	    phgPrintf("  Linear systems solved: [nits %d, res %0.2le ",
		      solver->nits, (double)solver->residual);
	}
	phgPrintf("%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
    }

    return et;
}

static void
Update(DOF *V, DOF *P, DOF *N, DOF *delta_v, FLOAT damp)
/* to keep the Quasi-Fermi potential constant when solving Poisson's equation */
{
    INT i;
   GRID *g = V->g;
   SIMPLEX *e;
    FLOAT dv;

    phgDofAXPBY(damp, delta_v, 1.0, &V);

    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED ||
	    !(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;

	dv = *(DofVertexData(delta_v, i)) * damp;
	*(P->data + i) *= Exp(-dv);
	*(N->data + i) *= Exp(dv);
    }
}

static void
Update_new(DOF *V, DOF *Log_P, DOF *Log_N, DOF *delta_v, FLOAT damp)
/* to keep the Quasi-Fermi potential constant when solving Poisson's equation */
{
    phgDofAXPBY(damp, delta_v, 1.0, &V);
    phgDofAXPBY(-1.0 * damp, delta_v, 1.0, &Log_P);
    phgDofAXPBY(damp, delta_v, 1.0, &Log_N);

}

static void
Current(DOF *V, DOF *P, DOF *N, CONTACT cont[], int ncont,
	FLOAT currents[][3])
{
    if (!strcmp(continuity_names[continuity_index], "zlamal"))
	current_zlamal(V, P, N, cont, ncont, currents);
    else if (!strcmp(continuity_names[continuity_index], "Average"))
	current_Average(V, P, N, cont, ncont, currents);
    else
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);
}

static void
Section_values(DOF *V, DOF *P, DOF *N)
{
  int i, j, Num = 64;
  FLOAT interval = 0.25, pos;
  FLOAT section_values[3];
  for(i = 1; i < Num; i++)
  {
     pos = i * interval;
     Calculate_section_values(V, P, N, pos, section_values);
     if(phgRank == 1 && fn_dofs != NULL)
     {
         FILE *f_in;
         f_in = fopen(fn_dofs, "a");
         if(f_in == NULL)
            phgError(1, "Can not open dofs output file %s.\n", fn_dofs);
         fprintf(f_in, "%lf  %lf  %lf  %lf\n", pos, section_values[0], section_values[1], section_values[2]);
         fclose(f_in);
     }
  }
}

static void
SolvePoisson(DOF *V, DOF *P, DOF *N, FLOAT damp)
{
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *delta;
    SOLVER *solver;
    int iter_newton;
    int maxit_newton;
    FLOAT r, rtol_v = 1.0e-4;
    if (damp < 0.5)
	maxit_newton = 1000;	/* it can take many iterations to get to the 
				   equilibrium state for MOS */
    else
	maxit_newton = 100;

    delta = phgDofNew(g, V->type, V->dim, "delta", DofNoAction);

    phgOptionsPush();
#if 0
    //poisson
    phgOptionsSetOptions("-solver  mumps -mumps_symmetry spd");
#elif 1				/* iterative method */
    phgOptionsSetOptions
	("-solver hypre -hypre_solver pcg -hypre_pc boomeramg");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-15 -solver_atol 1.0e-15 -solver_btol 1.0e-15");
    phgOptionsSetOptions("-solver_maxit 200");
#else /* __float128 */
   phgOptionsSetOptions
	("-solver pcg -pcg_pc_type solver -pcg_pc_opts \"-solver mumps -mumps_symmetry spd\"");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-30 -solver_atol 1.0e-30 -solver_btol 1.0e-30");
    phgOptionsSetOptions("-solver_maxit 5");
#endif

    for (iter_newton = 1; iter_newton <= maxit_newton; iter_newton++) {
	solver = phgSolverCreate(SOLVER_DEFAULT, delta, NULL);
	solver->mat->handle_bdry_eqns = FALSE;

	if (!strcmp(poisson_names[poisson_index], "fem"))
	    build_linear_system_poisson_fem(solver, V, P, N);
	else
	    phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

	phgDofSetDataByValue(delta, 0.0);
	elapsed_time(g, FALSE, NULL);
	phgSolverSolve(solver, TRUE, delta, NULL);
	elapsed_time(g, TRUE, solver);
	phgSolverDestroy(&solver);
	Update(V, P, N, delta, damp);

//	r = phgDofNormInftyVec(V);
	r = phgDofNormL2(V);
   r = r > 1.0 ? r : 1.0;
	//r = phgDofNormInftyVec(delta) / r;
	r = phgDofNormL2(delta)/r;
	phgPrintf("  rtol %0.6le\n", (double)r);
	if (r < rtol_v) {
      phgPrintf("\n||V|| = %.8e, ||P|| = %.8e, ||N|| = %.8e\n", phgDofNormL2(V), phgDofNormL2(P), phgDofNormL2(N));
	    phgPrintf("  Converged.\n");
	    break;
	}
    }

    if (iter_newton > maxit_newton)
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

    phgOptionsPop();
    phgDofFree(&delta);
}

static void
SolvePoisson_new(DOF *V, DOF *P, DOF *N, FLOAT damp)
{
    GRID *g = V->g;
    SIMPLEX *e;
    DOF *delta;
    SOLVER *solver;
    int iter_newton;
    int maxit_newton;
    FLOAT r, rtol_v = 1.0e-5;

    DOF *Log_P, *Log_N;
    Log_P = phgDofCopy(P, NULL, NULL, NULL);
    Log_N = phgDofCopy(N, NULL, NULL, NULL);

    INT i;
    for(i = 0; i < g->nvert; i++)
    {
       if(g->types_vert[i] == UNREFERENCED || !(verts_region[i] & MAT_SEMICONDUCTOR))
          continue;
       if(*(Log_P->data + i) > 0)
          *(Log_P->data + i) = Log(*(Log_P->data + i) / nie * SCALE);
       else
         *(Log_P->data + i) = 0.0;
       if(*(Log_N->data + i) > 0)
         *(Log_N->data + i) = Log(*(Log_N->data + i) / nie * SCALE);
       else
         *(Log_N->data + i) = 0.0;
    }

    if (damp < 0.5)
	maxit_newton = 1000000;	/* it can take many iterations to get to the 
				   equilibrium state for MOS */
    else
	maxit_newton = 1000000;

    delta = phgDofNew(g, V->type, V->dim, "delta", DofNoAction);

    phgOptionsPush();
#if 0
    //poisson
    phgOptionsSetOptions("-solver  mumps -mumps_symmetry spd");
#elif 1				/* iterative method */
    phgOptionsSetOptions
	("-solver hypre -hypre_solver pcg -hypre_pc boomeramg");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-15 -solver_atol 1.0e-15 -solver_btol 1.0e-15");
    phgOptionsSetOptions("-solver_maxit 1000");
#else /* __float128 */
   phgOptionsSetOptions
	("-solver pcg -pcg_pc_type solver -pcg_pc_opts \"-solver mumps -mumps_symmetry spd\"");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-30 -solver_atol 1.0e-30 -solver_btol 1.0e-30");
    phgOptionsSetOptions("-solver_maxit 5");
#endif

    for (iter_newton = 1; iter_newton <= maxit_newton; iter_newton++) {
	solver = phgSolverCreate(SOLVER_DEFAULT, delta, NULL);
	solver->mat->handle_bdry_eqns = FALSE;

	if (!strcmp(poisson_names[poisson_index], "fem"))
	    build_linear_system_poisson_fem_new(solver, V, Log_P, Log_N);
	else
	    phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

	phgDofSetDataByValue(delta, 0.0);
//	elapsed_time(g, FALSE, NULL);
	phgSolverSolve(solver, TRUE, delta, NULL);
//	elapsed_time(g, TRUE, solver);
	phgSolverDestroy(&solver);
	Update_new(V, Log_P, Log_N, delta, damp);

//	r = phgDofNormInftyVec(V);
	r = phgDofNormL2(V);
//   phgPrintf("||delta_V_mid|| = %.8e ||V_mid|| = %.8e\n", phgDofNormL2(delta), r);
   r = r > 1.0 ? r : 1.0;
	//r = phgDofNormInftyVec(delta) / r;
	r = phgDofNormL2(delta)/r;
	phgPrintf("  rtol %0.6le\n", (double)r);
	if (r < rtol_v) {
      phgPrintf("\n||V|| = %.8e, ||P|| = %.8e, ||N|| = %.8e\n", phgDofNormL2(V), phgDofNormL2(P), phgDofNormL2(N));
	    phgPrintf("  Converged.\n");
	    break;
	}
    }

    if (iter_newton > maxit_newton)
	phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);

    for(i = 0; i < g->nvert; i++)
    {
       if(g->types_vert[i] == UNREFERENCED || !(verts_region[i] & MAT_SEMICONDUCTOR))
          continue;
         *(P->data + i) = Exp(*(Log_P->data + i)) * nie / SCALE;
         *(N->data + i) = Exp(*(Log_N->data + i)) * nie / SCALE;
    }

    phgOptionsPop();
    phgDofFree(&delta);
    phgDofFree(&Log_P);
    phgDofFree(&Log_N);
}

static void
SolveContinuity(DOF *V, DOF *P, DOF *N, EQUATION eqn)
{
    GRID *g = V->g;
    SIMPLEX *e;
    SOLVER *solver;
    DOF *u;
    int iter_newton;
    int maxit_newton = 1;
    FLOAT r, rtol_v = 1.0e-6;

    assert(eqn == ELECTRON || eqn == HOLE);
    u = eqn == ELECTRON ? N : P;

    phgOptionsPush();
//#if 0

//    phgOptionsSetOptions("-solver mumps -mumps_symmetry unsym");
//#elif 1				/* iterative method */
  //  phgOptionsSetOptions("-solver petsc");
  //  phgOptionsSetOptions("-oem_options \"-ksp_type gmres\"");
 phgOptionsSetOptions
	("-petsc_pc_opts \"-solver hypre -hypre_solver boomeramg \
	    -hypre_amg_coarsen_type falgout -hypre_pc none\"");
    phgOptionsSetOptions("-solver_maxit 20000");
    phgOptionsSetOptions
	("-solver_rtol 1.0e-15 -solver_atol 1.0e-15 -solver_btol 1.0e-15");
//#elif 0
///    phgOptionsSetOptions("-solver petsc");
//   phgOptionsSetOptions("-oem_options \"-ksp_type gmres\"");
//   phgOptionsSetOptions("-oem_options \"-pc_asm_overlap 1\"");
//   phgOptionsSetOptions
//	("-oem_options \"-sub_pc_factor_mat_solver_package mumps\"");
//    phgOptionsSetOptions("-oem_options \"-mat_mumps_sym 0\"");
//    phgOptionsSetOptions("-oem_options \"-mat_mumps_icntl_4 0\"");
//   phgOptionsSetOptions("-solver_maxit 1000");
//   phgOptionsSetOptions
//	("-solver_rtol 1.0e-12 -solver_atol 1.0e-12 -solver_btol 1.0e-12");
//#else /* __float128 */
//   phgOptionsSetOptions
//	("-solver gmres -gmres_pc_type solver -gmres_pc_opts \"-solver mumps -mumps_symmetry unsym\"");
//    phgOptionsSetOptions
//	("-solver_rtol 1.0e-30 -solver_atol 1.0e-30 -solver_btol 1.0e-30");
//    phgOptionsSetOptions("-solver_maxit 5");
//#endif
    DOF *delta = phgDofNew(g, u->type, u->dim, "increment", DofNoAction);

    for (iter_newton = 1; iter_newton <= maxit_newton; iter_newton++) {
	solver = phgSolverCreate(SOLVER_DEFAULT, delta, NULL);
	solver->mat->handle_bdry_eqns = FALSE;
	if (!strcmp(continuity_names[continuity_index], "zlamal"))
	    build_linear_system_continuity_zlamal(solver, V, P, N, eqn);
   else if (!strcmp(continuity_names[continuity_index], "Average"))
	    build_linear_system_continuity_Average(solver, V, P, N, eqn);
	else
	    phgError(1, "%s:%d,unexpected error.\n", __FILE__, __LINE__);


	phgDofSetDataByValue(delta, 0.0);
//	elapsed_time(g, FALSE, NULL);
	phgSolverSolve(solver, TRUE, delta, NULL);
        phgPrintf("$$$$$$$");
//	elapsed_time(g, TRUE, solver);
	phgSolverDestroy(&solver);
	phgDofAXPBY(1.0,  delta, 1.0, &u);
   if(eqn == HOLE)
   phgPrintf("\n||P||_2 = %.8e ||delta_P|| = %.8e\n", phgDofNormL2(u), phgDofNormL2(delta));
    else
   phgPrintf("\n||N||_2 = %.8e ||delta_N|| = %.8e\n", phgDofNormL2(u), phgDofNormL2(delta));
    }

    phgDofFree(&delta);
    phgOptionsPop();
}

#undef r2
#undef B

static int
bc_map(int bctype)
{
    if (device == PN || device == NANO) {
	switch (bctype) {
	    case 1:
		return CATHODE;
	    case 2:
		return ANODE;
	    case 4:
		return JUNCTION_PN;
	    default:
		return UNDEFINED;
	}
    }

    return UNDEFINED;
}

static void Surface_phi(DOF *dof)
{
   GRID *g = dof->g;
   SIMPLEX *e;
   FLOAT area = 0.0, sum_area = 0.0, total_area = 0.0;
   const FLOAT *lambda = NULL;
   COORD center_vert;
   BTYPE interface;
   if(device == PN)
      interface = JUNCTION_PN;
   if(device == NANO)
      interface = ANODE;

   int face, vert, d;

   FLOAT e_value, dof_value, sum_value;

   dof_value = 0.0;
   ForAllElements(g, e)
   {
      for(face = 0; face < NFace; face++)
      {
          if(e->bound_type[face] & interface)
          {
               area = phgGeomGetFaceArea(g, e, face);
               for(d = 0; d < Dim; d++)
               {
                  center_vert[d] = 0.0;
               }
               for(vert = 0; vert < NVert; vert++)
               {
                  if(vert != face)
                  {
                     for(d = 0; d < Dim; d++)
                     {
                        center_vert[d] += g->verts[e->verts[vert]][d] / 3.0;
                     }
                  }
               }
               lambda = phgGeomXYZ2Lambda(g, e, center_vert[0], center_vert[1], center_vert[2]);
               phgDofEval(dof, e, lambda, &e_value);
               dof_value += area * e_value;
               sum_area += area;
          }
      }
   }
      
#if USE_MPI
   if(g->nprocs > 1)
   {
      MPI_Allreduce(&dof_value, &sum_value, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
      MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
   }
#endif

   sum_value /= total_area;
   sum_value *= KB * T;
   phgPrintf("analytic_density = %.2eM  phi_surface = %.8e\n", analytic_density, sum_value);
   //phgPrintf("total_area = %.8e\n",total_area);
   FLOAT y = 8 * EPSILONE0 * EPSILONE * KB * T * Q / 100 * analytic_density * NA / Pow(DM2CM, 3); 
   FLOAT x = surface_charge / Sqrt(y);
   FLOAT phi_analytic = 2 * KB * T * Log(x + Sqrt( x * x + 1));
   phgPrintf("phi_surface_analytic = %.8e\n", phi_analytic);
}

int
main(int argc, char *argv[])
{
    size_t mem, mem_peak;
    const char *fn = NULL;
    int pre_refines = 0;

    int iter_fp, maxit_fp = 100000;
    const FLOAT *err;
    FLOAT rtol = 1.0e-4, rtol_v = 1.0e-4;

    phgOptionsRegisterKeyword("poisson","discretization methods for Poisson's equations",poisson_names, &poisson_index);
    phgOptionsRegisterKeyword("continuity","discretization methods for continuity equations",continuity_names, &continuity_index);
    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("device", "device", &device);
    phgOptionsRegisterInt("surface", "surface", &surface);
    phgOptionsRegisterInt("analytic_test", "analytic_test", &analytic_test);
	phgOptionsRegisterFilename("fn_dofs", "dof section value output", &fn_dofs);

    /* PN */
    phgOptionsRegisterFloat("bias_cathode", "bias_cathode", &bias_cathode);
    phgOptionsRegisterFloat("bias_anode", "bias_anode", &bias_anode);
    phgOptionsRegisterFloat("anode_init", "anode_init", &anode_init);
    phgOptionsRegisterFloat("h", "bias increment", &h);
    phgOptionsRegisterFloat("surface_charge", "surface_charge", &surface_charge);
    phgOptionsRegisterFloat("bias_surface_charge", "bias surface charge", &bias_surface_charge);
    phgOptionsRegisterFloat("h_surface_charge", "bias surface charge increment", &h_surface_charge);
    phgOptionsRegisterFloat("analytic_density", "analytic_density", &analytic_density);
    phgOptionsRegisterFloat("bias_analytic_density", "bias analytic_density", &bias_analytic_density);
    phgOptionsRegisterFloat("h_analytic_density", "bias analytic_density increment", &h_analytic_density);
    phgOptionsRegisterFloat("analytic_period", "analytic_period", &analytic_period);
    phgOptionsRegisterFloat("analytic_potential", "analytic_potential", &analytic_potential);
    phgOptionsRegisterFloat("cylinder_r", "cylinder_r", &cylinder_r);
    phgOptionsRegisterFloat("cylinder_h", "cylinder_h", &cylinder_h);
    phgOptionsRegisterFloat("CM2MICRON", "CM2MICRON", &CM2MICRON);



    /** pre-set options **/
    phgOptionsPreset("-device 0");	/* 0 --> PN */
    phgOptionsPreset("-surface 1");
    phgOptionsPreset("-analytic_test 0");
    phgOptionsPreset("-mesh_file ./nano_mesh/nano_l1.mesh");
    phgOptionsPreset("-bias_anode 5.0 -bias_cathode 0.0 -h 0.1 -anode_init 0.0");
   phgOptionsPreset("-bias_surface_charge -2.0e-7 -h_surface_charge 2 -surface_charge -1.25e-8");
   phgOptionsPreset("-bias_analytic_density 1.0e-4 -h_analytic_density 2 -analytic_density 1.0e-4");
 //   phgOptionsPreset("-analytic_potential 0.0 -analytic_period 2.0");
    phgOptionsPreset("-cylinder_h 10.0 -cylinder_r 0.15");
    phgOptionsPreset("-CM2MICRON 1.0e5");
    phgOptionsPreset("-pre_refines 0");
    phgOptionsPreset("-poisson fem");
    phgOptionsPreset("-continuity Average");
    phgOptionsPreset("-fn_dofs ./dofs_l1_s0.txt");
    
phgOptionsPreset("-solver petsc -solver_maxit 100 -solver_rtol 1e-2 "
		     "-oem_options \"-ksp_type gmres -ksp_gmres_restart 50\" "
		     /*"-oem_options \"-ksp_monitor\""*/);
//  phgOptionsPreset("-solver petsc  -oem_options {-ksp_type gmres -ksp_gmres_restart 50 -ksp_pc_side right -pc_type asm -pc_asm_type restrict -pc_asm_overlap 2 -sub_ksp_type preonly -sub_pc_type lu }");
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);		/* WARNING: no idle processes are allowed */
    phgImportSetBdryMapFunc(bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines);
    phgBalanceGrid(g, 1.2, -1, NULL, 0.);

    InitConstant();

    phgPrintf
	("  %d vertices, %d elements, %d submeshes, load imbalance: %lg\n",
	 g->nvert_global, g->nleaf_global, g->nprocs, (double)g->lif);

    V = phgDofNew(g, DOF_P1, 1, "electric potential", DofInterpolation);
    P = phgDofNew(g, DOF_P1, 1, "hole concentration", DofInterpolation);
    N = phgDofNew(g, DOF_P1, 1, "electron concentration", DofInterpolation);

 //   Center(g, &tet_center, &tri_center, &volume);
    SetupVertsRegion(g, TRUE);
    /* initialize the device */
    if (device == PN) {
   anode = cathode = bias_cathode;	/* bias_cathode grounded */
 //  anode = anode_init;
  // cathode = bias_cathode;
	ohmic_contact = CATHODE | ANODE;
   if(analytic_test)
   {
      analytic_init();
   }
    }
    if(device == NANO)
    {
       anode = cathode = bias_cathode;
//        anode = anode_init;
  //      cathode = bias_cathode;
        ohmic_contact = CATHODE;
    }
	phgPrintf("\n *********** \n");

    InitialGuess(V, P, N);

   phgPrintf("\n  ================ EQUILIBRIUM STATE ================\n\n");
   FLOAT t0 = phgGetTime(NULL);
   for (iter_fp = 1; iter_fp <= maxit_fp; iter_fp++) {
  // for (iter_fp = 1; iter_fp <= 1; iter_fp++) {
	 RelativeError_new(V, P, N, FALSE);

	 phgPrintf("\n  ---- The %d-th Gummel iteration ----\n", iter_fp);

	 phgPrintf("\n  The Poisson equation:\n");
    SolvePoisson_new(V, P, N, 1.0);	/* to get thermal equilibrium state */
    //SolvePoisson(V, P, N, 1.0);	/* to get thermal equilibrium state */

	 phgPrintf("\n  The hole continuity equation:\n");
    SolveContinuity(V, P, N, HOLE);

    phgPrintf("\n  The electron continuity equation:\n");
    SolveContinuity(V, P, N, ELECTRON);


	    /* convergence criteria for gummel iterations */
	    err = RelativeError_new(V, P, N, TRUE);
	    phgPrintf("\n  rtol: V %0.6le, P %0.6le, N %0.6le\n", (double)err[0], (double)err[1], (double)err[2]);

	    if (err[0] < rtol_v && err[1] < rtol && err[2] < rtol) {
		phgPrintf("  Converged!\n");
		break;
	    }
   }
   FLOAT t1 = phgGetTime(NULL);
   phgPrintf("WallTime = %.4e s\n", t1 - t0);


	   if (iter_fp > maxit_fp)
	    phgError(1, "  Fail to converge.\n");

	/* compute terminal currents */
	if (device == PN) {
	    CONTACT cont[2] = { ANODE, CATHODE };
	    FLOAT currents[2][3];
	    Current(V, P, N, cont, 2, currents);
	    phgPrintf("\n  Currents (V_a: %0.4lfV, V_c: %0.4lfV): \n",
		      (double)anode, (double)cathode);
	    phgPrintf(" Total anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][0]),
		      (double)(Q * currents[1][0]));
	    phgPrintf(" Electron anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][1]),
		      (double)(Q * currents[1][1]));
	    phgPrintf("Hole anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][2]),
		      (double)(Q * currents[1][2]));
	}
   if(device == NANO)
   {
   if(surface)
      Surface_phi(V);
   }


    while (TRUE) {		/* bias ramp */
      if(BiasSweep())
	    break;

	TakeBoundary(V, P, N, 2);

	/* Gummel iteration for each bias */
   FLOAT t0 = phgGetTime(NULL);
	for (iter_fp = 1; iter_fp <= maxit_fp; iter_fp++) {
	    RelativeError(V, P, N, FALSE);

	    phgPrintf("\n  ---- The %d-th Gummel iteration ----\n", iter_fp);

	    phgPrintf("\n  The Poisson equation:\n");
	 //  SolvePoisson(V, P, N, 1.0);
	  SolvePoisson_new(V, P, N, 1.0);

     // phgPrintf("\n  The electron continuity equation:\n");
	 //  SolveContinuity(V, P, N, ELECTRON);

	    phgPrintf("\n  The hole continuity equation:\n");
	    SolveContinuity(V, P, N, HOLE);

            phgPrintf("\n  The electron continuity equation:\n");
            SolveContinuity(V, P, N, ELECTRON);

	    /* convergence criteria for gummel iterations */
	    err = RelativeError(V, P, N, TRUE);
	    phgPrintf("\n  rtol: V %0.6le, P %0.6le, N %0.6le\n",
		      (double)err[0], (double)err[1], (double)err[2]);

	    if (err[0] < rtol_v && err[1] < rtol && err[2] < rtol) {
		phgPrintf("  Converged!\n");
		break;
	    }
	}
   FLOAT t1 = phgGetTime(NULL);
   phgPrintf("WallTime = %.4e s\n", t1 - t0);

	if (iter_fp > maxit_fp)
	    phgError(1, "  Fail to converge.\n");

	/* compute terminal currents */
	if (device == PN) {
	    CONTACT cont[2] = { ANODE, CATHODE };
	    FLOAT currents[2][3];
	    Current(V, P, N, cont, 2, currents);
	    phgPrintf("\n  Currents (V_a: %0.4lfV, V_c: %0.4lfV): ",
		      (double)anode, (double)cathode);
	    phgPrintf("anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][0]),
		      (double)(Q * currents[1][0]));
	    phgPrintf("anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][1]),
		      (double)(Q * currents[1][1]));
	    phgPrintf("anode: %0.6le, cathode %0.6le\n",
		      (double)(Q * currents[0][2]),
		      (double)(Q * currents[1][2]));

#if 0
       FLOAT current_p = -1.0 * Dp * Q * analytic_density * NA / Pow(DM2CM,3) * anode * Beta * M_PI * cylinder_r * cylinder_r / cylinder_h/ CM2MICRON;
       FLOAT current_n = -1.0 * Dn * Q * analytic_density * NA / Pow(DM2CM,3) * anode * Beta * M_PI * cylinder_r * cylinder_r / cylinder_h / CM2MICRON;
    //   FLOAT current_p = -1.0 * Dp * Q * analytic_density * NA / Pow(DM2CM,3) * anode * Beta *  cylinder_r * cylinder_r / cylinder_h/ CM2MICRON;
    //   FLOAT current_n = -1.0 * Dn * Q * analytic_density * NA / Pow(DM2CM,3) * anode * Beta *  cylinder_r * cylinder_r / cylinder_h / CM2MICRON;
       phgPrintf("analytic current: %0.6le  %0.6le\n", (double)current_n, (double)current_p);
#endif
#if 0
if(Fabs(anode - 1.0) < 1.0e-4 || Fabs(anode - 2.0) < 1.0e-4 || Fabs(anode - 3.0) < 1.0e-4 || Fabs(anode - 4.0) < 1.0e-4 || Fabs(anode - 5.0) < 1.0e-4){
//if(Fabs(anode - 5.0) < 1.0e-4){
         if(phgRank == 1 && fn_dofs != NULL)
         {
            FILE *f_in;
            f_in = fopen(fn_dofs, "a");
            if(f_in == NULL)
               phgError(1, "Can not open dofs output file %s.\n", fn_dofs);
            fprintf(f_in, "\n anode = %0.4lf V\n", (double)anode);
            fclose(f_in);
         }
         Section_values(V, P, N);
       }
#endif
	}
   if(device == NANO)
   {
   if(surface)
      Surface_phi(V);
   }

#if 0
if(Fabs(anode - 1.0) < 1.0e-4 || Fabs(anode - 2.0) < 1.0e-4 || Fabs(anode - 3.0) < 1.0e-4 || Fabs(anode - 4.0) < 1.0e-4 || Fabs(anode - 5.0) < 1.0e-4){
//if(Fabs(anode - 5.0) < 1.0e-4){
    /* rescale the data for output */
    INT i;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) *= KB * T;
//	*(V->data + i) *= 1.0;
	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;
	*(P->data + i) *= SCALE / NA  * Pow(DM2CM, 3);
	*(N->data + i) *= SCALE / NA  * Pow(DM2CM, 3);
    }

    char vtk_name[40];
    sprintf(vtk_name, "nano_l1_s0_%0.2f.vtk", anode);
    phgPrintf("\n  \"%s\" created.\n", phgExportVTK(g, vtk_name, V, P, N, NULL));

    /* rescale the data for input */
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) /= KB * T;
//	*(V->data + i) *= 1.0;
	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;
	*(P->data + i) /= SCALE / NA  * Pow(DM2CM, 3);
	*(N->data + i) /= SCALE / NA  * Pow(DM2CM, 3);
    }
}

#endif
    }				/* end bias ramp */

    /* rescale the data for output */
    INT i;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	*(V->data + i) *= KB * T;
//	*(V->data + i) *= 1.0;
	if (!(verts_region[i] & MAT_SEMICONDUCTOR))
	    continue;
	*(P->data + i) *= SCALE / NA  * Pow(DM2CM, 3);
	*(N->data + i) *= SCALE / NA  * Pow(DM2CM, 3);
  /* *(P->data + i) *= SCALE;
   *(N->data + i) *= SCALE;*/
    }
    phgPrintf("\n||V|| = %.8e, ||P|| = %.8e, ||N|| = %.8e\n", phgDofNormL2(V), phgDofNormL2(P), phgDofNormL2(N));

#if 1
    char vtk_name[40];
    sprintf(vtk_name, "nano_l1.vtk");
    phgPrintf("\n  \"%s\" created.\n", phgExportVTK(g, vtk_name, V, P, N, NULL));
#endif
    if(analytic_test)
    {
       analytic_print();
       analytic_finalize();
    }
//    phgDofFree(&tet_center);
//    phgDofFree(&tri_center);
//    phgDofFree(&volume);
//    SetupVertsRegion(g, FALSE);

    phgPrintf
	("  %d vertices, %d elements, %d submeshes, load imbalance: %lg\n",
	 g->nvert_global, g->nleaf_global, g->nprocs, (double)g->lif);

    mem = phgMemoryUsage(g, &mem_peak);
    phgPrintf("\n  Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
	      (double)mem / (1024.0 * 1024.0),
	      (double)mem_peak / (1024.0 * 1024.0));

    phgDofFree(&V);
    phgDofFree(&P);
    phgDofFree(&N);

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
