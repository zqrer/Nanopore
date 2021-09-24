#ifndef PNP_COEF_H
#define PNP_COEF_H

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"phg.h"

static const char *poisson_names[] = { "fvm", "fem", NULL };
static const char *continuity_names[] =
    { "box", "mixed", "zlamal", "exp", "Average", NULL };
static int poisson_index = 0;
static int continuity_index = 0;

static enum { PN = 0, NANO = 1} device = PN;

typedef enum { POISSON = 0, ELECTRON = 1, HOLE = 2 } EQUATION;

typedef enum { REGION_N = 1, REGION_P = 2} REGION;

typedef enum { CATHODE = BDRY_USER1, ANODE = BDRY_USER2,
    JUNCTION_PN = BDRY_USER3, /* PN DIODE */
} CONTACT;

typedef enum { MAT_SEMICONDUCTOR = 1} MATERIAL;

/*---- constants from Sentaurus-Device manual ----*/
static const FLOAT KB = 8.617343e-5 /* eV/K */ ;	/* Boltzmann constant */
static const FLOAT Q = 1.602176487e-19 /* C */ ;	/* charge of a proton */
static const FLOAT EPSILONE0 = 8.8541878176e-12 /* F/m */ ;	/* vacumm permittivity */
static FLOAT CM2MICRON = 1.0e4;
static const FLOAT T = 300.0 /* K */ ;	/* absolute temperture */

/* silicon relative permittivity */
static const FLOAT EPSILONE = 80.0;
/* metal */
/* constant mobility model */
static FLOAT nie /* 1/cm^3 */ ;	/* intrinsic carrier concentration */
const FLOAT Dn = 2.03e-5 /* cm^2/s */ ;
const FLOAT Dp = 1.96e-5/* cm^2/s */ ;
const FLOAT NA = 6.022e+23 /*Avogadro constant*/; 
const FLOAT DM2CM = 10.0;

/*---- device specifications ----*/
static BTYPE ohmic_contact;
/** PN diode **/
static FLOAT bias_anode = 0.05, bias_cathode = 0. /* V */ ;
static FLOAT anode = 0., cathode = 0. /* V */ ;	/* runtime */
static FLOAT anode_init = 0.0;/*V*/
static FLOAT h = 0.05 /* V */ ;	/* bias increment */
static FLOAT anode_c = 0.5, cathode_c = 0.01 /*mol/L*/;
static FLOAT bias_anode_c = 0.5, bias_cathode_c = 0.01 /*mol/L*/;
static FLOAT anode_init_c = 0.01 /*mol/L*/;
static FLOAT h_c = 0.01 /*mol/L*/;

static SHORT *verts_region = NULL;
static const FLOAT SCALE = 1.0;

int surface  = 0;
static FLOAT surface_charge = 0.0 /* C/cm^2*/;
static FLOAT bias_surface_charge = -1.0e-7;
static FLOAT h_surface_charge = -5.0e-9; 

char *fn_dofs = NULL;

int analytic_test = 1; 
FLOAT analytic_period = 20.0;
FLOAT analytic_density = 0.1 /* mol/L */;
FLOAT analytic_potential = 0.0 /* V */;
FLOAT cylinder_r = 5.0;
FLOAT cylinder_h = 50.0;
FLOAT Kappa_2;
FLOAT Beta;
FLOAT bias_analytic_density = 0.1;
FLOAT h_analytic_density = 1.0;

DOF *analytic_p, *analytic_n, *analytic_u;
DOF *analytic_f_p, *analytic_f_n;
DOF *analytic_grad_u, *analytic_grad_n, *analytic_grad_p;
DOF *analytic_err_u, *analytic_err_p, *analytic_err_n;
DOF *analytic_err_grad_u, *analytic_err_grad_p, *analytic_err_grad_n;

DOF *V, *P, *N;
GRID *g;
SIMPLEX *e;

static void
InitConstant(void)
{
    /* the intrinsic carrier concentration,
     * bandgap narrowing effect neglected */
    nie =  0.1 * analytic_density * NA / Pow(DM2CM, 3) /* 1/cm^3 */ ;
    Kappa_2 = 100.0 * Q / (EPSILONE0 * KB * T); /* cm */
    Beta = 1 / (KB * T);
    if(analytic_test)
    {
      bias_anode = analytic_potential;
      anode_init = analytic_potential;
    }
}

#endif
