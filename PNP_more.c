#include"phg.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdarg.h>
#include<math.h>
#include"PNP_coefficient.h"


typedef void PNP_FUNC_1(FLOAT *);

/* For Calculate current */
//calculate the coordinate of the point c on the line a,b with c[d] = value
//lambda * a + (1 - lambda) * b = c
//lambda * a[d] + (1 - lambda) * b[2] = z
//lambda = (z - b[d]) / (a[d] - b[d])
//1 - lambda = (a[d] - z) / (a[d] - b[d])
//d short for direction d = 1, 2, ..., Dim;
static void linear_interpolation(COORD a, COORD b, COORD c, FLOAT value, int d) {
	if(a[d] == b[d]) {
		phgError(1, "\ninput points have same direction value");
	}
	else {
		int i;
//		c = (COORD *)malloc(sizeof(COORD));
		c[d] = value;
		for(i = 1; i < Dim; i++) {
			c[(d+i)%Dim] = (value - b[d]) / (a[d] - b[d]) * a[(d+i)%Dim]
				     + (value - a[d]) / (b[d] - a[d]) * b[(d+i)%Dim];
		}
	}
}

// (X, Y) to save all pairs of triangle vert (x[i], y[i]) i = 1, 2, 3, with X, Y belong to R^3
static void triangle_area(FLOAT *X, FLOAT *Y, FLOAT *area) {
	int i;
	*area = 0;
	for(i = 0; i < 3; i++) {
		*area += (X[i] * Y[(i+1)%3] - X[(i+1)%3] * Y[i]);
	}
	*area = Fabs(*area) * 0.5;
	free(X);
	free(Y);
}

// (X, Y) to save all pairs of tetrahetron vert (x[i], y[i]) i = 1, 2, 3, 4, with X, Y belong to R^4
static void tetrahetron_area(FLOAT *X, FLOAT *Y, FLOAT *area) {
	FLOAT *child_area = (FLOAT *)malloc(4 * sizeof(FLOAT));
	FLOAT *x, *y;
	int i, j;
	*area = 0;
	for(i = 0; i < 4; i++) {
		x = (FLOAT *)malloc(3 * sizeof(FLOAT));
		y = (FLOAT *)malloc(3 * sizeof(FLOAT));
		for(j = 0; j < 3; j++) {
			x[j] = X[(i+j)%4];
			y[j] = Y[(i+j)%4];
		}
		triangle_area(x, y, &child_area[i]);
		*area += child_area[i] * 0.5;
	}
}

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
#if 0
static void Calculate_current(int region, int d, DOF *u, DOF **p, FLOAT value) {
	FLOAT sum_current = 0;
	FLOAT p_current[2];
	FLOAT current[2];
	FLOAT current1[2], current2[2];
	FLOAT p_current1[2], p_current2[2];
	SIMPLEX *e;
	GRID *g = u->g;
	FLOAT area, *X, *Y;
	FLOAT sum_area = 0;
	FLOAT total_area = 0;
	
	FLOAT Z[2] = {1, -1};
	FLOAT D[2] = {Dp, Dn};

	COORD triangle_verts[3];
	COORD tetrahetron_verts[4];

	int i, j, k, l, vert;
	int type_flag;//number of p_tab
	int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on d direction larger than value
	int n_tab[NVert];//if smaller than value

	memset(current, 0, 2*sizeof(FLOAT));
	memset(current1, 0, 2*sizeof(FLOAT));
	memset(current2, 0, 2*sizeof(FLOAT));
	memset(p_current, 0, 2*sizeof(FLOAT));
	memset(p_current1, 0, 2*sizeof(FLOAT));
	memset(p_current2, 0, 2*sizeof(FLOAT));

	ForAllElements(g, e) {
//		if(e->region_mark == region) {
			memset(p_tab, 0, NVert*sizeof(int));
			memset(n_tab, 0, NVert*sizeof(int));
			type_flag = 0;
			i = 0;
			for(vert = 0; vert < NVert; vert++) {
				if(g->verts[e->verts[vert]][d] >= value) {
					p_tab[type_flag++] = vert;
				}
				else {
					n_tab[i++] = vert;
				}
			}
			if(type_flag == 1 || type_flag == 3) {
				for(i = 0; i < 3; i++) {
					if(type_flag == 1)
						linear_interpolation(g->verts[e->verts[p_tab[0]]],
								     g->verts[e->verts[n_tab[i]]], 
								     triangle_verts[i], value, d);
					if(type_flag == 3)
						linear_interpolation(g->verts[e->verts[n_tab[0]]],
								     g->verts[e->verts[p_tab[i]]], 
								     triangle_verts[i], value, d);
				}
				X = (FLOAT *)malloc(3 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
				for(i = 0; i < 3; i++) {
					X[i] = triangle_verts[i][(d+1)%3];
					Y[i] = triangle_verts[i][(d+2)%3];
				}
				triangle_area(X, Y, &area);
			}
			else if(type_flag == 2) {
				for(i = 0; i < 2; i++) {
					for(j = 0; j < 2; j++) {
						linear_interpolation(g->verts[e->verts[p_tab[i]]],
								     g->verts[e->verts[n_tab[j]]], 
								     tetrahetron_verts[2*i+j], value, d);
					}
				}
				X = (FLOAT *)malloc(4 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
				for(i = 0; i < 4; i++) {
					X[i] = tetrahetron_verts[i][(d+1)%3];
					Y[i] = tetrahetron_verts[i][(d+2)%3];
				}
				tetrahetron_area(X, Y, &area);
			}
			if(type_flag > 0 && type_flag < 4) {
				sum_area += area;
				for(i = 0; i < 2; i++) {
					if(type_flag != 2)  {
						//PNP_quad_channel_section(&value, d, e, triangle_verts, 0, NULL, 1, p[i], NULL);
						//current[i] += value;
						//PNP_quad_channel_section(&value, d, e, triangle_verts, 0, NULL, 1, u, p[i], NULL);
						//current[i] += ion[i].Z * value;
						current[i] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, p[i], NULL);
						current1[i] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, p[i], NULL);
						current[i] += Z[i] * PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, u, p[i], NULL);
						current2[i] += Z[i] * PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, u, p[i], NULL);
					}
					else {
						for(j = 0; j < 4; j++) {
							for(k = 0; k < 3; k++) {
								for(l = 0; l < 3; l++) {
									triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
								}
							}
							//PNP_quad_channel_section(&value, d, e, triangle_verts, 0, NULL, 1, p[i], NULL);
							//current[i] += value / 2.0;
							//PNP_quad_channel_section(&value, d, e, triangle_verts, 0, NULL, 1, u, p[i], NULL);
							//current[i] += ion[i].Z * value / 2.0;
							current[i] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, p[i], NULL) / 2.0;
							current1[i] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, p[i], NULL) / 2.0;
							current[i] += Z[i] * PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, u, p[i], NULL) / 2.0;
							current2[i] += Z[i] * PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 1, u, p[i], NULL) / 2.0;
						}
					}
				}
			}
//		}
	}
   FLOAT Ratio = 1.0;
   if(V_D)
   {
      func_D(0, 0, value, &Ratio);
   }
	for(i = 0; i < 2; i++) {
		//p: N/cm^3 | D: cm^2/s | I: A ->pA | grad and quad: mesh_unit -> cm
		current[i] *= -D[i] * Z[i] / CM2MICRON * 1.0e12 * Q * Ratio;
		current1[i] *= -D[i] * Z[i] / CM2MICRON * 1.0e12 * Q * Ratio;
		current2[i] *= -D[i] * Z[i] / CM2MICRON * 1.0e12 * Q * Ratio;
		MPI_Allreduce(&current[i], &p_current[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		MPI_Allreduce(&current1[i], &p_current1[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		MPI_Allreduce(&current2[i], &p_current2[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		sum_current += p_current[i];
		phgPrintf("\nion[%d]	current is %e pA", i, p_current[i]);
	}
	//print current section area
	phgPrintf("\n\ntotal	current is %e pA\n", sum_current);
	phgPrintf("\ntotal area is %e A^2\n", total_area);
	if(phgRank == 1 && fn_current != NULL) {
		FILE *f_in;
		f_in = fopen(fn_current, "a");
		if(f_in == NULL) {
			phgError(1, "Can not open current output file %s.\n", fn_current);
		}
		if(double_current)
			fprintf(f_in, "%le	%le	%le	%le	%le	%le	%le	%le	%le\n", anode, bias_surface_charge, volume_charge, value, 
				2 * p_current1[0], 2 * p_current1[1], 2 * p_current2[0], 2 * p_current2[1], 2 * SCALE * sum_current);
		else
			fprintf(f_in, "%le	%le	%le	%le	%le	%le	%le	%le	%le\n", anode, surface_charge, volume_charge, value, 
				p_current1[0], p_current1[1], p_current2[0], p_current2[1], SCALE * sum_current);
//			fprintf(f_in, "%le	%le	%le	%le	%le	%le	%le	%le\n", surface_charge, volume_charge, value, 
//				p_current1[0], p_current1[1], p_current2[0], p_current2[1], SCALE * sum_current);
		fclose(f_in);
	}
}
#endif
FLOAT Calculate_potential(int region, int d, DOF *u, FLOAT value) {
	FLOAT sum_potential = 0;
	FLOAT potential = 0;
	SIMPLEX *e;
	GRID *g = u->g;
	FLOAT area, *X, *Y;
	FLOAT sum_area = 0;
	FLOAT total_area = 0;

	COORD triangle_verts[3];
	COORD tetrahetron_verts[4];

	int i, j, k, l, vert;
	int type_flag;//number of p_tab
	int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on d direction larger than value
	int n_tab[NVert];//if smaller than value

	ForAllElements(g, e) {
//		if(e->region_mark == region) {
			memset(p_tab, 0, NVert*sizeof(int));
			memset(n_tab, 0, NVert*sizeof(int));
			type_flag = 0;
			i = 0;
			for(vert = 0; vert < NVert; vert++) {
				if(g->verts[e->verts[vert]][d] >= value) {
					p_tab[type_flag++] = vert;
				}
				else {
					n_tab[i++] = vert;
				}
			}
			if(type_flag == 1 || type_flag == 3) {
				for(i = 0; i < 3; i++) {
					if(type_flag == 1)
						linear_interpolation(g->verts[e->verts[p_tab[0]]],
								     g->verts[e->verts[n_tab[i]]], 
								     triangle_verts[i], value, d);
					if(type_flag == 3)
						linear_interpolation(g->verts[e->verts[n_tab[0]]],
								     g->verts[e->verts[p_tab[i]]], 
								     triangle_verts[i], value, d);
				}
				X = (FLOAT *)malloc(3 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
				for(i = 0; i < 3; i++) {
					X[i] = triangle_verts[i][(d+1)%3];
					Y[i] = triangle_verts[i][(d+2)%3];
				}
				triangle_area(X, Y, &area);
			}
			else if(type_flag == 2) {
				for(i = 0; i < 2; i++) {
					for(j = 0; j < 2; j++) {
						linear_interpolation(g->verts[e->verts[p_tab[i]]],
								     g->verts[e->verts[n_tab[j]]], 
								     tetrahetron_verts[2*i+j], value, d);
					}
				}
				X = (FLOAT *)malloc(4 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
				for(i = 0; i < 4; i++) {
					X[i] = tetrahetron_verts[i][(d+1)%3];
					Y[i] = tetrahetron_verts[i][(d+2)%3];
				}
				tetrahetron_area(X, Y, &area);
			}
			if(type_flag > 0 && type_flag < 4) {
				sum_area += area;
				if(type_flag != 2)  {
					potential += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, u, NULL);
				}
				else {
					for(j = 0; j < 4; j++) {
						for(k = 0; k < 3; k++) {
							for(l = 0; l < 3; l++) {
								triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
							}
						}
						potential += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, u, NULL) / 2.0;
					}
				}
			}
//		}
	}

	potential *= KB * T;
	MPI_Allreduce(&potential, &sum_potential, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
	sum_potential /= total_area;

	phgPrintf("\ntotal area is %e nm^2\n", total_area);
	
	return sum_potential;
}

static void
Calculate_section_values(DOF *V, DOF *P, DOF *N, FLOAT pos, FLOAT section_values[3])   
   /*section_values[0] : section_value for V,
     section_values[1] : section_value for P
     section_values[2] : section_value for N*/
{
	SIMPLEX *e;
	GRID *g = V->g;
	FLOAT area, *X, *Y;
	FLOAT sum_area = 0;
	FLOAT total_area = 0;

	COORD triangle_verts[3];
	COORD tetrahetron_verts[4];

   int d = 2;/* z axis */
	int i, j, k, l, vert;
	int type_flag;//number of p_tab
	int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on z direction larger than pos
	int n_tab[NVert];//if smaller than pos
  
   bzero(section_values, 3 * sizeof(section_values[0]));

	ForAllElements(g, e) {
//		if(e->region_mark == region) {
			memset(p_tab, 0, NVert * sizeof(int));
			memset(n_tab, 0, NVert * sizeof(int));
			type_flag = 0;
			i = 0;
			for(vert = 0; vert < NVert; vert++) {
				if(g->verts[e->verts[vert]][d] >= pos) {
					p_tab[type_flag++] = vert;
				}
				else {
					n_tab[i++] = vert;
				}
			}
			if(type_flag == 1 || type_flag == 3) {
				for(i = 0; i < 3; i++) {
					if(type_flag == 1)
						linear_interpolation(g->verts[e->verts[p_tab[0]]],
								     g->verts[e->verts[n_tab[i]]], 
								     triangle_verts[i], pos, d);
					if(type_flag == 3)
						linear_interpolation(g->verts[e->verts[n_tab[0]]],
								     g->verts[e->verts[p_tab[i]]], 
								     triangle_verts[i], pos, d);
				}
				X = (FLOAT *)malloc(3 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
				for(i = 0; i < 3; i++) {
					X[i] = triangle_verts[i][(d+1)%3];
					Y[i] = triangle_verts[i][(d+2)%3];
				}
				triangle_area(X, Y, &area);
			}
			else if(type_flag == 2) {
				for(i = 0; i < 2; i++) {
					for(j = 0; j < 2; j++) {
						linear_interpolation(g->verts[e->verts[p_tab[i]]],
								     g->verts[e->verts[n_tab[j]]], 
								     tetrahetron_verts[2*i+j], pos, d);
					}
				}
				X = (FLOAT *)malloc(4 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
				for(i = 0; i < 4; i++) {
					X[i] = tetrahetron_verts[i][(d+1)%3];
					Y[i] = tetrahetron_verts[i][(d+2)%3];
				}
				tetrahetron_area(X, Y, &area);
			}
			if(type_flag > 0 && type_flag < 4) {
				sum_area += area;
				if(type_flag != 2)  {
					section_values[0] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, V, NULL);
					section_values[1] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, P, NULL);
					section_values[2] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, N, NULL);
				}
				else {
					for(j = 0; j < 4; j++) {
						for(k = 0; k < 3; k++) {
							for(l = 0; l < 3; l++) {
								triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
							}
						}
						section_values[0] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, V, NULL) / 2.0;
						section_values[1] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, P, NULL) / 2.0;
						section_values[2] += PNP_quad_channel_section(NULL, d, e, triangle_verts, 0, NULL, 0, N, NULL) / 2.0;
					}
				}
			}
//		}
	}

	section_values[0] *= KB * T;
   section_values[1] *= SCALE / NA * Pow(DM2CM, 3);
   section_values[2] *= SCALE / NA * Pow(DM2CM, 3);
#if USE_MPI
      if (g->nprocs > 1) {
         FLOAT section_values0[3];
         memcpy(section_values0, section_values, 3 * sizeof(section_values0[0]));
         MPI_Allreduce(section_values0, section_values, 3, PHG_MPI_FLOAT, MPI_SUM, g->comm);
	      MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
      }
#endif
   for(i = 0; i < 2; i++)
      section_values[i] /= total_area;

//	phgPrintf("\ntotal area is %e nm^2\n", total_area);
}
static void line_integrate(FLOAT r, FLOAT z, int Num, FLOAT *GridValue, DOF *dof)
{
    int i;
    COORD *coord = phgAlloc(Num * sizeof(COORD));
    FLOAT theta = 2 * M_PI / Num;
    for(i = 0; i < Num; i++)
    {
       coord[i][0] = r * Cos(i * theta);
       coord[i][1] = r * Sin(i * theta);
       coord[i][2] = z;
    }
    phgInterGridDofEval(dof, Num, coord, GridValue, 1);
}
