#ifndef ANALYTIC_C
#define ANALYTIC_C

#define Cosx Cos(M_PI / analytic_period * x)
#define Cosy Cos(M_PI / analytic_period * y)
#define Cosz Cos(M_PI / analytic_period * z)
#define Sinx Sin(M_PI / analytic_period * x)
#define Siny Sin(M_PI / analytic_period * y)
#define Sinz Sin(M_PI / analytic_period * z)

static void func_analytic_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);
	*value = tmp_density + 0.5 * tmp_density * Cosx * Cosy * Cosz;
}

static void func_analytic_n(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);
	*value = tmp_density - 0.5 * tmp_density * Cosx * Cosy * Cosz;
}

static void func_grad_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *values) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);
	values[0] = -0.5 * tmp_density * M_PI / analytic_period * Sinx * Cosy * Cosz;
	values[1] = -0.5 * tmp_density * M_PI / analytic_period * Cosx * Siny * Cosz;
	values[2] = -0.5 * tmp_density * M_PI / analytic_period * Cosx * Cosy * Sinz;
}

static void func_grad_n(FLOAT x, FLOAT y, FLOAT z, FLOAT *values) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);
	values[0] = 0.5 * tmp_density * M_PI / analytic_period * Sinx * Cosy * Cosz;
	values[1] = 0.5 * tmp_density * M_PI / analytic_period * Cosx * Siny * Cosz;
	values[2] = 0.5 * tmp_density * M_PI / analytic_period * Cosx * Cosy * Sinz;
}

static void func_analytic_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);
	FLOAT tmp_potential = analytic_potential * Beta;	
	*value = Kappa_2 * tmp_density / 3.0 / EPSILONE / Pow(M_PI / analytic_period * CM2MICRON, 2) * Cosx * Cosy * Cosz + tmp_potential * (z + analytic_period / 2) / analytic_period;
}

static void func_grad_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *values) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);	
	FLOAT tmp_potential = analytic_potential * Beta;
	values[0] = -Kappa_2 * tmp_density / 3.0 / EPSILONE / (M_PI / analytic_period) / Pow(CM2MICRON, 2) * Sinx * Cosy * Cosz;
	values[1] = -Kappa_2 * tmp_density / 3.0 / EPSILONE / (M_PI / analytic_period) / Pow(CM2MICRON, 2) * Cosx * Siny * Cosz;
	values[2] = -Kappa_2 * tmp_density / 3.0 / EPSILONE / (M_PI / analytic_period) / Pow(CM2MICRON, 2) * Cosx * Cosy * Sinz + tmp_potential / analytic_period;
}

static void func_f_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);	
	FLOAT tmp_potential = analytic_potential * Beta;
/*	*value = -1.5 * tmp_density * Pow(M_PI / analytic_period * CM2MICRON, 2) * Cosx * Cosy * Cosz;   //laplace p
	*value += Kappa_2 * 0.5 * tmp_density * tmp_density / 3.0 / EPSILONE * (Pow(Sinx * Cosy * Cosz, 2) + Pow(Cosx * Siny * Cosz, 2) + Pow(Cosx * Cosy * Sinz, 2)) - tmp_density * tmp_potential / 2 / analytic_period * (M_PI / analytic_period) * Pow(CM2MICRON, 2) * Cosx * Cosy * Sinz;      //grad p grad u
	*value -= Kappa_2 * tmp_density * tmp_density / EPSILONE * (Cosx * Cosy * Cosz + 0.5 * Pow(Cosx * Cosy * Cosz, 2));                                                               // p laplace u*/
	*value = -1.5 * tmp_density * Pow(M_PI / analytic_period, 2) * Cosx * Cosy * Cosz;   //laplace p
	*value += Kappa_2 * 0.5 * tmp_density * tmp_density / 3.0 / EPSILONE / Pow(CM2MICRON, 2) * (Pow(Sinx * Cosy * Cosz, 2) + Pow(Cosx * Siny * Cosz, 2) + Pow(Cosx * Cosy * Sinz, 2)) - tmp_density * tmp_potential / 2 / analytic_period * (M_PI / analytic_period) * Cosx * Cosy * Sinz;      //grad p grad u
	*value -= Kappa_2 * tmp_density * tmp_density / EPSILONE / Pow(CM2MICRON, 2) * (Cosx * Cosy * Cosz + 0.5 * Pow(Cosx * Cosy * Cosz, 2));                                                               // p laplace u
	*value *= -1.0/SCALE;
}

static void func_f_n(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{
	FLOAT tmp_density = analytic_density * NA / Pow(DM2CM, 3);	
	FLOAT tmp_potential = analytic_potential * Beta;
/*	*value = 1.5 * tmp_density * Pow(M_PI / analytic_period * CM2MICRON, 2) * Cosx * Cosy * Cosz;   //laplace n
	*value += Kappa_2 * 0.5 * tmp_density * tmp_density / 3.0 / EPSILONE * (Pow(Sinx * Cosy * Cosz, 2) + Pow(Cosx * Siny * Cosz, 2) + Pow(Cosx * Cosy * Sinz, 2)) - tmp_density * tmp_potential / 2 / analytic_period * (M_PI / analytic_period) * Pow(CM2MICRON, 2) * Cosx * Cosy * Sinz;      //grad n grad u
	*value += Kappa_2 * tmp_density * tmp_density / EPSILONE * (Cosx * Cosy * Cosz - 0.5 * Pow(Cosx * Cosy * Cosz, 2));                                                               // n laplace u
   */
	*value = 1.5 * tmp_density * Pow(M_PI / analytic_period, 2) * Cosx * Cosy * Cosz;   //laplace n
	*value += Kappa_2 * 0.5 * tmp_density * tmp_density / 3.0 / EPSILONE / Pow(CM2MICRON, 2) * (Pow(Sinx * Cosy * Cosz, 2) + Pow(Cosx * Siny * Cosz, 2) + Pow(Cosx * Cosy * Sinz, 2)) - tmp_density * tmp_potential / 2 / analytic_period * (M_PI / analytic_period) * Cosx * Cosy * Sinz;      //grad n grad u
	*value += Kappa_2 * tmp_density * tmp_density / EPSILONE / Pow(CM2MICRON, 2) * (Cosx * Cosy * Cosz - 0.5 * Pow(Cosx * Cosy * Cosz, 2));                                                               // n laplace u
	*value *= -1.0/SCALE;
}

static void analytic_init() 
{
	FLOAT u_coef = KB * T;
	FLOAT p_coef = Pow(DM2CM, 3) / NA;
	const FLOAT zeros[3] = {0.0, 0.0, 0.0};

	if(analytic_test == 1) {
	//	analytic_p = phgDofNew(g, DOF_DEFAULT, 1, "analytic_p", func_analytic_p);
	//	analytic_n = phgDofNew(g, DOF_DEFAULT, 1, "analytic_n", func_analytic_n);
	//	analytic_u = phgDofNew(g, DOF_DEFAULT, 1, "analytic_u", func_analytic_u);

		analytic_p = phgDofNew(g, DOF_P1, 1, "analytic_p", func_analytic_p);
		analytic_n = phgDofNew(g, DOF_P1, 1, "analytic_n", func_analytic_n);
		analytic_u = phgDofNew(g, DOF_P1, 1, "analytic_u", func_analytic_u);

//		analytic_f_p = phgDofNew(g, DOF_DEFAULT, 1, "analytic_f_p", func_f_p);
//		analytic_f_n = phgDofNew(g, DOF_DEFAULT, 1, "analytic_f_n", func_f_n);

		analytic_f_p = phgDofNew(g, DOF_P1, 1, "analytic_f_p", func_f_p);
		analytic_f_n = phgDofNew(g, DOF_P1, 1, "analytic_f_n", func_f_n);

	//	analytic_grad_u = phgDofNew(g, DOF_DEFAULT, 3, "analytic_grad_u", func_grad_u);
	//	analytic_grad_p = phgDofNew(g, DOF_DEFAULT, 3, "analytic_grad_p", func_grad_p);
	//	analytic_grad_n = phgDofNew(g, DOF_DEFAULT, 3, "analytic_grad_n", func_grad_n);

		analytic_grad_u = phgDofNew(g, DOF_P1, 3, "analytic_grad_u", func_grad_u);
		analytic_grad_p = phgDofNew(g, DOF_P1, 3, "analytic_grad_p", func_grad_p);
		analytic_grad_n = phgDofNew(g, DOF_P1, 3, "analytic_grad_n", func_grad_n);

		//analytic_err_u = phgDofNew(g, DOF_DEFAULT, 1, "analytic_err_u", DofInterpolation);
	//	analytic_err_p = phgDofNew(g, DOF_DEFAULT, 1, "analytic_err_p", DofInterpolation);
	//	analytic_err_n = phgDofNew(g, DOF_DEFAULT, 1, "analytic_err_n", DofInterpolation);

		analytic_err_u = phgDofNew(g, DOF_P1, 1, "analytic_err_u", DofInterpolation);
		analytic_err_p = phgDofNew(g, DOF_P1, 1, "analytic_err_p", DofInterpolation);
		analytic_err_n = phgDofNew(g, DOF_P1, 1, "analytic_err_n", DofInterpolation);

//	analytic_err_grad_u = phgDofNew(g, DOF_DEFAULT, 3, "analytic_err_grad_u", DofInterpolation);
//		analytic_err_grad_p = phgDofNew(g, DOF_DEFAULT, 3, "analytic_err_grad_p", DofInterpolation);
//		analytic_err_grad_n = phgDofNew(g, DOF_DEFAULT, 3, "analytic_err_grad_n", DofInterpolation);

		analytic_err_grad_u = phgDofNew(g, DOF_P1, 3, "analytic_err_grad_u", DofInterpolation);
		analytic_err_grad_p = phgDofNew(g, DOF_P1, 3, "analytic_err_grad_p", DofInterpolation);
		analytic_err_grad_n = phgDofNew(g, DOF_P1, 3, "analytic_err_grad_n", DofInterpolation);

		phgDofSetDataByValue(analytic_err_u, 0.0);
		phgDofSetDataByValue(analytic_err_p, 0.0);
		phgDofSetDataByValue(analytic_err_n, 0.0);

		phgDofSetDataByValues(analytic_err_grad_u, zeros);
		phgDofSetDataByValues(analytic_err_grad_p, zeros);
		phgDofSetDataByValues(analytic_err_grad_n, zeros);
	}

}

static void analytic_print()
{	
	FLOAT u_coef = KB * T;
	FLOAT p_coef = Pow(DM2CM, 3) / NA;
	if(analytic_test == 1) {
		int i;
		SIMPLEX *e;
		FLOAT H1ReError[3], L2ReError[3];
		FLOAT sub_H1NormD[3] = {0.0, 0.0, 0.0}, sub_L2NormD[3] = {0.0, 0.0, 0.0};
		FLOAT sub_H1NormA[3] = {0.0, 0.0, 0.0}, sub_L2NormA[3] = {0.0, 0.0, 0.0};
#if USE_MPI 
		FLOAT H1NormD[3] = {0.0, 0.0, 0.0}, L2NormD[3] = {0.0, 0.0, 0.0};
		FLOAT H1NormA[3] = {0.0, 0.0, 0.0}, L2NormA[3] = {0.0, 0.0, 0.0};
#endif
		DOF *grad_u, *grad_p, *grad_n;
		grad_u = phgDofGradient(V, NULL, NULL, NULL);
		grad_p = phgDofGradient(P, NULL, NULL, NULL);
		grad_n = phgDofGradient(N, NULL, NULL, NULL);

      phgPrintf("\n");
		phgPrintf("||analytic_u||_2 = %lf\n", phgDofNormL2(analytic_u) * u_coef);
		phgPrintf("||numeric_u||_2 = %lf\n", phgDofNormL2(V) * u_coef);

		phgPrintf("||analytic_p||_2 = %lf\n", phgDofNormL2(analytic_p) * p_coef);
		phgPrintf("||numeric_p||_2 = %lf\n", phgDofNormL2(P) * p_coef);
		phgPrintf("||analytic_n||_2 = %lf\n", phgDofNormL2(analytic_n) * p_coef);
		phgPrintf("||numeric_n||_2 = %lf\n", phgDofNormL2(N) * p_coef);

		phgDofAXPBY(1.0, analytic_u, 0.0, &analytic_err_u);
		phgDofAXPY(-1.0, V, &analytic_err_u);
		phgDofAXPBY(1.0, analytic_p, 0.0, &analytic_err_p);
		phgDofAXPY(-1.0, P, &analytic_err_p);
		phgDofAXPBY(1.0, analytic_n, 0.0, &analytic_err_n);
		phgDofAXPY(-1.0, N, &analytic_err_n);

		phgDofAXPBY(1.0, analytic_grad_u, 0.0, &analytic_err_grad_u);
		phgDofAXPY(-1.0, grad_u, &analytic_err_grad_u);
		phgDofAXPBY(1.0, analytic_grad_p, 0.0, &analytic_err_grad_p);
		phgDofAXPY(-1.0, grad_p, &analytic_err_grad_p);
		phgDofAXPBY(1.0, analytic_grad_n, 0.0, &analytic_err_grad_n);
		phgDofAXPY(-1.0, grad_n, &analytic_err_grad_n);
	
		ForAllElements(g, e) {
			sub_L2NormD[0] += phgQuadDofDotDof(e, analytic_err_u, analytic_err_u, QUAD_DEFAULT);
			sub_L2NormD[1] += phgQuadDofDotDof(e, analytic_err_p, analytic_err_p, QUAD_DEFAULT);
			sub_L2NormD[2] += phgQuadDofDotDof(e, analytic_err_n, analytic_err_n, QUAD_DEFAULT);

			sub_L2NormA[0] += phgQuadDofDotDof(e, analytic_u, analytic_u, QUAD_DEFAULT);
			sub_L2NormA[1] += phgQuadDofDotDof(e, analytic_p, analytic_p, QUAD_DEFAULT);
			sub_L2NormA[2] += phgQuadDofDotDof(e, analytic_n, analytic_n, QUAD_DEFAULT);

			sub_H1NormD[0] += phgQuadDofDotDof(e, analytic_err_grad_u, analytic_err_grad_u, QUAD_DEFAULT);
			sub_H1NormD[1] += phgQuadDofDotDof(e, analytic_err_grad_p, analytic_err_grad_p, QUAD_DEFAULT);
			sub_H1NormD[2] += phgQuadDofDotDof(e, analytic_err_grad_n, analytic_err_grad_n, QUAD_DEFAULT);

			sub_H1NormA[0] += phgQuadDofDotDof(e, analytic_grad_u, analytic_grad_u, QUAD_DEFAULT);
			sub_H1NormA[1] += phgQuadDofDotDof(e, analytic_grad_p, analytic_grad_p, QUAD_DEFAULT);
			sub_H1NormA[2] += phgQuadDofDotDof(e, analytic_grad_n, analytic_grad_n, QUAD_DEFAULT);
		}

#if USE_MPI
		MPI_Allreduce(&sub_L2NormD, &L2NormD, 3, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		MPI_Allreduce(&sub_L2NormA, &L2NormA, 3, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		MPI_Allreduce(&sub_H1NormD, &H1NormD, 3, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		MPI_Allreduce(&sub_H1NormA, &H1NormA, 3, PHG_MPI_FLOAT, MPI_SUM, g->comm);
		for(i = 0; i < Dim; i++) {
			L2ReError[i] = Sqrt(L2NormD[i]) / Sqrt(L2NormA[i]);
			H1ReError[i] = L2ReError[i] + Sqrt(H1NormD[i]) / Sqrt(H1NormA[i]); 
		}
#else
		for(i = 0; i < Dim; i++) {
			L2ReError[i] = Sqrt(sub_L2NormD[i]) / Sqrt(sub_L2NormA[i]);
			H1ReError[i] = L2ReError[i] + Sqrt(sub_H1NormD[i]) / Sqrt(sub_H1NormA[i]); 
		}
#endif

		phgPrintf("||analytic_err_u||_L2 = %.6le\n", L2ReError[0]);
		phgPrintf("||analytic_err_p||_L2 = %.6le\n", L2ReError[1]);
		phgPrintf("||analytic_err_n||_L2 = %.6le\n", L2ReError[2]);	
		phgPrintf("||analytic_err_u||_H1 = %.6le\n", H1ReError[0]);
		phgPrintf("||analytic_err_p||_H1 = %.6le\n", H1ReError[1]);
		phgPrintf("||analytic_err_n||_H1 = %.6le\n", H1ReError[2]);

		phgDofFree(&grad_u);
		phgDofFree(&grad_p);
		phgDofFree(&grad_n);
	}	

}

static void analytic_finalize() 
{	
	if(analytic_test == 1) {
		phgDofFree(&analytic_u);
		phgDofFree(&analytic_p);
		phgDofFree(&analytic_n);

		phgDofFree(&analytic_f_p);
		phgDofFree(&analytic_f_n);

		phgDofFree(&analytic_grad_u);
		phgDofFree(&analytic_grad_p);
		phgDofFree(&analytic_grad_n);

		phgDofFree(&analytic_err_u);
		phgDofFree(&analytic_err_p);
		phgDofFree(&analytic_err_n);

		phgDofFree(&analytic_err_grad_u);
		phgDofFree(&analytic_err_grad_p);
		phgDofFree(&analytic_err_grad_n);
	}
}
	

#endif
