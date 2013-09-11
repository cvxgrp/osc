// A Splitting Method for Optimal Control  
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd 

// this file contains the code to perform a single cold-start

#include <stdio.h>
#include <stdlib.h>
#include <string.h>                                                                   
#include <sys/time.h>
#include <stdbool.h>
#include "math.h"
#include "../cholesky.h"
#include "../osc.h"
#include "prox.h"

#define MAX_ITERS 3000 // maximum number of iterations OSC will perform
#define EPS_ABS 0.001 // absolute tolerance (for testing convergence)
#define EPS_REL 0.001 // relative tolerance (for testing convergence)

// functions explained in osc.h

int main(int argc, char **argv)
{
	FILE * fp;
	if(open_file(argc, argv, 1, "data_KKT", &fp)==-1) return -1;
	all_data * data = malloc(sizeof(all_data));
	read_in_data(fp,data);
	fclose(fp);

	prob_vars * vars = malloc(sizeof(prob_vars));                                     
	init_vars(vars,data);
	
	FILE * fb;
	if(open_file(argc, argv, 2, "data_prox", &fb)==-1) return -1;
	prox_data * p_data = malloc(sizeof(prox_data));
	read_in_prox_data(fb,data,p_data);
	fclose(fb);

	bool relax_on = fabs(data->alpha-1) < 1e-3 ? false:true;
	printf("relaxation is %s\n",(relax_on)?"on":"off");

	struct timeval start,end;
	gettimeofday(&start,NULL);

	struct timeval lin_sys_st, lin_sys_end, prox_st,prox_end;
	double lin_sys_time=0, prox_time=0;

	double x_t_old[data->n*(data->T+1)];
	double u_t_old[data->m*(data->T+1)];
	int i;double rn, dn;
	for (i=0;i<MAX_ITERS;i++){
		//printf("iteration %i\n",i);
		memcpy(x_t_old,vars->x_t,sizeof(double)*data->n*(data->T+1));
		memcpy(u_t_old,vars->u_t,sizeof(double)*data->m*(data->T+1));
		
		gettimeofday(&lin_sys_st, NULL);  
		solve_lin_sys(vars,data);
		gettimeofday(&lin_sys_end,NULL);
		lin_sys_time += lin_sys_end.tv_sec*1e3 + lin_sys_end.tv_usec/1e3 - lin_sys_st.tv_sec*1e3 - lin_sys_st.tv_usec/1e3;
		
		//double x_hat[data->n*(data->T+1)];
		//memcpy(x_hat,vars->x,sizeof(double)*data->n*(data->T+1));
		//double u_hat[data->m*(data->T+1)];
		//memcpy(u_hat,vars->u,sizeof(double)*data->m*(data->T+1));

		if (relax_on) relax(vars,data);
		
		gettimeofday(&prox_st, NULL);  
		prox(vars,data,p_data);
		gettimeofday(&prox_end,NULL);
		prox_time += prox_end.tv_sec*1e3 + prox_end.tv_usec/1e3 - prox_st.tv_sec*1e3 - prox_st.tv_usec/1e3;
		
		update_dual_vars(vars,data);
		//print_all(data, vars);
		//memcpy(vars->x,x_hat,sizeof(double)*data->n*(data->T+1));
		//memcpy(vars->u,u_hat,sizeof(double)*data->m*(data->T+1));

		if(test_convergence(data,vars,x_t_old,u_t_old,&rn,&dn)) break;
	}
	gettimeofday(&end, NULL);  
	if (i==MAX_ITERS) printf("WARNING: OSC did not converge within %i iterations\n",i);
	else printf("OSC complete in: %4.4f seconds, %i iterations\n",end.tv_sec + end.tv_usec/1e6 - start.tv_sec - start.tv_usec/1e6,i+1);
	printf("primal residual norm is %1.4f, dual residual norm is %1.4f\n",rn,dn);
	printf("solving the linear system took an average of %4.2f ms\n",lin_sys_time/(i+1));
	printf("performing the prox step took an average of %4.4f ms\n",prox_time/(i+1));

	int n_cols = data->T+1<5 ? data->T+1 : 5;

	printf("\nfirst %i cols of u: \n",n_cols);
	for(int i=0;i<data->m;i++){
		for (int j=0;j<n_cols;j++){
			printf("%4.3f	",vars->u_t[j*data->m+i]);
		}
		printf("\n");
	}
	printf("\nfirst %i cols of x: \n",n_cols);
	for(int i=0;i<data->n;i++){
		for (int j=0;j<n_cols;j++){
			printf("%4.3f	",vars->x_t[j*data->n+i]);
		}
		printf("\n");
	}
	//print_all(data,vars);
	free_prox_data(data,p_data);
	free(p_data);
	free_data(data);
	free_vars(vars);
}

void print_all(all_data * data, prob_vars * vars){
	printf("\n u is \n");
	for(int i=0;i<data->m*(data->T+1);i++){
		printf("%f\n",vars->u[i]);
	}
	printf("\n x is \n");
	for(int i=0;i<data->n*(data->T+1);i++){
		printf("%f\n",vars->x[i]);
	}
	printf("\n u_t is \n");
	for(int i=0;i<data->m*(data->T+1);i++){
		printf("%f\n",vars->u_t[i]);
	}
	printf("\n x_t is \n");
	for(int i=0;i<data->n*(data->T+1);i++){
		printf("%f\n",vars->x_t[i]);
	}
	printf("\n y is \n");
	for(int i=0;i<data->m*(data->T+1);i++){
		printf("%f\n",vars->y[i]);
	}
	printf("\n z is \n");
	for(int i=0;i<data->n*(data->T+1);i++){
		printf("%f\n",vars->z[i]);
	}
}

void free_vars(prob_vars * vars){
	free(vars->x);free(vars->u);free(vars->x_t);free(vars->u_t);free(vars->z);free(vars->y);
	free(vars);
}

void free_data(all_data* data){
	free(data->Lp);free(data->Lr);free(data->Lx);free(data->D);free(data->RHS);free(data->P);free(data->Pinv);free(data->x_init);
	free(data);
}


void relax(prob_vars * vars, all_data* data){
//	#pragma omp parallel for
	for(int j=0;j<data->n*(data->T+1);j++){
		vars->x[j] = data->alpha*vars->x[j]+(1-data->alpha)*vars->x_t[j];
	}
//	#pragma omp parallel for
	for(int j=0;j<data->m*(data->T+1);j++){
		vars->u[j] = data->alpha*vars->u[j]+(1-data->alpha)*vars->u_t[j];
	}
}

int open_file(int argc, char ** argv, int idx, char * default_file, FILE ** fb) 
{
	if (argc<idx+1){
		printf("Not enough arguments supplied, using %s as default\n", default_file);
	}
	else{
		*fb = fopen(argv[idx], "r");
		if (*fb != NULL) return 0;
		else{
			printf("Couldn't open file %s, using %s as default\n", argv[idx],default_file);
			fclose(*fb);
		}
	}
	*fb = fopen(default_file, "r");
	if (*fb == NULL){
		printf("Couldn't open %s\n",default_file);
		return -1;
	}
	return 0;
}



int read_in_data(FILE * fp,all_data * data){
	struct timeval start, end;
	gettimeofday(&start,NULL);
	//MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT
	size_t NNZ;
	fscanf(fp, "%zu", &(data->n));
	fscanf(fp, "%zu", &(data->m));
	fscanf(fp, "%zu", &(data->T));
	fscanf(fp, "%zu", &NNZ);
	fscanf(fp, "%lf", &(data->rho));
	fscanf(fp, "%lf", &(data->alpha));
	data->nc = (2*data->n+data->m)*(data->T+1);
	data->x_init = malloc(sizeof(double)*data->n);
	for(int i = 0; i < data->n; i++)
	{
		fscanf(fp, "%lf", &(data->x_init[i]));
	}
	int *Ar = malloc(sizeof(int)*NNZ);
	for(int i = 0; i < NNZ; i++)
	{
		fscanf(fp, "%i", &Ar[i]);
	}
	int *Ap=malloc(sizeof(int)*(data->nc+1));
	for(int i = 0; i < data->nc+1; i++)
	{
		fscanf(fp, "%i", &Ap[i]);
	}
	double *Ax=malloc(sizeof(double)*NNZ);
	for(int i = 0; i < NNZ; i++)
	{
		fscanf(fp, "%lf", &Ax[i]);
	}
	data->RHS = malloc(sizeof(double)*data->nc);
	for(int i = 0; i < data->nc; i++)
	{   
		fscanf(fp, "%lf", &(data->RHS[i]));
	}

	gettimeofday(&end, NULL);  
	double time = end.tv_sec + end.tv_usec/1e6 - start.tv_sec - start.tv_usec/1e6;
	printf("KKT matrix factorization info:\n");
	gettimeofday(&start, NULL);  
	data->P = malloc(sizeof(int)*data->nc);
	data->Pinv = malloc(sizeof(int)*data->nc);
	double *info;
	cholesky_init(data->nc, Ap, Ar, data->P, &info);
	amd_info(info);
	
	//cholesky_factor(nc, Ap, Ar, Ax, NULL, NULL, &Lp, &Lr, &Lx, &diag, &LNZ); 
	cholesky_factor(data->nc, Ap, Ar, Ax, data->P, data->Pinv, &data->Lp, &data->Lr, &data->Lx, &data->D, &data->LNZ);
	gettimeofday(&end, NULL); 
	printf("reading in the data took %4.2lf seconds\n",time);
	printf("KKT matrix factorization took %4.2f ms\n",end.tv_sec*1e3 + end.tv_usec/1e3 - start.tv_sec*1e3 - start.tv_usec/1e3);
	free(Ar);free(Ap);free(Ax);free(info);
	return 0;
}

void init_vars(prob_vars * vars,all_data* data){
	vars->x = malloc(sizeof(double)*data->n*(data->T+1));
	vars->u = malloc(sizeof(double)*data->m*(data->T+1));
	vars->x_t = calloc(data->n*(data->T+1),sizeof(double));
	memcpy(vars->x_t,data->x_init,data->n*sizeof(double));
	vars->u_t = calloc(data->m*(data->T+1),sizeof(double));
	vars->z = calloc(data->n*(data->T+1),sizeof(double));
	vars->y = calloc(data->m*(data->T+1),sizeof(double));
}

double calc_norm_squared(double * A,size_t len){
	double norm2=0.0;
//	#pragma omp parallel for reduction(+: norm2)
	for(int i=0;i<len;i++){
		norm2 += pow(A[i],2);	
		}
	return norm2;
}


void scale_array(double * a,double b,int len){
//	#pragma omp parallel for
	for(int i=0;i<len;i++)
		a[i]*=b;
}

void add_array(double * a,double *b,int len){
//	#pragma omp parallel for
	for(int i=0;i<len;i++)
		a[i]+=b[i];
}

void sub_array(double * a,double *b,int len){
//	#pragma omp parallel for
	for(int i=0;i<len;i++)
		a[i]-=b[i];
}


bool test_convergence(all_data * data,prob_vars * vars,double *x_t_old,double *u_t_old,double*rn,double*dn){

	double u_tmp[data->m*(data->T+1)];
	memcpy(u_tmp,vars->u_t,sizeof(double)*data->m*(data->T+1));
	sub_array(u_tmp,vars->u,data->m*(data->T+1));
	double x_tmp[data->n*(data->T+1)];
	memcpy(x_tmp,vars->x_t,sizeof(double)*data->n*(data->T+1));
	sub_array(x_tmp,vars->x,data->n*(data->T+1));

	*rn = sqrt(calc_norm_squared(u_tmp,data->m*(data->T+1))+calc_norm_squared(x_tmp,data->n*(data->T+1)));
	double eps_primal = EPS_ABS*sqrt((data->T+1)*(data->n+data->m)) + EPS_REL*sqrt(fmaxl(calc_norm_squared(vars->x,data->n*(data->T+1))+calc_norm_squared(vars->u,data->m*(data->T+1)),calc_norm_squared(vars->x_t,data->n*(data->T+1))+ calc_norm_squared(vars->u_t,data->m*(data->T+1))));

	if(*rn>eps_primal){
		//printf("rn is %f\n",*rn);
		return false;
	}
	
	memcpy(u_tmp,vars->u_t,sizeof(double)*data->m*(data->T+1));
	sub_array(u_tmp,u_t_old,data->m*(data->T+1));
	memcpy(x_tmp,vars->x_t,sizeof(double)*data->n*(data->T+1));
	sub_array(x_tmp,x_t_old,data->n*(data->T+1));

	*dn = sqrt(calc_norm_squared(u_tmp,data->m*(data->T+1))+calc_norm_squared(x_tmp,data->n*(data->T+1)));
	(*dn)*=data->rho;
	double eps_dual = EPS_ABS*sqrt((data->T+1)*(data->n+data->m)) + EPS_REL*sqrt(calc_norm_squared(vars->z,data->n*(data->T+1))+ calc_norm_squared(vars->y,data->m*(data->T+1)));

	if(*dn>eps_dual){
		//printf("dn is %f\n",*dn);
		return false;
	}
	return true;
}

void update_dual_vars(prob_vars * vars,all_data *data){
	add_array(vars->z,vars->x_t,data->n*(data->T+1));
	sub_array(vars->z,vars->x,data->n*(data->T+1));
	add_array(vars->y,vars->u_t,data->m*(data->T+1));
	sub_array(vars->y,vars->u,data->m*(data->T+1));
}

void solve_lin_sys(prob_vars * vars,all_data *data){
	double rhs[data->nc];
	memcpy(rhs,data->RHS,sizeof(double)*data->nc);

//	#pragma omp parallel for
	for(int i=0;i<data->T+1;i++){
		for(int j=0;j<data->n;j++){
			rhs[i*(data->n+data->m)+j]+=data->rho*(vars->x_t[i*data->n+j]+vars->z[i*data->n+j]);
		}
		for(int j=0;j<data->m;j++){
			rhs[i*(data->n+data->m)+j+data->n]+=data->rho*(vars->u_t[i*data->m+j]+vars->y[i*data->m+j]);
		}
	}
	double *w = malloc(sizeof(double)*data->nc);
	cholesky_solve(data->nc, w, rhs, data->Lp, data->Lr, data->Lx, data->D, data->P);
//	#pragma omp parallel for
	for(int i=0;i<data->T+1;i++){
		for(int j=0;j<data->n;j++){
			vars->x[i*(data->n)+j]=w[i*(data->n+data->m)+j];
		}
		for(int j=0;j<data->m;j++){
			vars->u[i*(data->m)+j]=w[i*(data->n+data->m)+data->n+j];
		}
	}
	free(w);
}
