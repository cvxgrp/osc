// A Splitting Method for Optimal Control  
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd 

// this file contains the code to perform a single cold-start

#include <stdio.h>
#include <stdlib.h>
#include <string.h>                                                                   
#include <sys/time.h>
#include <stdbool.h>
#include "math.h"
#include "cholesky.h"
#include "osc.h"

#define MAX_ITERS 3000 // maximum number of iterations OSC will perform
#define EPS_ABS 0.001 // absolute tolerance (for testing convergence)
#define EPS_REL 0.001 // relative tolerance (for testing convergence)

// performs OSC algorithm and returns timing data
Timings osc(prob_vars * vars, all_data * data, prox_data * p_data){
	bool relax_on = fabs(data->alpha-1) < 1e-3 ? false:true;
	Timings tt={0};
	struct timeval start,end,lin_sys_st,lin_sys_end, prox_st,prox_end;

	double x_t_old[data->n*(data->T+1)];
	double u_t_old[data->m*(data->T+1)];
	int i;double rn, dn;

	gettimeofday(&start, NULL);  
	for (i=0;i<MAX_ITERS;i++){
		//printf("iteration %i\n",i);
		memcpy(x_t_old,vars->x_t,sizeof(double)*data->n*(data->T+1));
		memcpy(u_t_old,vars->u_t,sizeof(double)*data->m*(data->T+1));

		gettimeofday(&lin_sys_st, NULL);  
		solveLinSys(vars,data);
		gettimeofday(&lin_sys_end,NULL);
		tt.lin_sys_time += lin_sys_end.tv_sec*1e3 + lin_sys_end.tv_usec/1e3 - lin_sys_st.tv_sec*1e3 - lin_sys_st.tv_usec/1e3;
		if (relax_on) relax(vars,data);

		gettimeofday(&prox_st, NULL);  
		prox(vars,data,p_data);
		gettimeofday(&prox_end,NULL);
		tt.prox_time += prox_end.tv_sec*1e3 + prox_end.tv_usec/1e3 - prox_st.tv_sec*1e3 - prox_st.tv_usec/1e3;

		updateDualVars(vars,data);
		//printAll(data, vars);
		if(testConvergence(data,vars,x_t_old,u_t_old,&rn,&dn)) break;
	}
	gettimeofday(&end,NULL);
	tt.lin_sys_time /= (i+1);
	tt.prox_time /= (i+1);
	tt.total_time = end.tv_sec*1e3 + end.tv_usec/1e3 - start.tv_sec*1e3 - start.tv_usec/1e3;
	tt.itns = i+1;
	return tt;
}

void printAll(all_data * data, prob_vars * vars){
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

void freeVars(prob_vars * vars){
	free(vars->x);free(vars->u);free(vars->x_t);free(vars->u_t);free(vars->z);free(vars->y);
	free(vars);
}

void freeData(all_data* data){
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

double calcNormSquared(double * A,size_t len){
	double norm2=0.0;
//	#pragma omp parallel for reduction(+: norm2)
	for(int i=0;i<len;i++){
		norm2 += pow(A[i],2);	
		}
	return norm2;
}

void scaleArray(double * a,double b,int len){
//	#pragma omp parallel for
	for(int i=0;i<len;i++)
		a[i]*=b;
}

void addArray(double * a,double *b,int len){
//	#pragma omp parallel for
	for(int i=0;i<len;i++)
		a[i]+=b[i];
}

void subArray(double * a,double *b,int len){
//	#pragma omp parallel for
	for(int i=0;i<len;i++)
		a[i]-=b[i];
}

bool testConvergence(all_data * data,prob_vars * vars,double *x_t_old,double *u_t_old,double*rn,double*dn){

	double u_tmp[data->m*(data->T+1)];
	memcpy(u_tmp,vars->u_t,sizeof(double)*data->m*(data->T+1));
	subArray(u_tmp,vars->u,data->m*(data->T+1));
	double x_tmp[data->n*(data->T+1)];
	memcpy(x_tmp,vars->x_t,sizeof(double)*data->n*(data->T+1));
	subArray(x_tmp,vars->x,data->n*(data->T+1));

	*rn = sqrt(calcNormSquared(u_tmp,data->m*(data->T+1))+calcNormSquared(x_tmp,data->n*(data->T+1)));
	double eps_primal = EPS_ABS*sqrt((data->T+1)*(data->n+data->m)) + EPS_REL*sqrt(fmaxl(calcNormSquared(vars->x,data->n*(data->T+1))+calcNormSquared(vars->u,data->m*(data->T+1)),calcNormSquared(vars->x_t,data->n*(data->T+1))+ calcNormSquared(vars->u_t,data->m*(data->T+1))));

	if(*rn>eps_primal){
		//printf("rn is %f\n",*rn);
		return false;
	}
	
	memcpy(u_tmp,vars->u_t,sizeof(double)*data->m*(data->T+1));
	subArray(u_tmp,u_t_old,data->m*(data->T+1));
	memcpy(x_tmp,vars->x_t,sizeof(double)*data->n*(data->T+1));
	subArray(x_tmp,x_t_old,data->n*(data->T+1));

	*dn = sqrt(calcNormSquared(u_tmp,data->m*(data->T+1))+calcNormSquared(x_tmp,data->n*(data->T+1)));
	(*dn)*=data->rho;
	double eps_dual = EPS_ABS*sqrt((data->T+1)*(data->n+data->m)) + EPS_REL*sqrt(calcNormSquared(vars->z,data->n*(data->T+1))+ calcNormSquared(vars->y,data->m*(data->T+1)));

	if(*dn>eps_dual){
		//printf("dn is %f\n",*dn);
		return false;
	}
	return true;
}

void updateDualVars(prob_vars * vars,all_data *data){
	addArray(vars->z,vars->x_t,data->n*(data->T+1));
	subArray(vars->z,vars->x,data->n*(data->T+1));
	addArray(vars->y,vars->u_t,data->m*(data->T+1));
	subArray(vars->y,vars->u,data->m*(data->T+1));
}

void solveLinSys(prob_vars * vars,all_data *data){
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
	choleskySolve(data->nc, w, rhs, data->Lp, data->Lr, data->Lx, data->D, data->P);
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
