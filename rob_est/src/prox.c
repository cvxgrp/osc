// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// code to perform prox step for robust state estimation example 

#include <stdio.h> 
#include "../../src/osc.h"
#include "prox.h"
#include <stdlib.h>
#include <string.h>
#include "math.h"

#define SIG 0.1

void perturbData(double x_init[], all_data * data){
	for(int i=0;i<data->n;i++){
		int idx = i+(data->n+data->m)*(data->T+1);
		data->RHS[idx] = x_init[i]*(1+SIG*(2*((float)rand())/RAND_MAX-1));
	}
}

void readProxData(FILE * fb,all_data * data, prox_data ** p_data)
{
	*p_data = malloc(sizeof(prox_data)); 
	fscanf(fb, "%lf", &((*p_data)->M));
}

double nrm(double * x, int len){
	double nm=0;
	for(int i=0;i<len;i++) nm+=pow(x[i],2);
	return sqrt(nm);
}

void prox(prob_vars * vars, all_data * data, prox_data * p_data)
{
	double v[data->m*(data->T+1)];
	memcpy(v,vars->u,sizeof(double)*data->m*(data->T+1));
	subArray(v,vars->y,data->m*(data->T+1));

	memcpy(vars->x_t,vars->x,sizeof(double)*data->n*(data->T+1));
	subArray(vars->x_t,vars->z,data->n*(data->T+1));

	double nm,fac;
//	#pragma omp parallel for private(nm,fac)
	for(int i=0;i<data->T+1;i++){
		nm = nrm(&v[i*data->m],data->m);
		fac = 1-fminl(1/(1+data->rho),p_data->M/(data->rho*nm));
		for(int j=0;j<data->m;j++){
			vars->u_t[i*data->m+j]=fac*v[i*data->m+j];
		}
	}
}
void freeProxData(all_data * data,prox_data *p_data){
}
