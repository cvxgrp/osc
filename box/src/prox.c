// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// C file to perform prox step for box constrained example

#include <stdio.h> 
#include "prox.h"
#include <stdlib.h>
#include "../../src/osc.h"
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
	*(p_data) = malloc(sizeof(prox_data));
	fscanf(fb, "%lf", &((*p_data)->umax));
	fscanf(fb, "%lf", &((*p_data)->umin));
}

void prox(prob_vars * vars, all_data * data, prox_data * p_data)
{
	memcpy(vars->x_t,vars->x,sizeof(double)*data->n*(data->T+1));
	subArray(vars->x_t,vars->z,data->n*(data->T+1));
	memcpy(vars->u_t,vars->u,sizeof(double)*data->m*(data->T+1));
	subArray(vars->u_t,vars->y,data->m*(data->T+1));
//	#pragma omp parallel for
	for(int i=0;i<data->m*(data->T+1);i++){
		if(vars->u_t[i]>p_data->umax) vars->u_t[i]=p_data->umax;
		else if(vars->u_t[i]<p_data->umin) vars->u_t[i]=p_data->umin;
	}
}
void freeProxData(all_data * data,prox_data *p_data){
}
