// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// code for prox step for multi-period portfolio example

#include <stdio.h> 
#include "../../src/osc.h"
#include "prox.h"
#include <stdlib.h>
#include <string.h>
#include "math.h"

#define PI 3.141592654
#define SIG 0.1

double rand_gauss()
{
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;

	return Z;
}

void perturbData(double x_init[], all_data * data){
	for(int i=0;i<data->n;i++){
		int idx = i+(data->n+data->m)*(data->T+1);
		data->RHS[idx]=x_init[i]+SIG*rand_gauss();	
	}
}

void readProxData(FILE * fb,all_data * data, prox_data ** p_data)
{
	*p_data = malloc(sizeof(prox_data)); 
	(*p_data)->kappa=malloc(sizeof(double)*data->n);
	(*p_data)->x_term=malloc(sizeof(double)*data->n);
	for(int i=0;i<data->n;i++){
		fscanf(fb, "%lf", &((*p_data)->kappa[i]));
	}
	for(int i=0;i<data->n;i++){
		fscanf(fb, "%lf", &((*p_data)->x_term[i]));
	}
}

inline double abs_double(double x){
	return (x>0) ? x:-x;
}

inline double pos(double x){
	return (x>0) ? x:0;
}

inline double soft_thresh(double x, double gam){
	return x*pos(1-gam/abs_double(x));
}


void prox(prob_vars * vars, all_data * data, prox_data * p_data)
{
	double v[data->n*(data->T+1)];
	memcpy(v,vars->x,sizeof(double)*data->n*(data->T+1));
	subArray(v,vars->z,data->n*(data->T+1));
	
	double w[data->m*(data->T+1)];
	memcpy(w,vars->u,sizeof(double)*data->m*(data->T+1));
	subArray(w,vars->y,data->m*(data->T+1));
	
	int idx;double temp;
//	#pragma omp parallel for private(idx,temp)
	for(int i=0;i<data->T;i++){
		for(int j=0;j<data->n;j++){
			idx = i*data->n+j;
			temp = soft_thresh(w[idx],p_data->kappa[j]/data->rho);
			if (v[idx] + temp>=0){
				vars->x_t[idx] = v[idx];
				vars->u_t[idx] = temp;
			}
			else{
				vars->u_t[idx]=soft_thresh((w[idx]-v[idx])/2,p_data->kappa[j]/(2*data->rho));
				vars->x_t[idx]=-vars->u_t[idx];
			}
		}
	}
	for(int j=0;j<data->n;j++){
		idx = data->T*data->n+j;
			vars->u_t[idx]=soft_thresh((w[idx]-v[idx]+p_data->x_term[j])/2,p_data->kappa[j]/(2*data->rho));
			vars->x_t[idx]=p_data->x_term[j]-vars->u_t[idx];
		}
}
void freeProxData(all_data * data,prox_data *p_data){
	free(p_data->kappa);
	free(p_data->x_term);
}
