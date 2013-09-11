// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// code to perform prox step for supply chain example

/***********  Supply Chain - Proximal Step ***********/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>	
# include "math.h"
# include "prox.h"
# include "../osc.h"

#define SIG 0.1

void perturb_data(double x_init[], all_data * data){
	for(int i=0;i<data->n;i++){
		int idx = i+(data->n+data->m)*(data->T+1);
		data->RHS[idx] = x_init[i]*(1+SIG*(2*((float)rand())/RAND_MAX-1));
	}
}

double test_lambda(double * wi,double * vi, double * w, double v, prox_data * p_data, int i, double lambda){
	double sumu = 0;
	for ( int k = 1; k <= *(p_data->connections[i]); k++ ){ 
		wi[k-1] = -fminl(-fminl(w[p_data->connections[i][k]]-lambda,p_data->U),0);
		sumu += wi[k-1];
	}
	*vi =  -fminl(-fminl(v+lambda,p_data->C),0);
	return sumu;
}

void bisection(all_data * data,double * w, double* v, prox_data * p_data, int i)
{
	double lambda = 0,lambda_old = 0, vi;
	double * wi = malloc(sizeof(double)*p_data->connections[i][0]);
	
	/******** Compute lower and upper bound (search region) *********/
	double sumu = test_lambda(wi, &vi,w,*v,p_data,i,lambda);

	if ( sumu > vi ) {  /* Check condition */
		lambda = 1;
		sumu = test_lambda(wi, &vi,w,*v,p_data,i,lambda);
		while ( sumu > vi ) {
			lambda_old = lambda;
			lambda *= 2;
			sumu = test_lambda(wi,&vi,w,*v,p_data,i,lambda);
		}
	}
	double low = lambda_old;
	double up = lambda;

	/******** Bisection ********/
	double epsilon = 1e-3;

	while ( (up - low) > epsilon ) {
		lambda = (low + up) / 2;
		sumu = test_lambda(wi,&vi,w,*v,p_data,i,lambda);
		if ( sumu <= vi ) up = lambda; 
		else low = lambda; 
	}
	for ( int k = 1; k <= *(p_data->connections[i]); k++ ){ 
		w[p_data->connections[i][k]] = wi[k-1];
	}
	*v = vi;
	free(wi);
}

void read_in_prox_data(FILE * fb, all_data * data, prox_data * p_data)
{
	// get upper bound on capacity (C) and go to next line
	fscanf(fb,"%lf",&(p_data->C));
	// get upper bound on transportation load (U) and go to next line
	fscanf(fb,"%lf",&(p_data->U));
	// get nodes associated with sources
	fscanf(fb, "%u", &(p_data->numsource));

	p_data->idx_source = (unsigned int *) malloc(sizeof(unsigned int)*(p_data->numsource));
	for ( int i = 0; i < p_data->numsource; i++ ) {
		fscanf(fb, "%u", &(p_data->idx_source[i]));
	}

	// length[i] is the # of edges departing from node i
	unsigned int * length;
	length = (unsigned int *) malloc(sizeof(unsigned int)*data->n);

	p_data->connections = (unsigned int **) malloc(sizeof(unsigned int*)*data->n);
	for ( int i = 0; i < data->n; i++ ) {
		fscanf(fb,"%u", &length[i]);
		p_data->connections[i] = (unsigned int *) malloc(sizeof(unsigned int)*(length[i]+1));
		*(p_data->connections[i]) = length[i];
		for ( int j = 1; j <= length[i]; j++ ) {
			fscanf(fb, "%u", p_data->connections[i]+j);
		}
	}
	free(length);
}

void prox(prob_vars * vars, all_data * data, prox_data * p_data)
{
	memcpy(vars->x_t,vars->x,data->n*(data->T+1)*sizeof(double));
	sub_array(vars->x_t,vars->z,data->n*(data->T+1));

	memcpy(vars->u_t,vars->u,data->m*(data->T+1)*sizeof(double));
	sub_array(vars->u_t,vars->y,data->m*(data->T+1));

	/* The proximal step using bisection (x_t, u_t) */
	int idx;
	#pragma omp parallel for private(idx)
	for ( int t = 0; t < data->T+1; t++ ) {
		for ( int j = 0; j < p_data->numsource; j++ ) {	
			idx = t*data->m + p_data->idx_source[j];
			/* Saturate the inputs */
			vars->u_t[idx] =  -fminl(-fminl(vars->u_t[idx],p_data->U), 0);
		}

		for ( int i = 0; i < data->n; i++ ){  /* For every node... */
			idx = t*data->n+i;
			bisection(data,&(vars->u_t[t*data->m]), &(vars->x_t[idx]), p_data, i);
		}
	}
}


void free_prox_data(all_data * data, prox_data *p_data)
{
	free(p_data->idx_source);
	for ( int i = 0; i < data->n; i++ ) free(p_data->connections[i]);
	free(p_data->connections);
}
