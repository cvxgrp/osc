// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd

// this file contains the code to perform the warm-start simulations

#include <stdio.h>
#include <stdlib.h>
#include <string.h>                                                                   
#include <sys/time.h>
#include <stdbool.h>
#include "math.h"
#include "osc.h"
#include "run_osc.h"
#include "cholesky.h"

#define NUM_WARM 100 // num of warm start simulations to run
#define NUM_COLD 100 // num times to run cold start

void copyVars(prob_vars * new_vars,prob_vars *vars, all_data * data){
	memcpy(new_vars->x_t, vars->x_t,sizeof(double)*data->n*(data->T+1));
    memcpy(new_vars->u_t, vars->u_t,sizeof(double)*data->m*(data->T+1));
    memcpy(new_vars->z, vars->z,sizeof(double)*data->n*(data->T+1));
    memcpy(new_vars->y, vars->y,sizeof(double)*data->m*(data->T+1));
}

int main(int argc, char **argv)
{
	FILE * fp;
	if(openFile(argc, argv, 1, "data/data_KKT", &fp)==-1) return -1;
	all_data * data = malloc(sizeof(all_data));
	readInData(fp,data);
	fclose(fp);

	prob_vars * vars = malloc(sizeof(prob_vars));
	//initVars(vars,data);
	
	FILE * fb;
	if(openFile(argc, argv, 2, "data/data_prox", &fb)==-1) return -1;
	prox_data * p_data;
	readProxData(fb,data,&(p_data));
	fclose(fb);

	Timings cold_start;
	double av_t=0.0, av_lin_s=0.0, av_prox=0.0, worst_cd = 0.0;
	printf("Running %i cold starts\n", NUM_COLD);
	for(int i=0;i<NUM_COLD;i++){	
		initVars(vars,data);
		cold_start = osc(vars,data,p_data);
		av_t+=cold_start.total_time/NUM_COLD;
		av_lin_s+=cold_start.lin_sys_time/NUM_COLD;
		av_prox+=cold_start.prox_time/NUM_COLD;
		worst_cd = fmax(worst_cd,cold_start.total_time);
	}
	printf("COLD START:\n");
	printf("Iterations %i\n",(int)cold_start.itns);
	printf("Taking (on average) %4.2f ms\n",av_t);
	printf("Average time to solve linear system: %4.2f ms\n",av_lin_s);
	printf("Average time to take prox step: %4.4f ms\n",av_prox);
	printf("Worst case total solve time: %4.2f ms \n",worst_cd);

	Timings warm_start[NUM_WARM],average={0};
	double max_time = 0.0;

	printf("RUNNING %i WARM STARTS\n",NUM_WARM);
	prob_vars * new_vars = malloc(sizeof(prob_vars));
	initVars(new_vars,data);
	srand(0);	
	for (int j=0;j<NUM_WARM;j++){
		perturbData(data->x_init,data);
		copyVars(new_vars,vars,data);
		warm_start[j] = osc(new_vars,data,p_data);
		average.prox_time += warm_start[j].prox_time/NUM_WARM;
		average.lin_sys_time += warm_start[j].lin_sys_time/NUM_WARM;
		average.total_time += warm_start[j].total_time/NUM_WARM;
		average.itns += warm_start[j].itns/NUM_WARM;
   	max_time = fmax(warm_start[j].total_time, max_time); 
	}
	printf("COMPLETE\n");
	printf("Average num iterations %4.2f\n",average.itns);
	printf("Taking an average of %4.2f ms \n",average.total_time);
	printf("Average time to solve linear system: %4.2f ms\n",average.lin_sys_time);
	printf("Average time to take prox step: %4.4f ms\n",average.prox_time);
	printf("Worst case time to solve: %4.4f ms\n", max_time);


	freeProxData(data,p_data);
	free(p_data);
	freeVars(vars);
	freeVars(new_vars);
	freeData(data);
}

void initVars(prob_vars * vars,all_data* data){
  vars->x = malloc(sizeof(double)*data->n*(data->T+1));
  vars->u = malloc(sizeof(double)*data->m*(data->T+1));
  vars->x_t = calloc(data->n*(data->T+1),sizeof(double));
  memcpy(vars->x_t,data->x_init,data->n*sizeof(double));
  vars->u_t = calloc(data->m*(data->T+1),sizeof(double));
  vars->z = calloc(data->n*(data->T+1),sizeof(double));
  vars->y = calloc(data->m*(data->T+1),sizeof(double));
}

int openFile(int argc, char ** argv, int idx, char * default_file, FILE ** fb)
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


int readInData(FILE * fp,all_data * data){
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
	choleskyInit(data->nc, Ap, Ar, data->P, &info);
	amd_info(info);

	//choleskyFactor(nc, Ap, Ar, Ax, NULL, NULL, &Lp, &Lr, &Lx, &diag, &LNZ); 
	choleskyFactor(data->nc, Ap, Ar, Ax, data->P, data->Pinv, &data->Lp, &data->Lr, &data->Lx, &data->D, &data->LNZ);
	gettimeofday(&end, NULL); 
	printf("reading in the data took %4.2lf seconds\n",time);
	printf("KKT matrix factorization took %4.2f ms\n",end.tv_sec*1e3 + end.tv_usec/1e3 - start.tv_sec*1e3 - start.tv_usec/1e3);
	free(Ar);free(Ap);free(Ax);free(info);
	return 0;
}
