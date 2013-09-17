#ifndef RUN_OSC_H_GUARD                                                              
#define RUN_OSC_H_GUARD

// initialize the variables to zero (with correct memory allocation)
void initVars(prob_vars * vars,all_data* data);

// simple function to open supplied file
int openFile(int argc, char ** argv, int idx, char * default_file, FILE ** fb);

// read in the data from supplied data file into all_data struct
int readInData(FILE * fp,all_data * data);

void readProxData(FILE * fb,all_data * data,prox_data ** p_data);
void freeProxData(all_data * data,prox_data *p_data);
void perturbData(double x_init[], all_data * data);
#endif
