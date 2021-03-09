#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct 
{
	//input
	int nnodes;
	double *xcoord;
	double *ycoord;
	
	//params
	int model_type;
	double time_limit;
	char input_file[1000];
	int available_memory;

} instance;