#ifndef _PARSERS_H_
#define _PARSERS_H_

#include "tsp.h"

typedef struct run_options_t {
	char* model_name;
	int battery_test;
} *run_options;

void parse_command_line(int argc, char** argv, cplex_params params, run_options options);
instance parse_input_file(char* model_name, char* file_extension, enum model_folders folder);

#endif /* _PARSERS_H_*/
