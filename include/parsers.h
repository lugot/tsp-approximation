#ifndef INCLUDE_PARSERS_H_
#define INCLUDE_PARSERS_H_

#include "../include/tsp.h"

enum run_modes {
    NOT_SPECIFIED,
    SINGLE_MODEL,
    GENERATE,
    LOAD_DIR
};

typedef struct run_options_t {
    enum run_modes mode;
    char* model_name;
    int battery_test;
    enum model_folders folder;
} * run_options;

run_options create_options();
void free_options(run_options options);

void parse_command_line(int argc, char** argv, cplex_params params,
                        run_options options);
instance parse_input_file(char* model_name, char* file_extension,
                          enum model_folders folder);
instance* parse_input_dir(enum model_folders folder, char* file_extension,
                          int* ninstances, int nodes_lb, int nodes_ub);

#endif  // INCLUDE_PARSERS_H_
