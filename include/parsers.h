#ifndef INCLUDE_PARSERS_H_
#define INCLUDE_PARSERS_H_

#include "../include/tsp.h"

enum run_modes { NOT_SPECIFIED, SINGLE_INSTANCE, GENERATE, LOAD_DIR };

typedef struct run_options_t {
    enum run_modes mode;
    char* instance_name;
    char* instance_folder;
    int battery_test;
    int tests;
    int load_optimal;
} * run_options;

run_options create_options();
void free_options(run_options options);

void parse_command_line(int argc, char** argv, cplex_params params,
                        run_options options);
instance parse_input_file(char* instance_name, char* file_extension,
                          char* instance_folder);
instance* parse_input_dir_bounded(char* instance_folder, char* file_extension,
                                  int* ninstances, int load_optimal,
                                  int nodes_lb, int nodes_ub);
instance* parse_input_dir(char* instance_folder, char* file_extension,
                          int* ninstances, int load_optimal);

#endif  // INCLUDE_PARSERS_H_
