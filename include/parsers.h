#ifndef _PARSERS_H_
#define _PARSERS_H_

#include "tsp.h"

void parse_command_line(int argc, char** argv, instance inst);
void parse_input(instance inst);
void plot_instance_solution();

#endif /* _PARSERS_H_*/
