#ifndef _PARSERS_H_
#define _PARSERS_H_

#include "tsp.h"

void parse_command_line(int argc, char** argv, instance inst);
void parse_input_file(instance inst);
solution parse_optimal_tour(instance inst);

#endif /* _PARSERS_H_*/
