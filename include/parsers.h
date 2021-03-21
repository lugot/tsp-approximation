#ifndef _PARSERS_H_
#define _PARSERS_H_

#include "tsp.h"

void parse_command_line(int argc, char** argv, instance inst);
void parse_input_file(instance inst, const char* file_extension);

#endif /* _PARSERS_H_*/
