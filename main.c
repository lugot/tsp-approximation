#include "tsp.h"

void readinput(instance *inst);
void parse_command_line(int argc, char **argv, instance *inst);
void print_error(char *message);

int debug = 1;

int main(int argc, char **argv)
{

	if(debug)
		printf("Inizio\n");

	instance inst;

	parse_command_line(argc, argv, &inst);

	if(debug)
		printf("Mezzo\n");

	readinput(&inst);

	if(debug)
		printf("Fine\n");

	
	return 0;
}

void parse_command_line(int argc, char **argv, instance *inst)
{
	//default values
	inst->model_type = 0;
	inst->time_limit = 10000000;
	inst->available_memory = 1200;

	//read line
	for(int i=1; i<argc; i++)
	{
		if(strcmp(argv[i], "-f") == 0) {strcpy(inst->input_file, argv[++i]); continue;}
		if(strcmp(argv[i], "-file") == 0) {strcpy(inst->input_file, argv[++i]); continue;}
		if(strcmp(argv[i], "-filename") == 0) {strcpy(inst->input_file, argv[++i]); continue;}
		if(strcmp(argv[i], "-input_file") == 0) {strcpy(inst->input_file, argv[++i]); continue;}
		if(strcmp(argv[i], "-timelimit") == 0) {inst->time_limit = atof(argv[++i]); continue;}
		if(strcmp(argv[i], "-model_type") == 0) {inst->model_type = atoi(argv[++i]); continue;}
		if(strcmp(argv[i], "-model") == 0) {inst->model_type = atoi(argv[++i]); continue;}
		if(strcmp(argv[i], "-memory") == 0) {inst->available_memory = atoi(argv[++i]); continue;}
		if(strcmp(argv[i], "-available_memory") == 0) {inst->available_memory = atoi(argv[++i]); continue;}

	}


	if(1)
	{
		printf("Input file: %s\n", inst->input_file);
		printf("Model type: %d\n", inst->model_type);
		printf("Available memory: %d\n", inst->available_memory);
	}
	


}


void readinput(instance *inst)
{
	if(debug)
		printf("Apertura file\n");

	FILE *fin = fopen(inst->input_file, "r");
	if (fin == NULL) printf("File not found");
	
	if(debug)
		printf("File aperto\n");

	inst->nnodes = -1;

	char line[200];
	char *param;
	char *tok1;
	char *tok2;
	int section = 0;

	while(fgets(line, sizeof(line), fin) != NULL)
	{

		//if(debug)
		//	printf("Cycles\n");

		if (strlen(line) <= 1) continue;

		param = strtok(line, " :");

		if(strcmp(param, "NAME") == 0)
		{
			continue;
		}

		if(strcmp(param, "COMMENT") == 0)
		{
			tok1 = strtok(NULL, "");
		}

		if(strcmp(param, "TYPE") == 0)
		{
			tok1 = strtok(NULL, " : ");

			if(strncmp(tok1, "TSP", 3) != 0) print_error("Only handles TSP\n");
		}

		if(strcmp(param, "DIMENSION") == 0)
		{
			tok1 = strtok(NULL, " :");
			inst->nnodes = atoi(tok1);

			inst -> xcoord = (double *) calloc(inst->nnodes, sizeof(double));
			inst -> ycoord = (double *) calloc(inst->nnodes, sizeof(double));
			continue;
		}

		if(strcmp(param, "EDGE_WEIGHT_TYPE") == 0)
		{
			tok1 = strtok(NULL, " : ");
			
			if(strncmp(tok1, "ATT", 3) != 0) print_error("Only handles ATT");

		}

		if(strcmp(param, "NODE_COORD_SECTION") == 0)
		{
			tok1 = strtok(NULL, " :");
			
			if(inst->nnodes < 0) print_error("Define number of nodes before list of nodes");

			section++;
		}

		if(strcmp(param, "EOF") == 0)
		{
			tok1 = strtok(NULL, " :");
			break;
		}

		if(section == 1)
		{
			int i = atoi(param) - 1;
			if(i<0 || i>= inst->nnodes) print_error("Error in node indexing");
			tok1 = strtok(NULL, " :");
			tok2 = strtok(NULL, " :");
			inst->xcoord[i] = atof(tok1);
			inst->ycoord[i] = atof(tok2);
		}

	}

	fclose(fin);

}

void print_error(char *message)
{
	printf("%s\n", *message);
	exit(1);
}
