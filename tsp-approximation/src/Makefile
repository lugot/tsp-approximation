OBJS = main.o globals.o tsp.o parsers.o utils.o solvers.o union_find.o model_builder.o
HEADERS =
EXE = tsp_approx
all: $(EXE)
setting = -1
OS := $(shell uname)

setting = 1
CPLEX_HOME = /opt/ibm/ILOG/CPLEX_Studio201/cplex
CONCORDELIB = /opt/concorde
CCDIR = /opt/concorde
CC = gcc
AR = ar rc
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl ${CONCORDELIB}/libconcorde.a
INC = -I../include -I${CPLEX_HOME}/include/ilcplex -I${CCDIR}


# ---------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------
CFLAGS = -Wall -O0 -g
RM = rm -rf

.SUFFIXES:
.SUFFIXES: .o .c .cpp
.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
.cpp.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(EXE): $(OBJS) $(LIBUTILS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(OBJS) : $(HEADERS)

$(LIBUTILS): $(OBJS_LIBUTILS)
	$(AR) $(LIBUTILS) $(OBJS_LIBUTILS)

$(LIBUTILS) : $(HEADERS_LIBUTILS)

clean:
	$(RM) $(OBJS)
	$(RM) $(OBJS_LIBUTILS)
	$(RM) $(LIBUTILS)
	#$(RM) $(EXE)

again:
	make clean
	make

wow:
	@echo "                                      W O W W W W WWWWWWWWWWWWWWWWWWW"

who:
	@echo "you are user $(USER) with uname `uname` (OS = $(OS)) and you working with compiler setting $(setting)"

