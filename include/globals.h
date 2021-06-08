#ifndef INCLUDE_GLOBALS_H_
#define INCLUDE_GLOBALS_H_

#define XSMALL 1e5
#define EPSILON 1e-9
#define INF 1e9

#define HF_PERCENTAGE 80
#define HF_ITERATIONS 20
#define HF_INITIAL_PERC_TIME 0.1

#define SF_INITIAL_K 3
#define SF_MAX_K 15
#define SF_K_STEP 2
#define SF_INITIAL_PERC_TIME 0.2

#define GRASP_K 4
#define GRASP_VNS_PERC_TIME 0.3

#define VNS_K_START 5
#define VNS_K_MAX 20
#define VNS_K_STEP 1

#define TS_MAX_TENURE 100
#define TS_MIN_TENURE 20

extern int VERBOSE;
extern int EXTRA;
extern int SUPPRESS_CALLBACK;
extern int GEN_NNODES;

#endif  // INCLUDE_GLOBALS_H_
