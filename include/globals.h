#ifndef _GLOBALS_H_
#define _GLOBALS_H_

typedef struct{
    int debug;
    int printtree;
    double eps;
    int counts;
    double min_value;
    double max_value;
} globals;

extern globals GLOBALS;

#endif /* _GLOBALS_H_ */
