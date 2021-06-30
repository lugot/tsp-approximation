#ifndef INCLUDE_TRACKER_H_
#define INCLUDE_TRACKER_H_

typedef struct tracker_t {
    double* times;
    double* objs;

    int size;
    int capacity;
} * tracker;


tracker tracker_create();
void tracker_add(tracker t, double time, double obj);
double tracker_find(tracker t, double obj);
void tracker_free(tracker t);

#endif  // INCLUDE_TRACKER_H_

