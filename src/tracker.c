#include "../include/tracker.h"

#include <stdio.h>
#include <stdlib.h>

#include "../include/globals.h"

tracker tracker_create() {
    tracker t = (tracker)malloc(sizeof(struct tracker_t));

    t->size = 0;
    t->capacity = 10;
    t->times = (double*)malloc(t->capacity * sizeof(double));
    t->objs = (double*)malloc(t->capacity * sizeof(double));

    return t;
}

void tracker_add(tracker t, double time, double obj) {
    t->size++;

    if (t->size >= t->capacity) {
        t->capacity = 2 * t->capacity;
        t->times = (double*)realloc(t->times, t->capacity * sizeof(double));
        t->objs = (double*)realloc(t->objs, t->capacity * sizeof(double));
        /* TODO(lugot): SAFETY realloc is not safe! */
    }
    t->times[t->size - 1] = time;
    t->objs[t->size - 1] = obj;
}

double tracker_find(tracker t, double obj) {
    for (int i = 0; i < t->size; i++) {
        if (t->objs[i] < obj) return t->times[i];
    }

    return -1.0;
}

void tracker_print(tracker t) {
    for (int i = 0; i < t->size; i++) {
        printf("%3.2lf %lf\n", t->times[i], t->objs[i]);
    }
}

void tracker_free(tracker t) {
    free(t->times);
    free(t->objs);

    free(t);
}
