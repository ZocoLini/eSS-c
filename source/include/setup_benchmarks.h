#ifndef SETUP_H
#define SETUP_H
#include <AMIGO_problem.h>

extern double *min_dom;
extern double *max_dom;
extern int id_problem;
extern int dim;
extern AMIGO_problem *amigo;
extern int use_amigo;
extern int estime_init_cond;


void setup_benchmark(int );
double evaluate(double *);


#endif
