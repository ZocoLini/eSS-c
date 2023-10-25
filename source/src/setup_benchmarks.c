/**
 * @file setup_benchmarks.c
 * @author David R. Penas
 * @brief File containing functions about the load for different benchmarks.
 * When users typical want to put their benchmarks, they should modify this
 * file.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <benchmark_functions_SystemBiology.h>
#include <AMIGO_problem.h>
/**
 * @brief this the setup function where the usar.
 * @param exp experiment_total struct, the main struct of program.
 * @param idp parallel identification number.
 * @param first this variable is 1, when it is the first time that the program
 * enters in this function.
 * @return void function pointer of evaluation function of the selected.
 * benchmark.
 */

double *min_dom;
double *max_dom;
int id_problem=0;
int dim=0;
AMIGO_problem *amigo;
int use_amigo=0;
int estime_init_cond=0;


void setup_benchmark(int id) {
    const char *namealg, *custom, *noiseBBOB, *noiselessBBOB, *systemsBiology, *matlabproblem, *pythonNLP;
    char *name;
    
    // CONTROL ID 
    
    /// select a specific system biology problem:
    /// current_bench = 0 --> CIRCADIAN PROBLEM
    /// current_bench = 1 --> step PATHWAY PROBLEM
    /// current_bench = 2 --> NFKB PROBLEM
    /// current_bench = 3 --> B1 BIOPREDYN PROBLEM
    /// current_bench = 4 --> B2 BIOPREDYN PROBLEM
    /// current_bench = 5 --> B3 BIOPREDYN PROBLEM
    /// current_bench = 6 --> B4 BIOPREDYN PROBLEM
    /// current_bench = 7 --> B5 BIOPREDYN PROBLEM
    /// current_bench = 8 --> B6 BIOPREDYN PROBLEM
    load_benchmark_SystemBiology(id);
    id_problem = id;
    
}


double evaluate(double *solution) {
     
     // check bounds

     return evalSB_(solution);	

}



