#include <stdlib.h>
#include <string.h>
#include <AMIGO_problem.h>
#include <ggn.h>
#include <amigoRHS_CIRCADIAN.h>
#include <amigoRHS_3step.h>
#include <amigoRHS_NFKB.h>
#include <amigoRHS_B1.h>
#include <amigoRHS_B2.h>
#include <amigoRHS_B3.h>
#include <amigoRHS_B4.h>
#include <amigoRHS_B5.h>
#include <math.h>
#include <setup_benchmarks.h>

int load_benchmark_SystemBiology(int current_bench) {
    const char *path;
    int i,j;
    int counter, exito, init_cond;
    double point;
    int *index_non_obs;

    use_amigo=1;
    estime_init_cond=0; 
    if (current_bench == 0) {
	printf("LOAD CIRCADIAN PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/others/circadian/load_C.mat";
    } 
    else if (current_bench == 1) {
	printf("LOAD SSP PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/others/3-step_pathway/load_C.mat";
    } 
    else if (current_bench == 2) {
	printf("LOAD NFKB PROBLEM\n");    
        path = "./source/benchmarks/systemsBiology/others/Nfkb/load_C.mat";
    } 
    else if (current_bench == 3) {
	printf("LOAD B1 PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/BioPredyn/B1/load_C.mat";
    }  
    else if (current_bench == 4) {
	printf("LOAD B2 PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/BioPredyn/B2/load_C.mat";
    }  
    else if (current_bench == 5) {
	printf("LOAD B3 PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/BioPredyn/B3/load_C.mat";
    } 
    else if (current_bench == 6) {
	printf("LOAD B4 PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/BioPredyn/B4/load_C.mat";
    }  
    else if (current_bench == 7) {
	printf("LOAD B5 PROBLEM\n");
        path = "./source/benchmarks/systemsBiology/BioPredyn/B5/load_C.mat";
    } 
    else if (current_bench == 8) {
	printf("LOAD B6 PROBLEM\n");
	use_amigo=0;
    } else {
        exit(14);
    }
    int type=current_bench; 

// OPEN MAT CONDITIONAL 
    amigo =  openMatFileAMIGO(path); 
    if (type == 0) {
       	    set_AMIGO_problem_rhs(amigo, amigoRHS_CIRCADIAN, amigo_Y_at_tcon_CIRCADIAN);
      	    set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_CIRCADIAN, amigoRHS_get_sens_OBS_CIRCADIAN);
    } else if (type == 1) {
       	    set_AMIGO_problem_rhs(amigo, amigoRHS_MENDES, amigo_Y_at_tcon_MENDES);
       	    set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_MENDES, amigoRHS_get_sens_OBS_MENDES);
    } else if (type == 2) {
       	    set_AMIGO_problem_rhs(amigo, amigoRHS_NFKB, amigo_Y_at_tcon_NFKB);
       	    set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_NFKB, amigoRHS_get_sens_OBS_NFKB);    
    } else if (type == 3) {
            set_AMIGO_problem_rhs(amigo, amigoRHS_B1, amigo_Y_at_tcon_B1);
            set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_B1, amigoRHS_get_sens_OBS_B1);    
    } else if (type == 4) {
            set_AMIGO_problem_rhs(amigo, amigoRHS_B2, amigo_Y_at_tcon_B2);
            set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_B2, amigoRHS_get_sens_OBS_B2);   
    } else if (type == 5) {
            set_AMIGO_problem_rhs(amigo, amigoRHS_B3, amigo_Y_at_tcon_B3);
            set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_B3, amigoRHS_get_sens_OBS_B3);
    } else if (type == 6) {
            set_AMIGO_problem_rhs(amigo, amigoRHS_B4, amigo_Y_at_tcon_B4);
            set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_B4, amigoRHS_get_sens_OBS_B4);   
    } else if (type == 7) {
            set_AMIGO_problem_rhs(amigo, amigoRHS_B5, amigo_Y_at_tcon_B5);
            set_AMIGO_problem_obs_function(amigo, amigoRHS_get_OBS_B5, amigoRHS_get_sens_OBS_B5);
    }


    if ( estime_init_cond == 1) {
            init_cond = amigo->amigo_models[0]->n_states - amigo->amigo_models[0]->n_observables;
            dim = init_cond + amigo->nx;
            index_non_obs = (int *) malloc(  init_cond * sizeof(int) );
            counter = 0;
            for (i=0;i<amigo->amigo_models[0]->n_states;i++){
                exito=0;
                for (j=0;j<amigo->amigo_models[0]->n_observables;j++) {
                    if (amigo->amigo_models[0]->index_observables[j] == i){
                        exito = 1;
                        break;
                    }
                }
                if (exito == 0){
                    index_non_obs[counter]=i;
                    counter++;
                }
            }
            
            max_dom = (double *) malloc(dim * sizeof (double));
            min_dom = (double *) malloc(dim * sizeof (double));            
            /*for (i = 0; i < init_cond; i++) {
                exp->test.bench.max_dom[i] = 1.0;
                exp->test.bench.min_dom[i] = 0.0;
            }
            */
            counter = 0;
            for (i = init_cond; i < dim; i++) {

                max_dom[i] = amigo->UB[counter];
                min_dom[i] = amigo->LB[counter];
                counter++;
            }    
            
            free(index_non_obs);
            
     } 
    else {
	    if (use_amigo) {
           	dim = amigo->nx;
	   	max_dom = (double *) malloc(dim * sizeof (double));
           	min_dom = (double *) malloc(dim * sizeof (double));
           	for (i = 0; i < dim; i++) {
           	   max_dom[i] = amigo->UB[i];
           	   min_dom[i] = amigo->LB[i];
           	}
	    } else { 
    		if ( type == 8) {
           		dim = 37;
	                max_dom = (double *) malloc(dim * sizeof (double));
	                min_dom = (double *) malloc(dim * sizeof (double));
           		returnbounds( max_dom , min_dom );
    	    	}
    	    }
    }

    return 1;
    
}


void manage_init_cond(AMIGO_problem *amigo, double *U, double *U_aux) {
    int n_IC, counter;
    int i, j, exito, *index_non_obs;
    
    n_IC = amigo->amigo_models[0]->n_states - amigo->amigo_models[0]->n_observables;
    index_non_obs = (int *) malloc(n_IC * sizeof (int));
    counter = 0;
    for (i = 0; i < amigo->amigo_models[0]->n_states; i++) {
        exito = 0;
        for (j = 0; j < amigo->amigo_models[0]->n_observables; j++) {
            if (amigo->amigo_models[0]->index_observables[j] == i) {
                exito = 1;
                break;
            }
        }
        if (exito == 0) {
            index_non_obs[counter] = i;
            counter++;
        }
    }

    for (i=0;i<amigo->n_exp;i++){
        for (j=0;j<n_IC;j++) {
            amigo->amigo_models[i]->y0[index_non_obs[j]]=U[j];
        }
    }
    
    
    
    counter=0;
    for (i=n_IC;i<dim;i++){
        U_aux[counter]=U[i];
        counter++;
    }    
}


void repack_init_cond(AMIGO_problem *amigo, double *U, double *U_aux) {
    int n_IC,i,max,counter;
    
    n_IC = amigo->amigo_models[0]->n_states - amigo->amigo_models[0]->n_observables;
    max = dim;
    counter=0;
    for (i=n_IC; i<max ;i++){
        U[i]=U_aux[counter];
        counter++;
    }        
}

/*
int destroySystemBiology(experiment_total *exp) {

    if (exp->amigo != NULL) {
        free_AMIGO_problem(exp->amigo);
        exp->amigo=NULL;
    }
    free(exp->execution.transconst);
    free(exp->test.bench.logindex);
    free(exp->test.bench.log_max_dom );
    free(exp->test.bench.log_min_dom );


    return 1;

}
*/

double evalSB_(double *U) {
    double *U_aux;
    double cost;

    if (  use_amigo == 1 ) { 
        if (estime_init_cond == 1) {
            U_aux = (double *) malloc(amigo->nx*sizeof(double));
            manage_init_cond(amigo, U, U_aux);
            set_AMIGO_problem_pars(U_aux, amigo);
            cost = eval_AMIGO_problem_LSQ(amigo);
            amigo->nevals++;    
            repack_init_cond(amigo, U, U_aux);

            free(U_aux);
        } else {
            set_AMIGO_problem_pars(U, amigo);
            cost = eval_AMIGO_problem_LSQ(amigo);
            amigo->nevals++;
        }
        
        
    } else {
	if (id_problem == 8) cost= fitnessfunctionB6(U); 
    }

    return cost;

}




