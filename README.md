# saCeSS global optimization library  #

## Self-adaptive Cooperative enhanced Scatter Search (saCeSS)##

This is version 0.1 of a novel parallel global optimization code, self-adaptive 
cooperative enhanced scatter search (saCeSS). This code is distributed as a library with several parallel solvers based on the scatter search metaheuristic, incorporating several key new mechanisms: 

1. Asynchronous cooperation between parallel processes.
2. Coarse and fine-grained parallelism.
3. Self-tuning strategies.

The saCeSS code has been implemented using Fortran 90 and C. Parallelization has been implemented using MPI and openMP. It has been tested in Linux clusters running CentOS 6.7 and Debian 8.

The saCeSS library allows the solution of non-linear programming (NLP) and mixed-integer non-linear programming (MINLP) problems. It also provides efficient local solvers for nonlinear parameter estimation problems associated with complex models (e.g. those described by differential equations). The current distribution of saCeSS includes a set of optimization examples that can be used as benchmarks, taken from the BBOB and BioPreDyn testbeds.

## REFERENCES ##

### Main reference: ###

Penas, D.R., P. Gonzalez, J.A. Egea, R. Doallo and J.R. Banga (2017) Parameter estimation in large-scale systems biology models: a parallel and self-adaptive cooperative strategy. BMC Bioinformatics 18:52.

### Related previous papers: ###

Penas DR, P Gonzalez, JA Egea, JR Banga, R Doallo (2015) Parallel Metaheuristics in Computational Biology: An Asynchronous Cooperative Enhanced Scatter Search Method. Procedia Computer Science, 51:630-639.

Villaverde, A.F., J.A. Egea and J.R. Banga (2012) A cooperative strategy for parameter estimation in large scale systems biology models. BMC Systems Biology 6:75.

### Benchmark problems distributed with this library: ###

Villaverde AF, D Henriques, K Smallbone, S Bongard, J Schmid, D Cicin-Sain, A Crombach, J Saez-Rodriguez, K Mauch, E Balsa-Canto, P Mendes, J Jaeger and JR Banga (2015) BioPreDyn-bench: a suite of benchmark problems for dynamic modelling in systems biology. BMC Systems Biology 9:8.
