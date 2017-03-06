MODULE modsacess
    USE iso_c_binding
    USE scattersearchtypes
    USE scattersearchfunctions
    USE localsolver
    USE parallelscattersearchfunctions
   
#ifdef MPI2 
        
CONTAINS

    FUNCTION sacess(exp1, fitnessfunction, results1, maxfunevals, ftarget) RESULT (outresult)
    
        ! Declaración de variables
        USE common_functions
        USE qsort_module
        
        IMPLICIT NONE
        TYPE(local_solver_help_vars) :: local_solver_var
        TYPE(algorithm_common_vars) :: common_vars
        TYPE(time_ess) :: time
        TYPE(master_migration) :: migration_master
        TYPE(C_PTR), INTENT(INOUT) :: exp1, results1
        REAL (C_DOUBLE), INTENT(INOUT) :: ftarget
        INTEGER (C_LONG), INTENT(INOUT) :: maxfunevals
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj

        INTEGER :: stopOptimization, fin, nreset
        INTEGER(C_LONG) :: nfuneval, nfunevaltotal
        INTEGER :: i, j, nconst, dest, l
        REAL(C_DOUBLE) :: cputime1
        TYPE(opts) :: opts1, default1
        TYPE(problem) :: problem1
        TYPE(resultsf) :: results
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: tau
        INTEGER :: xbest_in_refset
        
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xbest,fbest, fbest_lastiter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xl_log, xu_log, randomV
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: fbest_new
        
        INTEGER, DIMENSION(:), ALLOCATABLE :: ncomb1, refset_change
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ncomb2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: MaxSubSet, MaxSubSet2
        
        TYPE(Refsettype) :: solutionset, refset, candidateset, childset
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: solutions
        INTEGER :: counter
      
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect
        INTEGER :: nrand
        INTEGER, DIMENSION(:), ALLOCATABLE :: index1, index2, index, diff_index
        INTEGER :: lb_p, auxsize
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: st_p
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: ppp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: hyper_x_L, hyper_x_U, &
                                                                                                factor, v1, v2, v3, new_comb
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: parents_index1, parents_index2
        INTEGER, DIMENSION(:), ALLOCATABLE :: members_update, candidate_update, members_to_update
        INTEGER(C_INT) :: outresult
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::  NAN, INF, threshold_slave
        INTEGER :: index_cooperative
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: aux_val_cooperative
        INTEGER :: tam, out_solver
        
        !COOPERATIVE VARIABLES
        INTEGER :: flag , cooperativempitestess, OMP_NUM_THREADS
        REAL(C_DOUBLE) :: initmpi, calctimempi
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: BKS_x
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: BKS_f
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: last_send
        INTEGER(C_LONG) :: lastevals
        INTEGER :: ncounter_slave_recv_sol, ncounter_slave_send_sol, pending_adaptation
        INTEGER(C_INT) :: getopenmpoption, openmp_pos

        CALL problem_specifications(exp1, problem1,opts1,common_vars%nvar,ftarget,maxfunevals) 
        CALL setnan(NAN)
        CALL setdblmax(INF)
        time%timeparallel = 0.0
        time%localsolvertime =  0.0        
        stopOptimization = 0        
        nfuneval = 0
        fin = 0
        nreset = 0
        common_vars%inititer = 0
        common_vars%init_comm  = 0
        last_send = INF
        lastevals = 0
        index_cooperative = -1
        ncounter_slave_recv_sol = 0
        ncounter_slave_send_sol = 0
        pending_adaptation = 0
        out_solver = 0
        
        CALL setNPROC(exp1, common_vars%NPROC)
        
        ! CHECKBOUNDS : Check if bounds have the same dimension
        CALL checkbounds(problem1,common_vars%status)
        !INIT BEST VARIABLES
        CALL initbestvars(problem1,xbest,fbest,common_vars%nvar)
        !INIT OLD BEST VARIABLES        
        CALL initoldbestvars(BKS_x,BKS_f,common_vars%nvar)
        !INIT INT BIN VARIABLES
        CALL initintbinvars(problem1)
        ! CHARGE PARALLEL ID
        CALL chargeid(exp1,common_vars%idp)
        ! INIT RANDOM GENERATOR        
        CALL initrngrandomparallel(exp1,common_vars%idp)
        ! CREATE MIGRATION MASTER: CREATE VARS TO USE IN SELF-ADAPTATION MODE
        CALL create_migration_master(migration_master,common_vars%NPROC) 
        ! CALC THE SIZE OF THE REFSET
        CALL calc_dim_refset(opts1%globaloptions%dim_refset, common_vars%nvar, opts1%useroptions%iterprint,common_vars%idp,&
        common_vars%NPROC) 
        ! IF THE CONFIGURATION IS HETEROGENEOUS, A DIFFERENT CONFIGURATION IS ASIGNED FOR EACH PROCESSROR
        CALL hete_param_eSS2(exp1, problem1, opts1, common_vars%NPROC, opts1%globaloptions%dim_refset, common_vars%idp,&
        common_vars%nvar)
        ! CALC THE SIZE OF THE NDIVERSE SET
        CALL calc_ndiverse(opts1%globaloptions%ndiverse, common_vars%nvar, opts1%useroptions%iterprint,common_vars%idp)
        ! INITIALIZATION OF LOCAL SOLVER VARS
        CALL initlocalsolvervarsess(exp1)
        ! INITIALIZATION OF OUTPUT VARS
        CALL initoutputvars(exp1)
        ! ENABLE DISTRIBUTED CRITERIA
        CALL setdistcriteria(exp1,1)
        
        if (opts1%globaloptions%ndiverse < opts1%globaloptions%dim_refset) then
            opts1%globaloptions%ndiverse = opts1%globaloptions%dim_refset
        end if
        
        common_vars%parallel = 1
        
        CALL chargecooperativeparametersfortran(exp1, opts1%globaloptions%dim_refset, tam, common_vars%idp, ftarget, maxfunevals)  
        CALL createcooperativetopologyess(exp1)
        
        CALL slave_information(common_vars%idp, opts1%globaloptions%dim_refset, opts1%globaloptions%ndiverse, 0)
        time%starttime = initmpi()
        CALL saveinittime(exp1,time%starttime)

        ! INIT LOCAL OPTIONS
        CALL initlocaloptions(opts1) 
        CALL init_inequality_constraints(problem1%neq, problem1%CL, problem1%CU)
   
        problem1%ineq = problem1%neq + problem1%ineq
       
        CALL calcnconst(problem1,nconst) 
       
        ! CHECK OUTPUT OBJECTIVE FUNCTION
        CALL check_output_obj_funct(problem1, opts1,common_vars%status, nconst)


        ! PASAMOS A LOGARITMO
        ALLOCATE(xl_log(size(problem1%XL)))
        ALLOCATE(xu_log(size(problem1%XU)))
        xl_log = problem1%XL
        xu_log = problem1%XU

        if (ALLOCATED(opts1%useroptions%log_var)) then    
            CALL converttolog2(xl_log, common_vars%nvar,opts1%useroptions%log_var)
            CALL converttolog2(xu_log, common_vars%nvar,opts1%useroptions%log_var)
        end if

        CALL build_aux_local_search(problem1,local_solver_var, opts1, exp1)
        CALL sizeseriallize(common_vars%nvar, nconst, common_vars%sizeser)
!-------------------------------------------------------------------------------------------------------!        
! OPEN MASTER ASYNCHRONOUS BUFFER && CREATE WINDOWS    
!-------------------------------------------------------------------------------------------------------!  
        CALL asynchinitmasterandwindows( exp1, common_vars%sizeser )     
!-------------------------------------------------------------------------------------------------------!  
        
!-------------------------------------------------------------------------------------------------------!        
! CALC INIT SOLUTION SET AND REFSET GLOBAL       
!-------------------------------------------------------------------------------------------------------!     
        if (common_vars%idp .NE. 0) then
            CALL create_init_solutions(exp1,opts1,solutions,common_vars%nvar,xl_log,xu_log);   
            if ((problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 1)) then
                CALL ssm_round_int(solutions, problem1%int_var + problem1%bin_var, problem1%XL, problem1%XU)
            end if
#ifdef OPENMP
            openmp_pos = getopenmpoption(exp1)
            if (openmp_pos .EQ. 1) then
                CALL evaluate_solutions_set_parallel(exp1,fitnessfunction,&
                solutionset,problem1,opts1, solutions, &
                common_vars%nvar, nfuneval, &
                nconst, xl_log, xu_log)
            else 
            CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, common_vars%nvar, &
                nfuneval, nconst, xl_log, xu_log)
            end if
#else
            CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, common_vars%nvar, &
                nfuneval, nconst, xl_log, xu_log)
#endif
            
            CALL create_refset (exp1, solutionset, refset, common_vars%nvar, opts1, nconst, opts1%globaloptions%dim_refset )

        else
            CALL create_refset_empty( refset, common_vars%nvar, opts1, nconst, opts1%globaloptions%dim_refset )
        end if
!-------------------------------------------------------------------------------------------------------!         


        ! Possible combinations among the opts1%globaloptions%dim_refset elements in Refset, taken in pairs
        !
        ALLOCATE(ncomb1(opts1%globaloptions%dim_refset))
        ncomb1 = (/(i, i = 1, opts1%globaloptions%dim_refset)/)
        call nchoosek_fpar(ncomb1, 2, ncomb2)
        
        MaxSubSet = (opts1%globaloptions%dim_refset**2 - opts1%globaloptions%dim_refset)/2
        MaxSubSet2 = 2 * MaxSubSet        

        !Check best value in refset
        if (common_vars%idp .NE. 0)   CALL check_best_value(refset, opts1, xbest, fbest )
        common_vars%iter = 0
        CALL optimization_begins(common_vars%idp)

        if ((opts1%useroptions%iterprint .eq. 1) .AND. (common_vars%idp .NE. 0)) then
                cputime1 = calctimeMPI(exp1,time%starttime)
                WRITE(*, '(A14, I3, A26, I10, A8, D25.10, A10, F12.2, A11, E10.1)') "[PROCESSOR ID=", &
                    common_vars%idp, "] Initial Pop: NFunEvals: ", nfuneval, &
                    " Bestf: ", fbest(1), " CPUTime: ", cputime1, " variance: ", variance_vector(refset%fpen)
        endif       
        CALL create_results(results, opts1, xbest, fbest, nfuneval, cputime1, common_vars%nvar)
        
!-------------------------------------------------------------------------------------------------------!
        
        ! Sort the Refset
        if (common_vars%idp .NE. 0 ) CALL sort_refset(refset,opts1, indvect)

        if (common_vars%idp .EQ. 0 ) then
            refset%x = NAN
            refset%fpen = NAN
            refset%f = NAN
            refset%penalty = NAN
            if (nconst .GT. 0) then
                refset%nlc = NAN         
            end if            
        end if

                
        if (opts1%globaloptions%combination .eq. 1) then         
            nrand = common_vars%nvar
        else if (opts1%globaloptions%combination .eq. 2) then
            nrand = 1
        end if

        ! CALC INDEXES
        auxsize = choose(ncomb1)
        ALLOCATE(index1(auxsize))
        ALLOCATE(index2(auxsize))
        ALLOCATE(index(auxsize))
        ALLOCATE(diff_index(auxsize))
        index1 = ncomb2(:, 1)
        index2 = ncomb2(:, 2)
        index = ncomb2(:, 2)
        
        CALL FUSION_VECTOR_INT(index1, index)

        ! CALC PPP
        diff_index = (index2 - index1)
        lb_p = 1 !Minimum ppp
        st_p = REAL(0.75d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        !Defines maximum ppp. max(ppp)= lb_p+ st_p
        !Example: lb_p=1.25 and st_p=0.5 varies ppp from
        !%1.25 to 1.75

        ALLOCATE (ppp(auxsize))
        ppp = st_p * REAL((diff_index - 1),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/ &
           REAL((opts1%globaloptions%dim_refset - 2),KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) + &
           REAL(lb_p,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        ! CALC HYPER
        ALLOCATE(hyper_x_L(common_vars%nvar,CEILING(MaxSubSet)))
        ALLOCATE(hyper_x_U(common_vars%nvar,CEILING(MaxSubSet)))
        
        i=1
        j=1
        hyper_x_L = reshape((/ ((problem1%XL(i), i = 1,common_vars%nvar), j = 1,  CEILING(MaxSubSet)) /), &
        (/ common_vars%nvar, CEILING(MaxSubSet) /))
        hyper_x_U = reshape((/ ((problem1%XU(i), i = 1,common_vars%nvar), j = 1,  CEILING(MaxSubSet)) /), &
        (/ common_vars%nvar, CEILING(MaxSubSet) /))



        ALLOCATE(refset_change(opts1%globaloptions%dim_refset))
        refset_change = 0


      
!-------------------------------------------------------------------------------------------------------!        
! INIT ASYNCHRONOUS STOPPING CRITERIA         
!-------------------------------------------------------------------------------------------------------!        
!        if (common_vars%idp .EQ. 0) then
!                CALL initcooperativestoppingcriteriaessmaster(exp1)
!        else
!                CALL initcooperativestoppingcriteriaessslave(exp1)
!        endif
!-------------------------------------------------------------------------------------------------------!
! INIT ASYNCHRONOUS STOPPING CRITERIA
!-------------------------------------------------------------------------------------------------------!
        CALL initcooperativestoppingcriteriaess( exp1 )
!-------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------!  
        

        cputime1 = calctimeMPI(exp1,time%starttime)
        CALL initprintfile(exp1, fbest, common_vars%parallel, common_vars%idp, cputime1, nfuneval, 1)    

        
        CALL printgant(exp1,cputime1,1) 
        lastevals = nfuneval
 
        
!-------------------------------------------------------------------------------------------------------!         
! INIT : Algorithm main loop
!-------------------------------------------------------------------------------------------------------! 
    do while (out_solver .EQ. 0)
!----------------------------------------------------------------------------------------------------------------------------------
! MIGRATION ASYNCHRONOUS REGION            
!----------------------------------------------------------------------------------------------------------------------------------   
!----------------------------------------------------------------------------------------------------------------------------------
!       MASTER CODE            
!----------------------------------------------------------------------------------------------------------------------------------   
          if ( common_vars%idp .EQ. 0 ) then
             CALL asynchronous_master_acess (exp1,opts1,problem1,refset,fin,nfuneval,fbest,stopOptimization,common_vars%sizeser, &
                                    nconst, common_vars%nvar, time%starttime, migration_master%vector_proc, &
                                    migration_master%ocurrences_send, common_vars  )   
             cputime1 = calctimeMPI(exp1,time%starttime)
             CALL printgant(exp1,cputime1,1)  
!----------------------------------------------------------------------------------------------------------------------------------
!       SLAVE CODE            
!----------------------------------------------------------------------------------------------------------------------------------        
          else
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printinititeration(exp1, common_vars%iter, cputime1)

                if ((common_vars%iter .NE. 0) ) then
                    cputime1 = calctimeMPI(exp1,time%starttime)

                    CALL asynchronous_slave_acess (exp1,opts1,problem1,refset,common_vars, nconst,BKS_x,BKS_f, &
                           time,refset_change,fitnessfunction,xl_log,xu_log, nfuneval, fbest, xbest, & 
                           ncounter_slave_recv_sol, ncounter_slave_send_sol,lastevals)
                           
                    CALL adapt_slave(exp1,opts1,problem1,fitnessfunction,common_vars,nfuneval,xl_log,xu_log, &
                           hyper_x_L,hyper_x_U,MaxSubSet,MaxSubSet2,ppp,index1, index2, index,nconst,refset,&
                           fbest,xbest,refset_change,local_solver_var,time, &
                           ncounter_slave_recv_sol, ncounter_slave_send_sol,pending_adaptation,lastevals)        

                    CALL select_better_solution(results,opts1,refset,xbest,fbest,local_solver_var%use_bestx,nfuneval,cputime1,fin)

                else 
                    BKS_f(1) = fbest(1)
                    BKS_x(:,1) = xbest
                end if

          end if
!----------------------------------------------------------------------------------------------------------------------------------   
            
          if ( fin .GT. 0 ) out_solver = 1 
            
          if ( (fin .LT. 1) .AND. (common_vars%idp .NE. 0 )) then
!---------------------------------------------------------------------------------------------------------             

            
            CALL create_childset ( childset, common_vars%nvar, CEILING(MaxSubSet2), nconst )   
            counter = 1
            if (.not.ALLOCATED(members_update)) then
                ALLOCATE(members_update(opts1%globaloptions%dim_refset))
                members_update =  1
            end if  

            CALL check_duplicated_replace (exp1, fitnessfunction, problem1,refset,opts1, parents_index1, parents_index2, &
                        xl_log, xu_log, refset_change, members_update, index1, index2, nfuneval, common_vars%nvar, nconst )
     
            CALL create_candidate(candidateset, refset)
            if (.not.ALLOCATED(candidate_update)) then
                ALLOCATE(candidate_update(opts1%globaloptions%dim_refset))
                candidate_update = 0
            end if

            ALLOCATE( factor(common_vars%nvar, size(ppp)))
            factor = reshape((/ ((ppp(j), i = 1, common_vars%nvar), j = 1, size(ppp)) /), (/ common_vars%nvar,size(ppp) /))
            factor = factor * (parents_index2 - parents_index1)/ REAL(1.5d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) 

            ALLOCATE( v1(common_vars%nvar, size(ppp)))
            ALLOCATE( v2(common_vars%nvar, size(ppp)))
            ALLOCATE( v3(common_vars%nvar, size(ppp)))

            v1 = parents_index1 - factor
            v2 = parents_index2 - factor
            v3 = REAL(2.0d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  * parents_index2  - parents_index1 - factor  

            CALL check_vector_in_hyper_bounds(exp1,opts1, v1, hyper_x_L, hyper_x_U )
            CALL check_vector_in_hyper_bounds(exp1,opts1, v3, hyper_x_L, hyper_x_U )
            CALL generate_new_comb_matrix(exp1, new_comb, ppp, MaxSubSet, nrand, v1,v2,v3)


#ifdef OPENMP
            openmp_pos = getopenmpoption(exp1)           
            if (openmp_pos .EQ. 1) then 
            CALL update_candidateset_with_new_comb_parallel( exp1,opts1, fitnessfunction,problem1,new_comb,& 
                candidateset,childset,candidate_update, members_update,MaxSubSet2,nrand,nconst,nfuneval,counter, &
                index)
            else
            CALL update_candidateset_with_new_comb( exp1,opts1,fitnessfunction,problem1,new_comb,candidateset,&
            childset,candidate_update,members_update,MaxSubSet2,nrand,nconst,nfuneval,counter, index )
            end if
#else              
            CALL update_candidateset_with_new_comb( exp1,opts1,fitnessfunction,problem1,new_comb,candidateset,&
            childset,candidate_update,members_update,MaxSubSet2,nrand,nconst,nfuneval,counter, index )
#endif

            
!----------------------------------------------------------------------------------------------------------------------------------
! CHECK PARALLEL STOPPING CRITERIA           
!----------------------------------------------------------------------------------------------------------------------------------              
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results, opts1, refset, xbest, fbest,local_solver_var%use_bestx,nfuneval, cputime1, fin)
            CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0)

            if (fin < 1) THEN
                do l = 1, common_vars%NPROC
                    if (fin .NE. 3) then
                        dest = l - 1
                        if (dest .NE. common_vars%idp) THEN
                            flag = cooperativempitestess(exp1, dest )
                            if (flag .EQ. 1) THEN
                                fin = 3
                            end if
                        end if
                    end if
                end do
            end if
            cputime1 = calctimeMPI(exp1,time%starttime)

           CALL asynchronousstoppingcriteriaess(exp1, fin, nfuneval, fbest, cputime1, opts1%useroptions%maxtime, stopOptimization)
           CALL printcheckvtr(exp1, fbest(1), nfuneval, cputime1,results%timevtr )

!----------------------------------------------------------------------------------------------------------------------------------
            
            
            ALLOCATE(members_to_update(size(candidate_update)))
            members_to_update = members_to_update * 0
            CALL index_of_ones(candidate_update, members_to_update)
            CALL reajust_index(members_to_update)

           CALL apply_beyond_to_members_to_update( exp1, opts1, fitnessfunction, problem1, members_to_update,  nfuneval, &
                  refset, candidateset, nconst, nrand, common_vars%nvar)                 
              
            if ( .not. ALLOCATED(fbest_lastiter) ) then
                ALLOCATE(fbest_lastiter(size(fbest)))
            end if
            fbest_lastiter=fbest
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results, opts1, refset, xbest, fbest, local_solver_var%use_bestx, &
                                                                                        nfuneval, cputime1, fin)    

            common_vars%iter=common_vars%iter+1
            local_solver_var%n_minimo=local_solver_var%n_minimo+1 
         
            if (common_vars%idp .NE. 0) then
                CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0) 
            end if
            
            CALL printiterationcesslog( exp1, fbest(1), nfuneval, cputime1, common_vars%iter)
            CALL update_refset_change(members_update, refset_change)  
            CALL remove_possible_stuck_members_in_refset(problem1, opts1, exp1, fitnessfunction, nconst, &
                        refset, refset_change,nfuneval, common_vars%nvar, xl_log, xu_log )               
            CALL check_the_improvement_of_fbest(results, opts1, fbest, xbest, fbest_lastiter, nreset, cputime1, fin, nfuneval)
  
            if (( opts1%localoptions%empty .eq. 0 )  .and. (fin .eq. 0)) then
                   
             if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (opts1%localoptions%bestx .gt. 0) .and. &
                       (local_solver_var%n_minimo .ge. local_solver_var%n_critico) ) then    
                  if ( local_solver_var%use_bestx .gt. 0 ) then

                    cputime1 = calctimeMPI(exp1,time%starttime)
                    CALL printgant(exp1,cputime1,1)          
                    CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                                childset,xbest,fbest,nfuneval,fbest_lastiter,refset,refset_change,1)
                    
                    cputime1 = calctimeMPI(exp1,time%starttime)
                    CALL printgant(exp1,cputime1,2)                       
                    local_solver_var%n_minimo=0
                  end if
              else 
                  if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (local_solver_var%n_minimo .ge. &
                  local_solver_var%n_critico) ) then
                    cputime1 = calctimeMPI(exp1,time%starttime)
                    CALL printgant(exp1,cputime1,1)                          
                    CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                                childset,xbest,fbest,nfuneval,fbest_lastiter,refset,refset_change,1)
                    
                    
                    cputime1 = calctimeMPI(exp1,time%starttime)
                    CALL printgant(exp1,cputime1,2)                       
                    local_solver_var%n_minimo=0
                  end if 
              end if
                
            end if
            
            

            fbest_lastiter=fbest
            
             
!----------------------------------------------------------------------------------------------------------------------------------
! CHECK PARALLEL STOPPING CRITERIA           
!----------------------------------------------------------------------------------------------------------------------------------              
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results, opts1, refset, xbest, fbest,local_solver_var%use_bestx,nfuneval, cputime1, fin)
            CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0)

            if (fin < 1) THEN
                do l = 1, common_vars%NPROC
                    if (fin .NE. 3) then
                        dest = l - 1
                        if (dest .NE. common_vars%idp) THEN
                            flag = cooperativempitestess(exp1, dest )
                            if (flag .EQ. 1) THEN
                                fin = 3
                            end if
                        end if
                    end if
                end do
            end if
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL asynchronousstoppingcriteriaess(exp1, fin, nfuneval, fbest, cputime1, opts1%useroptions%maxtime, stopOptimization)
            CALL printcheckvtr(exp1, fbest(1), nfuneval, cputime1,results%timevtr )
            
!----------------------------------------------------------------------------------------------------------------------------------              
           
            CALL finalize_algorithm(exp1,problem1,results,opts1,fitnessfunction,refset,xbest,fbest,fin,nfuneval,common_vars%NPROC, &
                    cputime1, nconst, common_vars%nvar,local_solver_var,common_vars%idp)            
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0) 
            CALL addlocalscounteress(exp1)
            if (fin > 0) then
                CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0)
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printgant(exp1,cputime1,1)                       
            end if

            
            CALL destroy_refsettype(candidateset)
            CALL destroy_refsettype(childset)
            
            
            if (ALLOCATED(parents_index1)) DEALLOCATE(parents_index1)
            if (ALLOCATED(parents_index2)) DEALLOCATE(parents_index2)

            if (ALLOCATED(factor))  DEALLOCATE(factor)
            if (ALLOCATED(v1))  DEALLOCATE(v1)
            if (ALLOCATED(v2))  DEALLOCATE(v2)
            if (ALLOCATED(v3))  DEALLOCATE(v3)
            if (ALLOCATED(candidate_update)) DEALLOCATE(candidate_update)
            if (ALLOCATED(members_update)) DEALLOCATE(members_update)
            if (ALLOCATED(new_comb))  DEALLOCATE(new_comb)
            if (ALLOCATED(members_to_update))  DEALLOCATE(members_to_update)

            if (opts1%useroptions%iterprint .eq. 1 ) then
                cputime1 = calctimeMPI(exp1,time%starttime)
            endif
            
            IF (ALLOCATED(indvect) ) DEALLOCATE(indvect)
           ! CALL sort_refset(refset,opts1, indvect)
           ! refset_change = refset_change(indvect)
            CALL printrefset(exp1, refset%fpen,  opts1%globaloptions%dim_refset, cputime1, common_vars%iter, refset_change ) 
            CALL printcooperative( exp1, refset%cooperative, opts1%globaloptions%dim_refset, cputime1,  common_vars%iter ) 
            
            CALL printenditeration(exp1)
            !index_cooperative, aux_val_cooperative


          end if
          
          
        end do

        CALL mpibarrieress()
        CALL ending_solver_message(common_vars%idp, fin)
        CALL mpibarrieress()    
        
!---------------------------------------------------------------------------------------------------------
! GATHER THE RESULTS IN MPI EXECUTION CASE
!---------------------------------------------------------------------------------------------------------
        cputime1 = calctimeMPI(exp1,time%starttime)
        CALL select_better_solution(results, opts1, refset, xbest, fbest, local_solver_var%use_bestx, nfuneval, cputime1, fin)  
        CALL returnminlocelement(exp1, fbest(1), fbest_new, i, common_vars%idp)
        CALL cooperativebcastelement(exp1,xbest,common_vars%nvar,i)
        CALL returnsumelementlong(exp1, nfuneval, nfunevaltotal)
        CALL returnavgelementintdouble(exp1, common_vars%iter, results%totaliter, common_vars%NPROC-1)
        CALL returnminelement(exp1, results%timevtr, results%totalvtr)
        cputime1 = calctimeMPI(exp1,time%starttime)
        outresult = 1
        results%timetotal = REAL(cputime1,KIND=C_DOUBLE)
        fbest=fbest_new
        CALL detranslationinterface( xbest, exp1 )
        CALL printmasterocurrencesend( exp1, migration_master%ocurrences_send, common_vars%NPROC, 1, &
        opts1%globaloptions%dim_refset, &
                        opts1%localoptions%n2,opts1%localoptions%balance)


!---------------------------------------------------------------------------------------------------------
        CALL mpibarrieress() 
        CALL seed_recount(exp1, common_vars%idp)
        CALL mpibarrieress()   
                
        
        CALL updateresultsess(exp1, results1, results%timetotal, nfunevaltotal, fbest(1), xbest, results%totaliter,results%totalvtr)

        
        CALL settotaltime(results1, results%timetotal)
        CALL setparalleltime(results1, time%timeparallel)
        CALL setlocalsolvertime(results1, time%localsolvertime)
        

        if (common_vars%idp .EQ. 0) then
                CALL  printresults(exp1, fbest, nfunevaltotal,results%timetotal, results%totalvtr,results%totaliter)
        end if
        
        
        ! Free memory
        CALL destroyasynchinitmasterandwindows(exp1)
        CALL destroy_opts(opts1)
        CALL destroy_opts(default1)
        CALL destroy_problem(problem1)

        
        if (ALLOCATED(xl_log)) DEALLOCATE(xl_log)
        if (ALLOCATED(xu_log)) DEALLOCATE(xu_log)   
        if (ALLOCATED(BKS_x)) DEALLOCATE(BKS_x)
        if (ALLOCATED(BKS_f)) DEALLOCATE(BKS_f)
        if (ALLOCATED(randomV)) DEALLOCATE(randomV)
        
        if (ALLOCATED(solutions)) DEALLOCATE(solutions)
        if (ALLOCATED(ncomb1)) DEALLOCATE(ncomb1)
        if (ALLOCATED(ncomb2)) DEALLOCATE(ncomb2)
                 
        if (ALLOCATED(indvect)) DEALLOCATE(indvect)
        if (ALLOCATED(index1)) DEALLOCATE(index1)
        if (ALLOCATED(index2)) DEALLOCATE(index2)
        if (ALLOCATED(index)) DEALLOCATE(index)
        if (ALLOCATED(diff_index)) DEALLOCATE(diff_index)
        if (ALLOCATED(factor)) DEALLOCATE(factor)
        if (ALLOCATED(ppp)) DEALLOCATE(ppp)
        if (ALLOCATED(hyper_x_L)) DEALLOCATE(hyper_x_L)
        if (ALLOCATED(hyper_x_U)) DEALLOCATE(hyper_x_U)
        if (ALLOCATED(refset_change)) DEALLOCATE(refset_change)
        if (ALLOCATED(opts1%localoptions%finish)) DEALLOCATE(opts1%localoptions%finish)
        if (ALLOCATED(xbest)) DEALLOCATE(xbest)
        if (ALLOCATED(fbest)) DEALLOCATE(fbest)
        if (ALLOCATED(fbest_lastiter)) DEALLOCATE(fbest_lastiter)
        CALL destroy_local_solver_vars(local_solver_var)     
        CALL destroy_migration_master(migration_master) 
        CALL mpibarrieress()
        CALL destroy_refsettype(solutionset)
        CALL mpibarrieress()
        CALL destroy_refsettype(refset)
        CALL mpibarrieress()
        CALL destroy_resultsf(results)

        CALL mpibarrieress()

    END FUNCTION sacess
#endif

    
END MODULE modsacess
