MODULE parallelscattersearchfunctions
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
    USE qsort_module
    USE scattersearchfunctions
    USE funcevalinterface
    USE localsolver
#ifdef OPENMP
    USE omp_lib
#endif
    
#ifdef MPI2
    
CONTAINS


! ----------------------------------------------------------------------
! SUBROUTINES asynchronous_master_acess
! ----------------------------------------------------------------------
SUBROUTINE asynchronous_master_acess (exp1,opts1,problem1,refset,fin,nfuneval,fbest,stopOptimization, sizeser, &
                                    nconst, nvar, starttime, vector_proc, ocurrences_send,common_vars)
    IMPLICIT NONE    
    TYPE(opts), INTENT(IN) :: opts1
    TYPE(C_PTR),  INTENT(INOUT) :: exp1
    TYPE(algorithm_common_vars) :: common_vars    
    TYPE(problem), INTENT(INOUT) :: problem1
    INTEGER, INTENT(INOUT) :: fin, stopOptimization
    INTEGER :: init
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) ::  ocurrences_send
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: vector_proc

    INTEGER(C_LONG), INTENT(IN) :: nfuneval
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: fbest
    REAL(C_DOUBLE), INTENT(INOUT) :: starttime
    INTEGER, INTENT(IN) :: sizeser, nconst, nvar
    TYPE(Refsettype), INTENT(INOUT) :: refset
    INTEGER :: flag, cooperativempitestess, dest, NPROC, i,resultflag,send_solution_to_slaves, checkcandidate
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: VTR
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::cputime3,cputime4
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE :: serializedata, &
                                                                               selectserializedata!, serializedatapending
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE :: BKS_M_x
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE :: BKS_M_f
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: last_send
    TYPE(Refsettype) :: newsol, candidatesol
    REAL(C_DOUBLE) :: calctimempi
    INTEGER :: idsender, first_communication
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: INF, lasttime, threshold_imp
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  ::  percentage    
    INTEGER :: diffusions, lastdiffusions, lastdiffusionsset, nrejects
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE :: slaveComm, slaveScore, slaveScore_aux
    INTEGER,DIMENSION(:),ALLOCATABLE :: orderindexscore, counter_adapt_request
    INTEGER :: resultflag_adapt, idslave_to_adapt, accept, MAX_ADAPT, counter_adapt, bestid, send_counter
    TYPE(adapt_vars), DIMENSION(:), ALLOCATABLE :: adapt_master, adapt_master_origin

    
    init = 0
    CALL setNPROC(exp1, NPROC)
    CALL setdblmax(INF)
    last_send = INF
    CALL init_adapt_vars(exp1,adapt_master, NPROC, common_vars%nvar)
    CALL init_adapt_vars(exp1,adapt_master_origin,NPROC, common_vars%nvar)
    ALLOCATE(slaveComm(NPROC-1))
    ALLOCATE(slaveScore(NPROC-1))
    ALLOCATE(counter_adapt_request(NPROC-1))
    slaveComm = 0d0
    slaveScore = 0d0
    counter_adapt_request = 0
    MAX_ADAPT = FLOOR(1d0 * REAL(NPROC-1,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/2d0)

    CALL setparallelsacessfieldsmaster(exp1, threshold_imp );

    nrejects = 0

    ALLOCATE(selectserializedata(sizeser,1))
    ALLOCATE(BKS_M_x(nvar,1))
    ALLOCATE(BKS_M_f(1))
    BKS_M_x=0
    BKS_M_f=INF
    lasttime=0d0
    first_communication = 1

    do while (fin < 1) 
        ! CHECK STOPPING CRITERIA
        cputime3 = calctimeMPI(exp1,starttime)
        if (fin < 1) THEN
            do i = 1, NPROC-1
                if (fin .NE. 3) then
                    dest = i - 1
                    if (dest .NE. 0) THEN
                        flag = cooperativempitestess(exp1, dest )
                        if (flag .EQ. 1) THEN
                            fin = 3
                        end if
                    end if
                end if
            end do
        end if
        CALL asynchronousstoppingcriteriaess(exp1, fin, nfuneval,BKS_M_f,cputime3,&
                                                    opts1%useroptions%maxtime,stopOptimization)
        CALL create_refset_empty ( newsol, nvar, opts1, nconst, 1 )
        CALL create_refset_empty ( candidatesol, nvar, opts1, nconst, 1 )
        candidatesol%fpen=BKS_M_f
        candidatesol%x=BKS_M_x
        ALLOCATE(serializedata(sizeser,1))
        send_solution_to_slaves = 0
        checkcandidate = 0
        
! MASTER RECEIVED FROM SLAVES              
        if (fin < 1) THEN
            do i = 1, NPROC-1
                    resultflag = 1
                    CALL receivesolutionsmaster(exp1, i, serializedata,sizeser,resultflag)

                    if (resultflag .EQ. 1) then
                        CALL deserializerefset(newsol, nvar, 1, serializedata)
                        if (newsol%fpen(1) .NE. last_send) then
                                cputime3 = calctimeMPI(exp1,starttime)
                                CALL printrecvmasterlog(exp1, newsol%fpen(1),candidatesol%fpen(1),cputime3, i)
                        end if
                        if ( newsol%fpen(1) .LT. candidatesol%fpen(1) ) then
                            candidatesol = newsol
                            checkcandidate=1
                            idsender = i
                            selectserializedata = serializedata
                        end if
                    end if
            end do
        end if

        
! MASTER SEND TO SLAVES      
        if ((checkcandidate .EQ. 1)) then
            if (first_communication .EQ. 0) then
                cputime3 = calctimeMPI(exp1,starttime)        
                percentage = ( (( BKS_M_f(1) -  candidatesol%fpen(1) )  /  BKS_M_f(1) ) * 100d0 )
                CALL calc_percentage( BKS_M_f(1) ,  candidatesol%fpen(1) , percentage)

                if (candidatesol%fpen(1) .NE. last_send) then
                        CALL printcomparenewsolutionmasterlog( exp1, candidatesol%fpen(1), BKS_M_f(1), percentage , &
                                       cputime3,threshold_imp ) 
                        slaveComm(idsender) = slaveComm(idsender) + 1d0
                        slaveScore(idsender) = slaveComm(idsender)*calctimeMPI(exp1,starttime)
                end if
                
                if ( (percentage .ge. threshold_imp ) .OR. ( candidatesol%fpen(1) .LT. problem1%vtr(1) )) then     
                    send_solution_to_slaves = 1   
                else
                    if (candidatesol%fpen(1) .NE. last_send) then 
                                nrejects = nrejects + 1
                    end if
                end if
                last_send = candidatesol%fpen(1)
            else
                first_communication = 0
                send_solution_to_slaves = 1  
                percentage = 0.0
            end if
        end if
        
        
        if (send_solution_to_slaves .EQ. 1) then
                CALL printstatemaster( candidatesol%fpen(1),  cputime3 )
                diffusions = diffusions + 1
                cputime3 = calctimeMPI(exp1,starttime)
                CALL printgant(exp1,cputime3,1)  
                ocurrences_send(idsender+1) =  ocurrences_send(idsender+1) + 1
                lasttime = cputime3 
                CALL masterupdatesolutionslaves(exp1 , selectserializedata, sizeser,idsender, init)        

                cputime4 = calctimeMPI(exp1,starttime)
                cputime3 = cputime4 - cputime3
                CALL printputmasterendlog(exp1 , cputime3, cputime4, candidatesol%fpen(1), idsender)
                CALL printpercentage(exp1, cputime4, candidatesol%fpen(1), percentage )
                ! UPDATE 
                BKS_M_f=candidatesol%fpen
                BKS_M_x=candidatesol%x
            
                cputime3 = calctimeMPI(exp1,starttime)
                CALL printgant(exp1,cputime3,3)  
        end if

        ! ADAPTATION THRESHOLD
        CALL adaptation_threshold2( threshold_imp, nrejects, NPROC )
        
        ! ADAPTATION SETTINGS
        resultflag_adapt = 0
        CALL adaptationcheck(exp1, resultflag_adapt, idslave_to_adapt)
        
        if ( resultflag_adapt .GT. 0 )  then 
                 ALLOCATE(orderindexscore(size(slaveScore)))
                 ALLOCATE(slaveScore_aux(size(slaveScore)))
                 slaveScore(idslave_to_adapt+1) = 0d0

                 slaveScore_aux = slaveScore
                 call QsortC_ind(slaveScore_aux, orderindexscore)
                 bestid = orderindexscore(size(slaveScore))

                if ((counter_adapt_request(idslave_to_adapt+1) .EQ. 0)) then
                    counter_adapt = 0
                    do i=1,NPROC-1
                        if (counter_adapt_request(i) .NE. 0) then
                            counter_adapt = counter_adapt + 1
                        end if
                    end do
                    
                    if (( counter_adapt .LT. MAX_ADAPT) .AND. (slaveScore(bestid) .NE. 0 )) then
                        accept = 1
                    else
                        accept = 0
                    end if
                    CALL adaptationsendmaster(exp1, idslave_to_adapt+1, adapt_master(bestid)%size_dim_pen, &
                        adapt_master(bestid)%ncounter_pen, adapt_master(bestid)%balance_pen, accept)

                    if (accept .EQ. 1) then
                            CALL printadaptationscores(exp1, slaveScore, NPROC-1)
                            counter_adapt_request(idslave_to_adapt+1) = 1   
                            slaveComm(idslave_to_adapt+1)  = 0d0
                            adapt_master(idslave_to_adapt+1)%size_dim_pen = adapt_master(bestid)%size_dim_pen
                            adapt_master(idslave_to_adapt+1)%ncounter_pen = adapt_master(bestid)%ncounter_pen
                            adapt_master(idslave_to_adapt+1)%balance_pen =  adapt_master(bestid)%balance_pen
                            cputime3 = calctimeMPI(exp1,starttime)
                            CALL printadaptationmaster(exp1, idslave_to_adapt+1, bestid,adapt_master(bestid)%balance_pen, &
                            adapt_master(bestid)%size_dim_pen,adapt_master(bestid)%ncounter_pen, accept,resultflag_adapt, &
                                                        cputime3 )
                    end if

                           
                else 
                    accept = 1
                    slaveComm(idslave_to_adapt+1)  = 0d0
                    adapt_master(idslave_to_adapt+1)%size_dim_pen = adapt_master_origin(idslave_to_adapt+1)%size_dim_pen
                    adapt_master(idslave_to_adapt+1)%ncounter_pen = adapt_master_origin(idslave_to_adapt+1)%ncounter_pen
                    adapt_master(idslave_to_adapt+1)%balance_pen =  adapt_master_origin(idslave_to_adapt+1)%balance_pen

                    CALL adaptationsendmaster(exp1, idslave_to_adapt+1, adapt_master(idslave_to_adapt+1)%size_dim_pen, &
                        adapt_master(idslave_to_adapt+1)%ncounter_pen, adapt_master(idslave_to_adapt+1)%balance_pen, accept)
                    counter_adapt_request(idslave_to_adapt+1) = 0    
                    cputime3 = calctimeMPI(exp1,starttime)                    
                    CALL printdesadaptationmaster(exp1, idslave_to_adapt+1, idslave_to_adapt+1, &
                            adapt_master(idslave_to_adapt+1)%balance_pen, &
                            adapt_master(idslave_to_adapt+1)%size_dim_pen, adapt_master(idslave_to_adapt+1)%ncounter_pen, &
                                                        cputime3)                    
                end if
                DEALLOCATE(orderindexscore)
                DEALLOCATE(slaveScore_aux)
        end if 
    


        CALL destroy_refsettype(newsol)
        CALL destroy_refsettype(candidatesol)
        if (ALLOCATED(serializedata) ) DEALLOCATE(serializedata)
    
    end do

    if (ALLOCATED(adapt_master)) DEALLOCATE(adapt_master)
    if (ALLOCATED(adapt_master_origin)) DEALLOCATE(adapt_master_origin)

    DEALLOCATE(slaveComm)
    DEALLOCATE(slaveScore)    
    DEALLOCATE(counter_adapt_request)
    if (ALLOCATED(selectserializedata) ) DEALLOCATE(selectserializedata)
    DEALLOCATE(BKS_M_x)
    DEALLOCATE(BKS_M_f)

END SUBROUTINE asynchronous_master_acess



SUBROUTINE adapt_slave(exp1,opts1,problem1,fitnessfunction,common_vars,nfuneval, xl_log,xu_log,hyper_x_L,hyper_x_U, &
            MaxSubSet,MaxSubSet2,ppp,index1, index2, index,nconst,refset,fbest,xbest,refset_change,local_solver_var,time, &
            ncounter_slave_recv_sol, ncounter_slave_send_sol, pending_adaptation, last_evals)
    IMPLICIT NONE
    TYPE(time_ess), INTENT(INOUT) :: time
    TYPE(local_solver_help_vars) :: local_solver_var    
    TYPE(C_PTR),  INTENT(INOUT) :: exp1
    TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
    TYPE(opts), INTENT(INOUT) :: opts1 
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
    INTEGER(C_LONG), INTENT(INOUT) :: last_evals, nfuneval
    INTEGER, INTENT(INOUT) :: pending_adaptation
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: hyper_x_L, hyper_x_U
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: MaxSubSet, MaxSubSet2
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ppp
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: index1, index2, index      
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE, INTENT(INOUT)  :: xl_log,xu_log
    INTEGER, INTENT(IN) :: nconst
    TYPE(Refsettype), INTENT(INOUT) :: refset
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) ::  fbest, xbest
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
    INTEGER, INTENT(INOUT) :: ncounter_slave_recv_sol, ncounter_slave_send_sol
    INTEGER :: adaptflag, ncounter, size_dim, accept
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: balance
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: cputime3,  calctimeMPI
    INTEGER(C_LONG) :: threshold_adapt_evals
    INTEGER :: code,  mult_num_sendSol, minimum_num_sendSol
    REAL(KIND =C_DOUBLE) :: evals_threshold

    CALL setparallelsacessfieldsslaves(exp1, evals_threshold, mult_num_sendSol, minimum_num_sendSol)

    threshold_adapt_evals = INT(REAL(common_vars%nvar)*evals_threshold,KIND=C_LONG)
 
    ! SEND OF ADAPTATION SIGNAL
    if (pending_adaptation .EQ. 0) then
            if (ncounter_slave_recv_sol .GT. (( ncounter_slave_send_sol * mult_num_sendSol )+ minimum_num_sendSol)) then
                code = 1
                pending_adaptation = 1 
            end if

            if ( (nfuneval - last_evals) .GT. (threshold_adapt_evals) )  then
                code = 3
                pending_adaptation = 1 
            end if

            if (pending_adaptation .EQ. 1) then
                CALL sendslaveadaptsignal(exp1,code)
                cputime3 = calctimeMPI(exp1,time%starttime)
                CALL printadaptationslave(exp1, ncounter_slave_recv_sol, ncounter_slave_send_sol, (nfuneval - last_evals), cputime3)
            end if
    end if
    
    
    ! RECEPTION OF ADAPTATION DATA
    CALL checkadaptsettings(exp1, size_dim, ncounter, balance, adaptflag, accept)
    if (adaptflag .EQ. 1) then
        if (accept .EQ. 1 ) then
            last_evals = nfuneval
            ncounter_slave_recv_sol = 0
            ncounter_slave_send_sol = 0
            CALL adaptation_configuration( exp1,opts1,problem1,refset,common_vars, fbest,xbest, &
                    refset_change, local_solver_var, size_dim, balance, ncounter, &
                    index1, index2, index,MaxSubSet, MaxSubSet2, hyper_x_L, hyper_x_U, ppp, nconst, xl_log, xu_log, &
                    fitnessfunction, nfuneval)
            local_solver_var%n_minimo = 0
            local_solver_var%n_critico = opts1%localoptions%n2 
            opts1%localoptions%balance = balance
            cputime3 = calctimeMPI(exp1,time%starttime)
            CALL printadaptation(exp1, balance, size_dim, ncounter, cputime3)   
        end if
        pending_adaptation = 0
    end if
    
    
END SUBROUTINE adapt_slave

        
! ----------------------------------------------------------------------
! SUBROUTINES asynchronous_slave_cess
! ----------------------------------------------------------------------
SUBROUTINE asynchronous_slave_acess (exp1,opts1,problem1,refset,common_vars, nconst,BKS_M_x,BKS_M_f, &
    time,refset_change,fitnessfunction,xl_log,xu_log, nfuneval,  fbest, xbest, &
    ncounter_slave_recv_sol, ncounter_slave_send_sol, last_evals )
    IMPLICIT NONE
    TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
    TYPE(time_ess), INTENT(INOUT) :: time
    TYPE(opts), INTENT(INOUT) :: opts1 
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(C_PTR),  INTENT(INOUT) :: exp1
    TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
    INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE, INTENT(INOUT)  :: xl_log,xu_log
    INTEGER(C_LONG), INTENT(INOUT) :: last_evals
    INTEGER, INTENT(IN) :: nconst
    TYPE(Refsettype), INTENT(INOUT) :: refset
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: BKS_M_x 
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: BKS_M_f
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) ::  fbest, xbest 
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: serializedata, serializedatawindow
    TYPE(Refsettype) :: bestsol, candidate_from_Master
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE ::  dist_one
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: percentage
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE ::  candidate_from_Master_norm 
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: refset_norm, dist
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: f11, trans1, trans2
    INTEGER, DIMENSION(:), ALLOCATABLE :: indvect, indvect2
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
    INTEGER :: DIMEN(2), DIMEN2(2), send, rest, flagreturn
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind11,ind12,    indmin
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: cputime3,  calctimeMPI
    INTEGER :: i, position, max_coop(1)
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: threshold_imp
    INTEGER, INTENT(INOUT) :: ncounter_slave_recv_sol, ncounter_slave_send_sol
    INTEGER :: max_loc2(1)

    threshold_imp = 0.5d0
    
    
    CALL create_refset_empty ( bestsol, common_vars%nvar, opts1, nconst, 1 )
    CALL create_refset_empty ( candidate_from_Master, common_vars%nvar, opts1, nconst, 1 )
    ALLOCATE(serializedata(common_vars%sizeser, 1))
    ALLOCATE(serializedatawindow(common_vars%sizeser,1))

    ! SLAVE RECEIVE SOLUTION FROM MASTER   
    CALL returnwindowvalueslave(exp1, serializedatawindow, common_vars%sizeser, flagreturn)
    
    if (flagreturn .EQ. 1) then
        CALL deserializerefset(candidate_from_Master, common_vars%nvar, 1, serializedatawindow)  
        cputime3 = calctimeMPI(exp1,time%starttime)
        CALL printreceivedslave(exp1, candidate_from_Master%fpen(1), cputime3)    
        
        if ((candidate_from_Master%fpen(1) .LT.  fbest(1) ) .OR. ( candidate_from_Master%fpen(1) .LT. problem1%vtr(1) )) then
            cputime3 = calctimeMPI(exp1,time%starttime)       
            CALL calc_percentage( fbest(1) , candidate_from_Master%fpen(1) , percentage)           
            CALL printreceptioncompare( exp1, candidate_from_Master%fpen(1), fbest(1), percentage , cputime3, threshold_imp ) 
            cputime3 = calctimeMPI(exp1,time%starttime) 
            CALL printiterationcess(exp1, candidate_from_Master%fpen(1), nfuneval, cputime3, 0, 0)   
 
            ncounter_slave_recv_sol = ncounter_slave_recv_sol + 1           

            if ( common_vars%init_comm .EQ. 0 ) then
                CALL sort_refset(refset,opts1, indvect)
                refset_change = refset_change(indvect)

                refset%x(:,opts1%globaloptions%dim_refset) = candidate_from_Master%x(:,1)
                refset%f(opts1%globaloptions%dim_refset) = candidate_from_Master%f(1)
                refset%fpen(opts1%globaloptions%dim_refset) = candidate_from_Master%fpen(1)
                refset%penalty(opts1%globaloptions%dim_refset) = candidate_from_Master%penalty(1)
                if (nconst .GT. 0 ) then
                         refset%nlc(:,opts1%globaloptions%dim_refset) = candidate_from_Master%nlc(:,1)
                end if
                refset%cooperative(opts1%globaloptions%dim_refset ) = 1
                refset_change( opts1%globaloptions%dim_refset ) = -100
                position = opts1%globaloptions%dim_refset
                common_vars%init_comm = 1
            else 
                do i=1,opts1%globaloptions%dim_refset
                        if ( refset%cooperative(i) .EQ. 1  ) then
                          refset%x(:,i) = candidate_from_Master%x(:,1)
                          refset%f(i) = candidate_from_Master%f(1)
                          refset%fpen(i) = candidate_from_Master%fpen(1)
                          refset%penalty(i) = candidate_from_Master%penalty(1)
                          if (nconst .GT. 0 ) then
                                refset%nlc(:,i) = candidate_from_Master%nlc(:,1)
                          end if
                          refset_change( i ) = -100
                          position = i
                          EXIT
                       end if
               end do
            end if
! END NOVO
            cputime3 = calctimeMPI(exp1,time%starttime)
            CALL printreplaceslavelog(exp1, candidate_from_Master%fpen(1), cputime3, position )                 
        else
            CALL printdiscardreceivedslave(exp1, candidate_from_Master%fpen(1), cputime3,  fbest(1)) 
        end if
        BKS_M_x = candidate_from_Master%x
        BKS_M_f = candidate_from_Master%fpen        
        
    end if
! END SLAVE RECEIVE SOLUTION FROM MASTER    
      
! SLAVE SEND SOLUTIONS TO MASTER
    send = 0
   
    CALL sort_refset(refset,opts1, indvect2)
    refset_change = refset_change(indvect2)
 
    if (refset%fpen(1) .LT.  BKS_M_f(1) ) then
        cputime3 = calctimeMPI(exp1,time%starttime)        
        CALL calc_percentage( BKS_M_f(1) , refset%fpen(1) , percentage) 
        if (percentage .ge. threshold_imp ) then
            send = 1
            ncounter_slave_send_sol = ncounter_slave_send_sol + 1
        end if
        CALL printcomparenewsolutionslavelog( exp1, refset%fpen(1), BKS_M_f(1), percentage , cputime3, threshold_imp ) 
    end if
    
    if (send .EQ. 1) then
        cputime3 = calctimeMPI(exp1,time%starttime)
        CALL printgant(exp1,cputime3,1) 
        last_evals = nfuneval
        bestsol%x(:,1) = refset%x(:,1)
        bestsol%fpen(1) = refset%fpen(1) 
        bestsol%f(1) = refset%f(1)
        bestsol%penalty(1) =  refset%penalty(1)
        if (nconst .GT. 0 ) then
            bestsol%nlc(:,1) = refset%nlc(:,1)
        end if 

        CALL serializerefset(bestsol, common_vars%nvar, 1, serializedata, rest ) 
        CALL sendbestsolutionslave (exp1, serializedata, common_vars%sizeser, common_vars%inititer)
        cputime3 = calctimeMPI(exp1,time%starttime)
        CALL printfinalsendslavelog( exp1 , bestsol%fpen(1) , cputime3)
        CALL printgant(exp1,cputime3,3)        

    end if
    
    
    
! END SLAVE SEND SOLUTIONS TO MASTER      
    
    if (ALLOCATED(indvect) )  DEALLOCATE(indvect)
    if (ALLOCATED(indvect2) ) DEALLOCATE(indvect2)
    CALL destroy_refsettype(bestsol)
    CALL destroy_refsettype(candidate_from_Master)
        
    IF (ALLOCATED(serializedatawindow) ) DEALLOCATE(serializedatawindow)
    IF (ALLOCATED(serializedata)) DEALLOCATE(serializedata)

    
END SUBROUTINE asynchronous_slave_acess


SUBROUTINE calc_percentage( fbest, candidate, percentage) 
    IMPLICIT NONE
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: fbest, candidate
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(OUT) :: percentage
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: trans1, trans2
 
    if (((fbest .GT. 0 ) .AND. ( candidate .LT. 0 )) .OR. ((fbest .LT. 0 ) .AND. &
        (candidate .LT. 0) .AND. (fbest .GT. candidate))) then
        trans1=fbest+ABS(candidate)+1
        trans2=1
    else if ((( fbest .LT. 0  ) .AND. ( candidate .GT. 0))  &
        .OR. ((fbest .LT. 0 ) .AND. ( candidate .LT. 0) &
        .AND. (fbest .LT. candidate))) then
        trans1=1
        trans2=candidate+ABS(fbest)+1
    else
        trans1=fbest
        trans2=candidate
    end if
    
    percentage =((( trans1 - trans2 )  /  trans1 ) * 100d0 )
END SUBROUTINE calc_percentage
        

SUBROUTINE update_index_cooperative (aux_val_cooperative, refset, index_cooperative)
    IMPLICIT NONE
    CHARACTER(len = 2) :: typesearch
    INTEGER, DIMENSION(:), ALLOCATABLE :: indexfind
    TYPE(Refsettype), INTENT(INOUT) :: refset
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: aux_val_cooperative
    INTEGER,  INTENT(INOUT) :: index_cooperative
    
    typesearch = "LT"
    CALL find_in_vector(refset%fpen, aux_val_cooperative, indexfind, typesearch)
    CALL reajust_index(indexfind)
    index_cooperative = indexfind(1)
            
            
END SUBROUTINE   update_index_cooperative           





SUBROUTINE adaptation_configuration( exp1,opts1,problem1,refset,common_vars, fbest,xbest, &
    refset_change, local_solver_var, new_dim_refset, new_balance, new_counter_local, index1, index2, index,&
    MaxSubSet, MaxSubSet2, hyper_x_L, hyper_x_U, ppp, nconst, xl_log, xu_log, fitnessfunction, nfuneval)
    IMPLICIT NONE                             
    INTEGER, INTENT(IN) :: nconst
    TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_var
    TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
    TYPE(opts), INTENT(INOUT) :: opts1 
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(C_PTR),  INTENT(INOUT) :: exp1
    TYPE(Refsettype), INTENT(INOUT) :: refset
    TYPE(Refsettype) :: refset_new 
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: fbest, xbest
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
    INTEGER, DIMENSION(:), ALLOCATABLE :: refset_change_aux
    INTEGER, INTENT(INOUT) :: new_dim_refset, new_counter_local
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: new_balance
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: hyper_x_L, hyper_x_U
    INTEGER :: auxsize, dim_refsetOLD
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: MaxSubSet, MaxSubSet2
    INTEGER, DIMENSION(:), ALLOCATABLE :: ncomb1
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ncomb2
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: index1, index2, index
    INTEGER, DIMENSION(:), ALLOCATABLE :: diff_index    
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ppp
    INTEGER :: i,j
    INTEGER, DIMENSION(:), ALLOCATABLE :: indvect, indvect2, indvect3      
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: new_refset_member
    TYPE(outfuction) :: outfunct
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: xl_log, xu_log
    TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
    INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
    INTEGER :: end1, includedata, exist_coop, minimum(1)

    CALL sort_refset(refset,opts1, indvect)
    refset_change = refset_change(indvect) 
        
    dim_refsetOLD = opts1%globaloptions%dim_refset
    opts1%globaloptions%dim_refset = new_dim_refset
    opts1%localoptions%n2 = new_counter_local
    opts1%localoptions%balance = new_balance
   
    ! RECALC ANY VARS:
    if (ALLOCATED(index1) ) DEALLOCATE(index1)     
    if (ALLOCATED(index2) ) DEALLOCATE(index2)     
    if (ALLOCATED(index) ) DEALLOCATE(index)     
    if (ALLOCATED(diff_index) ) DEALLOCATE(diff_index) 
    if (ALLOCATED(ppp) ) DEALLOCATE(ppp)
    ALLOCATE(refset_change_aux(new_dim_refset))
    refset_change_aux = 0   
 
    if (ALLOCATED(hyper_x_L) ) DEALLOCATE(hyper_x_L)
    if (ALLOCATED(hyper_x_U) ) DEALLOCATE(hyper_x_U)    

    ALLOCATE(ncomb1(opts1%globaloptions%dim_refset))
    ncomb1 = (/(i, i = 1, opts1%globaloptions%dim_refset)/)
    call nchoosek_fpar(ncomb1, 2, ncomb2)
    MaxSubSet = (REAL(opts1%globaloptions%dim_refset)**2d0 - REAL(opts1%globaloptions%dim_refset) )/2d0
    MaxSubSet2 = 2d0 * MaxSubSet  

    CALL calc_indexes(ncomb1,ncomb2,index1,index2,index,diff_index,auxsize)
    CALL calc_ppp(auxsize,ppp,opts1%globaloptions%dim_refset,diff_index)
    CALL calc_hyper(common_vars%nvar, hyper_x_L, hyper_x_U, MaxSubSet, problem1%XU, problem1%XL )    

    ! REBUILD REFSET
    if ( new_dim_refset .LT. dim_refsetOLD) then      
        end1 = new_dim_refset
    else if ( new_dim_refset .GT. dim_refsetOLD) then
        end1 = dim_refsetOLD
    else
        end1 = new_dim_refset
    end if

    CALL create_refset_empty( refset_new, common_vars%nvar, opts1, nconst, new_dim_refset )
    refset_new%x(:,1:end1) = refset%x(:,1:end1)
    refset_new%f(1:end1) = refset%f(1:end1)
    refset_new%fpen(1:end1) = refset%fpen(1:end1)
    refset_new%penalty(1:end1) = refset%penalty(1:end1)
    if (nconst .GT. 0) then
        refset_new%nlc(:,1:end1) = refset%nlc(:,1:end1)
    end if
    refset_new%cooperative(1:end1) = refset%cooperative(1:end1)
    refset_change_aux(1:end1) = refset_change(1:end1)

    do i=end1+1,new_dim_refset
                outfunct%include = 0
                do while ( outfunct%include .EQ. 0 )
                        ALLOCATE(new_refset_member(common_vars%nvar))
                        CALL random_Vector_SEED(exp1,new_refset_member, common_vars%nvar, xl_log, xu_log)
                        if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0)  &
                                CALL converttonormal(new_refset_member, common_vars%nvar)
                        outfunct = ssm_evalfc(exp1, fitnessfunction, new_refset_member, problem1, opts1, nconst, 0)
                        nfuneval = nfuneval + 1
                        if (outfunct%include .eq. 1) then
                                refset_new%x(:,i) = outfunct%x
                                refset_new%f(i) = outfunct%value
                                refset_new%fpen(i) = outfunct%value_penalty
                                if (nconst .GT. 0) then 
                                        refset_new%nlc(:,i)=outfunct%nlc
                                end if
                                refset_new%penalty(i) = outfunct%pena
                        end if
                        CALL destroy_outfuction(outfunct)
                        DEALLOCATE(new_refset_member)
                end do
    end do
    minimum = minloc(refset_new%fpen)
    if ( refset_new%fpen(minimum(1)) .GT. fbest(1)) then
         outfunct = ssm_evalfc(exp1, fitnessfunction,xbest, problem1, opts1, nconst, 0)
         nfuneval = nfuneval + 1
         refset_new%x(:,new_dim_refset) = outfunct%x
         refset_new%f(new_dim_refset) = outfunct%value
         refset_new%fpen(new_dim_refset) = outfunct%value_penalty
         if (nconst .GT. 0)  then
                refset_new%nlc(:,new_dim_refset)=outfunct%nlc
         end if
         refset_new%penalty(new_dim_refset) = outfunct%pena
         CALL destroy_outfuction(outfunct)
    end if
    
    CALL destroy_refsettype(refset)
    CALL create_refset_empty( refset, common_vars%nvar, opts1, nconst, new_dim_refset )
    refset%x = refset_new%x
    refset%f = refset_new%f
    refset%fpen = refset_new%fpen
    refset%penalty = refset_new%penalty
    if (nconst .GT. 0) then
            refset%nlc = refset_new%nlc
    end if
    refset%cooperative = refset_new%cooperative
    exist_coop = 1
    do i=1,new_dim_refset
        if (refset%cooperative(i) .EQ. 1) exist_coop = 0
    end do  
    if (exist_coop .EQ. 1) then 
            CALL sort_refset(refset,opts1, indvect3)
            refset%cooperative(new_dim_refset) = 1            
    end if
    DEALLOCATE(refset_change)
    ALLOCATE(refset_change(new_dim_refset)) 
    refset_change = refset_change_aux

    CALL destroy_refsettype(refset_new)
     
    if (ALLOCATED(refset_change_aux)) DEALLOCATE(refset_change_aux)
    if (ALLOCATED(indvect3) ) DEALLOCATE(indvect3)
    if (ALLOCATED(indvect)  ) DEALLOCATE(indvect)
END SUBROUTINE adaptation_configuration




SUBROUTINE adaptation_threshold2 ( threshold, nrejects, nslaves )
    IMPLICIT NONE
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: threshold
    INTEGER, INTENT(INOUT) :: nrejects
    INTEGER, INTENT(IN) :: nslaves
    INTEGER :: diff
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: threshold_min


    threshold_min = 0.5d0
    

    if ( nrejects .GE. (nslaves-1) ) then
        threshold = threshold/2d0    
        nrejects = 0
    end if
    
    if ( threshold .LE. threshold_min ) then
        threshold = threshold_min
    end if
    
END SUBROUTINE adaptation_threshold2





SUBROUTINE hete_param_eSS (exp1, problem1, opts1, NPROC, dim_refset, idp,nvar)
      TYPE(problem), INTENT(INOUT) :: problem1
      TYPE(opts), INTENT(INOUT) :: opts1
      TYPE(C_PTR),  INTENT(INOUT) :: exp1
      INTEGER, INTENT(IN) :: NPROC, nvar
      INTEGER, INTENT(INOUT) :: dim_refset
      INTEGER :: ishete, coopera, hete, iscoop, SMALL_DIM_FREC_MIG,MEDIUM_DIM_FREC_MIG, BIG_DIM_FREC_MIG
      INTEGER, INTENT(IN) :: idp
      INTEGER :: idb
      REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: nnn(2),nn(1),polim(3), dim, nvar2


      if ( idp .GE. 0  )  then

           if ( MOD(idp,10) .EQ. 0 )  then
                dim = 1
                opts1%localoptions%n1 = 1
                opts1%localoptions%n2 = 1
                opts1%localoptions%balance = 0.0
            else if ( MOD(idp,10) .EQ. 1 )  then
                dim = 3
                opts1%localoptions%n1 = 4
                opts1%localoptions%n2 = 4
                opts1%localoptions%balance =  0.0
            else if ( MOD(idp,10) .EQ. 2 )  then
                dim = 5
                opts1%localoptions%n1 = 10
                opts1%localoptions%n2 = 10
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 3 )  then
                dim = 10
                opts1%localoptions%n1 = 20
                opts1%localoptions%n2 = 20
                opts1%localoptions%balance = 0.5
            else if ( MOD(idp,10) .EQ. 4 )  then
                dim = 15
                opts1%localoptions%n1 = 100
                opts1%localoptions%n2 = 100
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 5 )  then
                dim = 12
                opts1%localoptions%n1 = 50
                opts1%localoptions%n2 = 50
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 6 )  then
                dim = 7.5
                opts1%localoptions%n1 = 15
                opts1%localoptions%n2 = 15
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 7 )  then
                dim = 5
                opts1%localoptions%n1 = 7
                opts1%localoptions%n2 = 7
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 8 )  then
                dim = 2
                opts1%localoptions%n1 = 2
                opts1%localoptions%n2 = 2
                opts1%localoptions%balance = 0.0
            else if ( MOD(idp,10) .EQ. 9 ) then
                dim = 0.5
                opts1%localoptions%n1 = 1
                opts1%localoptions%n2 = 1
                opts1%localoptions%balance =  0.0
            end if

         else
                dim = 10d0
                opts1%localoptions%n1 = 10
                opts1%localoptions%n2 = 10
                opts1%localoptions%balance = 0.5d0
         end if
         CALL getbench(exp1,idb)
         if ( idb .EQ. 3 ) then
             nvar2= REAL(nvar,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/50d0
         else
             nvar2=nvar
         end if

         polim = (/  REAL(1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                    REAL(-1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
            (-1d0) *dim * REAL( nvar2/1d0, KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) /)

         CALL SQUARE_ROOTS_POLINOM(polim, nnn)
         nn = maxval(nnn, MASK = nnn .GT. 0)
         dim_refset = CEILING(nn(1))
        if ( dim_refset .LT. 3 ) dim_refset=3
         opts1%globaloptions%dim_refset = dim_refset


         CALL updatenp(exp1, dim_refset)
 
    END SUBROUTINE

! ----------------------------------------------------------------------
! SUBROUTINES hete_param_eSS
! ----------------------------------------------------------------------
SUBROUTINE hete_param_eSS2 (exp1, problem1, opts1, NPROC, dim_refset, idp2, nvar )
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(opts), INTENT(INOUT) :: opts1 
    TYPE(C_PTR),  INTENT(INOUT) :: exp1
    INTEGER, INTENT(IN) :: NPROC, idp2, nvar
    INTEGER, INTENT(INOUT) :: dim_refset
    INTEGER :: ishete, coopera, hete, iscoop, SMALL_DIM_FREC_MIG, MEDIUM_DIM_FREC_MIG, BIG_DIM_FREC_MIG
    INTEGER :: idp, idb
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: nnn(2), nn(1),polim(3), dim, nvar2
    

            
       idp = idp2
       idp = idp - 1
       if ( idp .GE. 0  )  then     
           if ( MOD(idp,10) .EQ. 0 )  then
                dim = 1
                opts1%localoptions%n1 = 1
                opts1%localoptions%n2 = 1
                opts1%localoptions%balance = 0.0
            else if ( MOD(idp,10) .EQ. 1 )  then
                dim = 3
                opts1%localoptions%n1 = 4
                opts1%localoptions%n2 = 4
                opts1%localoptions%balance =  0.0
            else if ( MOD(idp,10) .EQ. 2 )  then
                dim = 5
                opts1%localoptions%n1 = 10
                opts1%localoptions%n2 = 10
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 3 )  then
                dim = 10
                opts1%localoptions%n1 = 20
                opts1%localoptions%n2 = 20
                opts1%localoptions%balance = 0.5
            else if ( MOD(idp,10) .EQ. 4 )  then
                dim = 15
                opts1%localoptions%n1 = 100
                opts1%localoptions%n2 = 100
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 5 )  then
                dim = 12
                opts1%localoptions%n1 = 50
                opts1%localoptions%n2 = 50
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 6 )  then
                dim = 7.5
                opts1%localoptions%n1 = 15
                opts1%localoptions%n2 = 15
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 7 )  then
                dim = 5
                opts1%localoptions%n1 = 7
                opts1%localoptions%n2 = 7
                opts1%localoptions%balance = 0.25
            else if ( MOD(idp,10) .EQ. 8 )  then
                dim = 2
                opts1%localoptions%n1 = 2
                opts1%localoptions%n2 = 2
                opts1%localoptions%balance = 0.0
            else if ( MOD(idp,10) .EQ. 9 ) then
                dim = 0.5
                opts1%localoptions%n1 = 1
                opts1%localoptions%n2 = 1
                opts1%localoptions%balance =  0.0
            end if
        else
                dim = 10d0
                opts1%localoptions%n2 = 10
                opts1%localoptions%balance = 0.5d0
        end if
        
         CALL getbench(exp1,idb)
         if ( idb .EQ. 3 ) then
             nvar2= REAL(nvar,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/50d0
         else
             nvar2=nvar
         end if

         polim = (/  REAL(1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                    REAL(-1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
           (-1d0) *dim * REAL( nvar2/1d0, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) /)                
           CALL SQUARE_ROOTS_POLINOM(polim, nnn)
           nn = maxval(nnn, MASK = nnn .GT. 0)
           dim_refset = CEILING(nn(1))          
        if ( dim_refset .LT. 3 ) dim_refset=3   
        opts1%globaloptions%dim_refset = dim_refset

           
            CALL updatenp(exp1, dim_refset)

    
    
END SUBROUTINE







#ifdef OPENMP
! ----------------------------------------------------------------------
! SUBROUTINES evaluate_solutions_set_parallel
! ----------------------------------------------------------------------   
   SUBROUTINE evaluate_solutions_set_parallel(exp1,fitnessfunction,solutionset,problem1,opts1,solutions, nvar, nfuneval,nconst, &
                     xl_log, xu_log)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        TYPE(Refsettype), INTENT(INOUT) :: solutionset
        TYPE(Refsettype) :: solutionset_aux
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(opts), INTENT(IN) :: opts1
        INTEGER, INTENT(IN) :: nconst
        REAL(C_DOUBLE) ::  cputime1, cputime2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: solutions
        INTEGER, INTENT(IN)  :: nvar
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER :: l_f_0, l_x_0, DIMEN(2), tempnum, counter, counter2, i, sizet, neval
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: newsolution, entervalues
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp2, temp3
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: NAN, calctimeMPI
        INTEGER(C_INT) :: error, setNAN, continueeval
        INTEGER :: OMP_NUM_THREADS, IDPROC, clocal
        TYPE(outfuction), DIMENSION(:), ALLOCATABLE ::  outfunct1_VECTOR
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        REAL :: SIZE_MIN
        IDPROC = 0
        error = setNAN(NAN)
        
        if (ALLOCATED(problem1%F0)) then
            l_f_0 = size(problem1%F0)
        else
            l_f_0 = 0
        end if

        if (ALLOCATED(problem1%X0)) then
            DIMEN = shape(problem1%X0)
            l_x_0 = DIMEN(2)
        else
            l_x_0 = 0
        end if

        ! Put x_0 without their corresponding problem1%F0 values at the
        ! beginning of solutions
        !print *, "l_x_0", l_x_0, "-","l_f_0", l_f_0
        if ( (l_x_0  - l_f_0) > 0) then
            DIMEN = shape(problem1%X0)
            if (ALLOCATED(temp2)) DEALLOCATE(temp2)
            ALLOCATE(temp2(DIMEN(1), DIMEN(2)- l_f_0))
            temp2 = problem1%X0(:,(l_f_0 + 1):DIMEN(2))
            CALL FUSION_MATRIX(temp2, solutions)
            tempnum = l_x_0 - l_f_0
            !DEALLOCATE(temp2)
        else
            tempnum = 0
        end if
         
        sizet = opts1%globaloptions%ndiverse + 5 + tempnum

        counter = 0
        neval = 0 
        IDPROC = -1
        clocal = 0
        ALLOCATE(outfunct1_VECTOR(sizet))
        
        OMP_NUM_THREADS = omp_get_max_threads()
           
        SIZE_MIN = CEILING(REAL(opts1%globaloptions%dim_refset)/REAL(OMP_NUM_THREADS))
 !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) FIRSTPRIVATE(IDPROC,clocal) PRIVATE(newsolution,i,continueeval)
 !  REDUCTION(+:neval,counter) DEFAULT(SHARED) 
        do i = 1, sizet
            continueeval = 0
            if (IDPROC .EQ. -1 ) IDPROC = OMP_GET_THREAD_NUM()
            if ( .not. ALLOCATED(newsolution)) ALLOCATE(newsolution(nvar))
            newsolution = solutions(:,i)
            do while (continueeval .EQ. 0)
                outfunct1_VECTOR(i) = ssm_evalfc(exp1, fitnessfunction, newsolution, problem1, opts1, nconst, IDPROC)
                neval = neval + 1

                if (outfunct1_VECTOR(i)%include .eq. 1) then
                        counter = counter + 1
                        clocal = clocal + 1
                        continueeval = 1
                else if ( (outfunct1_VECTOR(i)%include .eq. 0) .AND. &
                  ( REAL(clocal) .LT. SIZE_MIN)) then

                    CALL random_Vector_SEED(exp1,newsolution, nvar, xl_log, xu_log)   
                    if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0 )  &
                           CALL converttonormal2(newsolution, nvar, opts1%useroptions%log_var)
                    call destroy_outfuction(outfunct1_VECTOR(i))
                else
                    continueeval = 1
                end if
            end do
            if (ALLOCATED(newsolution)) DEALLOCATE(newsolution)
        end do
 !$OMP END PARALLEL DO
        
        ALLOCATE(solutionset%x(nvar, counter))
        solutionset%x =  REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        if (nconst .GT. 0 ) then 
            ALLOCATE(solutionset%nlc(nconst , counter))
            solutionset%nlc = REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) 
        end if
        ALLOCATE(solutionset%f(counter))
        solutionset%f =  REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(solutionset%fpen(counter))
        solutionset%fpen = REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(solutionset%penalty(counter))
        solutionset%penalty = REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        counter2 = 1
        do i = 1, counter
            if (outfunct1_VECTOR(i)%include .eq. 1) then
                solutionset%x(:,counter2) = outfunct1_VECTOR(i)%x
                solutionset%f(counter2) = outfunct1_VECTOR(i)%value
                solutionset%fpen(counter2) = outfunct1_VECTOR(i)%value_penalty
                solutionset%penalty(counter2) = outfunct1_VECTOR(i)%pena
                if (nconst .GT. 0) then 
                    solutionset%nlc(:,counter2) = outfunct1_VECTOR(i)%nlc
                end if
                counter2 = counter2 + 1
            end if
            call destroy_outfuction(outfunct1_VECTOR(i))
        end do
        DEALLOCATE(outfunct1_VECTOR)
        nfuneval = nfuneval + neval

        !resize solutionset
        counter = counter2-1
        if ( counter .LT. sizet )  then
                sizet = counter
                ALLOCATE(solutionset_aux%x(nvar, sizet))
                solutionset_aux%x(:,1:counter)=solutionset%x(:,1:counter)
                if (nconst .GT. 0 ) then
                        ALLOCATE(solutionset_aux%nlc(nconst , sizet))
                        solutionset_aux%nlc(:,1:counter)=solutionset%nlc(:,1:counter)
                end if
                ALLOCATE(solutionset_aux%f(sizet))
                solutionset_aux%f(1:counter)=solutionset%f(1:counter)
                ALLOCATE(solutionset_aux%fpen(sizet))
                solutionset_aux%fpen(1:counter)=solutionset%fpen(1:counter)
                ALLOCATE(solutionset_aux%penalty(sizet))
                solutionset_aux%penalty(1:counter)=solutionset%penalty(1:counter)

                DEALLOCATE(solutionset%x)
                ALLOCATE(solutionset%x(nvar, sizet))
                if (nconst .GT. 0 ) then
                        DEALLOCATE(solutionset%nlc)
                        ALLOCATE(solutionset%nlc(nconst , sizet))
                end if
                DEALLOCATE(solutionset%f)
                ALLOCATE(solutionset%f(sizet))
                DEALLOCATE(solutionset%fpen)
                ALLOCATE(solutionset%fpen(sizet))
                DEALLOCATE(solutionset%penalty)
                ALLOCATE(solutionset%penalty(sizet))

                solutionset%x = solutionset_aux%x
                if (nconst .GT. 0 ) then
                        solutionset%nlc=solutionset_aux%nlc
                end if
                solutionset%f=solutionset_aux%f
                solutionset%fpen=solutionset_aux%fpen
                solutionset%penalty=solutionset_aux%penalty

                DEALLOCATE(solutionset_aux%x)
                if (nconst .GT. 0 ) then
                        DEALLOCATE(solutionset_aux%nlc)
                end if
                DEALLOCATE(solutionset_aux%f)
                DEALLOCATE(solutionset_aux%fpen)
                DEALLOCATE(solutionset_aux%penalty)

        end if
            

        !Add points with f_0 values. We assume feasiblity
        if (l_f_0 .GT. 0) then
            DIMEN = shape(problem1%X0)
            if (ALLOCATED(temp2)) DEALLOCATE(temp2)
            ALLOCATE(temp2(DIMEN(1),l_f_0))
            temp2 = problem1%X0(:,1:l_f_0)
            call FUSION_MATRIX(temp2, solutionset%x)
            call FUSION_VECTOR(problem1%F0, solutionset%f)
            call FUSION_VECTOR(problem1%F0, solutionset%fpen)
            ALLOCATE(entervalues(l_f_0))
            entervalues = entervalues * 0
            call FUSION_VECTOR(entervalues, solutionset%penalty)
            DEALLOCATE(entervalues)
            
            if (ALLOCATED(solutionset%nlc)) then
                ALLOCATE(temp3(nconst,l_f_0))
                temp3=0
                CALL FUSION_MATRIX(temp3,solutionset%nlc)
            end if
        end if    
        
        
   END SUBROUTINE evaluate_solutions_set_parallel  
    
    
! ----------------------------------------------------------------------
! SUBROUTINES apply_beyond_to_members_to_update_parallel
! ----------------------------------------------------------------------   
    SUBROUTINE apply_beyond_to_members_to_update_parallel( exp1, opts1, fitnessfunction, problem1, members_to_update,  nfuneval, &
        refset, candidateset, nconst, nrand, nvar, timeparallel )
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1 
        TYPE(Refsettype), INTENT(INOUT) :: refset, candidateset
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj 
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER :: clock_rate, clock_start, clock_stop        
        INTEGER(C_LONG) :: neval
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  ::  members_to_update
        INTEGER, INTENT(IN) :: nconst, nvar , nrand
        REAL(C_DOUBLE), INTENT(INOUT) :: timeparallel
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp1, temp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: auxiliarvalue
        INTEGER :: OMP_NUM_THREADS, IDPROC
        TYPE(outfuction), DIMENSION(:), ALLOCATABLE :: outfunct1
        INTEGER :: i 
#ifdef TIMEPAR      
            CALL SYSTEM_CLOCK(count_rate=clock_rate)
            CALL SYSTEM_CLOCK(count=clock_start)
#endif              
        IDPROC = 0
!$OMP PARALLEL SHARED(OMP_NUM_THREADS)           
           OMP_NUM_THREADS = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL    
           

               
       ALLOCATE(outfunct1(OMP_NUM_THREADS)) 
       ALLOCATE(temp(nvar))
       ALLOCATE(temp1(nvar))     

            
       if (ALLOCATED(members_to_update)) then
                ! exploramos los miembros a actualizar (los que prometen más resultados)
                ! y aplicamos ls técnica go-beyond en ellos
            
 !$OMP PARALLEL PRIVATE(OMP_NUM_THREADS,IDPROC) DEFAULT(SHARED)
                IDPROC = OMP_GET_THREAD_NUM()
                neval = 0
 !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(temp,temp1,i,auxiliarvalue) REDUCTION(+:neval)
                do i = 1, size(members_to_update)

                    temp = refset%x(:,members_to_update(i))
                    temp1 = candidateset%x(:,members_to_update(i))
                    auxiliarvalue = candidateset%fpen(members_to_update(i))

                    ! ssm_beyond
                    CALL ssm_beyond(exp1, temp, temp1, auxiliarvalue, neval, &
                    fitnessfunction, nrand, nconst, outfunct1(IDPROC+1), opts1, problem1, IDPROC)
                          
                             
                    if (ALLOCATED(outfunct1(IDPROC+1)%x) ) then  
                        refset%x(:,members_to_update(i)) = outfunct1(IDPROC+1)%x
                        refset%f(members_to_update(i)) = outfunct1(IDPROC+1)%value
                        refset%fpen(members_to_update(i)) = outfunct1(IDPROC+1)%value_penalty
                        refset%penalty(members_to_update(i)) = outfunct1(IDPROC+1)%pena
                        if (nconst .GT. 0) refset%nlc(:,members_to_update(i)) = outfunct1(IDPROC+1) % nlc
                        !end if
                    else    
                        refset%x(:,members_to_update(i)) = candidateset%x(:,members_to_update(i))
                        refset%f(members_to_update(i)) = candidateset%f(members_to_update(i))
                        refset%fpen(members_to_update(i)) = candidateset%fpen(members_to_update(i))
                        refset%penalty(members_to_update(i)) = candidateset%penalty(members_to_update(i))
                        if (nconst .GT. 0) refset%nlc(:,members_to_update(i)) = candidateset%nlc(:,members_to_update(i))
                    end if
                    
                end do
 !$OMP END DO

 !$OMP END PARALLEL   
        
#ifdef TIMEPAR            
            CALL SYSTEM_CLOCK(count=clock_stop)
            timeparallel = timeparallel + REAL(clock_stop-clock_start)/REAL(clock_rate)
#endif                
       nfuneval = nfuneval + neval
       end if
            
            DEALLOCATE(temp)
            DEALLOCATE(temp1)
            call destroy_outfuction(outfunct1(IDPROC+1))
            if (ALLOCATED(outfunct1)) DEALLOCATE(outfunct1)    
            
    END SUBROUTINE apply_beyond_to_members_to_update_parallel
    
    
 
! ----------------------------------------------------------------------
! SUBROUTINES update_candidateset_with_new_comb_parallel
! ----------------------------------------------------------------------   
    SUBROUTINE update_candidateset_with_new_comb_parallel( exp1,opts1, fitnessfunction,problem1,new_comb,candidateset,childset, &
            candidate_update, members_update,MaxSubSet2,nvar,nconst,nfuneval,counter, index3 )
            IMPLICIT NONE
            TYPE(problem), INTENT(IN) :: problem1
            TYPE(C_PTR),  INTENT(INOUT) :: exp1
            TYPE(Refsettype), INTENT(INOUT) :: candidateset, childset
            TYPE(opts), INTENT(IN) :: opts1 
            TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj        
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: members_update, candidate_update
            INTEGER(C_LONG), INTENT(INOUT) ::  nfuneval
            INTEGER, INTENT(INOUT) :: counter
            INTEGER :: neval    
            INTEGER, INTENT(IN) :: nconst, nvar
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: MaxSubSet2
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  new_comb        
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: index3
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: newsolution
            INTEGER :: i, OMP_NUM_THREADS, IDPROC, ii
            TYPE(outfuction), DIMENSION(:), ALLOCATABLE :: outfunct1, outfunct1_VECTOR
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: cputime1,cputime2,calctimeMPI
            
            IDPROC = 0
            
       
            
            OMP_NUM_THREADS = omp_get_max_threads()

            ALLOCATE(outfunct1(OMP_NUM_THREADS)) 
            ALLOCATE(outfunct1_VECTOR(CEILING(MaxSubSet2)))

            neval = 0
 
 !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(newsolution,i,IDPROC)
 !REDUCTION(+:neval)
            do i = 1, CEILING(MaxSubSet2)
                ! Evaluate
                IDPROC = OMP_GET_THREAD_NUM()
                ALLOCATE(newsolution(nvar))
                newsolution = new_comb(:,i)
                
                
                outfunct1(IDPROC+1) = ssm_evalfc(exp1, fitnessfunction, newsolution, problem1, opts1, nconst,IDPROC)
                neval = neval + 1 
                outfunct1_VECTOR(i) = outfunct1(IDPROC+1)

                call destroy_outfuction(outfunct1(IDPROC+1))
                if (ALLOCATED(newsolution))  DEALLOCATE(newsolution)                
            end do
 !$OMP END PARALLEL DO
            
            if (ALLOCATED(outfunct1)) DEALLOCATE(outfunct1)
            counter = 1
            do i = 1, CEILING(MaxSubSet2)
                if (outfunct1_VECTOR(i)%include .eq. 1) then
                    childset%x(:,counter)=outfunct1_VECTOR(i)%x
                    childset%f(counter)= outfunct1_VECTOR(i)%value
                    childset%fpen(counter) = outfunct1_VECTOR(i)%value_penalty
                    childset%penalty(counter) = outfunct1_VECTOR(i)%pena
                    if (ALLOCATED(outfunct1_VECTOR(i)%nlc) ) childset%nlc(:,counter) = outfunct1_VECTOR(i)%nlc
                    childset%parent_index(counter) = index3(i)
                    counter = counter + 1
                    if (outfunct1_VECTOR(i)%value_penalty .LT. candidateset%fpen(index3(i))) then
                            candidateset%x(:,index3(i)) = outfunct1_VECTOR(i)%x
                            candidateset%fpen(index3(i)) = outfunct1_VECTOR(i)%value_penalty
                            candidateset%f(index3(i)) = outfunct1_VECTOR(i)%value
                            candidateset%penalty(index3(i)) = outfunct1_VECTOR(i)%pena
                            if (ALLOCATED(outfunct1_VECTOR(i)%nlc)) candidateset%nlc(:,index3(i)) = outfunct1_VECTOR(i)%nlc
                            candidate_update(index3(i)) = 1
                            members_update(index3(i)) = 0
                    end if

                end if
                call destroy_outfuction(outfunct1_VECTOR(i))              
            end do            
            
            CALL reajustchild(childset, counter, nvar, nconst)
            if (ALLOCATED(outfunct1_VECTOR)) DEALLOCATE(outfunct1_VECTOR)
            

            nfuneval = nfuneval + neval
            
            
                
    END SUBROUTINE update_candidateset_with_new_comb_parallel   

#endif



! ----------------------------------------------------------------------    
! SUBROUTINES create_refset_local
! ----------------------------------------------------------------------    
   SUBROUTINE create_refset_local ( exp1, problem1, opts1, refsetglobal, refset, nvar, nconst, tam )
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(Refsettype), INTENT(INOUT) :: refsetglobal, refset
        TYPE(Refsettype) :: refsetlocal
        INTEGER, INTENT(IN) :: nvar, nconst, tam
        INTEGER :: cond, iscoop, cond2, ishete
        TYPE(problem), INTENT(IN) :: problem1

        
        ! Initialize Refset
        ALLOCATE(refset%x(nvar,tam))
        refset%x = refset%x * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%f(tam))
        refset%f = refset%f * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%fpen(tam))
        refset%fpen = refset%fpen *REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%penalty(tam))
        refset%penalty = refset%penalty * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        if (nconst .GT. 0) then
            ALLOCATE(refset%nlc(nconst,tam))
            refset%nlc = refset%nlc * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))              
        end if
             
        
       
    

            CALL cooperativebcastelement(exp1, refsetglobal%x, nvar*opts1%globaloptions%dim_refset)
            CALL cooperativebcastelement(exp1, refsetglobal%f, opts1%globaloptions%dim_refset);
            CALL cooperativebcastelement(exp1, refsetglobal%fpen, opts1%globaloptions%dim_refset);
            CALL cooperativebcastelement(exp1, refsetglobal%penalty, opts1%globaloptions%dim_refset);     
            if (nconst .GT. 0) then
                CALL cooperativebcastelement(exp1, refsetglobal%nlc, opts1%globaloptions%dim_refset*nconst); 
            end if  
                
            refset%x =  refsetglobal%x(:,1:opts1%globaloptions%dim_refset)
            refset%f =  refsetglobal%f(1:opts1%globaloptions%dim_refset)
            refset%fpen =  refsetglobal%fpen(1:opts1%globaloptions%dim_refset)
            refset%penalty =  refsetglobal%penalty(1:opts1%globaloptions%dim_refset)
            if (nconst .GT. 0) then
                refset%nlc =  refsetglobal%nlc(:,1:opts1%globaloptions%dim_refset)
            end if              

        
        opts1%globaloptions%dim_refset = tam
        
        CALL destroy_refsettype(refsetlocal)
   END SUBROUTINE create_refset_local   
   
   


! ----------------------------------------------------------------------
! SUBROUTINES add_vect_proc
! ----------------------------------------------------------------------
SUBROUTINE add_vect_proc(vector_proc, id) 
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector_proc
        INTEGER, INTENT(INOUT) :: id
        INTEGER ,DIMENSION(:), ALLOCATABLE :: auxvect
        INTEGER :: DIM1, i, contador

        DIM1 = size( vector_proc )
        contador = 2
        
        ALLOCATE(auxvect(DIM1))

        do i=1, DIM1
              if  (vector_proc(i) .NE. id ) then 
                      auxvect(contador) =vector_proc(i)
                      contador = contador + 1
              else
                  auxvect(1) = id
              end if
        end do

        CALL MOVE_ALLOC( auxvect , vector_proc)
        

END SUBROUTINE


! ----------------------------------------------------------------------
! SUBROUTINES restart_acess
! ----------------------------------------------------------------------
SUBROUTINE restart_acess(exp1,problem1,opts1,fitnessfunction,refset,nvar,xl_log,xu_log,xbest,fbest,nconst,timeparallel,nfuneval,&
    refset_change, stage_1, stage_2, use_bestx) 
        TYPE(opts), INTENT(IN) :: opts1
        TYPE(problem), INTENT(INOUT) :: problem1
        TYPE(C_PTR),  INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        INTEGER, INTENT(IN) :: nvar,nconst
        INTEGER :: NPROC
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(Refsettype) :: solutionset
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
         
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:),ALLOCATABLE :: solutions
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: xl_log,xu_log
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: timeparallel
        REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: fbest, xbest       
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
        INTEGER, INTENT(INOUT) :: stage_1, stage_2, use_bestx
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect       
        

        CALL setNPROC(exp1, NPROC)
        CALL create_init_solutions(exp1,opts1,solutions,nvar,xl_log,xu_log)
        if (ALLOCATED(problem1%X0)) DEALLOCATE(problem1%X0)
        if (ALLOCATED(problem1%F0)) DEALLOCATE(problem1%F0)
        ALLOCATE(problem1%X0(nvar, opts1%globaloptions%dim_refset+1))
        ALLOCATE(problem1%F0( opts1%globaloptions%dim_refset+1  ))

        problem1%X0(:,1:opts1%globaloptions%dim_refset) = refset%x
        problem1%F0(1:opts1%globaloptions%dim_refset) = refset%fpen
        problem1%X0(:,opts1%globaloptions%dim_refset+1) = xbest
        problem1%F0(opts1%globaloptions%dim_refset+1) = fbest(1)


        if ((problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 1)) then
             CALL ssm_round_int(solutions, problem1%int_var + problem1%bin_var, problem1%XL, problem1%XU)
        end if

        CALL destroy_refsettype(refset)

        !CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1,solutions, nvar, &
        !            nfuneval, nconst)

       
        CALL create_refset ( exp1, solutionset, refset, nvar, opts1, nconst, opts1%globaloptions%dim_refset )
        refset_change = 0
         
        stage_1 = 1
        stage_2 = 0
        use_bestx = 0

        if (ALLOCATED(indvect)) DEALLOCATE(indvect)
        CALL sort_refset(refset,opts1, indvect)
        fbest=refset%fpen(1)
        xbest=refset%x(:,1)
        
        DEALLOCATE(indvect)

    END SUBROUTINE   
       
! ----------------------------------------------------------------------
! SUBROUTINES check_duplicated_replace
! ----------------------------------------------------------------------     
    SUBROUTINE unique_solution (problem1, XVAR, FVAR, nvar )
            TYPE(problem), INTENT(IN) :: problem1
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: XVAR
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: FVAR
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE:: XVAR_L
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: FVAR_L            
            INTEGER, INTENT(IN) :: nvar
            INTEGER, DIMENSION(:), ALLOCATABLE :: index_delete
            INTEGER :: i, j, contador, contadorlocal
            
            ALLOCATE(index_delete( size(FVAR) ))
            
            contador=0
            do i=1,size(FVAR)
                contadorlocal = 0
                if ( (i+1) .NE. size(FVAR) ) then
                    do j=i+1,size(FVAR)
                        if ( FVAR(i) .NE. FVAR(j) ) then
                            contadorlocal = contadorlocal + 1
                        end if
                    end do
                
                    if ( (size(FVAR) - i ) .EQ. contadorlocal) then
                        index_delete(i) = 0
                    else
                        index_delete(i) = 1
                        contador = contador + 1
                    end if
                
                else 
                    index_delete(i) = 0
                end if
            end do
            
            ALLOCATE(XVAR_L(nvar, size(index_delete) - contador))
            ALLOCATE(FVAR_L(size(index_delete) - contador))
            
            contador = 1
            
            do i=1,size(FVAR)
                if ( index_delete(i) .NE. 1 ) then
                    XVAR_L(:,contador) = XVAR(:,i)
                    FVAR_L(contador) = FVAR(i)
                    contador = contador + 1                    
                end if
            end do
            
            
            CALL MOVE_ALLOC(XVAR_L,XVAR)
            CALL MOVE_ALLOC(FVAR_L,FVAR)
            
            DEALLOCATE(index_delete)
            
    END SUBROUTINE unique_solution   
    
   
! ----------------------------------------------------------------------
! SUBROUTINES initoldbestvars
! ----------------------------------------------------------------------  
    SUBROUTINE initoldbestvars(oldbestx,oldfbest,nvar) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nvar
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: oldbestx
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: oldfbest     
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: INF

        CALL setdblmax(INF)

        ALLOCATE(oldbestx(nvar,1))
        ALLOCATE(oldfbest(1))
        oldbestx=0
        oldfbest=INF
    
    END SUBROUTINE initoldbestvars  
    
    
! ----------------------------------------------------------------------
! SUBROUTINES create_migration_master
! ----------------------------------------------------------------------  
    SUBROUTINE create_migration_master(migration_master,NPROC) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NPROC
        INTEGER :: i
        TYPE(master_migration), INTENT(INOUT) :: migration_master
        
        ALLOCATE(migration_master%vector_proc(NPROC-1))
        do i=1,NPROC-1
                migration_master%vector_proc(i) = i
        end do 
        
        ALLOCATE(migration_master%ocurrences_send(NPROC))
        migration_master%ocurrences_send=0
    
    END SUBROUTINE create_migration_master     
   
    
! ----------------------------------------------------------------------
! SUBROUTINES create_migration_master
! ----------------------------------------------------------------------  
    SUBROUTINE destroy_migration_master(migration_master) 
        IMPLICIT NONE
        TYPE(master_migration), INTENT(INOUT) :: migration_master
        
        IF (ALLOCATED(migration_master%ocurrences_send))  DEALLOCATE(migration_master%ocurrences_send)
        IF (ALLOCATED(migration_master%vector_proc)) DEALLOCATE(migration_master%vector_proc)
    
    END SUBROUTINE destroy_migration_master     
    
    
    
    

! ----------------------------------------------------------------------
! SUBROUTINES asynchronous_acess_dist
! ----------------------------------------------------------------------
SUBROUTINE asynchronous_acess_dist(exp1,opts1,problem1,refset,sizeser,nvar,nconst,oldbestx,oldfbest,&
                                                        fbest,xbest,init,idp, starttime, refset_change)
    TYPE(opts), INTENT(IN) :: opts1
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(C_PTR),  INTENT(INOUT) :: exp1
    INTEGER, INTENT(IN) :: sizeser, nconst, nvar, idp
    INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(INOUT) :: init
    TYPE(Refsettype), INTENT(INOUT) :: refset
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE, INTENT(INOUT) :: oldbestx
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: fbest, &
                                                                                          oldfbest,xbest
    REAL(C_DOUBLE), INTENT(INOUT) :: starttime
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE :: serializedata, serializedatawindow
    TYPE(Refsettype) :: bestsol, windowsol
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE :: auxbestx
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE ::  auxfbest, resultcomp
    INTEGER, DIMENSION(:), ALLOCATABLE :: refset_change
    INTEGER, DIMENSION(:), ALLOCATABLE :: indvect, indvect2
    INTEGER :: DIMEN(2), DIMEN2(2), send, rest, flagreturn, idsent, sent
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::cputime3,calctimeMPI


    CALL create_refset_empty ( bestsol, nvar, opts1, nconst, 1 )
    CALL create_refset_empty ( windowsol, nvar, opts1, nconst, 1 )
    ALLOCATE(serializedata(sizeser, 1))
    ALLOCATE(serializedatawindow(sizeser,1))
    ALLOCATE(auxbestx(nvar))
    ALLOCATE(auxfbest(1))


    ! miramos si existe novo valor, e si existe o recepcionamos
    CALL returnwindowvaluesdist(exp1, serializedatawindow, sizeser, flagreturn, idsent )
    if (flagreturn .EQ. 1) then
        CALL deserializerefset(windowsol, nvar, 1, serializedatawindow)
        cputime3 = calctimeMPI(exp1,starttime)
        CALL printreturnselectsolution(exp1, windowsol%fpen(1),oldfbest(1),cputime3, idsent)
        if (windowsol%fpen(1) .LT. oldfbest(1) ) then
          cputime3 = calctimeMPI(exp1,starttime)
          CALL printreplaceslavelog(exp1, windowsol%fpen(1), cputime3, opts1%globaloptions%dim_refset )
          CALL sort_refset(refset,opts1, indvect)
          refset_change = refset_change(indvect)
          refset%x(:,opts1%globaloptions%dim_refset) = windowsol%x(:,1)
          refset%f(opts1%globaloptions%dim_refset) = windowsol%f(1)
          refset%fpen(opts1%globaloptions%dim_refset) = windowsol%fpen(1)
          refset%penalty(opts1%globaloptions%dim_refset) = windowsol%penalty(1)
          if (nconst .GT. 0 ) then
            refset%nlc(:,opts1%globaloptions%dim_refset) =  windowsol%nlc(:,1)
          end if
          refset_change(opts1%globaloptions%dim_refset)=0
        end if
    end if

    CALL sort_refset(refset,opts1, indvect2)
    refset_change = refset_change(indvect2)
    send = 0
    ALLOCATE(resultcomp(1))

    ! Miramos se temos dato para enviar
    if (refset%f(1) .LT.  oldfbest(1) ) then
        cputime3 = calctimeMPI(exp1,starttime)
        resultcomp = (( ( oldfbest -  refset%f )  /  oldfbest ) * 100d0)
        CALL printcomparenewsolutionslavelog( exp1, auxfbest(1),oldfbest(1),resultcomp(1) , cputime3, resultcomp )

        if (resultcomp(1) .LT. 0.5d0 ) then
            send = 0
        else
            send = 1
            oldfbest(1) = refset%f(1)
            oldbestx(:,1) = refset%x(:,1)
        end if
    end if

    if (send .EQ. 1) then
        bestsol%x(:,1) = refset%x(:,1)
        bestsol%fpen(1) = refset%fpen(1)
        bestsol%f(1) = refset%f(1)
        bestsol%penalty(1) = refset%penalty(1)
        if (nconst .GT. 0 ) then
            bestsol%nlc(:,1) = refset%nlc(:,1)
        end if
        cputime3 = calctimeMPI(exp1,starttime)
        CALL serializerefset(bestsol, nvar, 1, serializedata, rest )
        CALL sendbestsolutiondist (exp1, serializedata, sizeser, init, sent)
        cputime3 = calctimeMPI(exp1,starttime)
        CALL printfinalsendslavelog( exp1 , bestsol%fpen(1) , cputime3)
    end if

    if (ALLOCATED(indvect) )  DEALLOCATE(indvect)
    if (ALLOCATED(indvect2) ) DEALLOCATE(indvect2)
    if (ALLOCATED(resultcomp)) DEALLOCATE(resultcomp)
    CALL destroy_refsettype(bestsol)
    CALL destroy_refsettype(windowsol)
    IF (ALLOCATED(auxfbest)) DEALLOCATE(auxfbest)
    IF (ALLOCATED(serializedatawindow) ) DEALLOCATE(serializedatawindow)
    IF (ALLOCATED(serializedata)) DEALLOCATE(serializedata)
    IF (ALLOCATED(auxbestx)) DEALLOCATE(auxbestx)


END SUBROUTINE




SUBROUTINE init_adapt_vars(exp1,adapt_master, NPROC, nvar)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: exp1
    TYPE(adapt_vars), DIMENSION(:), ALLOCATABLE :: adapt_master
    INTEGER, INTENT(IN) :: nvar
    INTEGER :: nvar2, idb
    REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: nnn(2),nn(1),polim(3), dim
    INTEGER :: i, idp, dim_refset
    INTEGER, INTENT(IN) :: NPROC
    
    ALLOCATE(adapt_master(NPROC-1))
    do i=1, NPROC-1
        idp = i-1
        if ( i .GE. 0  )  then
            
            if ( MOD(idp,10) .EQ. 0 )  then
                dim = 1
                adapt_master(i)%ncounter_pen = 1
                adapt_master(i)%balance_pen = 0.0
            else if ( MOD(idp,10) .EQ. 1 )  then
                dim = 3
                adapt_master(i)%ncounter_pen = 4
                adapt_master(i)%balance_pen =  0.0
            else if ( MOD(idp,10) .EQ. 2 )  then
                dim = 5
                adapt_master(i)%ncounter_pen = 10
                adapt_master(i)%balance_pen = 0.25
            else if ( MOD(idp,10) .EQ. 3 )  then
                dim = 10
                adapt_master(i)%ncounter_pen = 20
                adapt_master(i)%balance_pen = 0.5
            else if ( MOD(idp,10) .EQ. 4 )  then
                dim = 15
                adapt_master(i)%ncounter_pen = 100
                adapt_master(i)%balance_pen = 0.25
            else if ( MOD(idp,10) .EQ. 5 )  then
                dim = 12
                adapt_master(i)%ncounter_pen = 50
                adapt_master(i)%balance_pen = 0.25
            else if ( MOD(idp,10) .EQ. 6 )  then
                dim = 7.5
                adapt_master(i)%ncounter_pen = 15
                adapt_master(i)%balance_pen = 0.25
            else if ( MOD(idp,10) .EQ. 7 )  then
                dim = 5
                adapt_master(i)%ncounter_pen = 7
                adapt_master(i)%balance_pen = 0.25
            else if ( MOD(idp,10) .EQ. 8 )  then
                dim = 2
                adapt_master(i)%ncounter_pen = 2
                adapt_master(i)%balance_pen = 0.0
            else if ( MOD(idp,10) .EQ. 9 ) then
                dim = 0.5
                adapt_master(i)%ncounter_pen = 1
                adapt_master(i)%balance_pen =  0.0
            end if

         else
                dim = 10d0
                adapt_master(i)%ncounter_pen = 10
                adapt_master(i)%balance_pen = 0.5d0
         end if

         CALL getbench(exp1,idb)

         if ( idb .EQ. 3 ) then
             nvar2= REAL(nvar,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/50d0
         else
             nvar2=nvar
         end if

         polim = (/  REAL(1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                    REAL(-1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
           (-1d0) *dim * REAL( nvar2/1d0,KIND= SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) /)
         
         CALL SQUARE_ROOTS_POLINOM(polim, nnn)
         nn = maxval(nnn, MASK = nnn .GT. 0)
         dim_refset = CEILING(nn(1))
         if ( dim_refset .LT. 3 ) dim_refset=3 
         adapt_master(i)%size_dim_pen = dim_refset
         adapt_master(i)%pending = 0
    end do

    
END SUBROUTINE

#endif

END MODULE parallelscattersearchfunctions

