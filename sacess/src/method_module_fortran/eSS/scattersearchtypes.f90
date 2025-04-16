MODULE scattersearchtypes
    USE iso_c_binding
#ifdef MPI
        USE MPI
#endif

    INTEGER, PARAMETER :: PRECISION_D = 15
    INTEGER, PARAMETER :: RANGE_D = 307
   
    INTEGER, PARAMETER :: PRECISION_F = 7
    INTEGER, PARAMETER :: RANGE_F = 150  
    
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), PARAMETER :: TOL_PREC = &
            REAL(1d-18,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  
    
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), PARAMETER :: eps_m = &
            REAL(1d-18,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  
            
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), PARAMETER :: AUXPRECISIN = &
            REAL(10000d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  
       
    TYPE :: adapt_vars
        INTEGER :: size_dim_pen
        INTEGER :: ncounter_pen
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: balance_pen
        INTEGER :: pending
    END TYPE
            
    TYPE :: local_solver_help_vars
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: local_solutions
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: initial_points
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: local_solutions_values
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: initial_points_values
        INTEGER :: evals_per_iter
        INTEGER :: n_minimo
        INTEGER :: n_critico
        INTEGER :: stage_1
        INTEGER :: stage_2
        INTEGER :: use_bestx
        INTEGER, DIMENSION(:), ALLOCATABLE  :: N_UPPER
        INTEGER, DIMENSION(:), ALLOCATABLE  :: N_LOWER
        TYPE(C_FUNPTR) ::  fobj_mysqp
        INTEGER :: pass_misqp
    END TYPE local_solver_help_vars
    
    
    TYPE :: ess_common_vars
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xl_log
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xu_log

    END TYPE ess_common_vars
                
    TYPE :: algorithm_common_vars
        INTEGER :: NPROC
        INTEGER :: idp
        INTEGER :: parallel
        INTEGER :: cooperative
        INTEGER :: inititer
        INTEGER :: init_comm
        INTEGER :: sizeser ! SIZE OF SERIALIZATION DATA
        INTEGER :: iter
        INTEGER :: status
        INTEGER :: nvar
    END TYPE algorithm_common_vars    
            
    TYPE :: time_ess
        REAL(C_DOUBLE) :: timeparallel 
        REAL(C_DOUBLE) :: localsolvertime
        REAL(C_DOUBLE) :: starttime
        INTEGER :: clock_rate
        INTEGER :: clock_start
        INTEGER :: clock_stop
    END TYPE time_ess
            
    TYPE :: master_migration
        INTEGER(C_INT), DIMENSION(:), ALLOCATABLE :: vector_proc
        INTEGER, DIMENSION(:), ALLOCATABLE :: ocurrences_send
    END TYPE master_migration
    
    TYPE :: problem
        INTEGER :: empty 
        CHARACTER(len = 20) :: f
        INTEGER :: neq
        INTEGER :: int_var
        INTEGER :: bin_var
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE  :: vtr
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: XL
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: XU
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: CL
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: CU
        INTEGER :: ineq
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: X0
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: F0
    END TYPE problem

    TYPE :: uoptions
        INTEGER :: empty 
        INTEGER :: iterprint
        INTEGER :: plot
        INTEGER :: weight
        INTEGER :: nstuck_solution
        INTEGER :: strategy
        INTEGER :: inter_save        
        INTEGER(C_LONG_LONG) :: maxeval
        REAL(C_DOUBLE) :: maxtime
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: tolc
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: prob_bound
        INTEGER, DIMENSION(:), ALLOCATABLE :: log_var    
        INTEGER :: init_point
    END TYPE uoptions

    TYPE :: goptions
        INTEGER :: empty 
        INTEGER :: dim_refset ! -1 si es auto
        INTEGER :: ndiverse ! -1 si es auto
        INTEGER :: initiate
        INTEGER :: combination
        INTEGER :: regenerate
        INTEGER :: intens
        INTEGER :: diverse_criteria
        INTEGER :: n_stuck
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: tolf
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: tolx
        CHARACTER(len = 20,kind=C_CHAR) :: delete
    END TYPE goptions

    TYPE :: extraparametters
        INTEGER :: empty
        INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: texp ! nf2kb parametter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: yexp ! nf2kb parametter
        REAL(C_FLOAT) :: k1 ! solnp parametter
        REAL(C_FLOAT) :: k2 ! solnp parametter
        REAL(C_FLOAT) :: k3 ! solnp parametter
        REAL(C_FLOAT) :: k4 ! solnp parametter
    END TYPE extraparametters

    TYPE :: loptions
        INTEGER :: empty 
        INTEGER :: tol
        INTEGER :: iterprint
        INTEGER :: n1
        INTEGER :: n2
        INTEGER :: bestx
        INTEGER :: merit_filter
        INTEGER :: distance_filter
        INTEGER :: wait_maxdist_limit
        INTEGER :: wait_th_limit
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: threshold_local
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: thfactor
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: maxdistfactor
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: balance
        TYPE(extraparametters) :: extrap        
        CHARACTER(len = 20,kind=C_CHAR) :: solver
        CHARACTER(len = 20,kind=C_CHAR), DIMENSION(:), ALLOCATABLE :: finish
    END TYPE loptions

    TYPE :: opts
        INTEGER :: empty
        TYPE(uoptions) :: useroptions
        TYPE(goptions) :: globaloptions
        TYPE(loptions) :: localoptions
    END TYPE opts

    TYPE :: Refsettype
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: x
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: f
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: fpen
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: nlc
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: penalty
        INTEGER, DIMENSION(:), ALLOCATABLE :: parent_index
        INTEGER, DIMENSION(:), ALLOCATABLE :: cooperative
    END TYPE Refsettype
    
    TYPE :: resultsf
        TYPE(Refsettype) :: refsett
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: fbest
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xbest
        REAL(C_DOUBLE)  :: cputime
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: f
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: x
        REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: time
        INTEGER(C_LONG), DIMENSION(:), ALLOCATABLE :: neval
        INTEGER(C_LONG) :: numeval
        INTEGER(C_LONG_LONG) :: local_solutions
        INTEGER(C_LONG_LONG), DIMENSION(:), ALLOCATABLE :: local_solutions_values
        INTEGER(C_LONG_LONG) :: end_crit
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: totaliter 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: timevtr
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: totalvtr
        REAL(C_DOUBLE) :: timetotal
    END TYPE resultsf
    
        

    TYPE :: outfuction
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: value
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: value_penalty
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: pena
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: nlc
        INTEGER :: include
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: x
    END TYPE outfuction
    

    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), PARAMETER :: INFINITE = 1.7976931348623158D308
    
    
    TYPE :: topology_data_fortran
        INTEGER :: num_r
        INTEGER :: num_left
        INTEGER, DIMENSION(:), ALLOCATABLE :: rigth
        INTEGER, DIMENSION(:), ALLOCATABLE :: left
#ifdef MPI2
!        MPI_COMM :: comunicator
#endif
    END TYPE
    

END MODULE scattersearchtypes
