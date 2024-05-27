module Solver_MA87Module
    implicit none
    ! Lapack, BLAS, OpemMP and Metis are needed in some routines/functions
    interface SparseSystemSolver
        module procedure RealSparseSystemSolver_BVector
        module procedure RealSparseSystemSolver_BMatrix
        module procedure DPSparseSystemSolver_BMatrix
        module procedure DPSparseSystemSolver_BVector
    end interface SparseSystemSolver
    contains
    ! ----------------- ALL FUNCTIONS AND SUBROUTINES -----------------
    ! 2. Sparse System Solver
    ! 2.1. Real system (B Vector)
    function RealSparseSystemSolver_BVector(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma87_single
        use hsl_mc68_single
        use hsl_mc69_single
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma87_keep)                                                 :: keep
        type(ma87_control)                                              :: control
        type(ma87_info)                                                 :: info
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: RowVectorA,ColVectorA
        real, dimension(:), allocatable, intent(in)                     :: ValueVectorA, ValueVectorB
        real, dimension(:), allocatable                                 :: val,x
        ! Sparse Solver
        n = size(ValueVectorB)
        ne = size(ValueVectorA)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma87_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor
        call ma87_factor(n, ptr, row, val, order, keep, control, info)
        call stop_on_bad_flag("factor", info%flag)
        ! Solve
        call ma87_solve(x, order, keep, control, info)
        call stop_on_bad_flag("solve", info%flag)
        ! Finalize
        call ma87_finalise(keep, control)
    end function RealSparseSystemSolver_BVector
    ! 2.2. Real system (B Matrix)
    function RealSparseSystemSolver_BMatrix(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma87_single
        use hsl_mc68_single
        use hsl_mc69_single
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma87_keep)                                                 :: keep
        type(ma87_control)                                              :: control
        type(ma87_info)                                                 :: info
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: RowVectorA,ColVectorA 
        real, dimension(:), allocatable, intent(in)                     :: ValueVectorA
        real, dimension(:,:), allocatable, intent(in)                   :: ValueVectorB
        real, dimension(:), allocatable                                 :: val
        real, dimension(:,:), allocatable                               :: X
        ! Sparse Solver
        n = size(ValueVectorB,1)
        ne = size(ValueVectorA)
        nrhs = size(ValueVectorB,2)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma87_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor and solve
        call ma87_factor_solve(n, ptr, row, val, order, keep, control, info, nrhs, n, X)
        call stop_on_bad_flag("factor", info%flag)
        ! Finalize
        call ma87_finalise(keep, control)
    end function RealSparseSystemSolver_BMatrix
    ! 2.3. Double precision system (B Vector)
    function DPSparseSystemSolver_BVector(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma87_double
        use hsl_mc68_double
        use hsl_mc69_double
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma87_keep)                                                 :: keep
        type(ma87_control)                                              :: control
        type(ma87_info)                                                 :: info
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: RowVectorA,ColVectorA
        double precision, dimension(:), allocatable, intent(in)         :: ValueVectorA, ValueVectorB
        double precision, dimension(:), allocatable                     :: val,X
        ! Sparse Solver
        n = size(ValueVectorB)
        ne = size(ValueVectorA)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma87_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor
        call ma87_factor(n, ptr, row, val, order, keep, control, info)
        call stop_on_bad_flag("factor", info%flag)
        ! Solve
        call ma87_solve(X, order, keep, control, info)
        call stop_on_bad_flag("solve", info%flag)
        ! Finalize
        call ma87_finalise(keep, control)
    end function DPSparseSystemSolver_BVector
    ! 2.4. Real system (B Matrix)
    function DPSparseSystemSolver_BMatrix(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma87_double
        use hsl_mc68_double
        use hsl_mc69_double
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma87_keep)                                                 :: keep
        type(ma87_control)                                              :: control
        type(ma87_info)                                                 :: info
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: ColVectorA, RowVectorA
        double precision, dimension(:), allocatable, intent(in)         :: ValueVectorA
        double precision, dimension(:,:), allocatable, intent(in)       :: ValueVectorB
        double precision, dimension(:), allocatable                     :: val
        double precision, dimension(:,:), allocatable                   :: X
        ! Sparse Solver
        n = size(ValueVectorB,1)
        ne = size(ValueVectorA)
        nrhs = size(ValueVectorB,2)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma87_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor and solve
        call ma87_factor_solve(n, ptr, row, val, order, keep, control, info, nrhs, n, X)
        call stop_on_bad_flag("factor", info%flag)
        ! Finalize
        call ma87_finalise(keep, control)
    end function DPSparseSystemSolver_BMatrix
    ! Warring signal
    subroutine stop_on_bad_flag(context, flag)
        character(len=*), intent(in)                                 :: context
        integer, intent(in)                                          :: flag
        if(flag.eq.0) return
        write(*,*) "Failure during ", context, " with flag = ", flag
        stop
    end subroutine stop_on_bad_flag
end module Solver_MA87Module