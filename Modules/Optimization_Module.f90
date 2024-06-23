module Optimization_module
    use MMA_Module
    use FEA_Module
    implicit none
    type, extends(Structure)                                    :: Optimization
        character(len=3)                                        :: OptAlgorithm
        integer                                                 :: Iteration
        integer                                                 :: MaxTOIteration
        double precision                                        :: MinDensityValue
        double precision                                        :: PostProcesingFilter
        double precision                                        :: FilterRadius
        double precision                                        :: PenalFactor
        double precision                                        :: VolFraction
        double precision                                        :: MutationRate
        double precision                                        :: OptimizationTolerance
        double precision                                        :: Change
        double precision, dimension(:), allocatable             :: DensityVector    
        double precision, dimension(:), allocatable             :: FinalDensityVector
        double precision, dimension(:), allocatable             :: Compilance
        double precision, dimension(:), allocatable             :: DiffVolume
        double precision, dimension(:), allocatable             :: DiffCompilance
        double precision, dimension(:), allocatable             :: VolumePerElement
        ! Additional for OC
        double precision                                        :: L1
        double precision                                        :: L2
        ! Additional for MMA
        integer                                                 :: mm
        integer                                                 :: nn
        double precision, dimension(:,:), allocatable           :: xmin
        double precision, dimension(:,:), allocatable           :: xmax
        double precision, dimension(:,:), allocatable           :: XMMA
        double precision, dimension(:,:), allocatable           :: xold1
        double precision, dimension(:,:), allocatable           :: xold2
        double precision, dimension(:,:), allocatable           :: low
        double precision, dimension(:,:), allocatable           :: upp
        double precision                                        :: f0val
        double precision, dimension(:,:), allocatable           :: df0dx
        double precision, dimension(:,:), allocatable           :: fval
        double precision, dimension(:,:), allocatable           :: dfdx
        double precision                                        :: a0 
        double precision, dimension(:,:), allocatable           :: a 
        double precision, dimension(:,:), allocatable           :: c 
        double precision, dimension(:,:), allocatable           :: d 
        ! Post-processing
        integer, dimension(:), allocatable                      :: FinalGeometry
    contains
        procedure                                               :: MMAParameters
        procedure                                               :: SetMinDensityValue
        procedure                                               :: SetMaxIterations
        procedure                                               :: SetPostProcesingFilter
        procedure                                               :: SetFilterRadius
        procedure                                               :: SetPenalFactor
        procedure                                               :: SetVolFraction
        procedure                                               :: SetMutationRate
        procedure                                               :: UploadOptimizationParameters
        procedure                                               :: TopologyOptimizationProcess
    end type 
    contains
    ! ----------------------------------------------------------------- !
    !       subroutines to define the information required for TOP      !
    ! ----------------------------------------------------------------- !
    ! Additional routines for MMA Algorithm
    ! 1. setting MMA parameters
    subroutine MMAParameters(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision                                             :: xminval
        double precision                                             :: xmaxval
        xmaxval = 1.0d0
        xminval = 0.001d0
        Self%mm = 1
        Self%nn = Self%Ne
        allocate(Self%xmin(Self%nn,1))
        allocate(Self%xmax(Self%nn,1))
        allocate(Self%XMMA(Self%nn,1))
        allocate(Self%xold1(Self%nn,1))
        allocate(Self%xold2(Self%nn,1))
        allocate(Self%low(Self%nn,1))
        allocate(Self%upp(Self%nn,1))
        allocate(Self%fval(Self%mm,1))
        allocate(Self%dfdx(Self%mm,Self%nn))
        allocate(Self%df0dx(Self%nn,1))
        allocate(Self%a(Self%mm,1),Self%c(Self%mm,1),Self%d(Self%mm,1))
        Self%xmax = xmaxval
        Self%xmin = xminval
        Self%xold1(:,1) = Self%DensityVector
        Self%xold2(:,1) = Self%DensityVector
        Self%upp = 1.0d0
        Self%low = 1.0d0
        Self%a0 = 1.0d0
        Self%a = 0.0d0
        Self%c = 100.0d0    !100.0d0
        Self%d = 0.0d0
    end subroutine MMAParameters
    ! topology optimization subroutines
    ! 2. Optimization algoritm
    subroutine SetOptimizationAlgorithm(Self,OptAlgorithm)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        character(len=3)                                             :: OptAlgorithm
        Self%OptAlgorithm = OptAlgorithm
    end subroutine SetOptimizationAlgorithm
    ! 3. Minimun Density Value
    subroutine SetMinDensityValue(Self,MinDensityValue)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: MinDensityValue
        Self%MinDensityValue = MinDensityValue
    end subroutine SetMinDensityValue
    ! 4. Input max iterations
    subroutine SetMaxIterations(Self,MaxTOIteration)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        integer, intent(in)                                          :: MaxTOIteration
        Self%MaxTOIteration = MaxTOIteration
    end subroutine SetMaxIterations
    ! 5. Input Post-processing filter
    subroutine SetPostProcesingFilter(Self,PostProcesingFilter)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: PostProcesingFilter
        Self%PostProcesingFilter = PostProcesingFilter
    end subroutine SetPostProcesingFilter
    ! 6. Input filter radius
    subroutine SetFilterRadius(Self,FilterRadius)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: FilterRadius
        Self%FilterRadius = FilterRadius
    end subroutine SetFilterRadius
    ! 7. Input penal factor
    subroutine SetPenalFactor(Self,PenalFactor)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: PenalFactor
        Self%PenalFactor = PenalFactor
    end subroutine SetPenalFactor
    ! 8. Input max volumne fraction
    subroutine SetVolFraction(Self,VolFraction)
        implicit none
        double precision, intent(in)                                 :: VolFraction
        class(Optimization), intent(inout)                           :: Self
        Self%VolFraction = VolFraction
    end subroutine SetVolFraction
    ! 9. Input mutation/changing rate
    subroutine SetMutationRate(Self,Mutation)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: Mutation
        Self%MutationRate = Mutation
    end subroutine SetMutationRate
    ! 10. Input optimization tolerance
    subroutine SetOptimizationTolerance(Self,OptimizationTolerance)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: OptimizationTolerance
        Self%OptimizationTolerance = OptimizationTolerance
    end subroutine SetOptimizationTolerance
    ! 11. Updating Optimization parameters
    subroutine UploadOptimizationParameters(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        if (Self%DimAnalysis.eq.2) then
            Self%L1 = 0.0d0
            Self%L2 = 100000.0d0
        elseif(Self%DimAnalysis.eq.3) then
            Self%L1 = 0.0d0
            Self%L2 = 1000000000.0d0
        end if
        Self%Change = 1.0d0
        allocate(Self%DensityVector(Self%Ne));
        Self%DensityVector = 1.0d0
        allocate(Self%Compilance(Self%Ne));
        Self%Compilance = 0.0d0
        allocate(Self%DiffCompilance(Self%Ne));
        Self%DiffCompilance = 0.0d0
        allocate(Self%DiffVolume(Self%Ne));
        ! Volume of all elements
        call VolumeAllElement(Self)
        ! Assign the MMA Params.
        call MMAParameters(Self)
    end subroutine UploadOptimizationParameters
    ! 12. Printing TO Convergence
    Subroutine PrintingConvergence(Iteration,Self)
        implicit none
        integer, intent(in)                                          :: Iteration
        class(Optimization), intent(inout)                           :: Self
        ! Printing code
        write(unit=*, fmt=*) 'Ite,', Iteration,'PenalFactor',Self%PenalFactor,'ObjFunction,',Sum(Self%Compilance),'Change', &
                                Self%Change,'FracVolume,',Sum(Self%VolumePerElement*Self%DensityVector)/Sum(Self%VolumePerElement)
    end subroutine PrintingConvergence
    ! 13. Volume Calculation
    subroutine VolumeAllElement(Self)
        implicit none
        class(Optimization), intent(inout)                              :: Self
        integer                                                      :: i
        double precision                                             :: Area
        double precision, dimension(:), allocatable                  :: v1,v2,v3,v4,Centroid
        double precision, dimension(:,:), allocatable                :: Coordinates
        allocate(Self%VolumePerElement(Self%Ne))
        Self%VolumePerElement = 0.0d0
        do i = 1, Self%Ne, 1
            Area = 0.0d0
            if ((Self%ElementType.eq.'tria3').or.(Self%ElementType.eq.'tria6')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = [Coordinates(3,:) - Coordinates(2,:),0.0d0]
                v2 = [Coordinates(1,:) - Coordinates(2,:),0.0d0]
                Area = norm2(CrossProduct(v1,v2))/2.0d0
                Self%VolumePerElement(i) = Area*Self%Thickness
            elseif ((Self%ElementType.eq.'quad4').or.(Self%ElementType.eq.'quad8')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = [Coordinates(2,:) - Coordinates(1,:),0.0d0]
                v2 = [Coordinates(4,:) - Coordinates(1,:),0.0d0]
                v3 = [Coordinates(4,:) - Coordinates(3,:),0.0d0]
                v4 = [Coordinates(2,:) - Coordinates(3,:),0.0d0]
                Area = norm2(CrossProduct(v1,v2))/2.0d0 + norm2(CrossProduct(v3,v4))/2.0d0
                Self%VolumePerElement(i) = Area*Self%Thickness
            elseif ((Self%ElementType.eq.'tetra4').or.(Self%ElementType.eq.'tetra10')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = Coordinates(3,:) - Coordinates(1,:)
                v2 = Coordinates(2,:) - Coordinates(1,:)
                v3 = Coordinates(4,:) - Coordinates(1,:)
                Self%VolumePerElement(i) = abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
            elseif ((Self%ElementType.eq.'hexa8').or.(Self%ElementType.eq.'hexa20')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                Centroid = sum(Coordinates,1)/size(Coordinates(:,1))
                ! 1
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(2,:) - Centroid
                v3 = Coordinates(5,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 2
                v1 = Coordinates(2,:) - Centroid
                v2 = Coordinates(5,:) - Centroid
                v3 = Coordinates(6,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 3
                v1 = Coordinates(3,:) - Centroid
                v2 = Coordinates(7,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 4
                v1 = Coordinates(3,:) - Centroid
                v2 = Coordinates(4,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 5
                v1 = Coordinates(2,:) - Centroid
                v2 = Coordinates(6,:) - Centroid
                v3 = Coordinates(7,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 6
                v1 = Coordinates(2,:) - Centroid
                v2 = Coordinates(3,:) - Centroid
                v3 = Coordinates(7,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 7
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(5,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 8
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(4,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 9
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(2,:) - Centroid
                v3 = Coordinates(3,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 10
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(3,:) - Centroid
                v3 = Coordinates(4,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 11
                v1 = Coordinates(5,:) - Centroid
                v2 = Coordinates(6,:) - Centroid
                v3 = Coordinates(7,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 12
                v1 = Coordinates(5,:) - Centroid
                v2 = Coordinates(7,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
            end if
        end do
        !call FilePrinting(Self%VolumePerElement,'V','DataResults/.InternalData/VolumePerElement.txt')
    end subroutine VolumeAllElement
    ! ----------------------------------------------------------------- !
    !      subroutines required for the topology optimization proc.     !
    ! ----------------------------------------------------------------- !
    ! 1. Density filter
    subroutine DensityFilter(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        ! internal variables
        integer                                                      :: i,j,ne
        double precision                                             :: suma
        double precision                                             :: fac
        double precision, dimension(:), allocatable                  :: DiffCompilanceNew
        double precision, dimension(:), allocatable                  :: DiffVolumeNew
        double precision, dimension(:), allocatable                  :: Radius
        double precision, dimension(:,:), allocatable                :: RadPromE,PosPromE
        allocate(DiffCompilanceNew(self%Ne));      DiffCompilanceNew = 0.0d0
        allocate(DiffVolumeNew(Self%Ne));              DiffVolumeNew = 0.0d0
        allocate(RadPromE(self%Ne,Self%DimAnalysis));       RadPromE = 0.0d0
        allocate(PosPromE(self%Ne,Self%DimAnalysis));       PosPromE = 0.0d0
        ! Element position
        do i = 1, Self%Ne, 1
            PosPromE(i,:) = Sum(Self%Coordinates(Self%ConnectivityN(i,:),:),1)/Self%Npe
        end do
        
        !$omp parallel default(none) shared(Self, PosPromE, DiffCompilanceNew, DiffVolumeNew) private(i, j, suma, RadPromE, Radius, fac)
        !$omp for
        do i = 1, Self%Ne, 1
            suma = 0.0d0
            do j = 1, Self%DimAnalysis, 1
                RadPromE(:,j) = PosPromE(:,j) - PosPromE(i,j)
            end do
            Radius = sqrt(sum(RadPromE**2,2))
            do j = 1, Self%Ne, 1
                fac = Self%FilterRadius - Radius(j)
                suma = suma + max(0.0d0,fac)
                DiffCompilanceNew(i) = DiffCompilanceNew(i) + (max(0.0d0,fac))*Self%DensityVector(j)*Self%DiffCompilance(j)
                DiffVolumeNew(i) = DiffVolumeNew(i) + (max(0.0d0,fac))*Self%DensityVector(j)*Self%VolumePerElement(j)
            end do
            DiffCompilanceNew(i) = DiffCompilanceNew(i)/(Self%DensityVector(i)*suma)
            DiffVolumeNew(i) = DiffVolumeNew(i)/(Self%DensityVector(i)*suma)
        end do
        !$omp end for
        !$omp end parallel
        
        Self%DiffCompilance = DiffCompilanceNew
        Self%DiffVolume = DiffVolumeNew/(Self%DensityVector*Self%VolumePerElement)
        deallocate(DiffCompilanceNew,DiffVolumeNew,Radius,RadPromE)
    end subroutine DensityFilter
    ! 2. Optimization algorithm
    ! 2.1. Optimality criteria (OC)
    subroutine OptimalityCriteria(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        ! internal variables
        double precision                                             :: move
        double precision                                             :: xminval
        double precision                                             :: v1
        double precision                                             :: v2
        double precision                                             :: L1
        double precision                                             :: L2
        double precision                                             :: Lmid
        double precision, dimension(:), allocatable                  :: one
        double precision, dimension(:), allocatable                  :: DensityVectorNew
        double precision, dimension(:), allocatable                  :: DensityVectorOld
        write(unit=*, fmt=*) '-- Getting new solution (OC)'
        allocate(one(Self%ne));                               one = 1.0d0
        allocate(DensityVectorNew(Self%ne));     DensityVectorNew = 0.0d0
        DensityVectorOld = Self%DensityVector
        xminval = 0.01d0
        move = Self%MutationRate
        L1 = Self%L1
        L2 = Self%L2
        do while((L2-L1)/(L2+L1).gt.1e-3)
            Lmid = 0.5d0*(L2 + L1)
            DensityVectorNew = min(DensityVectorOld+move,DensityVectorOld*sqrt((-1.0d0*Self%DiffCompilance)/Lmid))
            DensityVectorNew = min(one,DensityVectorNew)
            DensityVectorNew = max(DensityVectorOld-move,DensityVectorNew)
            DensityVectorNew = max(xminval,DensityVectorNew)
            v1 = sum(DensityVectorNew*Self%VolumePerElement)
            v2 = sum(Self%VolumePerElement*Self%VolFraction)
            if ((v1-v2).gt.0.0d0) then
                L1 = Lmid
            else
                L2 = Lmid
            end if
        end do
        ! updating
        Self%DensityVector = DensityVectorNew
        Self%Change = maxval(abs(DensityVectorOld-DensityVectorNew))
        deallocate(DensityVectorNew,DensityVectorOld,one)
    end subroutine OptimalityCriteria
    ! 2.2. Method of Moving Asymptotes (MMA)
    subroutine MethodMovingAsympotes(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision                                             :: v1
        double precision                                             :: v2
        write(unit=*, fmt=*) '-- Getting new solution (MMA)'
        Self%xMMA(:,1) = Self%DensityVector
        v1 = sum(Self%DensityVector*Self%VolumePerElement)
        v2 = sum(Self%VolumePerElement*Self%VolFraction)
        Self%fval(1,1) = (v1 - v2)
        Self%dfdx(1,:) = Self%DiffVolume
        Self%f0val = sum(Self%Compilance)
        Self%df0dx(:,1) = Self%DiffCompilance
        ! applying the MMA solver interface of Svanberg
        call MMA_Solver_Interface(self%mm,Self%nn,Self%Iteration,Self%XMMA,Self%xmin,Self%xmax,Self%xold1,Self%xold2, &
                                Self%f0val,Self%df0dx,Self%fval,Self%dfdx,Self%low,Self%upp,Self%a0,Self%a,Self%c,Self%d)
        ! Replace solutions
        Self%xold2 = Self%xold1
        Self%xold1(:,1) = Self%DensityVector
        Self%DensityVector = Self%XMMA(:,1)
        Self%Change = maxval(abs(Self%DensityVector-Self%xold1(:,1)))
        !call FilePrinting(Self%DensityVector,'V','DataResults/DensityVector.txt')
    end subroutine MethodMovingAsympotes
    ! 3. Compilance and Diff Compilance
    subroutine GetCompilance(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, dimension(:), allocatable                  :: ObjFun
        ObjFun = Self%StrainEnergyE
        Self%Compilance = (Self%DensityVector**Self%PenalFactor)*ObjFun
        Self%DiffCompilance = -Self%PenalFactor*(Self%DensityVector**(Self%PenalFactor-1.0d0))*ObjFun
        deallocate(ObjFun)
    end subroutine GetCompilance
    ! 4. Final Topology
    subroutine FinalToplogy(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        integer                                                      :: i
        Self%FinalGeometry = Pack([(i,i=1,Self%Ne)],Self%DensityVector.ge.Self%PostProcesingFilter)
    end subroutine FinalToplogy
    ! 5. Topology Optimization process
    subroutine TopologyOptimizationProcess(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        call UploadOptimizationParameters(Self)
        call PreAssemblyRoutine(Self)
        Self%Iteration = 1
        write(unit=*, fmt=*) '---------- topology optimization process ----------'
        do
            write(unit=*, fmt=*) '-Updating mechanical parameters and solving the system'
            call UploadStructure(Self,Self%DensityVector,Self%PenalFactor)
            call SolveSystemMA87(Self)
            call ProcessingResults(Self)
            write(unit=*, fmt=*) '-Getting compilance and Diff'
            call GetCompilance(Self)
            write(unit=*, fmt=*) '-Applying Density Filtering'
            call DensityFilter(Self)
            write(unit=*, fmt=*) '-Getting new solution'
            if (Self%OptAlgorithm.eq.'OCM') then; call OptimalityCriteria(Self)    ; end if
            if (Self%OptAlgorithm.eq.'MMA') then; call MethodMovingAsympotes(Self) ; end if 
            ! Printing Convergence
            call PrintingConvergence(Self%Iteration,Self)
            if (Self%Change.lt.Self%OptimizationTolerance) then; write(unit=*, fmt=*) 'Precision achieved'        ; exit; end if
            if (Self%Iteration.ge.Self%MaxTOIteration)     then; write(unit=*, fmt=*) 'Max TO iterations reached' ; exit; end if 
            Self%Iteration = Self%Iteration + 1
        end do
        call FinalToplogy(Self)
    end subroutine TopologyOptimizationProcess
end module Optimization_module