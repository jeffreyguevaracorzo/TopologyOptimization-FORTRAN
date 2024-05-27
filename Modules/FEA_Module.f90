module FEA_Module
    use Base_Module
    use Solver_MA87Module
    implicit none
    ! Type variable for all the information of the structure
    type Structure
        character(len=20)                               :: AnalysisType     ! PlaneStress(2D), PlainStrain(2D), SolidIso(3D)
        character(len=20)                               :: ElementType      ! tri3 tri6, cuad4, cuad8, tetra4, tetra10, hexa8, hexa20
        integer                                         :: DimAnalysis      ! 2D or 3D
        integer                                         :: QuadGauss        ! Number of GaussPoints
        integer                                         :: N                ! Number of nodes
        integer                                         :: Ne               ! Number of elements
        integer                                         :: Npe              ! Number of nodes per element
        integer                                         :: NBc              ! Number of nodes with restrictions
        integer                                         :: Npl              ! Number of point loads
        integer                                         :: Ndl              ! Number of distributed loads
        integer, Dimension(:), allocatable              :: FreeD            ! Free degrees of freedom
        integer, Dimension(:), allocatable              :: FixedD           ! Fixed degrees of freedom
        integer, dimension(:,:), allocatable            :: ConnectivityN, ConnectivityD            ! Connectivity
        double precision                                :: YoungModulus, PoissonModulus, Thickness ! (in 3D thickness = 0.0)
        double precision, dimension(:), allocatable     :: GaussPoint       ! Gauss Approximation points
        double precision, dimension(:), allocatable     :: GaussWeights     ! Gauss Approximation weights
        double precision, dimension(:), allocatable     :: FGlobal_PL       ! Global vector of point loads
        double precision, dimension(:), allocatable     :: FGlobal_DL       ! Global vector of distributed loads
        double precision, dimension(:,:), allocatable   :: Coordinates      ! Coordinates of nodes
        double precision, dimension(:,:,:), allocatable :: BLocal           ! Matrix of derivatives of shape functions
        double precision, dimension(:,:,:), allocatable :: KLocal           ! Local stiffness matrix per element
        ! System of equations [K]{u}={f} in sparse format
        integer, dimension(:), allocatable              :: index_i_KGlobal  ! Indices for rows of the global stiffness matrix
        integer, dimension(:), allocatable              :: index_j_KGlobal  ! Indices for columns of the global stiffness matrix
        integer, dimension(:), allocatable              :: Rows_KGlobal     
        integer, dimension(:), allocatable              :: Cols_KGlobal     
        double precision, dimension(:), allocatable     :: value_KGlobal    ! Value of the global stiffness matrix at position i,j
        double precision, dimension(:), allocatable     :: value_FGlobal    ! Global load vector
        double precision, dimension(:), allocatable     :: value_UGlobal    ! Global displacement vector
        ! performing variables
        integer, dimension(:,:), allocatable            :: Node_Interaction
        integer, dimension(:,:), allocatable            :: Elem_Interaction 
        integer, dimension(:,:,:), allocatable          :: Location_KGlobal 
        double precision, dimension(:), allocatable     :: UGlobal          ! Global displacement vector
        ! Results
        double precision, dimension(:), allocatable     :: StrainEnergyE    ! Strain energy per element
        double precision, dimension(:,:), allocatable   :: Displacement     ! Global dispacement matrix
        double precision, dimension(:,:), allocatable   :: StrainE          ! Strain per element
        double precision, dimension(:,:), allocatable   :: StressE          ! Stress per element
    contains
        procedure                                       :: SetAnalysisType
        procedure                                       :: SetElementType
        procedure                                       :: SetYoungModulus
        procedure                                       :: SetPoissonModulus
        procedure                                       :: SetThickness
        procedure                                       :: SetGaussAprox
        procedure                                       :: SetCoordinates
        procedure                                       :: SetConnectivity
        procedure                                       :: SetBondaryConditions
        procedure                                       :: SetPointLoads
        procedure                                       :: SetDistributedLoads
        procedure                                       :: PreAssemblyRoutine
        procedure                                       :: GetKGlobalSparse
        procedure                                       :: GetFGlobalSparse
        procedure                                       :: UploadStructure
        procedure                                       :: SolveSystemMA87
        procedure                                       :: ProcessingResults
    end type Structure

    contains
    ! --------------- ADDITIONAL BASE FUNCTIONS AND SUBROUTINES ---------------
    ! 1. Reading files routines
    Subroutine ReadingfileInteger(Path,Nrow,Ncol,Matrix)
        implicit none
        character(len=*), intent(in)                         :: Path
        integer                                              :: ios,iounit,i,j
        integer, intent(inout)                               :: Nrow
        integer, intent(in)                                  :: Ncol
        integer, dimension(:,:), allocatable, intent(inout)  :: Matrix
        open(unit=iounit, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=iounit, fmt=*) !title
            read(unit=iounit, fmt=*) Nrow ; allocate(Matrix(Nrow,Ncol))
            read(unit=iounit, fmt=*) !references
            do i = 1, Nrow, 1
                read(unit=iounit, fmt=*) (Matrix(i,j), j = 1, Ncol, 1)
            end do
        close(iounit)
    end subroutine ReadingfileInteger
    Subroutine ReadingfileDP(Path,Nrow,Ncol,Matrix)
        implicit none
        character(len=*), intent(in)                                  :: Path
        integer                                                       :: ios,iounit,i,j
        integer, intent(inout)                                        :: Nrow
        integer, intent(in)                                           :: Ncol
        double precision, dimension(:,:), allocatable, intent(inout)  :: Matrix
        open(unit=iounit, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=iounit, fmt=*) !title
            read(unit=iounit, fmt=*) Nrow ; allocate(Matrix(Nrow,Ncol))
            read(unit=iounit, fmt=*) !references
            do i = 1, Nrow, 1
                read(unit=iounit, fmt=*) (Matrix(i,j), j = 1, Ncol, 1)
            end do
        close(iounit)
    end subroutine ReadingfileDP
    ! 2. Area
    function Area(Coordinates,Type) result(AreaAprox)
        character(len=*), intent(in)                                :: Type
        double precision                                            :: AreaAprox
        double precision, dimension(:), allocatable                 :: vector1, vector2, vector3, vector4
        double precision, dimension(:), allocatable                 :: Av1, Av2
        double precision, dimension(:,:), allocatable, intent(in)   :: Coordinates
        allocate(vector1(3),vector2(3),vector3(3),vector4(3))
        allocate(Av1(3), Av2(3))
        if ((Type.eq.'tetra4').or.(Type.eq.'tetra10')) then
            vector1 = Coordinates(2,:)-Coordinates(1,:)
            vector2 = Coordinates(3,:)-Coordinates(1,:)
            Av1 = CrossProduct(vector1,vector2)
            AreaAprox = abs(Norm(Av1))/2
        elseif ((Type.eq.'hexa8').or.(Type.eq.'hexa20')) then
            vector1 = Coordinates(1,:)-Coordinates(2,:)
            vector2 = Coordinates(3,:)-Coordinates(2,:)
            vector3 = Coordinates(1,:)-Coordinates(4,:)
            vector4 = Coordinates(3,:)-Coordinates(4,:)
            Av1 = CrossProduct(vector1,vector2)
            Av2 = CrossProduct(vector3,vector4)
            AreaAprox = abs(Norm(Av1))/2 + abs(Norm(Av2))/2
        end if
    end function Area
    ! 3. Sorting
    subroutine Sort(vector,n)
        implicit none
        integer, intent(in)                 :: n
        integer, allocatable, intent(inout) :: vector(:)
        ! internal variables
        integer                             :: i, j, temp
        integer, allocatable                :: dumb(:)
        ! Ordenado de valores
        do i = 1, n-1
            do j = 1, n-i
                if (vector(j) > vector(j+1)) then
                    temp = vector(j)
                    vector(j) = vector(j+1)
                    vector(j+1) = temp
                endif
            end do
        end do
        ! eliminar 0's y repetidos
        do i = 2, n, 1
            if (Vector(i-1).eq.vector(i)) then 
            vector(i-1) = 0
            end if
        end do
        vector = pack(vector,vector.gt.0)
    end subroutine Sort
    ! ----------------- SUBROUTINES FOR STRUCTURE INFORMATION ----------------- 
    ! 1. Input Analysis Type
    subroutine SetAnalysisType(Self,AnalysisType)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: AnalysisType
        Self%AnalysisType = AnalysisType
        if (Self%AnalysisType.eq.'PlaneStress'.or.Self%AnalysisType.eq.'PlaneStrain') then
            Self%DimAnalysis = 2
        elseif (Self%AnalysisType.eq.'SolidIso') then 
            Self%DimAnalysis = 3
        else
            stop "ERROR, Setting Analysis type"
        end if
    end subroutine SetAnalysisType
    ! 2. Input Element Type
    subroutine SetElementType(Self,ElementType)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: ElementType
        Self%ElementType = ElementType
        if ( Self%ElementType.eq.'tria3' ) then; Self%Npe = 3 ; end if
        if ( Self%ElementType.eq.'tria6' ) then; Self%Npe = 6 ; end if
        if ( Self%ElementType.eq.'quad4' ) then; Self%Npe = 4 ; end if
        if ( Self%ElementType.eq.'quad8' ) then; Self%Npe = 8 ; end if
        if ( Self%ElementType.eq.'tetra4' ) then; Self%Npe = 4 ; end if
        if ( Self%ElementType.eq.'tetra10' ) then; Self%Npe = 10 ; end if
        if ( Self%ElementType.eq.'hexa8' ) then; Self%Npe = 8 ; end if
        if ( Self%ElementType.eq.'hexa20' ) then; Self%Npe = 20 ; end if
    end subroutine SetElementType
    ! 3. Input Thickness (only for 3D cases)
    subroutine SetThickness(Self,Thickness)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(in)                                :: Thickness
        Self%Thickness = Thickness
    end subroutine SetThickness
    ! 4. Input Young Modulus
    subroutine SetYoungModulus(Self,YoungModulus)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(in)                                :: YoungModulus
        Self%YoungModulus = YoungModulus
    end subroutine SetYoungModulus
    ! 5. Input Poisson Modulus
    subroutine SetPoissonModulus(Self,PoissonModulus)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(in)                                :: PoissonModulus
        Self%PoissonModulus = PoissonModulus
    end subroutine SetPoissonModulus
    ! 6. Input Gauss Aproximation
    subroutine SetGaussAprox(Self,Gauss)
        implicit none
        class(Structure), intent(inout)                             :: Self 
        integer, intent(in)                                         :: Gauss
        double precision                                            :: LB, UB, Coef1, Coef2
        if (Gauss.le.0.or.Gauss.gt.5) stop "ERROR, GuassQuadrature (greater than 5 or <= to zero)" 
        Self%QuadGauss = Gauss
        select case (Gauss)
            case (1)
                allocate(Self%GaussPoint(1))
                allocate(Self%GaussWeights(1))
                Self%GaussPoint = [0.00000000d0]
                Self%GaussWeights = [2.00000000d0]
            case (2)
                allocate(Self%GaussPoint(2))
                allocate(Self%GaussWeights(2))
                Self%GaussPoint = [-0.57735026d0,0.57735026d0]
                Self%GaussWeights = [1.00000000d0,1.00000000d0]
            case (3)
                allocate(Self%GaussPoint(3))
                allocate(Self%GaussWeights(3))
                Self%GaussPoint = [-0.77459666d0,0.00000000d0,0.77459666d0]
                Self%GaussWeights = [0.55555555d0,0.88888888d0,0.55555555d0]
            case (4)
                allocate(Self%GaussPoint(4))
                allocate(Self%GaussWeights(4))
                Self%GaussPoint = [-0.86113631d0,-0.33998104d0,0.33998104d0,0.86113631d0]
                Self%GaussWeights = [0.34785484d0,0.65214515d0,0.65214515d0,0.34785484d0]
            case (5)
                allocate(Self%GaussPoint(5))
                allocate(Self%GaussWeights(5))
                Self%GaussPoint = [-0.90617984d0,-0.53846931d0,0.00000000d0,0.53846931d0,0.90617984d0]
                Self%GaussWeights = [0.23692688d0,0.47862867d0,0.56888888d0,0.47862867d0,0.23692688d0]
        end select
        LB = 0.0d0;  UB = 1.0d0
        Coef1 = (UB - LB)/2.0d0
        Coef2 = (UB + LB)/2.0d0
        ! changing coordinates and points due to boundary changes
        if (Self%ElementType.eq.'tria3'.or.Self%ElementType.eq.'tria6'.or. &
            Self%ElementType.eq.'tetra4'.or.Self%ElementType.eq.'tetra10') then
            Self%GaussPoint = Self%GaussPoint*Coef1 + Coef2
            Self%GaussWeights = Self%GaussWeights*Coef1
        end if
    end subroutine SetGaussAprox
    ! 7. Input Coordinates
    subroutine SetCoordinates(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        call ReadingfileDP(Path,Self%N,Self%DimAnalysis,Self%Coordinates)
        write(unit=*, fmt=*) '- Coordinates'
        !call FilePrinting(Self%Coordinates,'DataResults/.InternalData/CheckingCoordinatesLecture.txt')
    end subroutine SetCoordinates
    ! 8. Input Connectivity
    subroutine SetConnectivity(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,k
        call ReadingfileInteger(Path,Self%Ne,Self%Npe,Self%ConnectivityN)
        allocate(Self%ConnectivityD(Self%Ne,Self%Npe*Self%DimAnalysis))
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = Self%DimAnalysis-1, 0, -1
                    Self%ConnectivityD(i,j*Self%DimAnalysis - k) = Self%ConnectivityN(i,j)*Self%DimAnalysis - k
                end do
            end do
        end do
        write(unit=*, fmt=*) '- Connectivity'
        !call FilePrinting(Self%ConnectivityN,'DataResults/.InternalData/CheckingConnectivityNLecture.txt')
        !call FilePrinting(Self%ConnectivityD,'DataResults/.InternalData/CheckingConnectivityDLecture.txt')
    end subroutine SetConnectivity
    ! 9. Input Boundary Conditions
    subroutine SetBondaryConditions(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer, dimension(:,:), allocatable                        :: BoundaryC
        integer                                                     :: i,j,k
        call ReadingfileInteger(Path,Self%NBc,(Self%DimAnalysis+1),BoundaryC)
        i = (Self%N)*(Self%DimAnalysis) - count(BoundaryC(:,2:).eq.1)
        j = count(BoundaryC(:,2:).eq.1)
        k = 1
        allocate(Self%FreeD(i))
        allocate(Self%FixedD(j))
        Self%FreeD = 0
        Self%FixedD = 0
        ! constrained?
        do i = 1, Self%NBc, 1
            do j = 1, Self%DimAnalysis, 1
                if (BoundaryC(i,j+1).eq.1) then
                    Self%FixedD(k) = BoundaryC(i,1)*Self%DimAnalysis - (Self%DimAnalysis-j)
                    k = k + 1
                else
                    cycle
                end if
            end do
        end do
        ! free!
        k = 1
        do i = 1, Self%N*Self%DimAnalysis, 1
            if (any(Self%FixedD.eq.i)) then
                cycle
            else
                Self%FreeD(k) = i
                k = k + 1
            end if
        end do
        write(unit=*, fmt=*) '- Boundary conditions/Constrains'
        !call FilePrinting(Self%FreeD,'V','DataResults/.InternalData/CheckingFreeDofLecture.txt')
        !call FilePrinting(Self%FixedD,'V','DataResults/.InternalData/CheckingFixedDofLecture.txt')
    end subroutine SetBondaryConditions
    ! 10. Input Point Loads
    subroutine SetPointLoads(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,node,dim
        double precision, dimension(:,:), allocatable               :: PointLoads
        dim = Self%DimAnalysis
        call ReadingfileDP(Path,Self%Npl,(dim+1),PointLoads)
        allocate(Self%FGlobal_PL(Self%N*dim)) ; Self%FGlobal_PL = 0.0d0
        ! assembly load vector
        do i = 1, Self%Npl, 1
            node = int(PointLoads(i,1))
            do j = 1, dim, 1
                Self%FGlobal_PL(node*dim-(dim-j)) = Self%FGlobal_PL(node*dim-(j-dim)) + PointLoads(i,j+1)
            end do
        end do
        write(unit=*, fmt=*) '- Point Loads'
        !call FilePrinting(Self%FGlobal_PL,'V','DataResults/.InternalData/CheckingPointLoadsLecture.txt')
    end subroutine SetPointLoads
    ! 11. Input Distributed Loads
    subroutine SetDistributedLoads(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: n,i,j,k,Npf,dim,node
        double precision                                            :: FaceArea, length
        double precision, dimension(:), allocatable                 :: Vector, ForcePerNode
        double precision, dimension(:,:), allocatable               :: DistributedLoads
        double precision, dimension(:,:), allocatable               :: Coordinates
        ! calculation and assembly of load vectors
        dim = Self%DimAnalysis
        allocate(Self%FGlobal_DL(Self%N*dim)); Self%FGlobal_DL = 0.0d0;
        allocate(ForcePerNode(dim));              ForcePerNode = 0.0d0;
        allocate(Vector(dim));                          vector = 0.0d0;
        ! nodes per face
        if (Self%ElementType.eq.'tria3'.or.Self%ElementType.eq.'quad4') then; Npf = 2; end if
        if (Self%ElementType.eq.'tria6'.or.Self%ElementType.eq.'quad8') then; Npf = 3; end if
        if (Self%ElementType.eq.'tetra4')                               then; Npf = 3; n=3; allocate(Coordinates(3,3)); end if
        if (Self%ElementType.eq.'tetra10')                              then; Npf = 6; n=3; allocate(Coordinates(3,3)); end if
        if (Self%ElementType.eq.'hexa8')                                then; Npf = 4; n=4; allocate(Coordinates(4,3)); end if
        if (Self%ElementType.eq.'hexa20')                               then; Npf = 8; n=4; allocate(Coordinates(4,3)); end if
        ! reading file
        call ReadingfileDP(Path,Self%Ndl,(Npf+dim),DistributedLoads)
        select case (dim)
            case (2)
                do i = 1, Self%Ndl, 1
                    Vector = Self%Coordinates(DistributedLoads(i,2),1:dim) - Self%Coordinates(DistributedLoads(i,1),1:dim)
                    length = Norm(Vector)
                    FaceArea = length*Self%Thickness
                    ForcePerNode = (FaceArea/Npf)*(DistributedLoads(i,(Npf+1):))
                    do j = 1, Npf, 1
                        node = int(DistributedLoads(i,j))
                        do k = 1, dim, 1
                            Self%FGlobal_DL(node*dim-(dim-k)) = ForcePerNode(k)
                        end do
                    end do
                end do
            case (3)
                do i = 1, Self%Ndl, 1
                    Coordinates = Self%Coordinates(int(DistributedLoads(i,1:n)),:)
                    FaceArea = Area(Coordinates,Self%ElementType)
                    ForcePerNode = (FaceArea/Npf)*DistributedLoads(i,(Npf+1):)
                    do j = 1, Npf, 1
                        node = int(DistributedLoads(i,j))
                        do k = 1, dim, 1
                            Self%FGlobal_DL(node*dim-(dim-k)) = ForcePerNode(k)
                        end do
                    end do
                end do
        end select
        write(unit=*, fmt=*) '- Distributed Loads'
        !call FilePrinting(Self%FGlobal_PL,'V','DataResults/.InternalData/CheckingDistributedLoadsLecture.txt')
    end subroutine SetDistributedLoads
    ! 12. Read Files (Connectivity,Coordinates,BoundaryConditions,etc...)
    subroutine ReadFiles(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        call SetCoordinates(Self,'DataStructure/Coordinates.txt')
        call SetConnectivity(Self,'DataStructure/Connectivity.txt')
        call SetBondaryConditions(Self,'DataStructure/Constrains.txt')
        call SetPointLoads(Self,'DataStructure/PointLoads.txt')
        call SetDistributedLoads(Self,'DataStructure/DistributedLoads.txt')
    end subroutine ReadFiles
    ! ------------ FINITE ELEMENT ANALYSISS FUNCTIONS AND SUBROUTINES ------------
    ! 1. Derivative of the shape functions
    subroutine DiffFormFunction(Self,DiffFunction,e,n,z)
        ! dN1/de     dN2/de     dN3/de      ....     dNn/de
        ! dN1/dn     dN2/dn     dN3/dn      ....     dNn/dn
        ! dN1/dz     dN2/dz     dN3/dz      ....     dNn/dz (3D)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                             :: e, n
        double precision, intent(inout), optional                   :: z
        double precision, dimension(:,:), allocatable, intent(out)  :: DiffFunction
        if (Self%ElementType.eq.'tria3') then
            allocate(DiffFunction(2,3))
            !  line 1
            DiffFunction(1,1) = 1.0d0
            DiffFunction(1,2) = 0.0d0
            DiffFunction(1,3) = - 1.0d0
            !  line 2
            DiffFunction(2,1) = 0.0d0
            DiffFunction(2,2) = 1.0d0
            DiffFunction(2,3) = - 1.0d0
        elseif (Self%ElementType.eq.'tria6') then
            allocate(DiffFunction(2,6))
            !  line 1
            DiffFunction(1,1) = 4.0d0*e + 4.0d0*n - 3.0d0 
            DiffFunction(1,2) = 4.0d0*e - 1.0d0 
            DiffFunction(1,3) = 0.0d0 
            DiffFunction(1,4) = 4.0d0 - 4.0d0*n - 8.0d0*e 
            DiffFunction(1,5) = 4.0d0*n 
            DiffFunction(1,6) = - 4.0d0*n 
            !  line 2
            DiffFunction(2,1) = 4.0d0*e + 4.0d0*n - 3.0d0 
            DiffFunction(2,2) = 0.0d0 
            DiffFunction(2,3) = 4.0d0*n - 1.0d0 
            DiffFunction(2,4) = - 4.0d0*e 
            DiffFunction(2,5) = 4.0d0*e 
            DiffFunction(2,6) = 4.0d0 - 8.0d0*n - 4.0d0*e 
        elseif (Self%ElementType.eq.'quad4') then
            allocate(DiffFunction(2,4))
            !  line 1
            DiffFunction(1,1) = n/4.0d0 - 1.0d0/4.0d0 
            DiffFunction(1,2) = 1.0d0/4.0d0 - n/4.0d0 
            DiffFunction(1,3) = n/4.0d0 + 1.0d0/4.0d0 
            DiffFunction(1,4) = - n/4.0d0 - 1.0d0/4.0d0 
            !  line 2
            DiffFunction(2,1) = e/4.0d0 - 1.0d0/4.0d0 
            DiffFunction(2,2) = - e/4.0d0 - 1.0d0/4.0d0 
            DiffFunction(2,3) = e/4.0d0 + 1.0d0/4.0d0 
            DiffFunction(2,4) = 1.0d0/4.0d0 - e/4.0d0 
        elseif (Self%ElementType.eq.'quad8') then
            allocate(DiffFunction(2,8))
            !  line 1
            DiffFunction(1,1) = - (e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) - ((n - 1.0d0)*(e + n + 1.0d0))/4.0d0 
            DiffFunction(1,2) = ((n - 1.0d0)*(n - e + 1.0d0))/4.0d0 - (e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(1,3) = (e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) + ((n + 1.0d0)*(e + n - 1.0d0))/4.0d0 
            DiffFunction(1,4) = (e/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) + ((n + 1.0d0)*(e - n + 1.0d0))/4.0d0 
            DiffFunction(1,5) = e*(n - 1.0d0) 
            DiffFunction(1,6) = 1.0d0/2.0d0 - n**2.0d0/2.0d0 
            DiffFunction(1,7) = - e*(n + 1.0d0) 
            DiffFunction(1,8) = n**2.0d0/2.0d0 - 1.0d0/2.0d0 
            !  line 2
            DiffFunction(2,1) = - (e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) - (e/4.0d0 - 1.0d0/4.0d0)*(e + n + 1.0d0) 
            DiffFunction(2,2) = (e/4.0d0 + 1.0d0/4.0d0)*(n - e + 1.0d0) + (e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(2,3) = (e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) + (e/4.0d0 + 1.0d0/4.0d0)*(e + n - 1.0d0) 
            DiffFunction(2,4) = (e/4.0d0 - 1.0d0/4.0d0)*(e - n + 1.0d0) - (e/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(2,5) = e**2.0d0/2.0d0 - 1.0d0/2.0d0 
            DiffFunction(2,6) = - 2.0d0*n*(e/2.0d0 + 1.0d0/2.0d0) 
            DiffFunction(2,7) = 1.0d0/2.0d0 - e**2.0d0/2.0d0 
            DiffFunction(2,8) = 2.0d0*n*(e/2.0d0 - 1.0d0/2.0d0) 
        elseif (Self%ElementType.eq.'tetra4') then
            allocate(DiffFunction(3,4))
            !  line 1
            DiffFunction(1,1) = - 1.0d0 
            DiffFunction(1,2) = 1.0d0
            DiffFunction(1,3) = 0.0d0
            DiffFunction(1,4) = 0.0d0
            !  line 2
            DiffFunction(2,1) = - 1.0d0
            DiffFunction(2,2) = 0.0d0
            DiffFunction(2,3) = 1.0d0
            DiffFunction(2,4) = 0.0d0
            !  line 3
            DiffFunction(3,1) = - 1.0d0
            DiffFunction(3,2) = 0.0d0
            DiffFunction(3,3) = 0.0d0
            DiffFunction(3,4) = 1.0d0
        elseif (Self%ElementType.eq.'tetra10') then
            allocate(DiffFunction(3,10))
            !  line 1
            DiffFunction(1,1) = 4.0d0*e + 4.0d0*n + 4.0d0*z - 3.0d0 
            DiffFunction(1,2) = 4.0d0*e - 1.0d0 
            DiffFunction(1,3) = 0.0d0 
            DiffFunction(1,4) = 0.0d0 
            DiffFunction(1,5) = 4.0d0 - 4.0d0*n - 4.0d0*z - 8.0d0*e 
            DiffFunction(1,6) = 4.0d0*n 
            DiffFunction(1,7) = - 4.0d0*n 
            DiffFunction(1,8) = 4.0d0*z 
            DiffFunction(1,9) = 0 
            DiffFunction(1,10) = - 4.0d0*z 
            !  line 2
            DiffFunction(2,1) = 4.0d0*e + 4.0d0*n + 4.0d0*z - 3.0d0 
            DiffFunction(2,2) = 0.0d0 
            DiffFunction(2,3) = 4.0d0*n - 1.0d0 
            DiffFunction(2,4) = 0.0d0 
            DiffFunction(2,5) = - 4.0d0*e 
            DiffFunction(2,6) = 4.0d0*e 
            DiffFunction(2,7) = 4.0d0 - 8.0d0*n - 4.0d0*z - 4.0d0*e 
            DiffFunction(2,8) = 0.0d0 
            DiffFunction(2,9) = 4.0d0*z 
            DiffFunction(2,10) = - 4.0d0*z 
            !  line 3
            DiffFunction(3,1) = 4.0d0*e + 4.0d0*n + 4.0d0*z - 3.0d0 
            DiffFunction(3,2) = 0.0d0 
            DiffFunction(3,3) = 0.0d0 
            DiffFunction(3,4) = 4.0d0*z - 1.0d0 
            DiffFunction(3,5) = - 4.0d0*e 
            DiffFunction(3,6) = 0.0d0 
            DiffFunction(3,7) = - 4.0d0*n 
            DiffFunction(3,8) = 4.0d0*e 
            DiffFunction(3,9) = 4.0d0*n 
            DiffFunction(3,10) = 4.0d0 - 4.0d0*n - 8.0d0*z - 4.0d0*e 
        elseif (Self%ElementType.eq.'hexa8') then
            allocate(DiffFunction(3,8))
            !  line 1
            DiffFunction(1,1) = - ((n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,2) = ((n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,3) = - ((n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,4) = ((n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,5) = ((n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,6) = - ((n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,7) = ((n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,8) = - ((n + 1.0d0)*(z + 1.0d0))/8.0d0 
            !  line 2
            DiffFunction(2,1) = - (e/8.0d0 - 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,2) = (e/8.0d0 + 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,3) = - (e/8.0d0 + 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,4) = (e/8.0d0 - 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,5) = (e/8.0d0 - 1.0d0/8.0d0)*(z + 1.0d0) 
            DiffFunction(2,6) = - (e/8.0d0 + 1.0d0/8.0d0)*(z + 1.0d0) 
            DiffFunction(2,7) = (e/8.0d0 + 1.0d0/8.0d0)*(z + 1.0d0) 
            DiffFunction(2,8) = - (e/8.0d0 - 1.0d0/8.0d0)*(z + 1.0d0) 
            !  line 3
            DiffFunction(3,1) = - (e/8.0d0 - 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,2) = (e/8.0d0 + 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,3) = - (e/8.0d0 + 1.0d0/8.0d0)*(n + 1.0d0) 
            DiffFunction(3,4) = (e/8.0d0 - 1.0d0/8.0d0)*(n + 1.0d0) 
            DiffFunction(3,5) = (e/8.0d0 - 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,6) = - (e/8.0d0 + 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,7) = (e/8.0d0 + 1.0d0/8.0d0)*(n + 1.0d0) 
            DiffFunction(3,8) = - (e/8.0d0 - 1.0d0/8.0d0)*(n + 1.0d0) 
        elseif (Self%ElementType.eq.'hexa20') then
            allocate(DiffFunction(3,20))
            !  line 1
            DiffFunction(1,1) = - (e*n*z*(n - 1.0d0)*(z - 1.0d0))/8.0d0 - (n*z*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,2) = (e*n*z*(n - 1.0d0)*(z - 1.0d0))/8.0d0 + (n*z*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,3) = - (e*n*z*(n + 1.0d0)*(z - 1.0d0))/8.0d0 - (n*z*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,4) = (e*n*z*(n + 1.0d0)*(z - 1.0d0))/8.0d0 + (n*z*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,5) = (e*n*z*(n - 1.0d0)*(z + 1.0d0))/8.0d0 + (n*z*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,6) = - (e*n*z*(n - 1.0d0)*(z + 1.0d0))/8.0d0 - (n*z*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,7) = (e*n*z*(n + 1.0d0)*(z + 1.0d0))/8.0d0 + (n*z*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,8) = - (e*n*z*(n + 1.0d0)*(z + 1.0d0))/8.0d0 - (n*z*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,9) = - (e*(n - 1.0d0)*(z - 1.0d0))/2.0d0 
            DiffFunction(1,10) = (e*(n + 1.0d0)*(z - 1.0d0))/2.0d0 
            DiffFunction(1,11) = - (e*(n + 1.0d0)*(z + 1.0d0))/2.0d0 
            DiffFunction(1,12) = (e*(n - 1.0d0)*(z + 1.0d0))/2.0d0 
            DiffFunction(1,13) = - ((n**2.0d0 - 1.0d0)*(z - 1.0d0))/4.0d0 
            DiffFunction(1,14) = ((n**2.0d0 - 1.0d0)*(z - 1.0d0))/4.0d0 
            DiffFunction(1,15) = - ((n**2.0d0 - 1.0d0)*(z + 1.0d0))/4.0d0 
            DiffFunction(1,16) = ((n**2.0d0 - 1.0d0)*(z + 1.0d0))/4.0d0 
            DiffFunction(1,17) = - ((z**2.0d0 - 1.0d0)*(n - 1.0d0))/4.0d0 
            DiffFunction(1,18) = ((z**2.0d0 - 1.0d0)*(n - 1.0d0))/4.0d0 
            DiffFunction(1,19) = - ((z**2.0d0 - 1.0d0)*(n + 1.0d0))/4.0d0 
            DiffFunction(1,20) = ((z**2.0d0 - 1.0d0)*(n + 1.0d0))/4.0d0 
            !  line 2
            DiffFunction(2,1) = - (e*n*z*(e - 1.0d0)*(z - 1.0d0))/8.0d0 - (e*z*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,2) = (e*n*z*(e + 1.0d0)*(z - 1.0d0))/8.0d0 + (e*z*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,3) = - (e*n*z*(e + 1.0d0)*(z - 1.0d0))/8.0d0 - (e*z*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,4) = (e*n*z*(e - 1.0d0)*(z - 1.0d0))/8.0d0 + (e*z*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,5) = (e*n*z*(e - 1.0d0)*(z + 1.0d0))/8.0d0 + (e*z*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,6) = - (e*n*z*(e + 1.0d0)*(z + 1.0d0))/8.0d0 - (e*z*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,7) = (e*n*z*(e + 1.0d0)*(z + 1.0d0))/8.0d0 + (e*z*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,8) = - (e*n*z*(e - 1.0d0)*(z + 1.0d0))/8.0d0 - (e*z*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,9) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,10) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,11) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,12) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,13) = - 2.0d0*n*(e/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,14) = 2.0d0*n*(e/4.0d0 + 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,15) = - 2.0d0*n*(e/4.0d0 + 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,16) = 2.0d0*n*(e/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,17) = - (e/4.0d0 - 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            DiffFunction(2,18) = (e/4.0d0 + 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            DiffFunction(2,19) = - (e/4.0d0 + 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            DiffFunction(2,20) = (e/4.0d0 - 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            !  line 3
            DiffFunction(3,1) = - (e*n*z*(e - 1.0d0)*(n - 1.0d0))/8.0d0 - (e*n*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,2) = (e*n*z*(e + 1.0d0)*(n - 1.0d0))/8.0d0 + (e*n*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,3) = - (e*n*z*(e + 1.0d0)*(n + 1.0d0))/8.0d0 - (e*n*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,4) = (e*n*z*(e - 1.0d0)*(n + 1.0d0))/8.0d0 + (e*n*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,5) = (e*n*z*(e - 1.0d0)*(n - 1.0d0))/8.0d0 + (e*n*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,6) = - (e*n*z*(e + 1.0d0)*(n - 1.0d0))/8.0d0 - (e*n*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,7) = (e*n*z*(e + 1.0d0)*(n + 1.0d0))/8.0d0 + (e*n*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,8) = - (e*n*z*(e - 1.0d0)*(n + 1.0d0))/8.0d0 - (e*n*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,9) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,10) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,11) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,12) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,13) = - (e/4.0d0 - 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,14) = (e/4.0d0 + 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,15) = - (e/4.0d0 + 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,16) = (e/4.0d0 - 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,17) = - 2.0d0*z*(e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,18) = 2.0d0*z*(e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,19) = - 2.0d0*z*(e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,20) = 2.0d0*z*(e/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
        else
            stop "ERROR, in DiffFormFunction, problem with ElementType"
        end if
    end subroutine DiffFormFunction
    ! 2. Elascitity Tensor
    subroutine ElasticityTensor(Self,ETensor)
        implicit none
        class(Structure), intent(inout)                                   :: Self    
        double precision, dimension(:,:), allocatable, intent(inout)      :: ETensor
        double precision                                                  :: E,V,Constant
        E = Self%YoungModulus
        V = Self%PoissonModulus    
        if (Self%AnalysisType.eq.'PlaneStress') then
            allocate(ETensor(3,3))
            Constant = E/(1.0d0 - V**2.0d0)
            ETensor(1,:) = [1.0d0,V,0.0d0]
            ETensor(2,:) = [V,1.0d0,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0-V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'PlaneStrain') then
            allocate(ETensor(3,3))
            Constant = E/((1.0d0 + V)*(1.0d0 - 2.0d0*V))
            ETensor(1,:) = [1.0d0-V,V,0.0d0]
            ETensor(2,:) = [V,1.0d0-V,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'SolidIso') then
            allocate(ETensor(6,6))
            Constant = E*(1.0d0-v)/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0,V/(1.0d0-v),V/(1.0d0-v),0.0d0,0.0d0,0.0d0]
            ETensor(2,:) = [V/(1.0d0-v),1.0d0,V/(1.0d0-v),0.0d0,0.0d0,0.0d0]
            ETensor(3,:) = [V/(1.0d0-v),V/(1.0d0-v),1.0d0,0.0d0,0.0d0,0.0d0]
            ETensor(4,:) = [0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v)),0.0d0,0.0d0]
            ETensor(5,:) = [0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v)),0.0d0]
            ETensor(6,:) = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v))]
            ETensor = Constant*ETensor
        end if
    end subroutine ElasticityTensor
    ! 3. Stiffness matrices (Local)
    subroutine GetKlocal(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                             :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)  :: DensityVector 
        ! internal variables
        integer                                                     :: el,i,j,k,l,m
        double precision                                            :: e,n,z,w1,w2,w3,DetJacobian,Fac
        double precision, dimension(:,:), allocatable               :: Be,Jacobian,InvJacobian,D
        double precision, dimension(:,:), allocatable               :: DiffN,DiffNXY,ElementCoordinates
        ! ---- 1st part ----
        ! in this section the local stiffness matrix (dense) of each of the elements is calculated using gauss 
        ! quadrature and each of them is stored in the array KLocal[element,[K_x,K_y]].
        if (Self%DimAnalysis.eq.2) then
            if (allocated(Self%KLocal).and.allocated(Self%BLocal)) then
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            else
                allocate(Self%KLocal(Self%Ne,Self%Npe*2,Self%Npe*2))
                allocate(Self%BLocal(Self%Ne,3,2*Self%Npe))
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            end if
            allocate(Be(3,2*Self%Npe));                                    Be = 0.0d0;
            allocate(Jacobian(2,2));                                 Jacobian = 0.0d0;
            allocate(InvJacobian(2,2));                           InvJacobian = 0.0d0;
            allocate(ElementCoordinates(Self%Npe,2));      ElementCoordinates = 0.0d0;
            do el = 1, Self%Ne, 1
                ! SIMP PenalFactorization
                call ElasticityTensor(Self,D)
                D = (DensityVector(el)**PenalFactor)*D
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        Fac = DetJacobian*w1*w2*(Self%Thickness)
                        ! Be
                        Be = 0.0d0
                        do k = 1, size(DiffN,2), 1
                            Be(1,2*k-1) = DiffNxy(1,k)
                            Be(2,2*k) = DiffNxy(2,k)
                            Be(3,2*k-1) = DiffNxy(2,k)
                            Be(3,2*k) = DiffNxy(1,k)
                        end do
                        ! Be Local
                        Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + Fac*Be
                        ! K Local
                        Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + Fac*(matmul(transpose(Be),matmul(D,Be)))*&
                                                                    w1*w2*(Self%Thickness)
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.3) then
            if (allocated(Self%KLocal).and.allocated(Self%BLocal)) then
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            else
                allocate(Self%KLocal(Self%Ne,Self%Npe*3,Self%Npe*3))
                allocate(Self%BLocal(Self%Ne,6,3*Self%Npe))
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            end if
            allocate(Be(6,3*Self%Npe));                                    Be = 0.0d0;
            allocate(Jacobian(3,3));                                 Jacobian = 0.0d0;
            allocate(InvJacobian(3,3));                           InvJacobian = 0.0d0;
            allocate(ElementCoordinates(Self%Npe,3));      ElementCoordinates = 0.0d0;
            do el = 1, Self%Ne, 1
                ! SIMP PenalFactorization
                call ElasticityTensor(Self,D)
                D = (DensityVector(el)**PenalFactor)*D
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k)
                            w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            Fac = DetJacobian*w1*w2*w3
                            ! Be
                            do l = 1, size(DiffN,2), 1
                                Be(1,3*l-2) = DiffNxy(1,l)
                                Be(2,3*l-1) = DiffNxy(2,l)
                                Be(3,3*l) = DiffNxy(3,l)
                                Be(4,3*l-2) = DiffNxy(2,l)
                                Be(4,3*l-1) = DiffNxy(1,l)
                                Be(5,3*l-1) = DiffNxy(3,l)
                                Be(5,3*l) = DiffNxy(2,l)
                                Be(6,3*l-2) = DiffNxy(3,l)
                                Be(6,3*l) = DiffNxy(1,l)
                            end do
                            ! Be Local
                            Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + Fac*Be
                            ! K Local
                            Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + Fac*(matmul(transpose(Be),matmul(D,Be)))
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if
        deallocate(Be,Jacobian,InvJacobian,DiffNXY,ElementCoordinates)
    end subroutine GetKlocal
    ! 4. Node-Element interaction
    subroutine PreAssemblyRoutine(self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        integer                                                     :: i,j,k,i1,IndexRow,IndexCol
        integer, dimension(:), allocatable                          :: InPosition
        logical, dimension(:), allocatable                          :: InLogical
        integer, dimension(:,:), allocatable                        :: Node_Interaction,Elem_Interaction
        ! -------------------------------------------------------------------------------------- !
        ! note: In this part, a preliminary assembly of the stiffness matrix (in sparse form)    !
        !       is made. The idea is to locate the position of each element of the matrix in     !
        !       the sparse vector so that the assembly of the global matrix in each iteration    !
        !       can be faster.                                                                   !
        ! -------------------------------------------------------------------------------------- !
        ! get de element and DoF incidence for each DoF (only the free degres)
        j = size(Self%FreeD)
        ! maximum of 50 elements per node
        allocate(Elem_Interaction(j,50))
        Elem_Interaction = 0
        ! considering the quad20 20*20 
        allocate(Node_Interaction(j,600))
        Node_Interaction = 0

        !$omp parallel do private(InLogical,InPosition)
        do i = 1, size(Self%FreeD), 1
            ! element interaction
            InLogical = any(Self%ConnectivityD.eq.Self%FreeD(i),2)
            InPosition = pack([(k,k=1,Self%Ne)],InLogical)
            j = count(InLogical)
            Elem_Interaction(i,1:j) = InPosition 
            ! DoF interaction
            InPosition = reshape(Self%ConnectivityD(Elem_Interaction(i,1:j),:),[j*(Self%Npe)*Self%DimAnalysis])
            ! Only free degres
            do k = 1, size(InPosition), 1
                if(any(InPosition(k).eq.Self%FixedD)) InPosition(k) = 0
            end do
            ! Sorting ang eliminating 0s
            call sort(InPosition,j*Self%Npe*Self%DimAnalysis)
            InLogical = InPosition.ge.Self%FreeD(i)
            InPosition = pack(InPosition,InLogical)
            j = size(InPosition)
            Node_Interaction(i,1:j) = InPosition
        end do
        !$omp end parallel do 

        ! this print the interaction of the Free DoF with the others
        !call FilePrinting(Elem_Interaction,'DataResults/.InternalData/Elem_Interaction.txt')
        !call FilePrinting(Node_Interaction,'DataResults/.InternalData/Node_Interaction.txt')

        ! pre-assembly
        i = count(Node_Interaction.ne.0)
        allocate(Self%index_i_KGlobal(i)) ; Self%index_i_KGlobal = 0
        allocate(Self%index_j_KGlobal(i)) ; Self%index_j_KGlobal = 0
        allocate(Self%Location_KGlobal(Self%Ne,Self%Npe*Self%DimAnalysis,Self%Npe*Self%DimAnalysis))
        Self%Location_KGlobal = 0
        k = 1
        do i = 1, size(Self%FreeD), 1                               ! Col
            do j = 1, count(Node_Interaction(i,:).ne.0), 1          ! Row
                Self%index_i_KGlobal(k) = findloc(Self%FreeD,Node_Interaction(i,j),1) ! Row-indx
                Self%index_j_KGlobal(k) = i                                           ! Col-indx
                do i1 = 1, count(Elem_Interaction(i,:).ne.0), 1
                    IndexRow = findloc(Self%ConnectivityD(Elem_Interaction(i,i1),:),Node_Interaction(i,j),1)
                    IndexCol = findloc(Self%ConnectivityD(Elem_Interaction(i,i1),:),Node_Interaction(i,1),1)
                    if ((IndexCol.eq.0).or.(IndexRow.eq.0)) cycle
                    Self%Location_KGlobal(Elem_Interaction(i,i1),IndexRow,IndexCol) = k
                end do
                k = k+1
            end do
        end do
    end subroutine PreAssemblyRoutine
    ! 4. Global Stifness Matrix (global sparse form)
    subroutine GetKGlobalSparse(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        integer                                                     :: n,i,j,k
        logical, dimension(:), allocatable                          :: InLogical
        n = size(Self%index_i_KGlobal)
        if(allocated(Self%value_KGlobal)) then
            Self%value_KGlobal = 0.0d0
        else
            allocate(Self%value_KGlobal(n)); Self%value_KGlobal = 0.0d0
        end if
        do i = 1, Self%Ne, 1
            do j = 1, Self%DimAnalysis*Self%Npe, 1
                do k = 1, Self%DimAnalysis*Self%Npe, 1
                    n = Self%Location_KGlobal(i,j,k)    
                    if (n.eq.0) cycle   ! doesnt exist 
                    Self%value_KGlobal(n) = Self%value_KGlobal(n) + Self%KLocal(i,j,k)
                end do
            end do
        end do
        ! eliminating null or zero elements (StiffnessMatrix=0.0)
        InLogical = Self%value_KGlobal.ne.0.0d0
        Self%Rows_KGlobal = pack(Self%index_i_KGlobal,InLogical)
        Self%Cols_KGlobal = pack(Self%index_j_KGlobal,InLogical)
        Self%value_KGlobal = pack(Self%value_KGlobal,InLogical)
        !call FilePrinting(Self%Rows_KGlobal,'V','DataResults/.InternalData/Rows_KGlobal.txt')
        !call FilePrinting(Self%Cols_KGlobal,'V','DataResults/.InternalData/Cols_KGlobal.txt') 
        !call FilePrinting(Self%value_KGlobal,'V','DataResults/.InternalData/Vals_KGlobal.txt') 
    end subroutine GetKGlobalSparse
    ! 5. Global Load Vector (global sparse form)
    subroutine GetFGlobalSparse(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        ! Assembly
        allocate(Self%value_FGlobal(size(Self%FreeD)))
        Self%value_FGlobal = 0.0d0
        Self%value_FGlobal = Self%FGlobal_PL(Self%FreeD) + Self%FGlobal_DL(Self%FreeD) 
    end subroutine GetFGlobalSparse
    ! 5. Updating Data Structure
    subroutine UploadStructure(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                             :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)  :: DensityVector 
        ! Getting Stiffness Matrix
        call GetKlocal(Self,DensityVector,PenalFactor)
        call GetKGlobalSparse(Self)
        ! Getting Load Vector
        call GetFGlobalSparse(Self)
    end subroutine UploadStructure
    ! 6. Solve Structural System
    subroutine SolveSystemMA87(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        if (allocated(Self%UGlobal)) then
            Self%UGlobal = 0.0d0;
        else
            allocate(Self%UGlobal(Self%n*Self%DimAnalysis));          Self%UGlobal = 0.0d0;
        end if
        ! Apply HSL-MA87 Solver
        Self%value_UGlobal = SparseSystemSolver(Self%Rows_KGlobal,Self%Cols_KGlobal,Self%value_KGlobal,Self%value_FGlobal)
        ! Linking Solution
        Self%UGlobal(Self%FreeD) = Self%value_UGlobal
        ! mechanical variables
        deallocate(Self%value_FGlobal,Self%value_UGlobal,Self%value_KGlobal)
    end subroutine SolveSystemMA87
    ! 7. Processing Results
    subroutine ProcessingResults(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        ! internal variables
        integer                                                     :: i,j
        integer, dimension(:), allocatable                          :: index
        double precision                                            :: energy
        double precision, dimension(:), allocatable                 :: Strain, Stress
        double precision, dimension(:,:), allocatable               :: D
        logical                                                     :: T1,T2,T3,T4
        call ElasticityTensor(Self,D)
        ! allocating
        T1 = allocated(Self%Displacement)
        T2 = allocated(Self%StrainE)
        T3 = allocated(Self%StressE)
        T4 = allocated(Self%StrainEnergyE)
        if(T1.and.T2.and.T3.and.T4) then
            Self%StrainE = 0.0d0
            Self%StressE = 0.0d0
            Self%StrainEnergyE = 0.0d0
        else
            allocate(Self%Displacement(Self%N,Self%DimAnalysis))
            if (Self%DimAnalysis.eq.2) then
                allocate(Strain(3),Stress(3))
                Strain = 0.0d0; Stress = 0.0d0;
                allocate(Self%StrainE(Self%Ne,3))
                Self%StrainE = 0.0d0
                allocate(Self%StressE(Self%Ne,3))
                Self%StressE = 0.0d0
            elseif (Self%DimAnalysis.eq.3) then
                allocate(Strain(6),Stress(6))
                Strain = 0.0d0; Stress = 0.0d0
                allocate(Self%StrainE(Self%Ne,6))
                Self%StrainE = 0.0d0
                allocate(Self%StressE(Self%Ne,6))
                Self%StressE = 0.0d0
            end if
            allocate(Self%StrainEnergyE(Self%Ne)) 
            Self%StrainEnergyE = 0.0d0
        end if
        ! star processing
        ! Displacement
        Self%Displacement = Transpose(reshape(Self%UGlobal,[Self%DimAnalysis,Self%N]))
        ! (calculation per element)
        do i = 1, Self%Ne, 1
            ! 1. Strain     Strain_e = Be*Ue        
            ! exx - eyy - ezz(3D) - exy - eyz(3D) - exz(3D)
            Strain = matmul(Self%BLocal(i,:,:),Self%UGlobal(Self%ConnectivityD(i,:)))
            Self%StrainE(i,:) = Strain
            ! 2. Stress      Sigma_e = D*Strain_e   
            ! Sxx - Syy - Szz(3D) - Sxy - Syz(3D) - Sxz(3D)
            Stress = matmul(D,Strain)
            Self%StressE(i,:) = Stress
            ! 3. Energy    Ee = 0.5*Sigma_e*Strain_e
            energy = (0.5d0)*dot_product(Self%StrainE(i,:),Self%StressE(i,:)) 
            Self%StrainEnergyE(i) = energy
        end do
        !call FilePrinting(Self%Displacement,'DataResults/Displacement.txt')
        !call FilePrinting(Self%StrainEnergyE,'V','DataResults/StrainEnergyE.txt')
        !call FilePrinting(Self%StrainE,'DataResults/StrainE.txt')
        !call FilePrinting(Self%StressE,'DataResults/StressE.txt')
    end subroutine ProcessingResults
end module FEA_Module