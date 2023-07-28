module FEAModule
    implicit none
    ! ----------------------------------------------------------------- !
    !     this variable groups all the information of the structure     !
    ! ----------------------------------------------------------------- !
    type                                            :: Structure        ! Name of the type of derived structure
        character(len=20)                           :: DimAnalysis      ! 2D or 3D
        character(len=20)                           :: AnalysisType     ! PlaneStress(2D), PlainStrain(2D), HookeLaw(3D)
        character(len=20)                           :: ElementType      ! tri3 tri6, cuad4, cuad8, tetra4, tetra10, hexa8, hexa20
        integer                                     :: QuadGauss        !
        integer                                     :: N                ! Number of nodes
        integer                                     :: Ne               ! Number of elements
        integer                                     :: Npe              ! Number of nodes per element
        integer                                     :: NBc              ! Number of nodes with restrictions
        integer                                     :: Npl              ! Number of point loads
        integer                                     :: Ndl              ! Number of distributed loads
        integer, Dimension(:), allocatable          :: FreeD            ! Free degrees of freedom
        integer, Dimension(:), allocatable          :: FixedD           ! Fixed degrees of freedom
        integer, dimension(:,:), allocatable        :: ConnectivityN    ! Node connectivity per element
        integer, dimension(:,:), allocatable        :: ConnectivityD    ! Degrees of freedom Connectivity per element
        integer, dimension(:,:), allocatable        :: BoundaryC        ! Boundary Conditions
        double precision                                        :: YoungModulus     ! 
        double precision                                        :: PoissonModulus   !
        double precision                                        :: Thickness        ! Only for 2D analysis (in 3D thickness = 0.0)
        double precision, dimension(:), allocatable             :: GaussPoint       ! Gauss Approximation points
        double precision, dimension(:), allocatable             :: GaussWeights     ! Gauss Approximation weights
        double precision, dimension(:), allocatable             :: FGlobal_PL       ! Global vector of point loads
        double precision, dimension(:), allocatable             :: FGlobal_DL       ! Global vector of distributed loads
        double precision, dimension(:,:), allocatable           :: Coordinates      ! Coordinates of nodes
        double precision, dimension(:,:), allocatable           :: PointLoads       ! 
        double precision, dimension(:,:), allocatable           :: DistributedLoads ! 
        double precision, dimension(:,:,:), allocatable         :: BLocal           ! Matrix of derivatives of shape functions
        double precision, dimension(:,:,:), allocatable         :: KLocal           ! Local stiffness matrix per element
        ! System of equations [K]{u}={f} in sparse format               
        integer, dimension(:), allocatable          :: index_i_KGlobal  ! Indices for rows of the global stiffness matrix
        integer, dimension(:), allocatable          :: index_j_KGlobal  ! Indices for columns of the global stiffness matrix
        integer, dimension(:,:), allocatable        :: Node_Interaction
        integer, dimension(:,:), allocatable        :: Element_Interaction 
        double precision, dimension(:), allocatable             :: value_KGlobal    ! Value of the global stiffness matrix at position i,j
        double precision, dimension(:), allocatable             :: value_FGlobal    ! Global load vector
        double precision, dimension(:), allocatable             :: value_UGlobal    ! Global displacement vector
        ! Results
        double precision, dimension(:), allocatable             :: UGlobal          ! Global displacement vector
        double precision, dimension(:), allocatable             :: StrainEnergyE    ! Strain energy per element
        double precision, dimension(:), allocatable             :: StrainEnergyN    ! Strain energy per node
        double precision, dimension(:,:), allocatable           :: StrainE         ! Strain per element
        double precision, dimension(:,:), allocatable           :: StrainN         ! Strain per node
        double precision, dimension(:,:), allocatable           :: StressE          ! Stress per element
        double precision, dimension(:,:), allocatable           :: StressN          ! Stress per node
    contains
        procedure                                   :: SetDimAnalysis
        procedure                                   :: SetAnalysisType
        procedure                                   :: SetElementType
        procedure                                   :: SetYoungModulus
        procedure                                   :: SetPoissonModulus
        procedure                                   :: SetThickness
        procedure                                   :: SetGaussAprox
        procedure                                   :: SetCoordinates
        procedure                                   :: SetConnectivity
        procedure                                   :: SetBondaryConditions
        procedure                                   :: SetPointLoads
        procedure                                   :: SetDistributedLoads
        procedure                                   :: UploadStructure
        procedure                                   :: HSL87Solver
        procedure                                   :: ProcessingResults
    end type
contains
    ! ----------------------------------------------------------------- !
    !      base subroutines for reading, writting  and operations       !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine ReadDoublePrecisionFile(Path,nr,nc,Array)
        implicit none
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,ios
        integer, intent(out)                                        :: nr
        integer, intent(in)                                         :: nc
        double precision, dimension(:,:), allocatable, intent(out)              :: Array
        !read file
        open(unit=1, file=Path, iostat=ios, status='old', Action='read')
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=1, fmt=*)
            read(unit=1, fmt=*) nr
            allocate(Array(nr,nc))
            read(unit=1, fmt=*)
            do i = 1, nr, 1
                read(unit=1, fmt=*) (Array(i,j), j=1, nc, 1)
            end do 
        close(1)
    end subroutine ReadDoublePrecisionFile
    ! listo
    subroutine ReadIntegerFile(Path,nr,nc,Array)
        implicit none
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,ios
        integer, intent(out)                                        :: nr
        integer, intent(in)                                         :: nc
        integer, dimension(:,:), allocatable, intent(out)           :: Array
        !read file
        open(unit=1, file=Path, iostat=ios, status='old', Action='read')
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=1, fmt=*)
            read(unit=1, fmt=*) nr
            allocate(Array(nr,nc))
            read(unit=1, fmt=*)
            do i = 1, nr, 1
                read(unit=1, fmt=*) (Array(i,j), j=1, nc, 1)
            end do 
        close(1)
    end subroutine ReadIntegerFile
    ! listo
    function Area(Coordinates,Type) result(AreaAprox)
        character(len=*), intent(in)                                :: Type
        double precision, dimension(:), allocatable                 :: vector1, vector2, vector3, vector4
        double precision, dimension(:), allocatable                 :: Av1, Av2
        double precision, dimension(:,:), allocatable, intent(in)   :: Coordinates
        double precision                                            :: AreaAprox
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
    ! listo
    function CrossProduct(v1, v2) result(v3)
        double precision, dimension(3), intent(in)    :: v1, v2
        double precision, dimension(3)                :: v3
        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
    end function CrossProduct
    ! listo
    function Determinant(A) result(Det)
        implicit none
        double precision                              :: Det
        double precision, intent(in)                  :: A(:,:)
        if (size(A(1,:)).eq.2) then
            det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        elseif (size(A(1,:)).eq.3) then
            det = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)&
                 -A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)&
                 +A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
        else
            write(*,*) 'ERROR, revisar subrutina para determinante de matriz'
            call exit
        end if
    end function Determinant
    ! listo
    function Inverse(A) result(Ainv)
        ! this function use lapack library
        implicit none
        double precision,intent(in)         :: A(:,:)
        double precision                    :: Ainv(size(A,1),size(A,2))
        ! internal variables
        double precision                    :: work(size(A,1))
        integer                             :: n,info,ipiv(size(A,1))
        Ainv = A
        n = size(A,1)
        call DGETRF(n,n,Ainv,n,ipiv,info)
        !call SGETRF(n,n,Ainv,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        call DGETRI(n,Ainv,n,ipiv,work,n,info)
        !call SGETRI(n,Ainv,n,ipiv,work,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function Inverse
    ! listo
    function Norm(vec) result(Ans)
        implicit none
        double precision, dimension(:), allocatable, intent(in)   :: vec
        ! internal variables
        double precision                                          :: sum_of_squares, Ans
        integer                                                   :: i
        sum_of_squares = 0.0d0
        do i = 1, size(vec)
        sum_of_squares = sum_of_squares + vec(i)**2
        end do
        Ans = sqrt(sum_of_squares)
    end function Norm
    ! listo
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
    ! escritura
    subroutine WriteIntegerArray(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                        :: path
        integer, dimension(:,:), allocatable, intent(in)    :: Array
        ! internal value
        integer                                             :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array,1), 1
                write(unit=1, fmt=*) (Array(i,j),j=1,size(Array,2),1)
            end do    
        close(1)
    end subroutine WriteIntegerArray
    ! escritura
    subroutine WriteDoublePrecisionArray(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                                :: path
        double precision, dimension(:,:), allocatable, intent(in)   :: Array
        ! internal value
        integer                                                     :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array,1), 1
                write(unit=1, fmt=*) (Array(i,j),j=1,size(Array,2),1)
            end do    
        close(1)
    end subroutine WriteDoublePrecisionArray
    ! escritura
    subroutine WriteIntegerVector(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                        :: path
        integer, dimension(:), allocatable, intent(in)      :: Array
        ! internal value
        integer                                             :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array), 1
                write(unit=1, fmt=*) Array(i)
            end do    
        close(1)
    end subroutine WriteIntegerVector
    ! escritura
    subroutine WriteDoublePrecisionVector(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                                :: path
        double precision, dimension(:), allocatable, intent(in)     :: Array
        ! internal value
        integer                                                     :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array), 1
                write(unit=1, fmt=*) Array(i)
            end do    
        close(1)
    end subroutine WriteDoublePrecisionVector
    ! ----------------------------------------------------------------- !
    !       subroutines to define the information required for FEA      !
    ! ----------------------------------------------------------------- !
    ! input parameter
    subroutine SetDimAnalysis(Self,DimAnalysis)
        implicit none
        class(Structure), intent(inout)                              :: Self
        character(len=*), intent(in)                                 :: DimAnalysis
        Self%DimAnalysis = DimAnalysis
    end subroutine
    ! input parameter
    subroutine SetAnalysisType(Self,AnalysisType)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: AnalysisType
        Self%AnalysisType = AnalysisType
    end subroutine SetAnalysisType
    ! input parameter
    subroutine SetElementType(Self,ElementType)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: ElementType
        Self%ElementType = ElementType
        if ( Self%ElementType.eq.'tria3' )   then; Self%Npe = 3 ;  end if
        if ( Self%ElementType.eq.'tria6' )   then; Self%Npe = 6 ;  end if
        if ( Self%ElementType.eq.'quad4' )   then; Self%Npe = 4 ;  end if
        if ( Self%ElementType.eq.'quad8' )   then; Self%Npe = 8 ;  end if
        if ( Self%ElementType.eq.'tetra4' )  then; Self%Npe = 4 ;  end if
        if ( Self%ElementType.eq.'tetra10' ) then; Self%Npe = 10 ; end if
        if ( Self%ElementType.eq.'hexa8' )   then; Self%Npe = 8 ;  end if
        if ( Self%ElementType.eq.'hexa20' )  then; Self%Npe = 20 ; end if
    end subroutine SetElementType
    ! input parameter
    subroutine SetThickness(Self,Thickness)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(in)                                :: Thickness
        Self%Thickness = Thickness
    end subroutine SetThickness
    ! input parameter
    subroutine SetYoungModulus(Self,YoungModulus)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(in)                                :: YoungModulus
        Self%YoungModulus = YoungModulus
    end subroutine SetYoungModulus
    ! input parameter
    subroutine SetPoissonModulus(Self,PoissonModulus)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(in)                                :: PoissonModulus
        Self%PoissonModulus = PoissonModulus
    end subroutine SetPoissonModulus
    ! input parameter
    subroutine SetGaussAprox(Self,Gauss)
        implicit none
        class(Structure), intent(inout)                             :: Self 
        integer, intent(in)                                         :: Gauss
        Self%QuadGauss = Gauss
        select case (Gauss)
            case (1)
                allocate(Self%GaussPoint(1));       Self%GaussPoint = [0.00000000]
                allocate(Self%GaussWeights(1));   Self%GaussWeights = [2.00000000]
            case (2)
                allocate(Self%GaussPoint(2));       Self%GaussPoint = [-0.57735026,0.57735026]
                allocate(Self%GaussWeights(2));   Self%GaussWeights = [1.00000000,1.00000000]
            case (3)
                allocate(Self%GaussPoint(3));       Self%GaussPoint = [-0.77459666,0.00000000,0.77459666]
                allocate(Self%GaussWeights(3));   Self%GaussWeights = [0.55555555,0.88888888,0.55555555]
            case (4)
                allocate(Self%GaussPoint(4));       Self%GaussPoint = [-0.86113631,-0.33998104,0.33998104,0.86113631]
                allocate(Self%GaussWeights(4));   Self%GaussWeights = [0.34785484,0.65214515,0.65214515,0.34785484]
            case (5)
                allocate(Self%GaussPoint(5));       Self%GaussPoint = [-0.90617984,-0.53846931,0.00000000,0.53846931,0.90617984]
                allocate(Self%GaussWeights(5));   Self%GaussWeights = [0.23692688,0.47862867,0.56888888,0.47862867,0.23692688]
        end select
    end subroutine SetGaussAprox
    ! input parameter
    subroutine SetCoordinates(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        ! internal variables
        integer                                                     :: dim
        if ( Self%DimAnalysis.eq.'2D' ) then; dim = 2 ; end if
        if ( Self%DimAnalysis.eq.'3D' ) then; dim = 3 ; end if
        call ReadDoublePrecisionFile(Path,Self%N,3,Self%Coordinates)
    end subroutine SetCoordinates
    ! input parameter
    subroutine SetConnectivity(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        ! internal variables
        integer                                                     :: i,j,k,dim
        ! calculate ConnectivityD
        if ( Self%DimAnalysis.eq.'2D' ) then; dim = 2 ; end if
        if ( Self%DimAnalysis.eq.'3D' ) then; dim = 3 ; end if
        call ReadIntegerFile(Path,Self%Ne,Self%Npe,Self%ConnectivityN)
        allocate(Self%ConnectivityD(Self%Ne,Self%Npe*dim))
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = dim-1, 0, -1
                    Self%ConnectivityD(i,j*dim - k) = Self%ConnectivityN(i,j)*dim - k
                end do
            end do
        end do
        call WriteIntegerArray('DataStructure/AdditionalData/ConnectivityD.txt',Self%ConnectivityD)
    end subroutine SetConnectivity
    ! input parameter
    subroutine SetBondaryConditions(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        ! internal variables
        integer                                                     :: dumb
        integer                                                     :: i,j,k,l,dim 
        integer, dimension(:), allocatable                          :: InPosition
        logical, dimension(:), allocatable                          :: InLogical
        if ( Self%DimAnalysis.eq.'2D' ) then; dim = 2 ; end if
        if ( Self%DimAnalysis.eq.'3D' ) then; dim = 3 ; end if
        call ReadIntegerFile(Path,Self%NBc,(dim+1),Self%BoundaryC)
        ! Get the degress of freedom
        i = (Self%N)*(dim) - count(Self%BoundaryC(:,2:).eq.1)
        j = count(Self%BoundaryC(:,2:).eq.1)
        k = 1
        allocate(Self%FreeD(i))
        allocate(Self%FixedD(j))
        Self%FreeD = 0
        Self%FixedD = 0
        ! constrained?
        do i = 1, Self%NBc, 1
            do j = 1, Dim, 1
                if (Self%BoundaryC(i,j+1).eq.1) then
                    Self%FixedD(k) = Self%BoundaryC(i,1)*dim - (dim-j)
                    k = k + 1
                else
                    cycle
                end if
            end do
        end do
        ! free!
        k = 1
        do i = 1, Self%N*dim, 1
            if (any(Self%FixedD.eq.i)) then
                cycle
            else
                Self%FreeD(k) = i
                k = k + 1
            end if
        end do
        ! printing
        call WriteIntegerVector('DataStructure/AdditionalData/FreeDof.txt',Self%FreeD)
        call WriteIntegerVector('DataStructure/AdditionalData/FixedDof.txt',Self%FixedD)
        ! get de element and DoF incidence for each DoF (only the free degres)
        j = size(Self%FreeD)
        ! maximum of 20 elements per node
        allocate(Self%Element_Interaction(j,20))
        Self%Element_Interaction = 0
        ! considering the quad20 20*20 
        allocate(Self%Node_Interaction(j,400))
        Self%Node_Interaction = 0
        do i = 1, size(Self%FreeD), 1
            ! element
            InLogical = any(Self%ConnectivityD.eq.Self%FreeD(i),2)
            InPosition = pack([(k,k=1,Self%Ne)],InLogical)
            j = count(InLogical)
            Self%Element_Interaction(i,1:j) = InPosition 
            ! other DoF
            InPosition = reshape(Self%ConnectivityD(Self%Element_Interaction(i,1:j),:),[j*(Self%Npe)*dim])
            call sort(InPosition,j*Self%Npe*dim)
            j = size(InPosition)
            Self%Node_Interaction(i,1:j) = InPosition
        end do
        call WriteIntegerArray('DataStructure/AdditionalData/FreeDofInteraction_Dof.txt',Self%Node_Interaction)
        call WriteIntegerArray('DataStructure/AdditionalData/FreeDofInteraction_Element.txt',Self%Element_Interaction)
    end subroutine SetBondaryConditions
    ! input parameter
    subroutine SetPointLoads(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        ! internal variables
        integer                                                     :: i,j,dim,node 
        if ( Self%DimAnalysis.eq.'2D' ) then; dim = 2 ; end if
        if ( Self%DimAnalysis.eq.'3D' ) then; dim = 3 ; end if
        call ReadDoublePrecisionFile(Path,Self%Npl,(dim+1),Self%PointLoads)
        allocate(Self%FGlobal_PL(Self%N*dim)) ; Self%FGlobal_PL = 0.0d0
        ! assembly load vector
        do i = 1, Self%Npl, 1
            node = int(Self%PointLoads(i,1))
            do j = 1, dim, 1
                Self%FGlobal_PL(node*dim-(dim-j)) = Self%FGlobal_PL(node*dim-(j-dim)) + Self%PointLoads(i,j+1)
            end do
        end do
        ! writting
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/GlobalPointLoads.txt',Self%FGlobal_PL)
    end subroutine SetPointLoads
    ! input parameter
    subroutine SetDistributedLoads(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        ! internal variables
        integer                                                     :: n,i,j,k,Npf,dim,node
        double precision                                            :: FaceArea, length
        double precision, dimension(:), allocatable                 :: Vector, ForcePerNode
        double precision, dimension(:,:), allocatable               :: Coordinates
        if ( Self%DimAnalysis.eq.'2D' ) then; dim = 2 ; end if
        if ( Self%DimAnalysis.eq.'3D' ) then; dim = 3 ; end if
        ! calculation and assembly of load vectors
        allocate(Self%FGlobal_DL(Self%N*dim)); Self%FGlobal_DL = 0.0d0;
        allocate(ForcePerNode(dim));              ForcePerNode = 0.0d0;
        allocate(Vector(dim));                          vector = 0.0d0;
        ! nodes per face
        if ( Self%ElementType.eq.'tria3'.or.Self%ElementType.eq.'quad4' ) then; Npf = 2; end if
        if ( Self%ElementType.eq.'tria6'.or.Self%ElementType.eq.'quad8' ) then; Npf = 3; end if
        if ( Self%ElementType.eq.'tetra4' ) then;  Npf = 3; n=3; allocate(Coordinates(3,3)); end if
        if ( Self%ElementType.eq.'tetra10' ) then; Npf = 6; n=3; allocate(Coordinates(3,3)); end if
        if ( Self%ElementType.eq.'hexa8' ) then;   Npf = 4; n=4; allocate(Coordinates(4,3)); end if
        if ( Self%ElementType.eq.'hexa20' ) then;  Npf = 8; n=4; allocate(Coordinates(4,3)); end if
        ! reading file
        call ReadDoublePrecisionFile(Path,Self%Ndl,(Npf+dim),Self%DistributedLoads)
        select case (dim)
            case (2)
                do i = 1, Self%Ndl, 1
                    Vector = Self%Coordinates(Self%DistributedLoads(i,2),1:dim) - Self%Coordinates(Self%DistributedLoads(i,1),1:dim)
                    length = Norm(Vector)
                    FaceArea = length*Self%Thickness
                    ForcePerNode = (FaceArea/Npf)*(Self%DistributedLoads(i,(Npf+1):))
                    do j = 1, Npf, 1
                        node = int(Self%DistributedLoads(i,j))
                        do k = 1, dim, 1
                            Self%FGlobal_DL(node*dim-(dim-k)) = ForcePerNode(k)
                        end do
                    end do
                end do
            case (3)
                do i = 1, Self%Ndl, 1
                    Coordinates = Self%Coordinates(int(Self%DistributedLoads(i,1:n)),:)
                    FaceArea = Area(Coordinates,Self%ElementType)
                    ForcePerNode = (FaceArea/Npf)*Self%DistributedLoads(i,(Npf+1):)
                    do j = 1, Npf, 1
                        node = int(Self%DistributedLoads(i,j))
                        do k = 1, dim, 1
                            Self%FGlobal_DL(node*dim-(dim-k)) = ForcePerNode(k)
                        end do
                    end do
                end do
        end select
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/GlobalDistributedLoads.txt',Self%FGlobal_DL)
    end subroutine SetDistributedLoads
    ! ----------------------------------------------------------------- !
    !              subroutines that define the FEM procedure            !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine DiffFormFunction(Self,DiffFunction,e,n,z)
        ! dN1/de     dN2/de     dN3/de      ....     dNn/de
        ! dN1/dn     dN2/dn     dN3/dn      ....     dNn/dn
        ! dN1/dz     dN2/dz     dN3/dz      ....     dNn/dz (caso 3D)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                             :: e, n
        double precision, intent(inout), optional                   :: z
        double precision, dimension(:,:), allocatable, intent(out)  :: DiffFunction
        if (Self%ElementType.eq.'tria3') then
            allocate(DiffFunction(2,3))
            DiffFunction(1,:) = [-1.0d0,1.0d0,0.0d0]
            DiffFunction(2,:) = [-1.0d0,0.0d0,1.0d0]
        elseif (Self%ElementType.eq.'tria6') then
            allocate(DiffFunction(2,6))
            DiffFunction(1,:) = [4.0d0*e+4.0d0*n-3.0d0,4.0d0*e-1.0d0,0.0d0,4.0d0-4.0d0*n-8.0d0*e,4.0d0*n,-4.0d0*n]
            DiffFunction(2,:) = [4.0d0*e+4.0d0*n-3.0d0,0.0d0,4.0d0*n-1.0d0,-4.0d0*e,4.0d0*e,4.0d0-8.0d0*n-4.0d0*e]
        elseif (Self%ElementType.eq.'quad4') then
            allocate(DiffFunction(2,4))
            DiffFunction(1,:) = [n/4.0d0-1.0d0/4.0d0,1.0d0/4.0d0-n/4.0d0,n/4.0d0+1.0d0/4.0d0,-n/4.0d0-1.0d0/4.0d0]
            DiffFunction(2,:) = [e/4.0d0-1.0d0/4.0d0,-e/4.0d0-1.0d0/4.0d0,e/4.0d0+1.0d0/4.0d0,1.0d0/4.0d0-e/4.0d0]
        elseif (Self%ElementType.eq.'quad8') then
            allocate(DiffFunction(2,8))
            DiffFunction(1,:) = [-(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)-((n-1.0d0)*(e+n+1.0d0))/4.0d0, &
                            ((n-1.0d0)*(n-e+1.0d0))/4.0d0-(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)+((n+1.0d0)*(e+n-1.0d0))/4.0d0, &
                            (e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)+((n+1.0d0)*(e-n+1.0d0))/4.0d0,  &
                            e*(n-1.0d0),1.0d0/2.0d0-n**2.0d0/2.0d0,-e*(n+1.0d0),n**2.0d0/2.0d0-1.0d0/2.0d0]
            DiffFunction(2,:) = [-(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)-(e/4.0d0-1.0d0/4.0d0)*(e+n+1.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n-e+1.0d0)+(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)+(e/4.0d0+1.0d0/4.0d0)*(e+n-1.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(e-n+1.0d0)-(e/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            e**2.0d0/2.0d0-1.0d0/2.0d0,-n*(e+1.0d0),1.0d0/2.0d0-e**2.0d0/2.0d0,n*(e-1.0d0)]
        elseif (Self%ElementType.eq.'tetra4') then
            allocate(DiffFunction(3,4))
            DiffFunction(1,:) = [1.0d0,0.0d0,0.0d0,-1.0d0]
            DiffFunction(2,:) = [0.0d0,1.0d0,0.0d0,-1.0d0]
            DiffFunction(3,:) = [0.0d0,0.0d0,1.0d0,-1.0d0]
        elseif (Self%ElementType.eq.'tetra10') then
            allocate(DiffFunction(3,10))
            DiffFunction(1,:) = [4.0d0*e-1.0d0,0.0d0,0.0d0,4.0d0*e+4.0d0*n+4.0d0*z-3.0d0,4.0d0*n,0.0d0,4.0d0*z, &
                                4.0d0-4.0d0*n-4.0d0*z-8.0d0*e,-4.0d0*n,0.0d0]
            DiffFunction(2,:) = [0.0d0,4.0d0*n-1.0d0,0.0d0,4.0d0*e+4.0d0*n+4.0d0*z-3.0d0,4.0d0*e,4.0d0*z,0.0d0,-4.0d0*e, &
                                4.0d0-8.0d0*n-4.0d0*z-4.0d0*e,0.0d0]
            DiffFunction(3,:) = [0.0d0,0.0d0,4.0d0*z-1.0d0,4.0d0*e+4.0d0*n+4.0d0*z-3.0d0,0.0d0,4.0d0*n,4.0d0*e,-4.0d0*e, &
                                -4.0d0*n,0.0d0]
        elseif (Self%ElementType.eq.'hexa8') then
            allocate(DiffFunction(3,8))
            DiffFunction(1,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)+((n-1.0d0)*(e+n+z))/8.0d0, &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)-((n-1.0d0)*(e+n+z))/8.0d0, &
                             (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)+((n+1.0d0)*(e+n+z))/8.0d0, &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)-((n+1.0d0)*(e+n+z))/8.0d0, &
                            -(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)-((n-1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                             (e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)+((n-1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                            -(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)-((n+1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                             (e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)+((n+1.0d0)*(e+n+z-2.0d0))/8.0d0]
            DiffFunction(2,:) = [(e/8.0d0-1.0d0/8.0d0)*(e+n+z)+(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(e+n+z)-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                             (e/8.0d0+1.0d0/8.0d0)*(e+n+z)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(e+n+z)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(e+n+z-2.0d0)-(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0), &
                             (e/8.0d0+1.0d0/8.0d0)*(e+n+z-2.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(e+n+z-2.0d0)-(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0), &
                             (e/8.0d0-1.0d0/8.0d0)*(e+n+z-2.0d0)+(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)]
            DiffFunction(3,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0),-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                             (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0),-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0),(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0),(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)]
        elseif (Self%ElementType.eq.'hexa20') then
            allocate(DiffFunction(3,20))
            DiffFunction(1,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(z-3.0d0)+((n-1.0d0)*(z-3.0d0)*(e+n+z))/8.0d0, &
                            -((n-1.0d0)*(2.0d0*e+z-3.0d0)*(e+n+z))/8.0d0-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+z-3.0d0) &
                            -2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z), &
                            ((n+1.0d0)*(e+n+z)*(2.0d0*e+2.0d0*n+z-3.0d0))/8.0d0 &
                            +(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+z-3.0d0) &
                            +2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            -((n+1.0d0)*(2.0d0*n+z-3.0d0)*(e+n+z))/8.0d0-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*n+z-3.0d0), &
                            ((n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0)*(e+n+z-2.0d0))/8.0d0 &
                            +(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0) &
                            +2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            -((n-1.0d0)*(2.0d0*n+z+1.0d0)*(e+n+z-2.0d0))/8.0d0-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*n+z+1.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(z+1.0d0)+((n+1.0d0)*(z+1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                            -((n+1.0d0)*(2.0d0*e+z+1.0d0)*(e+n+z-2.0d0))/8.0d0-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+z+1.0d0) &
                            -2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z-2.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n-1.0d0)+(e*(n-1.0d0)*(e+n+z))/2.0d0, &
                            -((n**2.0d0-1.0d0)*(e+n+z))/4.0d0-(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0)-(e*(n+1.0d0)*(e+n+z))/2.0d0, &
                            ((n**2.0d0-1.0d0)*(e+n+z))/4.0d0+(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            0.0d0, &
                            ((n**2.0d0-1.0d0)*(e+n+z-2.0d0))/4.0d0+(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0)+(e*(n+1.0d0)*(e+n+z-2.0d0))/2.0d0, &
                            -((n**2.0d0-1.0d0)*(e+n+z-2.0d0))/4.0d0-(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(((e+n+z-1.0d0)**2.0d0-1.0d0)*(n-1.0d0))/4.0d0 &
                            -(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (((e+n+z-1.0d0)**2.0d0-1.0d0)*(n-1.0d0))/4.0d0 &
                            +(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            -(((e+n+z-1.0d0)**2.0d0-1.0d0)*(n+1.0d0))/4.0d0 &
                            -(e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (((e+n+z-1.0d0)**2.0d0-1.0d0)*(n+1.0d0))/4.0d0 &
                            +(e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0)]
            DiffFunction(2,:) = [(e/8.0d0-1.0d0/8.0d0)*(z-3.0d0)*(e+n+z)+(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(z-3.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(2.0d0*e+z-3.0d0)*(e+n+z) &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+z-3.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(e+n+z)*(2.0d0*e+2.0d0*n+z-3.0d0) &
                            +(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+z-3.0d0) &
                            +2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            -(e/8.0d0-1.0d0/8.0d0)*(2.0d0*n+z-3.0d0)*(e+n+z)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*n+z-3.0d0) &
                            -2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            (e/8.0d0-1.0d0/8.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0)*(e+n+z-2.0d0) &
                            +(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0) &
                            +2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(2.0d0*n+z+1.0d0)*(e+n+z-2.0d0) &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*n+z+1.0d0) &
                            -2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(z+1.0d0)*(e+n+z-2.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(z+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(2.0d0*e+z+1.0d0)*(e+n+z-2.0d0) &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+z+1.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(e+n+z)+(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n-1.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0)-2.0d0*n*(e/4.0d0+1.0d0/4.0d0)*(e+n+z), &
                            -(e**2.0d0/4.0d0-1.0d0/4.0d0)*(e+n+z)-(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0)+2.0d0*n*(e/4.0d0-1.0d0/4.0d0)*(e+n+z), &
                            0.0d0, &
                            (e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0)+2.0d0*n*(e/4.0d0+1.0d0/4.0d0)*(e+n+z-2.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(e+n+z-2.0d0)+(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0)-2.0d0*n*(e/4.0d0-1.0d0/4.0d0)*(e+n+z-2.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            -(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            +(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            -(e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            +(e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0)]
            DiffFunction(3,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(z-3.0d0)+(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z), &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+z-3.0d0)-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z), &
                            (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+z-3.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*n+z-3.0d0)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            (e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0) &
                            +(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*n+z+1.0d0)-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(z+1.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z-2.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+z+1.0d0)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z-2.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n-1.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            0.0d0, &
                            (e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0)]
        end if
    end subroutine DiffFormFunction
    ! listo
    subroutine ElasticityTensor(Self,ETensor)
        implicit none
        class(Structure), intent(inout)                             :: Self    
        double precision, dimension(:,:), allocatable, intent(inout)            :: ETensor
        double precision                                                        :: E,V,Constant
        E = Self%YoungModulus
        V = Self%PoissonModulus    
        if (Self%AnalysisType.eq.'PlaneStress') then
            allocate(ETensor(3,3))
            Constant = E/(1.0d0-V**2d0)
            ETensor(1,:) = [1.0d0,V,0.0d0]
            ETensor(2,:) = [V,1.0d0,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0-V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'PlaneStrain') then
            allocate(ETensor(3,3))
            Constant = E/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0-V,V,0.0d0]
            ETensor(2,:) = [V,1.0d0-V,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'HookeLaw') then
            allocate(ETensor(6,6))
            Constant = E/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0-V,V,V,0.0d0,0.0d0,0.0d0]
            ETensor(2,:) = [V,1.0d0-V,V,0.0d0,0.0d0,0.0d0]
            ETensor(3,:) = [V,V,1.0d0-V,0.0d0,0.0d0,0.0d0]
            ETensor(4,:) = [0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0,0.0d0,0.0d0]
            ETensor(5,:) = [0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0,0.0d0]
            ETensor(6,:) = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0]
            ETensor = Constant*ETensor
        end if
    end subroutine ElasticityTensor
    ! listo
    subroutine SetKlocal(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(inout)                                 :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)      :: DensityVector 
        ! internal variables
        integer                                                         :: el,i,j,k,l,m,ios
        integer                                                         :: row,col,loc1,loc2
        double precision                                                :: e,n,z,w1,w2,w3,DetJacobian,Kl
        double precision, dimension(:,:), allocatable                   :: Be,Jacobian,InvJacobian,D
        double precision, dimension(:,:), allocatable                   :: DiffN,DiffNXY,ElementCoordinates
        ! ---- 1st part ----
        ! in this section the local stiffness matrix (dense) of each of the elements is calculated using gauss 
        ! quadrature and each of them is stored in the array KLocal[element,[K_x,K_y]].
        open(unit=1,file='DataStructure/AdditionalData/StiffnessMatrixArray.txt',iostat=ios,status="replace",action="write")
        ! this is just to check
        write(unit=*, fmt=*) '(Getting Local Stiffness matrix)'
        if (Self%DimAnalysis.eq.'2D') then
            if (allocated(Self%KLocal).and.allocated(Self%BLocal)) then
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            else
                allocate(Self%KLocal(Self%Ne,Self%Npe*2,Self%Npe*2)); Self%KLocal = 0.0d0;
                allocate(Self%BLocal(Self%Ne,3,2*Self%Npe));          Self%BLocal = 0.0d0;
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
                        ! Jacobian
                        do k = 1, 2, 1
                            do l = 1, 2, 1
                                Jacobian(k,l) = dot_product(DiffN(k,:),ElementCoordinates(:,l))
                            end do
                        end do
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        ! Be
                        Be = 0.0d0
                        do k = 1, size(DiffN,2), 1
                            Be(1,2*k-1) = DiffNxy(1,k)
                            !Be(1,2*k) = 0.0d0
                            !Be(2,2*k-1) = 0.0d0
                            Be(2,2*k) = DiffNxy(2,k)
                            Be(3,2*k-1) = DiffNxy(2,k)
                            Be(3,2*k) = DiffNxy(1,k)
                        end do
                        ! Be Local
                        Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + Be*w1*w2
                        ! K Local
                        Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + DetJacobian*(matmul(transpose(Be),matmul(D,Be)))*&
                                                                    w1*w2*(Self%Thickness)
                        deallocate(DiffN)
                    end do
                end do
                !write(unit=1, fmt=*) reshape(Self%KLocal(el,:,:),[(Self%Npe*2)*(Self%Npe*2)])
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.'3D') then
            if (allocated(Self%KLocal).and.allocated(Self%BLocal)) then
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            else
                allocate(Self%KLocal(Self%Ne,Self%Npe*3,Self%Npe*3)); Self%KLocal = 0.0d0;
                allocate(Self%BLocal(Self%Ne,6,3*Self%Npe));          Self%BLocal = 0.0d0;
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
                            ! Jacobian
                            do l = 1, 2, 1
                                do m = 1, 2, 1
                                    Jacobian(l,m) = dot_product(DiffN(l,:),ElementCoordinates(:,m))
                                end do
                            end do
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            ! Be
                            Be = 0.0d0
                            do l = 1, size(DiffN,2), 1
                                Be(1,3*l-2) = DiffNxy(1,l)
                                !Be(1,3*l-1) = 0.0d0
                                !Be(1,3*l) = 0.0d0
                                !Be(2,3*l-2) = 0.0d0
                                Be(2,3*l-1) = DiffNxy(2,l)
                                !Be(2,3*l) = 0.0d0
                                !Be(3,3*l-2) = 0.0d0
                                !Be(3,3*l-1) = 0.0d0
                                Be(3,3*l) = DiffNxy(3,l)
                                Be(4,3*l-2) = DiffNxy(2,l)
                                Be(4,3*l-1) = DiffNxy(1,l)
                                !Be(4,3*l) = 0.0d0
                                !Be(5,3*l-2) = 0.0d0
                                Be(5,3*l-1) = DiffNxy(3,l)
                                Be(5,3*l) = DiffNxy(2,l)
                                Be(6,3*l-2) = DiffNxy(3,l)
                                !Be(6,3*l-1) = 0.0d0
                                Be(6,3*l) = DiffNxy(1,l)
                            end do
                            ! Be Local
                            Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + Be*w1*w2*w3
                            ! K Local
                            Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + (-1)*DetJacobian*(matmul(transpose(Be),matmul(D,Be)))*&
                                                                            w1*w2*w3
                            deallocate(DiffN)
                        end do
                    end do
                end do
                !write(unit=1, fmt=*) reshape(Self%KLocal(el,:,:),[Self%Npe*3**2])
                deallocate(D)
            end do
        end if
        close(1)        
        deallocate(Be,Jacobian,InvJacobian,DiffNXY,ElementCoordinates)
        ! ---- 2nd part ----
        ! in this section the sparse stiffness matrix (free degrees only) is assembled for further solution
        write(unit=*, fmt=*) '(Assembling stiffness sparse matrix)'
        Kl = 0.0d0
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/KSparseIndex.txt',iostat=ios,status="replace",action="write")
        open(unit=2,file='DataStructure/AdditionalData/SparseSystem/KSparseValue.txt',iostat=ios,status="replace",action="write")
        m = 0
        do j = 1, size(Self%FreeD), 1                                   ! column
            col = Self%FreeD(j)
            do i = 1, size(Self%Node_Interaction,2), 1                  ! row
                row = Self%Node_Interaction(j,i)
                if (row.eq.0) exit                                      ! exit (doenst interact)
                if (any(Self%FixedD.eq.row)) cycle                      ! discard (fixed DoF)
                if (col.gt.row) cycle                                   ! discard (bottom diagonal only)
                m = m + 1
                do k = 1,size(Self%Element_Interaction,2), 1            ! element
                    el = Self%Element_Interaction(j,k)
                    if (el.eq.0) exit                                   ! exit (no more elements)
                    loc1 = findloc(Self%ConnectivityD(el,:),row,1)
                    if (loc1.eq.0) cycle                                ! discard (doesnt exist in element)
                    loc2 = findloc(Self%ConnectivityD(el,:),col,1)
                    Kl = Kl + Self%KLocal(el,loc1,loc2)
                end do
                ! adding the elements to the vectors
                write(unit=1, fmt=*) findloc(Self%FreeD,row,dim=1), findloc(Self%FreeD,col,dim=1)
                write(unit=2, fmt=*) Kl
                Kl = 0d0                                                 ! reset
            end do
        end do
        close(2)
        close(1)
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/Resume.txt',iostat=ios,status="replace",action="write")
            write(unit=1,fmt=*) size(Self%FreeD), m 
        close(1)
    end subroutine SetKlocal
    ! listo
    subroutine UploadStructure(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                             :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)  :: DensityVector 
        write(unit=*, fmt=*) '1. Uploading structure'
        ! 1. load vector
        allocate(Self%value_FGlobal(size(Self%FreeD)));     Self%value_FGlobal = 0.0d0
        Self%value_FGlobal = Self%FGlobal_PL(Self%FreeD) + Self%FGlobal_DL(Self%FreeD)
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/SparseSystem/LoadVectorSparse.txt',Self%value_FGlobal)
        ! 2. Stifsness matrix
        call SetKlocal(Self,DensityVector,PenalFactor)
    end subroutine UploadStructure
    ! listo
    subroutine SparseSolver(Self)
        use hsl_ma87_double
        use hsl_mc68_double
        use hsl_mc69_double
        implicit none
        class(Structure), intent(inout)                             :: Self
        ! HSL-MA87 variables
        !integer, parameter                                          :: wp = kind(0d0)
        type(mc68_control)                                          :: control68
        type(mc68_info)                                             :: info68
        type(ma87_keep)                                             :: keep
        type(ma87_control)                                          :: control
        type(ma87_info)                                             :: info
        integer                                                     :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                          :: crow, ccol, ptr, row, order, map
        double precision, dimension(:), allocatable                 :: cval, val, x
        ! Organizing variables
        n = size(Self%value_FGlobal);         ne = size(Self%index_j_KGlobal);
        allocate(crow(ne));                 crow = Self%index_i_KGlobal;
        allocate(ccol(ne));                 ccol = Self%index_j_KGlobal;
        allocate(cval(ne));                 cval = Self%value_KGlobal;
        allocate(x(n));                        x = Self%value_FGlobal;
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, n, n, ne, crow, ccol, &
            ptr, row, flag, val_in=cval, val_out=val, lmap=lmap, map=map)
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
        Self%value_UGlobal = X
        ! Finalize
        call ma87_finalise(keep, control)
    end subroutine SparseSolver
    ! listo
    subroutine stop_on_bad_flag(context, flag)
        character(len=*), intent(in)                                 :: context
        integer, intent(in)                                          :: flag
        if(flag.eq.0) return
        write(*,*) "Failure during ", context, " with flag = ", flag
        stop
    end subroutine stop_on_bad_flag
    ! listo
    subroutine HSL87Solver(Self)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, dimension(:,:), allocatable                   :: Displacement
        integer                                                         :: i,dim,n,ne,ios
        if (Self%DimAnalysis.eq.'2D') then;     dim=2;
        elseif (Self%DimAnalysis.eq.'3D') then; dim=3;
        end if
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/Resume.txt',iostat=ios, status="old", action="read")
            read(1,*) n, ne
        close(1)
        allocate(Self%index_i_KGlobal(ne));              Self%index_i_KGlobal = 0;
        allocate(Self%index_j_KGlobal(ne));              Self%index_j_KGlobal = 0;
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/KSparseIndex.txt',iostat=ios, status="old", action="read")
            do i = 1, ne, 1
                read(1,*) Self%index_i_KGlobal(i), Self%index_j_KGlobal(i)
            end do
        close(1)
        allocate(Self%value_KGlobal(ne));                Self%value_KGlobal = 0.0d0;
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/KSparseValue.txt',iostat=ios, status="old", action="read")
            do i = 1, ne, 1
                read(1,*) Self%value_KGlobal(i)
            end do
        close(1)
        if (allocated(Self%UGlobal)) then
            Self%UGlobal = 0.0d0;
            Self%value_UGlobal = 0.0d0;
        else
            allocate(Self%UGlobal(Self%n*dim));                Self%UGlobal = 0.0d0;
            allocate(Self%value_UGlobal(n));             Self%value_UGlobal = 0.0d0;
        end if
        ! Sparse solver HSL-Ma87
        write(unit=*, fmt=*) '2. Applying Solver'
        call SparseSolver(Self)
        Self%UGlobal(Self%FreeD) = Self%value_UGlobal
        Displacement = transpose(reshape(Self%UGlobal,[dim,Self%n]))
        call WriteDoublePrecisionArray('NumericalResults/Displacement.txt',Displacement)
        ! deallocating the sparse system, it is not longer necessary
        deallocate(Self%value_FGlobal,Self%value_KGlobal,Self%index_i_KGlobal,Self%index_j_KGlobal)
    end subroutine HSL87Solver
    ! listo
    subroutine ProcessingResults(Self)
        ! this subroutine calculates stresses per element, per node, deformations and strain energy.
        implicit none
        class(Structure), intent(inout)                                     :: Self
        integer                                                             :: i,j
        integer, dimension(:), allocatable                                  :: index
        double precision                                                    :: energy
        double precision, dimension(:), allocatable                         :: Strain, Stress
        double precision, dimension(:,:), allocatable                       :: D
        write(unit=*, fmt=*) '3. Processing Results'
        call ElasticityTensor(Self,D)
        ! checking
        if (Self%DimAnalysis.eq.'2D') then
            allocate(Strain(3),Stress(3));                                          Strain = 0.0d0;         Stress = 0.0d0;
            if(allocated(Self%StrainE).and.allocated(Self%StrainN)) then;     Self%StrainE = 0.0d0;   Self%StrainN = 0.0d0;
            else ; allocate(Self%StrainE(Self%Ne,3),Self%StrainN(Self%N,3));  Self%StrainE = 0.0d0;   Self%StrainN = 0.0d0;
            end if
            if(allocated(Self%StressE).and.allocated(Self%StressE)) then;     Self%StressE = 0.0d0;   Self%StressN = 0.0d0;
            else; allocate(Self%StressE(Self%Ne,3),Self%StressN(Self%N,3));   Self%StressE = 0.0d0;   Self%StressN = 0.0d0;
            end if 
        elseif (Self%DimAnalysis.eq.'3D') then
            allocate(Strain(6),Stress(6));                                             Strain = 0.0;          Stress = 0.0;
            if(allocated(Self%StrainE).and.allocated(Self%StrainN)) then;     Self%StrainE = 0.0;   Self%StrainN = 0.0;
            else ; allocate(Self%StrainE(Self%Ne,6),Self%StrainN(Self%N,6));  Self%StrainE = 0.0;   Self%StrainN = 0.0;
            end if
            if(allocated(Self%StressE).and.allocated(Self%StressN)) then;        Self%StressE = 0.0;    Self%StressN = 0.0;
            else; allocate(Self%StressE(Self%Ne,6),Self%StressN(Self%N,6));      Self%StressE = 0.0;    Self%StressN = 0.0;
            end if 
        end if
        if (allocated(Self%StrainEnergyE).and.allocated(Self%StrainEnergyN)) then
            Self%StrainEnergyE = 0.0d0; Self%StrainEnergyN = 0.0d0;
        else 
            allocate(Self%StrainEnergyE(Self%Ne),Self%StrainEnergyN(Self%N))
            Self%StrainEnergyE = 0.0d0; Self%StrainEnergyN = 0.0d0;
        end if
        ! star processing
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
            energy = 0.5*dot_product(Self%StrainE(i,:),Self%StressE(i,:)) 
            Self%StrainEnergyE(i) = energy
        end do
        ! (node-weighted calculation)
        do i = 1, Self%N, 1
            ! locating the nodes 
            index = pack([(j,j=1,Self%Ne)],any(Self%ConnectivityN.eq.i,2))
            ! 1. Strain
            Strain = sum(Self%StrainE(index,:),1)
            Self%StrainN(i,:) = Strain/size(index)
            ! 2. Stress
            Stress = sum(Self%StressE(index,:),1)
            Self%StressN(i,:) = Stress/size(index)
            ! 3. Energy
            Energy = sum(Self%StrainEnergyE(index))
            Self%StrainEnergyN(i) = Energy/size(index)
        end do
        call WriteDoublePrecisionArray('NumericalResults/StrainE.txt',Self%StrainE)
        call WriteDoublePrecisionArray('NumericalResults/StrainN.txt',Self%StrainN)
        call WriteDoublePrecisionArray('NumericalResults/StressE.txt',Self%StressE)
        call WriteDoublePrecisionArray('NumericalResults/StressN.txt',Self%StressN)
        call WriteDoublePrecisionVector('NumericalResults/EnergyE.txt',Self%StrainEnergyE)
        call WriteDoublePrecisionVector('NumericalResults/EnergyN.txt',Self%StrainEnergyN)
    end subroutine ProcessingResults
end module FEAModule
