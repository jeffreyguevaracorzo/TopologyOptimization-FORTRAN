module Paraview_Module
    implicit none
    type PostprocessingInfo
        character(len=15)                           :: PostprocessingPath
        character(len=15)                           :: Selection
    contains 
        procedure                                   :: PlottingSelection
    end type PostprocessingInfo
    ! interfaces
    interface DensityParaviewPostProcessing
        module procedure DensityParaviewPostProcessingReal
        module procedure DensityParaviewPostProcessingDP
    end interface DensityParaviewPostProcessing
    interface DisplacementParaviewPostProcessing
        module procedure DisplacementParaviewPostProcessingReal
        module procedure DisplacementParaviewPostProcessingDP
    end interface DisplacementParaviewPostProcessing
    interface StrainParaviewPostProcessing
        module procedure StrainParaviewPostProcessingReal
        module procedure StrainParaviewPostProcessingDP
    end interface StrainParaviewPostProcessing
    interface StressParaviewPostProcessing
        module procedure StressParaviewPostProcessingReal
        module procedure StressParaviewPostProcessingDP
    end interface StressParaviewPostProcessing
    interface StrainEnergyParaviewPostProcessing
        module procedure StrainEnergyParaviewPostProcessingReal
        module procedure StrainEnergyParaviewPostProcessingDP
    end interface StrainEnergyParaviewPostProcessing
    contains
    ! ** base procedures
    ! Input postprocessing type
    subroutine PlottingSelection(self,info)
        implicit none
        character(len=*), intent(in)                :: info
        class(PostprocessingInfo), intent(inout)    :: self 
        self%Selection = info
    end subroutine PlottingSelection
    ! Generate postprocessing files
    subroutine GetParaviewFiles(Self1,Self2)
        use Optimization_module
        implicit none
        class(Optimization), intent(inout)          :: Self1
        class(PostprocessingInfo), intent(inout)    :: Self2
        Self1%ConnectivityN = Self1%ConnectivityN(Self1%FinalGeometry,:)
        if(Self2%Selection.eq.'DensityDist') then
            Self1%FinalDensityVector = Self1%DensityVector(Self1%FinalGeometry)
            call DensityParaviewPostProcessing('Paraview','Density',Self1%Coordinates,Self1%ConnectivityN,Self1%ElementType,&
                Self1%FinalDensityVector)
        elseif(Self2%Selection.eq.'Displacement') then
            call DisplacementParaviewPostProcessing('Paraview','Displacement',Self1%Coordinates,Self1%ConnectivityN, &
                Self1%ElementType,Self1%Displacement)
        elseif(Self2%Selection.eq.'Strain') then
            Self1%StrainE = Self1%StrainE(Self1%FinalGeometry,:)
            call StrainParaviewPostProcessing('Paraview','Strain',Self1%Coordinates,Self1%ConnectivityN,Self1%ElementType, &
                Self1%StrainE)
        elseif(Self2%Selection.eq.'Stress') then
            Self1%StressE = Self1%StressE(Self1%FinalGeometry,:)
            call StressParaviewPostProcessing('Paraview','Stress',Self1%Coordinates,Self1%ConnectivityN,Self1%ElementType, &
                Self1%StressE)
        elseif(Self2%Selection.eq.'StrainEnergy') then
            Self1%StrainEnergyE = Self1%StrainEnergyE(Self1%FinalGeometry)
            call StrainEnergyParaviewPostProcessing('Paraview','StrainEnergy',Self1%Coordinates,Self1%ConnectivityN, &
                Self1%ElementType,Self1%StrainEnergyE)
        end if
    end subroutine GetParaviewFiles
    ! --------------- ADDITIONAL BASE FUNCTIONS AND SUBROUTINES ---------------
    Subroutine GenerateGeometryParaviewFileReal(Path,Element,DimAnalysis,Connect,Coord)
        implicit none
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        integer, intent(in)                                        :: DimAnalysis
        integer, dimension(:,:), allocatable, intent(in)           :: Connect
        real, dimension(:,:), allocatable, intent(in)              :: Coord
        ! internal variables
        integer                                                    :: i,ios,iounit
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_geometryfile"
            write(unit=iounit,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=iounit,fmt='(a)') 'This is the 2st description line of the EnSight Gold geometry example'
            write(unit=iounit,fmt='(a)') 'node id given'
            write(unit=iounit,fmt='(a)') 'element id given'
            write(unit=iounit,fmt='(a)') 'extents'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') 'part'
            write(unit=iounit,fmt='(a)') '1'
            if (DimAnalysis.eq.2) then 
                write(unit=iounit,fmt='(a)') '2D uns-elements (description line for part 1)'
            elseif (DimAnalysis.eq.3) then
                write(unit=iounit,fmt='(a)') '3D uns-elements (description line for part 1)'
            else
                stop "ERROR DimAnalysis generating paraview geometry file"
            end if
            write(unit=iounit,fmt='(a)') 'coordinates'
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt=*) size(Coord,1)
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) i
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) Coord(i,1)
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) Coord(i,2)
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                if (DimAnalysis.eq.2) then 
                    write(unit=iounit,fmt=*) 0.0d0
                elseif (DimAnalysis.eq.3) then
                    write(unit=iounit,fmt=*) Coord(i,3)
                else
                    stop "ERROR DimAnalysis generating paraview geometry file"
                end if
            end do
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt='(a)') Element
            write(unit=iounit,fmt=*) size(Connect,1)
            do i = 1, size(Connect,1), 1
                write(unit=iounit,fmt=*) i
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Connect,1), 1
                write(unit=iounit,fmt=*) Connect(i,:)
            end do
        close(iounit)
    end subroutine GenerateGeometryParaviewFileReal
    subroutine GenerateGeometryParaviewFileDP(Path,Element,DimAnalysis,Connect,Coord)
        implicit none
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        integer, intent(in)                                        :: DimAnalysis
        integer, dimension(:,:), allocatable, intent(in)           :: Connect
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord
        ! internal variables
        integer                                                    :: i,ios,iounit
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_geometryfile"
            write(unit=iounit,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=iounit,fmt='(a)') 'This is the 2st description line of the EnSight Gold geometry example'
            write(unit=iounit,fmt='(a)') 'node id given'
            write(unit=iounit,fmt='(a)') 'element id given'
            write(unit=iounit,fmt='(a)') 'extents'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') 'part'
            write(unit=iounit,fmt='(a)') '1'
            if (DimAnalysis.eq.2) then 
                write(unit=iounit,fmt='(a)') '2D uns-elements (description line for part 1)'
            elseif (DimAnalysis.eq.3) then
                write(unit=iounit,fmt='(a)') '3D uns-elements (description line for part 1)'
            else
                stop "ERROR DimAnalysis generating paraview geometry file"
            end if
            write(unit=iounit,fmt='(a)') 'coordinates'
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt=*) size(Coord,1)
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) i
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) Coord(i,1)
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) Coord(i,2)
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                if (DimAnalysis.eq.2) then 
                    write(unit=iounit,fmt=*) 0.0d0
                elseif (DimAnalysis.eq.3) then
                    write(unit=iounit,fmt=*) Coord(i,3)
                else
                    stop "ERROR DimAnalysis generating paraview geometry file"
                end if
            end do
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt='(a)') Element
            write(unit=iounit,fmt=*) size(Connect,1)
            do i = 1, size(Connect,1), 1
                write(unit=iounit,fmt=*) i
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Connect,1), 1
                write(unit=iounit,fmt=*) Connect(i,:)
            end do
        close(iounit)
    end subroutine GenerateGeometryParaviewFileDP
    subroutine GenerateEscalarParaviewFileReal(Path,Element,Result)
        implicit none
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        real, dimension(:), allocatable, intent(in)                :: Result
        ! internal variables
        integer                                                    :: i,ios,iounit
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_Coordfile"
            write(unit=iounit,fmt='(a)') 'Scalar File'
            write(unit=iounit,fmt='(a)') 'part'
            write(unit=iounit,fmt='(a)') '1'
            write(unit=iounit,fmt='(a)') Element
            do i = 1, size(Result), 1
                write(unit=iounit,fmt=*) Result(i)
            end do
        close(iounit)
    end subroutine GenerateEscalarParaviewFileReal
    subroutine GenerateEscalarParaviewFileDP(Path,Element,Result)
        implicit none
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        double precision, dimension(:), allocatable, intent(in)    :: Result
        ! internal variables
        integer                                                    :: i,ios,iounit
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_Coordfile"
            write(unit=iounit,fmt='(a)') 'Scalar File'
            write(unit=iounit,fmt='(a)') 'part'
            write(unit=iounit,fmt='(a)') '1'
            write(unit=iounit,fmt='(a)') Element
            do i = 1, size(Result), 1
                write(unit=iounit,fmt=*) Result(i)
            end do
        close(iounit)
    end subroutine GenerateEscalarParaviewFileDP
    ! ------------ FINITE ELEMENT ANALYSISS FUNCTIONS AND SUBROUTINES ------------
    subroutine DensityParaviewPostProcessingReal(Path,ResultName,Coord,Connect,Element,Density)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        real, dimension(:), allocatable, intent(in)                :: Density           ! Result variable
        real, dimension(:,:), allocatable, intent(in)              :: Coord             ! Coordinates
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        ! internal variables
        real, dimension(:), allocatable                            :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,2,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '.esca'
            Result = Density(:)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,3,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '.esca'
            Result = Density(:)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
        else
            stop "ERROR Paraview - Geometry -"
        end if
    end subroutine DensityParaviewPostProcessingReal
    subroutine DensityParaviewPostProcessingDP(Path,ResultName,Coord,Connect,Element,Density)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        double precision, dimension(:), allocatable, intent(in)    :: Density           ! Result variable
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord             ! Coordinates
        ! internal variables
        double precision, dimension(:), allocatable                :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,2,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '.esca'
            Result = Density(:)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,3,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '.esca'
            Result = Density(:)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
        else
            stop "ERROR Paraview - Geometry -"
        end if
    end subroutine DensityParaviewPostProcessingDP
    subroutine DisplacementParaviewPostProcessingReal(Path,ResultName,Coord,Connect,Element,Displacement)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        real, dimension(:,:), allocatable, intent(in)              :: Displacement      ! Result variable
        real, dimension(:,:), allocatable, intent(in)              :: Coord             ! Coordinates
        ! internal variables
        real, dimension(:), allocatable                            :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '_X' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_X' // ' ' // Path1
                Path1 = ResultName // '_Y' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Y' // ' ' // Path1
                Path1 = ResultName // '_Eq' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Eq' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,2,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '_X' //  '.esca'
            Result = Displacement(:,1)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Y' //  '.esca'
            Result = Displacement(:,2)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Eq' //  '.esca'
            Result = sqrt(Displacement(:,1)**2+Displacement(:,2)**2)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '_X' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_X' // ' ' // Path1
                Path1 = ResultName // '_Y' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Y' // ' ' // Path1
                Path1 = ResultName // '_Z' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Z' // ' ' // Path1
                Path1 = ResultName // '_Eq' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Eq' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,3,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '_X' //  '.esca'
            Result = Displacement(:,1)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Y' //  '.esca'
            Result = Displacement(:,2)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Z' //  '.esca'
            Result = Displacement(:,3)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Eq' //  '.esca'
            Result = sqrt(Displacement(:,1)**2 + Displacement(:,2)**2 + Displacement(:,3)**2)
            call GenerateEscalarParaviewFileReal(Path1,'coordinates',Result)
        else
            stop "ERROR Paraview - Geometry -"
        end if
    end subroutine DisplacementParaviewPostProcessingReal
    subroutine DisplacementParaviewPostProcessingDP(Path,ResultName,Coord,Connect,Element,Displacement)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        double precision, dimension(:,:), allocatable, intent(in)  :: Displacement      ! Result variable
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord             ! Coordinates
        ! internal variables
        double precision, dimension(:), allocatable                :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '_X' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_X' // ' ' // Path1
                Path1 = ResultName // '_Y' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Y' // ' ' // Path1
                Path1 = ResultName // '_Eq' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Eq' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,2,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '_X' //  '.esca'
            Result = Displacement(:,1)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Y' //  '.esca'
            Result = Displacement(:,2)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Eq' //  '.esca'
            Result = sqrt(Displacement(:,1)**2+Displacement(:,2)**2)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // '_X' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_X' // ' ' // Path1
                Path1 = ResultName // '_Y' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Y' // ' ' // Path1
                Path1 = ResultName // '_Z' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Z' // ' ' // Path1
                Path1 = ResultName // '_Eq' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per node:   ' // ResultName // '_Eq' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,3,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // '_X' //  '.esca'
            Result = Displacement(:,1)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Y' //  '.esca'
            Result = Displacement(:,2)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Z' //  '.esca'
            Result = Displacement(:,3)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
            Path1 = Path // '/' // ResultName // '_Eq' //  '.esca'
            Result = sqrt(Displacement(:,1)**2 + Displacement(:,2)**2 + Displacement(:,3)**2)
            call GenerateEscalarParaviewFileDP(Path1,'coordinates',Result)
        else
            stop "ERROR Paraview - Geometry -"
        end if
    end subroutine DisplacementParaviewPostProcessingDP
    subroutine StrainParaviewPostProcessingReal(Path,ResultName,Coord,Connect,Element,StrainElement)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        real, dimension(:,:), allocatable, intent(in)              :: StrainElement     ! Result variable
        real, dimension(:,:), allocatable, intent(in)              :: Coord             ! Coordinates
        ! internal variables
        real, dimension(:), allocatable                            :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,2,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StrainElement(:,1)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StrainElement(:,2)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StrainElement(:,3)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt(StrainElement(:,1)**2+StrainElement(:,2)**2-StrainElement(:,1)*StrainElement(:,2) &
                    +3*StrainElement(:,3)**2)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_zz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_zz' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_yz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yz' // ' ' // Path1
                Path1 = ResultName // 'E_xz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xz' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,3,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StrainElement(:,1)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StrainElement(:,2)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_zz' //  '.esca'
            Result = StrainElement(:,3)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StrainElement(:,4)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yz' //  '.esca'
            Result = StrainElement(:,5)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xz' //  '.esca'
            Result = StrainElement(:,6)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt((StrainElement(:,1)-StrainElement(:,2))**2+(StrainElement(:,2)-StrainElement(:,3))**2 &
                        +(StrainElement(:,3)-StrainElement(:,1))**2+3*(StrainElement(:,4)**2+StrainElement(:,5)**2 &
                        +StrainElement(:,6)*2))
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
        else
            stop "ERROR Paraview - Strain -"
        end if
    end subroutine StrainParaviewPostProcessingReal
    subroutine StrainParaviewPostProcessingDP(Path,ResultName,Coord,Connect,Element,StrainElement)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        double precision, dimension(:,:), allocatable, intent(in)  :: StrainElement     ! Result variable
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord             ! Coordinates
        ! internal variables
        double precision, dimension(:), allocatable                :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,2,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StrainElement(:,1)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StrainElement(:,2)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StrainElement(:,3)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt(StrainElement(:,1)**2+StrainElement(:,2)**2-StrainElement(:,1)*StrainElement(:,2)&
                    +3*StrainElement(:,3)**2)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_zz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_zz' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_yz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yz' // ' ' // Path1
                Path1 = ResultName // 'E_xz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xz' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,3,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StrainElement(:,1)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StrainElement(:,2)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_zz' //  '.esca'
            Result = StrainElement(:,3)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StrainElement(:,4)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yz' //  '.esca'
            Result = StrainElement(:,5)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xz' //  '.esca'
            Result = StrainElement(:,6)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt((StrainElement(:,1)-StrainElement(:,2))**2+(StrainElement(:,2)-StrainElement(:,3))**2 &
                    +(StrainElement(:,3)-StrainElement(:,1))**2+3*(StrainElement(:,4)**2+StrainElement(:,5)**2+ &
                    StrainElement(:,6)*2))
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
        else
            stop "ERROR Paraview - Strain -"
        end if
    end subroutine StrainParaviewPostProcessingDP
    subroutine StressParaviewPostProcessingReal(Path,ResultName,Coord,Connect,Element,StressElement)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        real, dimension(:,:), allocatable, intent(in)              :: StressElement     ! Result variable
        real, dimension(:,:), allocatable, intent(in)              :: Coord             ! Coordinates
        ! internal variables
        real, dimension(:), allocatable                            :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,2,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StressElement(:,1)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StressElement(:,2)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StressElement(:,3)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt(StressElement(:,1)**2+StressElement(:,2)**2-StressElement(:,1)*StressElement(:,2) &
                    +3*StressElement(:,3)**2)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_zz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_zz' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_yz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yz' // ' ' // Path1
                Path1 = ResultName // 'E_xz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xz' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,3,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StressElement(:,1)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StressElement(:,2)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_zz' //  '.esca'
            Result = StressElement(:,3)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StressElement(:,4)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yz' //  '.esca'
            Result = StressElement(:,5)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xz' //  '.esca'
            Result = StressElement(:,6)
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt((StressElement(:,1)-StressElement(:,2))**2+(StressElement(:,2)-StressElement(:,3))**2 &
                +(StressElement(:,3)-StressElement(:,1))**2+3*(StressElement(:,4)**2+StressElement(:,5)**2+StressElement(:,6)*2))
            call GenerateEscalarParaviewFileReal(Path1,Element,Result)
        else
            stop "ERROR Paraview - Stress -"
        end if
    end subroutine StressParaviewPostProcessingReal
    subroutine StressParaviewPostProcessingDP(Path,ResultName,Coord,Connect,Element,StressElement)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        double precision, dimension(:,:), allocatable, intent(in)  :: StressElement     ! Result variable
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord             ! Coordinates
        ! internal variables
        double precision, dimension(:), allocatable                :: Result
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,2,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StressElement(:,1)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StressElement(:,2)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StressElement(:,3)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt(StressElement(:,1)**2+StressElement(:,2)**2-StressElement(:,1)*StressElement(:,2) &
                    +3*StressElement(:,3)**2)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                ! per element
                Path1 = ResultName // 'E_xx' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xx' // ' ' // Path1
                Path1 = ResultName // 'E_yy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yy' // ' ' // Path1
                Path1 = ResultName // 'E_zz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_zz' // ' ' // Path1
                Path1 = ResultName // 'E_xy' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xy' // ' ' // Path1
                Path1 = ResultName // 'E_yz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_yz' // ' ' // Path1
                Path1 = ResultName // 'E_xz' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_xz' // ' ' // Path1
                Path1 = ResultName // 'E_VonMises' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E_VonMises' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,3,Connect,Coord)
            ! Escalar files
            ! - per element
            Path1 = Path // '/' // ResultName // 'E_xx' //  '.esca'
            Result = StressElement(:,1)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yy' //  '.esca'
            Result = StressElement(:,2)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_zz' //  '.esca'
            Result = StressElement(:,3)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xy' //  '.esca'
            Result = StressElement(:,4)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_yz' //  '.esca'
            Result = StressElement(:,5)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_xz' //  '.esca'
            Result = StressElement(:,6)
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
            Path1 = Path // '/' // ResultName // 'E_VonMises' //  '.esca'
            Result = sqrt((StressElement(:,1)-StressElement(:,2))**2+(StressElement(:,2)-StressElement(:,3))**2 &
                +(StressElement(:,3)-StressElement(:,1))**2+3*(StressElement(:,4)**2+StressElement(:,5)**2+StressElement(:,6)*2))
            call GenerateEscalarParaviewFileDP(Path1,Element,Result)
        else
            stop "ERROR Paraview - Stress -"
        end if
    end subroutine StressParaviewPostProcessingDP
    subroutine StrainEnergyParaviewPostProcessingReal(Path,ResultName,Coord,Connect,Element,StrainEnergyE)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        real, dimension(:), allocatable, intent(in)                :: StrainEnergyE
        real, dimension(:,:), allocatable, intent(in)              :: Coord             ! Coordinates
        ! internal variables
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // 'E' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,2,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // 'E' //  '.esca'
            call GenerateEscalarParaviewFileReal(Path1,Element,StrainEnergyE)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // 'E' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileReal(Path1,Element,3,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // 'E' //  '.esca'
            call GenerateEscalarParaviewFileReal(Path1,Element,StrainEnergyE)
        else
            stop "ERROR Paraview - Geometry -"
        end if
    end subroutine StrainEnergyParaviewPostProcessingReal
    subroutine StrainEnergyParaviewPostProcessingDP(Path,ResultName,Coord,Connect,Element,StrainEnergyE)
        implicit none
        character(len=*), intent(in)                               :: Path              ! Folder name
        character(len=*), intent(in)                               :: ResultName        ! Result name
        character(len=*), intent(in)                               :: Element           ! Element type
        integer, dimension(:,:), allocatable, intent(in)           :: Connect           ! conectivity
        double precision, dimension(:), allocatable, intent(in)    :: StrainEnergyE
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord             ! Coordinates
        ! internal variables
        character(len=100)                                         :: Path1,CaseName
        integer                                                    :: i,ios,iounit
        if (size(Coord,2).eq.2) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // 'E' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,2,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // 'E' //  '.esca'
            call GenerateEscalarParaviewFileDP(Path1,Element,StrainEnergyE)
        elseif (size(Coord,2).eq.3) then
            ! case file (principal file)
            CaseName =  Path // '/' // ResultName //'.case'
            open(unit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=iounit,fmt='(a)') 'FORMAT'
                write(unit=iounit,fmt='(a)') 'type: ensight gold'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'GEOMETRY'
                write(unit=iounit,fmt='(a)') 'model:     ' // ResultName //'.geom'
                write(unit=iounit,fmt='(a)') ''
                write(unit=iounit,fmt='(a)') 'VARIABLE'
                Path1 = ResultName // 'E' //  '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // ResultName // 'E' // ' ' // Path1
            close(iounit)
            ! Geometry File
            Path1 = Path // '/' // ResultName //'.geom'
            call GenerateGeometryParaviewFileDP(Path1,Element,3,Connect,Coord)
            ! Escalar files
            Path1 = Path // '/' // ResultName // 'E' //  '.esca'
            call GenerateEscalarParaviewFileDP(Path1,Element,StrainEnergyE)
        else
            stop "ERROR Paraview - Geometry -"
        end if
    end subroutine StrainEnergyParaviewPostProcessingDP
end module Paraview_Module