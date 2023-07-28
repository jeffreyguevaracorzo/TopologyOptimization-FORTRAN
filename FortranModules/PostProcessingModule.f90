module ParaviewModule
    implicit none
    type                                            :: PostprocessingInfo
        character(len=50)                           :: PostprocessingPath
        character(len=50)                           :: Selection
        double precision                            :: FilterValue
    contains
        procedure                                   :: SetPostprocessingPath
        procedure                                   :: SetFilterValue
        procedure                                   :: PlottingSelection
        procedure                                   :: GenerateParaviewFiles
    end type
contains
    ! ----------------------------------------------------------------- !
    !                           base subroutines                        !
    ! ----------------------------------------------------------------- !
    ! this routine generates all files to display the paraview 
    ! results (node information).
    subroutine ParaviewPostProcessingPerNode(Path,ResultName,Coordinates,Connectivity,Element,Result)
        implicit none
        ! input
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: ResultName
        character(len=*), intent(in)                               :: Element
        integer, dimension(:,:), allocatable, intent(in)           :: Connectivity
        double precision, dimension(:), allocatable, intent(in)    :: Result
        double precision, dimension(:,:), allocatable, intent(in)  :: Coordinates
        ! internal 
        integer                                                    :: i, ios
        character(len=100)                                         :: Path1
        character(len=100)                                         :: Path2
        character(len=100)                                         :: Path3
        ! ----------- Process -----------
        Path1 = trim(Path) // '/' // trim(ResultName) // '.case'
        Path2 = trim(Path) // '/' // trim(ResultName) // '.geom'
        Path3 = trim(Path) // '/' // trim(ResultName) // '.esca'
        write(unit=*, fmt=*) Path1
        write(unit=*, fmt=*) Path2
        write(unit=*, fmt=*) Path3
        ! Dimensional check
        if((size(Result)).ne.(size(Coordinates,1))) then; stop "Error dimension Postprocessing"; end if 
        ! Case File        
        open(unit=1, file = Path1, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_casefile"
            write(unit=1,fmt='(a)') 'FORMAT'
            write(unit=1,fmt='(a)') 'type: ensight gold'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') 'GEOMETRY'
            write(unit=1,fmt='(a)') 'model:     ' // ResultName //'.geom'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') 'VARIABLE'
            write(unit=1,fmt='(a)') 'scalar per node:     ' // ResultName // ' ' // ResultName // '.esca'
        close(1)
        ! Geometry File
        open(unit=1, file=Path2, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_geometryfile"
            write(unit=1,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'node id given'
            write(unit=1,fmt='(a)') 'element id given'
            write(unit=1,fmt='(a)') 'extents'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') 'part'
            write(unit=1,fmt='(a)') '1'
            write(unit=1,fmt='(a)') '2D uns-elements (description line for part 1)'
            write(unit=1,fmt='(a)') 'coordinates'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt=*) size(Coordinates,1)
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) i
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,1)
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,2)
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,3)
            end do
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') Element
            write(unit=1,fmt=*) size(Connectivity,1)
            do i = 1, size(Connectivity,1), 1
                write(unit=1,fmt=*) i
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Connectivity,1), 1
                write(unit=1,fmt=*) Connectivity(i,:)
            end do
        close(1)
        ! Scalar File
        open(unit=1, file=Path3, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_coordinatesfile"
            write(unit=1,fmt='(a)') 'Scalar File'
            write(unit=1,fmt='(a)') 'part'
            write(unit=1,fmt='(a)') '1'
            write(unit=1,fmt='(a)') 'coordinates'
            do i = 1, size(Result), 1
                write(unit=1,fmt=*) Result(i)
            end do
        close(1)
    end subroutine ParaviewPostProcessingPerNode
    ! this routine generates all files to display the paraview 
    ! results (element information).
    subroutine ParaviewPostProcessingPerElement(Path,ResultName,Coordinates,Connectivity,Element,Result)
        implicit none
        ! input
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: ResultName
        character(len=*), intent(in)                               :: Element
        integer, dimension(:,:), allocatable, intent(in)           :: Connectivity
        double precision, dimension(:), allocatable, intent(in)    :: Result
        double precision, dimension(:,:), allocatable, intent(in)  :: Coordinates
        ! internal 
        integer                                                    :: i, ios
        character(len=100)                                         :: Path1
        character(len=100)                                         :: Path2
        character(len=100)                                         :: Path3
        ! ----------- Process -----------
        Path1 = trim(Path) // '/' // trim(ResultName) // '.case'
        Path2 = trim(Path) // '/' // trim(ResultName) // '.geom'
        Path3 = trim(Path) // '/' // trim(ResultName) // '.esca'
        write(unit=*, fmt=*) Path1
        write(unit=*, fmt=*) Path2
        write(unit=*, fmt=*) Path3
        ! Dimensional check
        if((size(Result)).ne.(size(Connectivity,1))) then; stop "Error dimension Postprocessing"; end if 
        ! Case File
        open(unit=1, file=Path1, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_casefile"
            write(unit=1,fmt='(a)') 'FORMAT'
            write(unit=1,fmt='(a)') 'type: ensight gold'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') 'GEOMETRY'
            write(unit=1,fmt='(a)') 'model:     ' // ResultName // '.geom'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') 'VARIABLE'
            write(unit=1,fmt='(a)') 'scalar per element:     ' // ResultName // ' ' // ResultName // '.esca'
        close(1)
        ! Geometry File
        open(unit=1, file=Path2, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_geometryfile"
            write(unit=1,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'This is the 2nd description line of the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'node id given'
            write(unit=1,fmt='(a)') 'element id given'
            write(unit=1,fmt='(a)') 'extents'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') 'part'
            write(unit=1,fmt='(a)') '1'
            write(unit=1,fmt='(a)') '2D uns-elements (description line for part 1)'
            write(unit=1,fmt='(a)') 'coordinates'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt=*) size(Coordinates,1)
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) i
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,1)
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,2)
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,3)
            end do
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') Element
            write(unit=1,fmt=*) size(Connectivity,1)
            do i = 1, size(Connectivity,1), 1
                write(unit=1,fmt=*) i
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Connectivity,1), 1
                write(unit=1,fmt=*) Connectivity(i,:)
            end do
        close(1)
        ! Scalar File
        open(unit=1, file=Path3, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_coordinatesfile"
            write(unit=1,fmt='(a)') 'Per_element scalar values for the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'part'
            write(unit=1,fmt='(a)') '1'
            write(unit=1,fmt='(a)') Element
            do i = 1, size(Result), 1
                write(unit=1,fmt=*) Result(i)
            end do
        close(1)
    end subroutine ParaviewPostProcessingPerElement
    ! ----------------------------------------------------------------- !
    !          subroutines to define the information required           !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine SetPostprocessingPath(self,path)
        implicit none
        character(len=*), intent(in)                :: path
        class(PostprocessingInfo), intent(inout)    :: self
        self%PostprocessingPath = path
    end subroutine
    ! listo
    subroutine PlottingSelection(self,info)
        implicit none
        character(len=*), intent(in)                :: info
        class(PostprocessingInfo), intent(inout)    :: self 
        self%Selection = info
    end subroutine PlottingSelection
    ! listo
    subroutine SetFilterValue(self,FilterValue)
        implicit none
        double precision, intent(in)                :: FilterValue
        class(PostprocessingInfo), intent(inout)    :: self 
        Self%FilterValue = FilterValue
    end subroutine SetFilterValue
    ! ----------------------------------------------------------------- !
    !             subroutines to define the Postprocessing              !
    ! ----------------------------------------------------------------- !
    ! pendiente
    subroutine GenerateParaviewFiles(self1,self2)
        use TopologyOptimizationModule
        implicit none
        class(PostprocessingInfo), intent(inout)        :: self1
        class(Optimization), intent(inout)              :: Self2
        ! internal values
        integer                                         :: i,j,dim
        integer, dimension(:), allocatable              :: Index
        integer, dimension(:,:), allocatable            :: Connectivity
        double precision, dimension(:,:), allocatable   :: Coordinates
        double precision, dimension(:), allocatable     :: ResultVector
        logical, dimension(:), allocatable              :: InLogical
        if (self2%DimAnalysis.eq.'2D') then; dim = 2; end if
        if (self2%DimAnalysis.eq.'3D') then; dim = 3; end if
        ! ---------------------------------------------------------------!
        ! Applying filter (deleting elements with low density)
        Coordinates = Self2%Coordinates
        j = size(Self2%DensityVector)
        InLogical = Self2%DensityVector.gt.self1%FilterValue
        Index = Pack([(i,i=1,j)],InLogical)
        Connectivity = Self2%ConnectivityN(Index,:)
        ! ---------------------------------------------------------------!
        ! star postprocessing
        if (self1%Selection.eq.'Stress') then !listo
            select case (dim)
                case (2)
                    allocate(ResultVector(Self2%Ne))
                    ResultVector = self2%StressE(Index,1)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressXX',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressYY',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,3)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressXY',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                case (3)
                    allocate(ResultVector(Self2%Ne))
                    ResultVector = self2%StressE(Index,1)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressXX',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressYY',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,3)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressZZ',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,4)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressXY',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,5)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressYZ',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StressE(Index,6)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStressXZ',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
            end select
        elseif (self1%Selection.eq.'StressVonMises') then !listo
            select case (dim)
                case (2)
                    ResultVector = sqrt(self2%StressE(Index,1)**2 - self2%StressE(Index,1)*self2%StressE(:,2) &
                                        + self2%StressE(Index,2)**2 + 3*self2%StressE(Index,3)**2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'StressVonMises',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
                case (3)
                    ResultVector = sqrt(self2%StressE(Index,1)**2 + self2%StressE(Index,2)**2 + self2%StressE(Index,3)**2 &
                                     - (self2%StressE(Index,1)*self2%StressE(Index,2) &
                                      + self2%StressE(Index,2)*self2%StressE(Index,3) &
                                      + self2%StressE(Index,1)*self2%StressE(Index,3)) + 3*(self2%StressE(Index,4)**2 &
                                      + self2%StressE(Index,5)**2 + self2%StressE(Index,6)**2))
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'StressVonMises',Coordinates, &
                                                        Connectivity,self2%ElementType,ResultVector)
            end select
        elseif (self1%Selection.eq.'Strain') then !listo
            select case (dim)
                case (2)
                    allocate(ResultVector(Self2%Ne))
                    ResultVector = self2%StrainE(Index,1)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainXX',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainYY',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,3)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainXY',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                case (3)
                    allocate(ResultVector(Self2%Ne))
                    ResultVector = self2%StrainE(Index,1)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainXX',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainYY',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,3)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainZZ',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,4)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainXY',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,5)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainYZ',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                    ResultVector = self2%StrainE(Index,6)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainXZ',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
            end select
        elseif (self1%selection.eq.'StrainEqv') then !list
            select case (dim)
                case (2)
                    ResultVector = sqrt(self2%StrainE(Index,1)**2 + self2%StrainE(Index,2)**2 + 2*self2%StrainE(Index,3)**2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainEqv',Coordinates, &
                                                          Connectivity,self2%ElementType,ResultVector)
                case (3)
                    ResultVector = sqrt(self2%StrainE(Index,1)**2 + self2%StrainE(Index,2)**2 + self2%StrainE(Index,3)**2 &
                                    + 2*self2%StrainE(Index,4)**2 + 2*self2%StrainE(Index,5)**2 + 2*self2%StrainE(Index,6)**2)
                    call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EStrainEqv',Coordinates, &
                                                           Connectivity,self2%ElementType,ResultVector)
            end select
        elseif (self1%Selection.eq.'DensEnergy') then !listo
            ResultVector = self2%StrainEnergyE(Index)
            call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'EEnergy',Coordinates, &
                                                  Connectivity,self2%ElementType,ResultVector)
        elseif (Self1%Selection.eq.'TopOpt') then
            ResultVector = Self2%DensityVector(Index)
            call ParaviewPostProcessingPerElement(self1%PostprocessingPath,'DensityDistribution',Coordinates, &
                                                  Connectivity,self2%ElementType,ResultVector)
        end if
    end subroutine GenerateParaviewFiles
end module ParaviewModule