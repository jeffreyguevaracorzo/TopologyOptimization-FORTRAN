program MainOptimization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                      UNIVERSIDAD INDUSTRIAL DE SANTANDER                                   !
    !                                    Facultad de ingenierias fisico-mecanicas                                !
    !                                 Programa de doctorado en ingenieria mecanica                               !
    !                                                                                                            !
    !                                                                                                            !
    !         PROGRAMA DE OPTIMIZATION TOPOLOGICA USANDO EL METODO DE ELEMENTOS FINITOS - LINEAL ELASTICO        !
    !                                                   (OPTOP-L)                                                !
    !                                                                                                            !                                           
    !                                       Desarrollado por: Jeffrey Guevara - Universidad Industrial de Sant.  !
    !                                          Asesorado por: Oscar Begambre - Universidad Industrial de Sant.   !
    !                                                                                                            !
    ! Nota: Este software es de código abierto y se proporciona únicamente con fines académicos para ensenanza e !
    ! investigacion. Puede ser usado, modificado y distribuido libremente. No se ofrece ningún tipo de garantía  !
    ! ni soporte técnico para este. Se espera que los usuarios que comprendan y asuman la responsabilidad total  !
    ! de su uso. En caso de ser usado para una investigacion citar el repositorio y el siguiente articulo.       !
    !                                                                                                            !
    ! Importante: El programa OPTOP-L se usa la libreria HSL-LA87, de modo que es necesario tener instaladas     !
    !             las librerias LaPack y BLAS (preferiblemente las versiones optimizadas) junto con la libreria  !
    !             Metis.                                                                                         !                                  
    !                                                                                                            !
    ! Cita:                                                                                                      !
    ! Repositorio:                                                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 1. Librerias usadas
    use Optimization_module
    use Paraview_Module
    implicit none
    type (Optimization)                                  :: OptimizationModel
    ! 2. Ingreso de parametros y propiedades de estructura
    write(unit=*, fmt=*) '1. Structure Properties'
    ! 2.1. Mechanical parameters
    call SetAnalysisType(OptimizationModel,'HookeLaw')               ! Set PlaneStress(2D), PlainStrain(2D), HookeLaw(3D)
    call SetElementType(OptimizationModel,'hexa8')                   ! tria3 tria6, quad4, quad8, tetra4, tetra10, hexa8, hexa20
    call SetThickness(OptimizationModel,0.0d0)                       ! Only 2D
    call SetYoungModulus(OptimizationModel,1000.0d0)                 ! Young modulus
    call SetPoissonModulus(OptimizationModel,0.3d0)                  ! Poisson modulus
    call SetGaussAprox(OptimizationModel,3)                          ! can use up to 5 gauss points
    ! 2.2. Optimization parameters
    call SetVolFraction(OptimizationModel,0.5d0)                     ! limiting volume fraction (lower limit of optimization)
    call SetMutationRate(OptimizationModel,0.2d0)                    ! rate of change/mutation of optimization
    call SetFilterRadius(OptimizationModel,10.0d0*1.5)               ! smoothing filter radius (1.5 times de FE size)
    call SetMaxIterations(OptimizationModel,150)                     ! maximum number of iterations
    call SetPenalFactor(OptimizationModel,5.0d0)                     ! SIMP method penalty factor
    ! 2.3 PostProcessing parameters
    call SetPostProcesingFilter(OptimizationModel,0.2d0)             ! PostProcessing filter (to eliminate intermediate densities)
    ! 3. Lectura de archivos .txt 
    write(unit=*, fmt=*) '2. Reading .txt files'
    call ReadFiles(OptimizationModel)
    write(unit=*, fmt=*) '2.1 Pre Assembly Routine'
    call PreAssemblyRoutine(OptimizationModel)
    ! 4. Applying the topology optimization meth.
    call TopologyOptimizationProcess(OptimizationModel)
    ! 5. Postprocesado con paraview
    write(unit=*, fmt=*) '3. PostProcessing (PARAVIEW)'
    OptimizationModel%ConnectivityN = OptimizationModel%ConnectivityN(OptimizationModel%FinalGeometry,:)
    ! note: leave whatever you want to postprocess and generate the file for viewing with paraview. 
    !       What you don't want to postprocess just comment the line.
    ! 5.1. Densities
    OptimizationModel%FinalDensityVector = OptimizationModel%DensityVector(OptimizationModel%FinalGeometry)
    call DensityParaviewPostProcessing('Paraview','Density',OptimizationModel%Coordinates,OptimizationModel%ConnectivityN, &
    OptimizationModel%ElementType,OptimizationModel%FinalDensityVector)
    ! 5.2 Displacement
    !OptimizationModel%ConnectivityN = OptimizationModel%ConnectivityN(OptimizationModel%FinalGeometry,:)
    !call DisplacementParaviewPostProcessing('Paraview','Displacement',OptimizationModel%Coordinates, &
    !OptimizationModel%ConnectivityN,OptimizationModel%ElementType,OptimizationModel%Displacement)
    ! 5.3 Strain
    !OptimizationModel%ConnectivityN = OptimizationModel%ConnectivityN(OptimizationModel%FinalGeometry,:)
    !OptimizationModel%StrainE = OptimizationModel%StrainE(OptimizationModel%FinalGeometry,:)
    !call StrainParaviewPostProcessing('Paraview','Strain',OptimizationModel%Coordinates,OptimizationModel%ConnectivityN,&
    !OptimizationModel%ElementType,OptimizationModel%StrainN,OptimizationModel%StrainE)
    ! 5.4 Stress
    !OptimizationModel%ConnectivityN = OptimizationModel%ConnectivityN(OptimizationModel%FinalGeometry,:)
    !OptimizationModel%StressE = OptimizationModel%StrainE(OptimizationModel%FinalGeometry,:)
    !call StressParaviewPostProcessing('Paraview','Stress',OptimizationModel%Coordinates,OptimizationModel%ConnectivityN, &
    !OptimizationModel%ElementType,OptimizationModel%StressN,OptimizationModel%StressE)
    ! 5.5 StrainEnergy
    !OptimizationModel%ConnectivityN = OptimizationModel%ConnectivityN(OptimizationModel%FinalGeometry,:)
    !OptimizationModel%StrainEnergyE = OptimizationModel%StrainEnergyE(OptimizationModel%FinalGeometry,:)
    !call StrainEnergyParaviewPostProcessing('Paraview','StrainEnergy',OptimizationModel%Coordinates, &
    !OptimizationModel%ConnectivityN,OptimizationModel%ElementType,OptimizationModel%StrainEnergyN,OptimizationModel%StrainEnergyE)
end program MainOptimization